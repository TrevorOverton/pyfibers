"""Scan drive voltage and refine the first observed activation boundary per tau.

Rectangular control is always included. All coarse trials are retained so
nonmonotonic responses (including loss of detection at high drive) are visible.
A finite scan cannot exclude activation windows between sampled amplitudes.
"""
import argparse
import csv
import multiprocessing as mp
import os
import time
import traceback
import uuid
from pathlib import Path

import numpy as np
from ansys_mrg import Experiment, Pulse, add_common_arguments, config_from_args, write_json


def accepted(summary, criterion="distal", include_offset=True, allow_end=False):
    spike = summary["any_active_node_spike"] if criterion == "any" else summary["detector_spike_count"] > 0
    return bool(spike and (allow_end or not summary["end_initiation_flag"])
                and (include_offset or summary["first_spike_phase"] == "during_first_phase"))


def scan_boundary(evaluate, magnitudes, relative_tolerance=0.01, max_bisections=30):
    """evaluate returns bool. Refine first sampled false->true bracket only."""
    if evaluate(0.0):
        raise RuntimeError("Zero-drive trial meets activation criterion. Inspect baseline/model.")
    states = [bool(evaluate(float(a))) for a in magnitudes]
    reversals = any(a and not b for a, b in zip(states[:-1], states[1:]))
    if not any(states):
        return {"status": "no_activation_observed_on_grid", "lower_v": None, "upper_v": None,
                "coarse_nonmonotonic": reversals}
    i = states.index(True)
    low, high = (0.0 if i == 0 else float(magnitudes[i - 1])), float(magnitudes[i])
    for _ in range(max_bisections):
        if high - low <= relative_tolerance * high:
            break
        mid = (low + high) / 2
        if evaluate(mid):
            high = mid
        else:
            low = mid
    return {"status": "bracketed" if high-low <= relative_tolerance*high else "iteration_limit",
            "lower_v": low, "upper_v": high, "coarse_nonmonotonic": reversals}



TRIAL_FIELDS = ["waveform", "tau_ms", "polarity", "magnitude_v", "accepted",
                "first_active_spike_ms", "first_detector_spike_ms", "detector_spike_count",
                "first_spike_phase", "end_initiation_flag", "max_vm_mv"]
THRESHOLD_FIELDS = ["waveform", "tau_ms", "polarity", "normalization", "criterion", "status",
                    "lower_v", "upper_v", "coarse_nonmonotonic", "signed_upper_drive_v",
                    "upper_first_phase_peak_abs_v"]


def run_case(task):
    """Top-level spawn-safe worker. No NEURON objects cross process boundaries.

    Each process handles only one complete amplitude search before exiting,
    so sections, recorders, and solver state cannot leak into another case.
    """
    index, options, polarity, tau, case_path = task
    a = argparse.Namespace(**options)
    directory = Path(case_path)
    directory.mkdir(parents=True, exist_ok=True)
    kind = "rectangular" if tau is None else "exponential"
    label = f"case {index + 1}: {kind}, tau={tau}, polarity={polarity:+d}"
    try:
        exp = Experiment(config_from_args(a))
        write_json(directory / "configuration.json", {**exp.metadata, "arguments": options,
                   "worker_pid": os.getpid(), "case_index": index, "tau_ms": tau, "polarity": polarity})
        pulse = Pulse(kind, a.onset, a.duration, 0.1 if tau is None else tau,
                      a.normalization, a.fall_tau)
        pulse.validate()
        grid = np.geomspace(a.min_voltage, a.max_voltage, a.scan_points)
        cache = {}
        with (directory / "trials.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=TRIAL_FIELDS)
            writer.writeheader(); f.flush()

            def evaluate(magnitude):
                if magnitude in cache:
                    return cache[magnitude]
                result = exp.run(polarity * magnitude, pulse, a.dt, a.tstop)
                summary = result["summary"]
                ok = accepted(summary, a.criterion, not a.first_phase_only, a.allow_end)
                cache[magnitude] = ok
                row = {"waveform": kind, "tau_ms": tau, "polarity": polarity,
                       "magnitude_v": magnitude, "accepted": ok}
                row.update({k: summary[k] for k in TRIAL_FIELDS if k in summary})
                writer.writerow(row); f.flush()
                print(f"[{label}] {magnitude:.6g} V: accepted={ok}, "
                      f"distal={summary['detector_spike_count']}, "
                      f"end={summary['end_initiation_flag']}", flush=True)
                return ok

            boundary = scan_boundary(evaluate, grid, a.rtol)
        row = {"waveform": kind, "tau_ms": tau, "polarity": polarity,
               "normalization": a.normalization, "criterion": a.criterion, **boundary}
        high = boundary["upper_v"]
        row["signed_upper_drive_v"] = None if high is None else polarity * high
        row["upper_first_phase_peak_abs_v"] = None if high is None else high * pulse.peak_fraction
        write_json(directory / "threshold.json", row)
        if high is not None:
            result = exp.run(polarity * high, pulse, a.dt, a.tstop)
            exp.save(result, directory)
        return index, row, exp.metadata
    except Exception:
        (directory / "error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        raise RuntimeError(f"{label} failed. Details: {directory / 'error.txt'}")


def write_combined(output, tasks, completed, include_partial=False):
    """Only the parent writes shared summaries; stable order independent of timing."""
    temp = output / "thresholds.csv.tmp"
    with temp.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=THRESHOLD_FIELDS)
        writer.writeheader()
        writer.writerows(completed[i] for i in sorted(completed))
    temp.replace(output / "thresholds.csv")
    temp = output / "trials.csv.tmp"
    with temp.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=TRIAL_FIELDS)
        writer.writeheader()
        for index, _, _, _, case_path in tasks:
            source = Path(case_path) / "trials.csv"
            if (index in completed or include_partial) and source.exists():
                with source.open(newline="") as cf:
                    for row in csv.DictReader(cf):
                        # A hard interruption can leave an incomplete last line.
                        if None not in row and all(v is not None for v in row.values()):
                            writer.writerow(row)
    temp.replace(output / "trials.csv")


def plot_thresholds(rows, a, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    for polarity in a.polarities:
        subset = [r for r in rows if r["polarity"] == polarity and r["waveform"] == "exponential"]
        subset.sort(key=lambda r: r["tau_ms"])
        line, = ax.plot([r["tau_ms"] for r in subset],
                        [r["upper_v"] if r["upper_v"] is not None else np.nan for r in subset],
                        "o-", label=f"Exponential, polarity {polarity:+d}")
        control = next(r for r in rows if r["polarity"] == polarity and r["waveform"] == "rectangular")
        if control["upper_v"] is not None:
            ax.axhline(control["upper_v"], color=line.get_color(), ls="--", label=f"Rectangular, polarity {polarity:+d}")
    ax.set(xscale="log", xlabel="Rise time constant (ms)", ylabel="Bracket upper drive magnitude (V)",
           title=f"Activation boundary; {a.duration:g} ms first phase, {a.normalization} normalization")
    ax.legend(fontsize=8); ax.grid(alpha=0.25)
    fig.savefig(output / "threshold_vs_tau.png", dpi=160); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(p)
    p.add_argument("--taus", type=float, nargs="+", default=[0.01, 0.03, 0.1, 0.3, 1.0])
    p.add_argument("--polarities", type=int, choices=[-1, 1], nargs="+", default=[1, -1])
    p.add_argument("--min-voltage", type=float, default=1.0)
    p.add_argument("--max-voltage", type=float, default=1000.0)
    p.add_argument("--scan-points", type=int, default=12)
    p.add_argument("--rtol", type=float, default=0.01)
    p.add_argument("--criterion", choices=["distal", "any"], default="distal")
    p.add_argument("--first-phase-only", action="store_true",
                   help="Require first crossing before shutoff; may exclude delayed onset responses")
    p.add_argument("--allow-end", action="store_true", help="Include end-initiated responses (diagnostic)")
    p.add_argument("--workers", type=int, default=0,
                   help="Processes; 0=auto (up to 12, half logical CPUs, and case count); 1=sequential")
    p.add_argument("--output", default="results/thresholds")
    a = p.parse_args()
    if (not np.all(np.isfinite([a.min_voltage, a.max_voltage, a.rtol, *a.taus])) or
            not 0 < a.min_voltage < a.max_voltage or a.scan_points < 2 or not 0 < a.rtol < 1 or min(a.taus) <= 0):
        p.error("Use positive finite taus, 0 < min-voltage < max-voltage, scan-points >= 2, and 0 < rtol < 1.")
    if a.workers < 0:
        p.error("--workers must be >= 0")
    # Remove duplicate cases to avoid redundant simulations; preserve input order.
    a.taus = list(dict.fromkeys(a.taus))
    a.polarities = list(dict.fromkeys(a.polarities))
    a.csv = str(Path(a.csv).resolve())
    if not Path(a.csv).is_file():
        p.error(f"CSV does not exist: {a.csv}")
    output = Path(a.output).resolve(); output.mkdir(parents=True, exist_ok=True)
    # Unique case tree prevents an interrupted/repeated run from mixing old trial logs.
    run_id = uuid.uuid4().hex[:12]
    tasks = []
    for polarity in a.polarities:
        for tau in [None, *a.taus]:
            kind = "rectangular" if tau is None else "exponential"
            directory = output / "cases" / run_id / f"{kind}_tau_{tau}_polarity_{polarity}"
            tasks.append((len(tasks), vars(a).copy(), polarity, tau, str(directory)))
    cpu_count = os.cpu_count() or 1
    workers = min(len(tasks), a.workers or min(12, max(1, cpu_count // 2)))
    # Set BEFORE spawning: child imports of NumPy/SciPy then use one native
    # math-library thread each instead of oversubscribing the machine.
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[key] = "1"
    completed = {}
    metadata = {"arguments": vars(a), "workers": workers, "logical_cpus": cpu_count,
                "start_method": "spawn", "status": "running", "run_id": run_id,
                "case_directories": [t[4] for t in tasks]}
    write_json(output / "configuration.json", metadata)
    write_combined(output, tasks, completed)
    # Prevent an old successful plot from appearing to describe a new failed run.
    (output / "threshold_vs_tau.png").unlink(missing_ok=True)
    print(f"Running {len(tasks)} independent cases with {workers} worker processes. "
          "Each case performs its amplitude scan and bisection sequentially.", flush=True)
    started = time.perf_counter()
    try:
        # spawn also works on Windows. Recycling after one case is intentional:
        # NEURON's global model state is discarded with the worker interpreter.
        with mp.get_context("spawn").Pool(workers, maxtasksperchild=1) as pool:
            for index, row, model_metadata in pool.imap_unordered(run_case, tasks, chunksize=1):
                completed[index] = row
                metadata.update(model_metadata)
                metadata["completed_cases"] = len(completed)
                write_json(output / "configuration.json", metadata)
                write_combined(output, tasks, completed)
                print(f"Completed {len(completed)}/{len(tasks)} cases "
                      f"({time.perf_counter() - started:.1f} s elapsed)", flush=True)
            pool.close()
            pool.join()
        rows = [completed[i] for i in sorted(completed)]
        plot_thresholds(rows, a, output)
        metadata["status"] = "completed"
    except BaseException as exc:
        # Pool context has terminated workers before their partial logs are read.
        metadata["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        metadata["error"] = str(exc)
        write_combined(output, tasks, completed, include_partial=True)
        raise
    finally:
        metadata["elapsed_seconds"] = time.perf_counter() - started
        metadata["completed_cases"] = len(completed)
        write_json(output / "configuration.json", metadata)
    print(f"Saved {output}. Missing thresholds mean none observed on the sampled grid, not proven absence.")


if __name__ == "__main__":
    mp.freeze_support()
    main()
