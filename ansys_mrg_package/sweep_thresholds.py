"""Scan drive voltage and refine the first observed activation boundary per tau.

Rectangular control is always included. All coarse trials are retained so
nonmonotonic responses (including loss of detection at high drive) are visible.
A finite scan cannot exclude activation windows between sampled amplitudes.
"""
import argparse
import csv
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
    p.add_argument("--first-phase-only", action="store_true", help="Require first active-node crossing before first-phase shutoff; can exclude delayed onset responses")
    p.add_argument("--allow-end", action="store_true", help="Accept flagged end-initiated responses (diagnostic only)")
    p.add_argument("--output", default="results/thresholds")
    a = p.parse_args()
    if (not np.all(np.isfinite([a.min_voltage, a.max_voltage, a.rtol, *a.taus])) or
            not 0 < a.min_voltage < a.max_voltage or a.scan_points < 2 or not 0 < a.rtol < 1 or min(a.taus) <= 0):
        p.error("Use positive finite taus, 0 < min-voltage < max-voltage, scan-points >= 2, and 0 < rtol < 1.")
    exp = Experiment(config_from_args(a))
    output = Path(a.output); output.mkdir(parents=True, exist_ok=True)
    write_json(output / "configuration.json", {**exp.metadata, "arguments": vars(a)})
    grid = np.geomspace(a.min_voltage, a.max_voltage, a.scan_points)
    rows = []
    with (output / "trials.csv").open("w", newline="") as f:
        fields = ["waveform", "tau_ms", "polarity", "magnitude_v", "accepted", "first_active_spike_ms",
                  "first_detector_spike_ms", "detector_spike_count", "first_spike_phase", "end_initiation_flag", "max_vm_mv"]
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
        for polarity in a.polarities:
            for tau in [None, *a.taus]:
                kind = "rectangular" if tau is None else "exponential"
                pulse = Pulse(kind, a.onset, a.duration, 0.1 if tau is None else tau,
                              a.normalization, a.fall_tau)
                pulse.validate()
                cache = {}

                def evaluate(magnitude):
                    if magnitude in cache:
                        return cache[magnitude]
                    result = exp.run(polarity * magnitude, pulse, a.dt, a.tstop)
                    s = result["summary"]
                    ok = accepted(s, a.criterion, not a.first_phase_only, a.allow_end)
                    cache[magnitude] = ok
                    row = {"waveform": kind, "tau_ms": tau, "polarity": polarity,
                           "magnitude_v": magnitude, "accepted": ok}
                    row.update({k: s[k] for k in fields if k in s})
                    writer.writerow(row); f.flush()
                    print(f"{kind:11s} tau={str(tau):>5} polarity={polarity:+d} "
                          f"drive={magnitude:.6g} V: accepted={ok}, "
                          f"distal={s['detector_spike_count']}, end={s['end_initiation_flag']}", flush=True)
                    return ok

                boundary = scan_boundary(evaluate, grid, a.rtol)
                row = {"waveform": kind, "tau_ms": tau, "polarity": polarity,
                       "normalization": a.normalization, "criterion": a.criterion, **boundary}
                high = boundary["upper_v"]
                row["signed_upper_drive_v"] = None if high is None else polarity * high
                row["upper_first_phase_peak_abs_v"] = None if high is None else high * pulse.peak_fraction
                rows.append(row)
                # Save progress and representative near-threshold traces, independently per case.
                with (output / "thresholds.csv").open("w", newline="") as tf:
                    tw = csv.DictWriter(tf, fieldnames=list(row)); tw.writeheader(); tw.writerows(rows)
                if high is not None:
                    result = exp.run(polarity * high, pulse, a.dt, a.tstop)
                    exp.save(result, output / f"{kind}_tau_{tau}_polarity_{polarity}")
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
    print(f"Saved {output}. Missing thresholds mean none observed on the sampled grid, not proven absence.")


if __name__ == "__main__":
    main()
