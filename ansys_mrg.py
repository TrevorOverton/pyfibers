"""Ansys spatial potentials -> PyFibers MRG extracellular stimulation.

Units: CSV distance mm, CSV potential V, model distance um, model potential mV,
all times ms. This module does not modify PyFibers or its stimulation.py.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import warnings
from dataclasses import asdict, dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator


def package_version(package):
    """Best-effort diagnostic: standalone installs may lack pip metadata."""
    try:
        return version(package)
    except PackageNotFoundError:
        # NEURON's standalone installer can provide an importable module
        # without a corresponding Python distribution metadata directory.
        module_version = getattr(sys.modules.get(package), "__version__", None)
        return str(module_version) if module_version is not None else "unknown (no package metadata)"


@dataclass(frozen=True)
class ModelConfig:
    csv_path: str = str(Path(__file__).with_name("nerve_voltage_example.csv"))
    reference_voltage_v: float = 100.0
    diameter_um: float = 10.0
    n_nodes: int = 25
    center_mm: float | None = None  # None = middle of CSV path
    temperature_c: float = 37.0
    passive_end_nodes: int = 2
    interpolation: str = "linear"
    remove_common_offset: bool = True
    detect_location: float = 0.9


@dataclass(frozen=True)
class Pulse:
    kind: str = "rectangular"  # rectangular, exponential, ramp, biphasic
    onset_ms: float = 1.0
    duration_ms: float = 2.0  # first phase; held fixed across tau sweep
    tau_ms: float = 0.1
    normalization: str = "asymptote"  # exponential: asymptote or peak
    fall_tau_ms: float = 0.0  # zero = abrupt off; >0 = exponential decay
    gap_ms: float = 0.05  # biphasic interphase gap
    second_duration_ms: float = 2.0

    def validate(self):
        vals = [self.onset_ms, self.duration_ms, self.tau_ms, self.fall_tau_ms,
                self.gap_ms, self.second_duration_ms]
        if not np.all(np.isfinite(vals)):
            raise ValueError("Pulse parameters must be finite.")
        if self.kind not in {"rectangular", "exponential", "ramp", "biphasic"}:
            raise ValueError("Unknown pulse kind.")
        if self.onset_ms <= 0 or self.duration_ms <= 0:
            raise ValueError("Use positive onset and duration (allow a baseline period).")
        if self.tau_ms <= 0 or min(self.fall_tau_ms, self.gap_ms) < 0 or self.second_duration_ms <= 0:
            raise ValueError("Invalid time constant, gap, or second-phase duration.")
        if self.normalization not in {"asymptote", "peak"}:
            raise ValueError("normalization must be asymptote or peak")
        if self.kind == "biphasic" and self.fall_tau_ms:
            raise ValueError("fall_tau_ms is only supported for monophasic waveforms.")

    @property
    def off_ms(self):
        return self.onset_ms + self.duration_ms

    @property
    def end_ms(self):
        if self.kind == "biphasic":
            return self.off_ms + self.gap_ms + self.second_duration_ms
        return self.off_ms + 8 * self.fall_tau_ms

    @property
    def peak_fraction(self):
        if self.kind == "exponential" and self.normalization == "asymptote":
            return float(-np.expm1(-self.duration_ms / self.tau_ms))
        return 1.0

    def __call__(self, time_ms):
        """Callable accepted directly by ScaledStim; scalar and array support."""
        # Round at sub-picosecond precision so floating arithmetic does not
        # add or remove one sample at mathematically grid-aligned phase edges.
        t = np.round(np.asarray(time_ms, dtype=float), 12)
        onset, off = round(self.onset_ms, 12), round(self.off_ms, 12)
        u = t - onset
        inside = (t >= onset) & (t < off)
        if self.kind == "exponential":
            value = -np.expm1(-np.maximum(u, 0) / self.tau_ms)
            if self.normalization == "peak":
                value /= -np.expm1(-self.duration_ms / self.tau_ms)
        elif self.kind == "ramp":
            value = np.maximum(u, 0) / self.duration_ms
        else:
            value = np.ones_like(t)
        out = np.where(inside, value, 0.0)
        if self.kind == "biphasic":
            start2 = round(self.off_ms + self.gap_ms, 12)
            out = np.where((t >= start2) & (t < round(self.end_ms, 12)),
                           -self.duration_ms / self.second_duration_ms, out)
        elif self.fall_tau_ms > 0:
            tail = self.peak_fraction * np.exp(-np.maximum(t - self.off_ms, 0) / self.fall_tau_ms)
            out = np.where(t >= off, tail, out)
        return float(out) if out.ndim == 0 else out


def read_profile(path):
    """Read the supplied explicit-unit CSV. Sort rows; reject ambiguous duplicates."""
    with open(path, encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or not {"Distance [mm]", "Voltage [V]"}.issubset(reader.fieldnames):
            raise ValueError('Expected CSV headers "Distance [mm]","Voltage [V]". Convert units explicitly.')
        rows = [(float(r["Distance [mm]"]), float(r["Voltage [V]"])) for r in reader]
    a = np.asarray(rows, dtype=float)
    if a.ndim != 2 or len(a) < 3 or not np.isfinite(a).all():
        raise ValueError("CSV must contain at least three finite distance/voltage pairs.")
    a = a[np.argsort(a[:, 0])]
    if np.any(np.diff(a[:, 0]) <= 0):
        raise ValueError("CSV distances must be unique; duplicate path positions are ambiguous.")
    return a[:, 0], a[:, 1]


def sample_profile(x_mm, v_v, targets_mm, method="linear", remove_offset=True):
    if np.min(targets_mm) < x_mm[0] or np.max(targets_mm) > x_mm[-1]:
        raise ValueError(f"Axon spans {min(targets_mm):.4f}..{max(targets_mm):.4f} mm; "
                         f"CSV covers {x_mm[0]:.4f}..{x_mm[-1]:.4f} mm. "
                         "Change center/nodes or export a longer path. Extrapolation is disabled.")
    # Change the spatially constant gauge BEFORE converting/interpolating.
    offset_v = float(np.mean(v_v)) if remove_offset else 0.0
    centered_mv = (v_v - offset_v) * 1000.0
    if method == "linear":
        values = np.interp(targets_mm, x_mm, centered_mv)
    elif method == "pchip":
        values = PchipInterpolator(x_mm, centered_mv, extrapolate=False)(targets_mm)
    else:
        raise ValueError("interpolation must be linear or pchip")
    return values, offset_v


def crossing_times(t, voltage, level=0.0):
    inds = np.flatnonzero((voltage[:-1] < level) & (voltage[1:] >= level))
    return t[inds] + (level - voltage[inds]) * (t[inds + 1] - t[inds]) / (voltage[inds + 1] - voltage[inds])


def summarize_spikes(t, vm, passive, detect_index, pulse):
    """0 mV upward crossings are candidates, not a proof of regenerative propagation.

    A conservative end flag checks the earliest active-node crossing, including
    ties within one integration step. Full traces are saved for inspection.
    """
    crossings = [crossing_times(t, row) for row in vm]
    if any(np.any(c < pulse.onset_ms) for c in crossings):
        raise RuntimeError("Spontaneous/baseline spike: the trial cannot define a stimulus threshold.")
    first = [float(c[0]) if len(c) else None for c in crossings]
    active = range(passive, len(vm) - passive)
    fired = [i for i in active if first[i] is not None]
    earliest = min((first[i] for i in fired), default=None)
    ties = [i for i in fired if first[i] <= earliest + 1.01 * (t[1] - t[0])] if fired else []
    edge = any(i <= passive + 1 or i >= len(vm) - passive - 2 for i in ties)
    phase = ("none" if earliest is None else
             "during_first_phase" if earliest < pulse.off_ms else "after_first_phase")
    return {
        "any_active_node_spike": bool(fired),
        "detector_spike_count": int(len(crossings[detect_index])),
        "first_detector_spike_ms": first[detect_index],
        "first_active_spike_ms": earliest,
        "earliest_active_nodes": ties,
        "end_initiation_flag": edge,
        "first_spike_phase": phase,
        "first_crossing_ms_by_node": first,
        "crossings_ms_by_node": [c.tolist() for c in crossings],
        "max_vm_mv": float(np.max(vm)),
    }


class Experiment:
    def __init__(self, config: ModelConfig):
        from pyfibers import FiberModel, build_fiber
        self.config = config
        if (not np.isfinite(config.reference_voltage_v) or config.reference_voltage_v <= 0
                or not np.isfinite(config.diameter_um) or config.diameter_um <= 0
                or not np.isfinite(config.temperature_c)):
            raise ValueError("Reference voltage and diameter must be positive; temperature must be finite.")
        if config.n_nodes < 9 or config.n_nodes % 2 != 1:
            raise ValueError("Use an odd n_nodes >= 9.")
        if config.passive_end_nodes < 1 or 2 * config.passive_end_nodes + 5 > config.n_nodes:
            raise ValueError("Keep at least one passive node at each end and five active nodes.")
        self.x_mm, self.v_v = read_profile(config.csv_path)
        self.fiber = build_fiber(FiberModel.MRG_INTERPOLATION, diameter=config.diameter_um,
                                n_nodes=config.n_nodes, temperature=config.temperature_c,
                                passive_end_nodes=config.passive_end_nodes)
        center = (self.x_mm[0] + self.x_mm[-1]) / 2 if config.center_mm is None else config.center_mm
        if not np.isfinite(center):
            raise ValueError("center_mm must be finite")
        coords = np.asarray(self.fiber.longitudinal_coordinates)
        self.section_mm = center + (coords - (coords[0] + coords[-1]) / 2) / 1000
        self.profile_mv, self.offset_v = sample_profile(
            self.x_mm, self.v_v, self.section_mm, config.interpolation, config.remove_common_offset)
        self.fiber.potentials = self.profile_mv.copy()
        self.node_indices = [i for i, sec in enumerate(self.fiber.sections) if sec in self.fiber.nodes]
        self.node_mm = self.section_mm[self.node_indices]
        self.detect_index = self.fiber.loc_index(config.detect_location)
        if not config.passive_end_nodes <= self.detect_index < config.n_nodes - config.passive_end_nodes:
            raise ValueError("Detection location maps to a passive end node; move it inward.")
        self.fiber.record_vm()
        self.metadata = {
            "model": asdict(config), "resolved_center_mm": float(center),
            "csv_sha256": hashlib.sha256(Path(config.csv_path).read_bytes()).hexdigest(),
            "pyfibers_version": package_version("pyfibers"), "neuron_version": package_version("neuron"),
            "fiber_length_mm": float(self.fiber.length / 1000),
            "section_count": len(self.fiber.sections), "detector_node": self.detect_index,
            "csv_potential_ptp_mv": float(np.ptp(self.v_v) * 1000),
            "mapped_potential_ptp_mv": float(np.ptp(self.profile_mv)),
            "removed_common_offset_v": self.offset_v,
        }
        print(f"CSV {self.x_mm[0]:g}..{self.x_mm[-1]:g} mm; axon "
              f"{self.section_mm[0]:.4f}..{self.section_mm[-1]:.4f} mm; "
              f"mapped potential range {np.ptp(self.profile_mv):.6g} mV "
              f"at {config.reference_voltage_v:g} V reference.")

    def run(self, drive_v, pulse: Pulse, dt_ms=0.001, tstop_ms=None):
        from pyfibers import ScaledStim
        from neuron import h
        pulse.validate()
        if not np.isfinite(drive_v) or not np.isfinite(dt_ms) or dt_ms <= 0:
            raise ValueError("Drive must be finite and dt must be positive and finite.")
        features = [pulse.duration_ms]
        if pulse.kind == "exponential":
            features.append(pulse.tau_ms)
        if pulse.fall_tau_ms:
            features.append(pulse.fall_tau_ms)
        if pulse.kind == "biphasic":
            features.append(pulse.second_duration_ms)
            if pulse.gap_ms:
                features.append(pulse.gap_ms)
        if dt_ms > min(features) / 10 + 1e-12:
            raise ValueError("dt must be <= 1/10 of the shortest waveform timescale; then check dt convergence.")
        requested_stop = pulse.end_ms + 5.0 if tstop_ms is None else tstop_ms
        if not np.isfinite(requested_stop) or requested_stop <= pulse.end_ms:
            raise ValueError("tstop must exceed waveform end; allow time for propagation and offset responses.")
        n = int(np.ceil(requested_stop / dt_ms))
        stop = n * dt_ms
        h.CVode().active(0)
        h.secondorder = 0
        stim = ScaledStim(waveform=pulse, dt=dt_ms, tstop=stop)
        # PyFibers warns if callable peak is not exactly 1; expected for an
        # unnormalized truncated exponential. Do not renormalize it implicitly.
        count, last = stim.run_sim(float(drive_v / self.config.reference_voltage_v), self.fiber,
                                  ap_detect_location=self.config.detect_location,
                                  ap_detect_threshold=0, fail_on_end_excitation=False)
        t = np.asarray(stim.time).copy()
        vm = np.asarray([np.asarray(v).copy() for v in self.fiber.vm])
        if not np.isfinite(vm).all() or vm.shape[1] != len(t):
            raise RuntimeError("Invalid recorded membrane voltage or time dimensions.")
        summary = summarize_spikes(t, vm, self.config.passive_end_nodes, self.detect_index, pulse)
        sample_times = np.arange(len(t) - 1) * dt_ms
        applied = np.asarray(pulse(sample_times)) * drive_v
        summary.update({"drive_v": float(drive_v), "scale_factor": float(drive_v / self.config.reference_voltage_v),
                        "dt_ms": dt_ms, "tstop_ms": float(t[-1]), "pulse": asdict(pulse),
                        "nominal_first_phase_peak_v": float(drive_v * pulse.peak_fraction),
                        "sampled_peak_abs_v": float(np.max(np.abs(applied))),
                        "signed_voltage_area_v_ms": float(applied.sum() * dt_ms),
                        "absolute_voltage_area_v_ms": float(np.abs(applied).sum() * dt_ms),
                        "pyfibers_detector_count": int(count),
                        "pyfibers_last_detector_ms": None if not count else float(last)})
        return {"summary": summary, "time_ms": t, "vm_mv": vm,
                "applied_time_ms": sample_times, "applied_drive_v": applied}

    def save(self, result, output):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        summary = result["summary"]
        write_json(output / "run.json", {**self.metadata, **summary})
        np.savez_compressed(output / "traces.npz", time_ms=result["time_ms"],
                            vm_mv=result["vm_mv"], node_position_mm=self.node_mm,
                            section_position_mm=self.section_mm, reference_potential_mv=self.profile_mv,
                            applied_time_ms=result["applied_time_ms"], applied_drive_v=result["applied_drive_v"])
        np.savetxt(output / "node_traces.csv", np.column_stack([result["time_ms"], result["vm_mv"].T]),
                   delimiter=",", header="time_ms," + ",".join(f"node_{i}_vm_mv" for i in range(len(self.node_mm))), comments="")
        np.savetxt(output / "mapped_profile.csv", np.column_stack([self.section_mm, self.profile_mv]),
                   delimiter=",", header="position_mm,potential_mv_at_reference_drive", comments="")
        with (output / "spike_events.csv").open("w", newline="") as f:
            w = csv.writer(f); w.writerow(["node", "position_mm", "crossing_time_ms", "active_node"])
            for i, times in enumerate(summary["crossings_ms_by_node"]):
                for t in times:
                    w.writerow([i, self.node_mm[i], t,
                                self.config.passive_end_nodes <= i < self.config.n_nodes - self.config.passive_end_nodes])
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        ax = axes[0, 0]
        ax.plot(self.x_mm, 1000 * (self.v_v - self.offset_v), color="0.65", label="CSV")
        ax.plot(self.section_mm, self.profile_mv, color="tab:blue", label="Mapped sections")
        ax.set(xlabel="Path distance (mm)", ylabel="Extracellular potential (mV)",
               title=f"Spatial profile at {self.config.reference_voltage_v:g} V drive")
        ax.legend()
        axes[0, 1].step(result["applied_time_ms"], result["applied_drive_v"], where="post")
        axes[0, 1].set(xlabel="Time (ms)", ylabel="Drive voltage (V)", title="Applied waveform")
        for i in sorted(set([self.fiber.loc_index(0.1), self.fiber.loc_index(0.5), self.detect_index])):
            axes[1, 0].plot(result["time_ms"], result["vm_mv"][i], label=f"Node {i}: {self.node_mm[i]:.2f} mm")
        axes[1, 0].set(xlabel="Time (ms)", ylabel="Membrane potential (mV)", title="Nodal responses")
        axes[1, 0].legend(fontsize=8)
        im = axes[1, 1].pcolormesh(result["time_ms"], self.node_mm, result["vm_mv"],
                                  shading="auto", cmap="viridis", vmin=-90, vmax=40, rasterized=True)
        axes[1, 1].set(xlabel="Time (ms)", ylabel="Node position (mm)", title="Inspect initiation and propagation")
        fig.colorbar(im, ax=axes[1, 1], label="Membrane potential (mV)")
        fig.savefig(output / "simulation.png", dpi=160)
        plt.close(fig)


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False), encoding="utf-8")


def add_common_arguments(parser):
    parser.add_argument("--csv", default=ModelConfig().csv_path)
    parser.add_argument("--reference-voltage", type=float, default=100.0)
    parser.add_argument("--diameter", type=float, default=10.0, help="Fiber diameter in um")
    parser.add_argument("--nodes", type=int, default=25)
    parser.add_argument("--center-mm", type=float, default=None, help="Center on CSV path (default midpoint)")
    parser.add_argument("--temperature", type=float, default=37.0)
    parser.add_argument("--passive-end-nodes", type=int, default=2)
    parser.add_argument("--interpolation", choices=["linear", "pchip"], default="linear")
    parser.add_argument("--keep-offset", action="store_true")
    parser.add_argument("--detect-location", type=float, default=0.9)
    parser.add_argument("--onset", type=float, default=1.0)
    parser.add_argument("--duration", type=float, default=2.0, help="First-phase duration, ms")
    parser.add_argument("--normalization", choices=["asymptote", "peak"], default="asymptote")
    parser.add_argument("--fall-tau", type=float, default=0.0)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--tstop", type=float, default=None)


def config_from_args(a):
    return ModelConfig(a.csv, a.reference_voltage, a.diameter, a.nodes, a.center_mm,
                       a.temperature, a.passive_end_nodes, a.interpolation, not a.keep_offset, a.detect_location)
