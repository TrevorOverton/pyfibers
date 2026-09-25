"""Numerical checks; add --integration to run actual NEURON/MRG controls."""
import argparse
from dataclasses import replace
import numpy as np
from ansys_mrg import ModelConfig, Pulse, Experiment, read_profile, sample_profile, summarize_spikes
from sweep_thresholds import scan_boundary


def numerical_checks():
    x, v = read_profile(ModelConfig().csv_path)
    assert len(x) == 1001 and np.isclose(x[-1] - x[0], 200)
    mapped, offset = sample_profile(np.array([0., 1., 2.]), np.array([50., 51., 52.]), np.array([.5, 1.5]))
    np.testing.assert_allclose(mapped, [-500, 500])
    assert offset == 51
    for method in ["linear", "pchip"]:
        a, _ = sample_profile(x, v, np.array([90., 100., 110.]), method)
        b, _ = sample_profile(x, v + 25, np.array([90., 100., 110.]), method)
        np.testing.assert_allclose(a, b, atol=1e-9)
    try:
        sample_profile(x, v, [-1.])
    except ValueError:
        pass
    else:
        raise AssertionError("Out-of-bounds mapping must fail")
    pulse = Pulse("exponential", duration_ms=2, tau_ms=.4)
    assert pulse(0) == 0 and pulse(pulse.off_ms) == 0
    assert np.isclose(pulse(1.4), 1 - np.exp(-1))
    peak = replace(pulse, normalization="peak")
    assert np.isclose(peak(peak.off_ms - 1e-9), 1)
    smooth = replace(pulse, fall_tau_ms=.2)
    assert np.isclose(smooth(smooth.off_ms), smooth.peak_fraction)
    assert np.isclose(smooth(smooth.off_ms + .2), smooth.peak_fraction / np.e)
    biph = Pulse("biphasic", duration_ms=.2, gap_ms=.1, second_duration_ms=.4)
    time = np.arange(0, 3, .001)
    assert abs(biph(time).sum() * .001) < 1e-10
    boundary = scan_boundary(lambda a: 3 <= a < 8, [1, 2, 4, 10])
    assert boundary["lower_v"] < 3 <= boundary["upper_v"]
    assert boundary["coarse_nonmonotonic"]
    absent = scan_boundary(lambda a: False, [1, 10])
    assert absent["upper_v"] is None
    # Synthetic crossing sequences verify onset and edge classification.
    t = np.arange(0, 6, .001)
    vm = np.full((9, len(t)), -80.)
    for i, onset in [(3, 1.2), (4, 1.3), (5, 1.4), (6, 1.5)]:
        vm[i, (t >= onset) & (t < onset + .1)] = 30
    s = summarize_spikes(t, vm, 1, 6, Pulse())
    assert s["detector_spike_count"] == 1 and not s["end_initiation_flag"]
    print("PASS: units, gauge, interpolation bounds, waveforms, spike classification, and nonmonotonic scan.")


def integration_checks():
    exp = Experiment(ModelConfig())
    p = Pulse("rectangular", duration_ms=.1)
    result = exp.run(100, p, .001, 3)
    assert not result["summary"]["any_active_node_spike"]
    baseline = exp.run(0, p, .001, 3)
    assert not baseline["summary"]["any_active_node_spike"]
    # Positive control uses a synthetic point source, NOT the Ansys CSV.
    potentials = exp.fiber.point_source_potentials(0, 250, exp.fiber.length / 2, 1, .3)
    exp.fiber.potentials = potentials.copy()
    positive = exp.run(-10, p, .001, 3)  # -0.1 mA point-source current
    assert positive["summary"]["detector_spike_count"] > 0
    assert not positive["summary"]["end_initiation_flag"]
    exp.fiber.potentials = potentials.copy() + 50000  # common offset 50 V
    shifted = exp.run(-10, p, .001, 3)
    np.testing.assert_allclose(positive["vm_mv"], shifted["vm_mv"], atol=1e-3, rtol=0)
    exp.fiber.potentials = potentials.copy()
    repeated = exp.run(-10, p, .001, 3)
    np.testing.assert_allclose(positive["vm_mv"], repeated["vm_mv"], atol=1e-6, rtol=0)
    fine = exp.run(-10, p, .0005, 3)
    assert fine["summary"]["detector_spike_count"] == positive["summary"]["detector_spike_count"]
    assert abs(fine["summary"]["first_detector_spike_ms"] - positive["summary"]["first_detector_spike_ms"]) < .05
    print("PASS: actual Ansys trial, zero drive, positive point-source control, common-offset invariance,")
    print("      repeated-trial reset, and positive-control timing under dt halving.")
    print("CSV 100 V, 0.1 ms:", result["summary"]["detector_spike_count"], "distal spikes")
    print("Synthetic control:", positive["summary"]["detector_spike_count"], "distal spikes at",
          positive["summary"]["first_detector_spike_ms"], "ms")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integration", action="store_true")
    args = parser.parse_args()
    numerical_checks()
    if args.integration:
        integration_checks()
