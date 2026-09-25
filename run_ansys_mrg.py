"""Run one Ansys -> MRG simulation. See README.md for setup and interpretation."""
import argparse
import json
from ansys_mrg import Experiment, Pulse, add_common_arguments, config_from_args


def main():
    p = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(p)
    p.add_argument("--voltage", type=float, default=100.0, help="Signed electrode drive V; positive retains CSV polarity")
    p.add_argument("--waveform", choices=["rectangular", "exponential", "ramp", "biphasic"], default="rectangular")
    p.add_argument("--tau", type=float, default=0.1, help="Exponential rise constant, ms")
    p.add_argument("--gap", type=float, default=0.05)
    p.add_argument("--second-duration", type=float, default=2.0)
    p.add_argument("--output", default="results/single")
    a = p.parse_args()
    pulse = Pulse(a.waveform, a.onset, a.duration, a.tau, a.normalization,
                  a.fall_tau, a.gap, a.second_duration)
    exp = Experiment(config_from_args(a))
    result = exp.run(a.voltage, pulse, a.dt, a.tstop)
    exp.save(result, a.output)
    print(json.dumps({k: v for k, v in result["summary"].items() if "by_node" not in k}, indent=2))
    print(f"Saved traces and plots to {a.output}")


if __name__ == "__main__":
    main()
