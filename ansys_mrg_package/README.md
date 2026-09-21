# Ansys extracellular voltages -> PyFibers MRG axon

This standalone package replaces the point-source stimulus in your example with
your Ansys CSV. It uses PyFibers' installed `ScaledStim`; do not overwrite the
package's `stimulation.py` with these files. Your original files are unchanged.

## Files

- `run_ansys_mrg.py`: run one pulse and save traces, plots, and spike events.
- `sweep_thresholds.py`: rectangular control plus exponential-rise threshold
  searches for each requested time constant and polarity.
- `ansys_mrg.py`: CSV loading, physical coordinate mapping, waveforms, model,
  recording, and analysis. Edit this file to add custom waveforms/protocols.
- `check_setup.py`: numerical checks, with optional actual NEURON controls.
- `nerve_voltage_example.csv`: your original example data.
- `requirements.txt`: dependencies; `README.md`: this guide.

## Install and run

Use your existing environment if it already runs the uploaded example. These
files were tested using Python 3.12, PyFibers 0.10.1, and NEURON 9.0.2 on Linux.
For a new environment:

```bash
python -m venv .venv
```

Activate with `.venv\Scripts\activate` in Windows Command Prompt,
`.venv\Scripts\Activate.ps1` in PowerShell, or `source .venv/bin/activate` on
Linux/macOS. Then:

```bash
python -m pip install -r requirements.txt
pyfibers_compile
python check_setup.py --integration
python run_ansys_mrg.py
```

NEURON must be installed for the SAME Python environment, with its mechanism
compiler available on PATH. On Windows, if `pip` cannot install a compatible
NEURON wheel, install NEURON using its official Windows instructions, verify
`python -c "import neuron"`, and install the remaining requirements into that
compatible environment. PyFibers documents Windows mechanism compilation and
finding `pyfibers_compile.exe` in the environment's Scripts folder. Compilation
must report success. If your original example already works, no recompilation
is normally needed. Do not ignore “mechanisms not found” messages.

Official setup: https://wmglab-duke.github.io/pyfibers/
NEURON: https://nrn.readthedocs.io/en/latest/

The default single run is +100 V, rectangular, 2 ms duration, onset 1 ms,
10 um fiber diameter, 25 nodes, 37 C, dt = 0.001 ms. It saves files to
`results/single`; plots are saved without requiring an interactive window.
The default CSV path is relative to these scripts, so it works from other
working directories too. Output paths are relative to your working directory.

## Example commands

Match the original example's 0.1 ms pulse duration (with a longer baseline):

```bash
python run_ansys_mrg.py --voltage 100 --duration 0.1 --output results/rectangle
```

First-order rising response; tau = 0.1 ms, truncated at 2 ms:

```bash
python run_ansys_mrg.py --waveform exponential --voltage 100 --tau 0.1 --duration 2 --output results/exp_tau_0p1
```

Reverse the CSV's electrode polarity and add a smooth fall:

```bash
python run_ansys_mrg.py --waveform exponential --voltage -100 --tau 0.1 --fall-tau 0.3 --output results/smooth
```

Biphasic pulse with equal-and-opposite voltage area:

```bash
python run_ansys_mrg.py --waveform biphasic --duration 0.1 --second-duration 0.1 --gap 0.05 --voltage -100 --output results/biphasic
```

A voltage-area-balanced pulse is charge-balanced only if the electrode/tissue
load makes current proportional to that voltage with time-independent
conductance. This CSV alone does not supply electrode current or charge.

Compare thresholds across rise constants and both polarities:

```bash
python sweep_thresholds.py --taus 0.01 0.03 0.1 0.3 1 --duration 2 --min-voltage 1 --max-voltage 1000 --scan-points 12 --output results/tau_asymptote
```

Repeat at equal first-phase peak instead of equal asymptotic drive:

```bash
python sweep_thresholds.py --taus 0.01 0.03 0.1 0.3 1 --duration 2 --normalization peak --output results/tau_peak
```

These sweeps may take many minutes. Progress prints after each simulation;
`trials.csv` is flushed after each trial. The default search range is an
illustrative computational range, not a prediction of physiological threshold
or a recommendation for experimental hardware. No result on that grid does
not establish that a spike is impossible at other amplitudes or placements.

Move or lengthen the axon without stretching the field:

```bash
python run_ansys_mrg.py --center-mm 110 --nodes 51 --output results/longer
```

Change diameter, temperature, polarity, duration, position, node count, and
interpolation with the corresponding flags. Use `--help` for all options.
Only longitudinal placement along the exported path can be changed with this
CSV. A different transverse electrode distance, axon orientation/path, electrode
geometry, or tissue conductivity generally requires another Ansys field export.

## What is actually applied?

For a spatial profile calculated at reference drive Vref = 100 V, the code uses

```
Ve(section, t) [mV] = 1000 * (Vcsv(section) [V] - C [V])
                     * (signed_drive [V] / Vref [V]) * w(t)
```

Positive drive preserves the CSV polarity; negative drive reverses it. “Positive”
is not automatically anodic at every location, nor “negative” always cathodic:
interpret the signs relative to your Ansys electrode and return definitions.
The scalar passed into `ScaledStim.run_sim` is drive/reference, not a current
in mA. Extracellular potentials are applied to every section (nodes, paranodes,
and internodes), rather than directly replacing membrane voltage.

C is the spatially constant mean of the CSV voltage, removed by default to
reduce numerical cancellation. It changes the voltage reference while retaining
all spatial differences. Use `--keep-offset` to verify invariance. A common
spatial offset is not the axon's transmembrane voltage. The integration checks
confirm common-offset invariance for this isolated MRG setup.

This separable spatial-times-temporal model assumes a linear, quasi-static
volume conductor with a fixed field shape and electrode boundary conditions
that scale together. If only one electrode is scaled while other nonzero
boundary voltages stay fixed, use separate source fields and superposition.
If electrode interfaces, dispersive tissue, or displacement currents cause
frequency-dependent field shape or delay, one static CSV is insufficient:
use a transient field export or an appropriate transfer-function model. An
imposed exponential waveform does not infer tissue/electrode dynamics from a
static Ansys solution.

## Physical placement and your specific CSV

The supplied CSV contains 1,001 samples from 0 to 200 mm at 0.2 mm spacing.
At a 100 V reference drive its potentials span approximately 50.01858 to
50.05476 V: a total spatial variation of 36.17908 mV on top of a ~50 V offset.

A 25-node, 10 um MRG-interpolation axon is approximately 26.94 mm long.
The default placement centers it at 100 mm, sampling section centers from
86.5324 to 113.4676 mm. Across this segment the reference field varies by
9.53846 mV. This centered placement is an explicit starting assumption, not
an inferred anatomical position.

The code converts section coordinates from um to mm, translates them to the
requested center, and interpolates at the physical section centers. It does
not compress the full 200 mm field into the 26.94 mm axon. Out-of-range
coordinates raise an error; extrapolation/clamping is disabled. CSV distance
must be arc length along the intended axon path, not a Cartesian coordinate
of an arbitrarily curved path. Headers explicitly identify the required units:
`Distance [mm]` and `Voltage [V]`.

Linear interpolation is the default, matching the usual PyFibers resampling
approach. PCHIP (`--interpolation pchip`) is available for sensitivity checks;
it cannot recover field detail that was not exported. Validate Ansys mesh
resolution and path-sampling convergence, especially near electrodes and where
field gradients change. Do not interpret extra interpolated samples as extra
physical information.

## Does rise time affect threshold?

Yes: amplitude alone does not determine excitation. The membrane integrates
input over time, while voltage-dependent channel states also evolve. Slow
depolarization can permit sodium-channel inactivation and potassium-channel
activation before a regenerative spike develops (accommodation). A faster
rise can therefore require a lower peak/asymptotic drive in some conditions.
The direction and size of the effect must be measured for this MRG geometry,
waveform family, polarity, and activation criterion; it is not a universal
law that every faster waveform is more efficient.

The exponential waveform during the first phase is

```
w(t) = 1 - exp(-(t - onset)/tau),  onset <= t < onset + duration
```

With `--normalization asymptote`, drive is the eventual plateau amplitude,
which a finite pulse may never reach. Its continuous first-phase peak is
`drive * (1 - exp(-duration/tau))`. Thus at fixed duration, changing tau changes
rise rate, achieved peak, and voltage-time area simultaneously. At the onset,
the slope is drive/tau. This stimulus tau is not the membrane time constant.

With `--normalization peak`, divide that waveform by
`1 - exp(-duration/tau)`. Every continuous waveform then has the same peak
immediately before shutoff, but voltage-time area still varies. The saved
sampled peak reports what the finite timestep actually applies. Both comparisons
are useful, but neither alone isolates channel accommodation. To investigate
mechanism, also examine nodal channel gates (PyFibers offers `record_gating`),
compare equal-area protocols, and vary pulse duration independently.

The default fall is abrupt. A slow rise followed by abrupt shutoff may produce
an offset response. `--fall-tau` adds exponential decay from the actual
first-phase endpoint, with automatic simulation time extending eight decay
constants plus 5 ms. The decay continues to the simulation end (it is not
cut to zero at eight constants). A late spike does not, by timing alone, prove
that the falling edge caused it; compare abrupt and smooth shutoff and inspect
all node traces.

## Spike and threshold definitions

- A spike candidate is an upward crossing of 0 mV in nodal membrane voltage.
  PyFibers APCount also uses 0 mV here. This avoids labelling small subthreshold
  depolarizations as spikes, but threshold crossings alone cannot prove active
  propagation in a strongly driven field.
- All nodes are recorded. Counts at the selected detector (default location
  0.9) are separate from whether any active node crosses the threshold.
- `first_spike_phase` reports whether the earliest active-node crossing occurs
  during the first phase or after its end. It is a timing classification.
- `end_initiation_flag` flags earliest crossings in the first or last two
  active nodes next to passive ends, including ties within one timestep.
  PyFibers also emits its own end-excitation warnings. Either flag is a reason
  to inspect the space-time plot and repeat with a longer axon. These are
  diagnostics, not an exhaustive detector of every boundary artifact.
- The default sweep accepts a detector spike anywhere in the observation
  window, provided earliest active-node initiation is not flagged near an end.
  It is a **distal detection threshold with an end filter**, not a rigorous
  proof of propagated activation or of a unique initiation mechanism. Both
  delayed onset responses and true offset responses can contribute.
- `--criterion any` instead searches for local activation anywhere in active
  nodes. `--first-phase-only` additionally requires the earliest active-node
  crossing before first-phase shutoff. This conservative timing filter can
  exclude legitimate delayed onset spikes, especially with short pulses; it
  does not distinguish onset/offset mechanisms conclusively. `--allow-end`
  includes flagged end initiation for diagnostic comparisons. These choices
  change the scientific question; keep them fixed across comparisons.
- The sweep runs a zero-drive control, scans the full logarithmic grid, and
  bisects the first observed false-to-true interval to 1% relative width by
  default. It records both bounds and flags coarse nonmonotonicity. Strong
  stimulation can cause conduction block or other response changes, so do not
  assume activation persists at every higher voltage. A finite grid can miss
  narrow activation windows; densify the grid where needed. A numerical solver
  failure aborts instead of being treated as “no spike.”

For propagation claims, inspect the space-time map for a traveling spike with
sequential latencies across multiple nodes, extending away from the initiation
site. Repeat near the reported threshold with half dt, longer observation
windows, longer fibers, and finer FEM/path sampling. In a long axon, default
5 ms post-stimulus observation may need to be increased with `--tstop`.

## Output files

Single runs create `simulation.png`, `run.json` (parameters, file hash, versions,
spike summary, actual sampled peak, signed/absolute voltage-time area),
`traces.npz` (all nodal Vm, times, mapped spatial field, sampled drive),
`node_traces.csv`, `mapped_profile.csv`, and `spike_events.csv`.
The sampled drive is held on intervals `[applied_time_ms[i], +dt)`; recorded
Vm is the solver response. Voltage-time area is not charge, and neither area
nor squared voltage alone establishes energy without current/impedance data.

Sweeps create `configuration.json`, `trials.csv`, `thresholds.csv`,
`threshold_vs_tau.png`, and a complete single-run output for each bracket's
upper endpoint. Missing graph points mean no accepted activation was observed
on the sampled amplitude grid; inspect trial logs for local, end, or late
responses that the chosen criterion excluded. A PyFibers warning that the
waveform peak is not exactly 1 is expected for truncated exponentials, ramps,
and asymmetric biphasic pulses; the actual sampled amplitude is saved and no
implicit renormalization is performed.

## Other factors to test

| Factor | Why it matters | Suggested comparison |
| --- | --- | --- |
| Pulse duration | Membrane charging and strength-duration behavior | Threshold versus duration at fixed field and polarity |
| Polarity and phase order | Depolarizing/hyperpolarizing regions differ | Both signs, monophasic versus biphasic |
| Fall time and interphase gap | Offset responses and recovery between phases | Abrupt versus smooth off; vary biphasic gap |
| Spatial variation along the axon | Neighboring sections are coupled by axial current | Axon placement and new electrode/path field exports |
| Diameter, node spacing, and node alignment | Changes spatial coupling and channel geometry | Diameter sweep; longitudinal shifts smaller than an internode |
| Repetition rate and prior pulses | Refractory period, afterpotentials, accommodation, and possible block | Extend Pulse to paired pulses/trains; allow adequate observation |
| Temperature and membrane properties | Change channel kinetics and excitability | Controlled temperature/model comparisons |
| Fiber length and numerical resolution | Can change apparent thresholds through numerical/boundary artifacts | Longer fibers; dt and FEM convergence |

In a simplified uniform cable, the spatial second derivative of extracellular
potential is a useful activating-function approximation. The heterogeneous
double-cable MRG model uses the full section-level potentials and conductances;
a second-derivative plot alone is not a quantitative activation threshold.
Large absolute extracellular voltage is not sufficient if it is nearly uniform
along an isolated fiber. Finite ends, bends, and changing geometry also matter.

## Validation performed for this package

`check_setup.py` verifies unit conversion, gauge invariance in interpolation,
no extrapolation, waveform values and biphasic discrete area, onset/end
classification, and first-boundary search with a synthetic nonmonotonic response.

With `--integration`, it also runs your actual CSV, a zero-drive control, and
a synthetic point-source positive control on a real compiled MRG model. It
checks common-offset invariance of Vm, repeated-trial reproducibility, and
positive-control latency/count with dt halved. The synthetic positive control
is only a software check, not an Ansys result. The supplied CSV at +100 V,
0.1 ms rectangular duration, centered 25-node axon, dt=0.001 ms, tstop=3 ms
produced no detected spikes. This does not establish a threshold or the effect
of tau in your Ansys field; run and validate the sweeps for those conclusions.

## Sources

- PyFibers installation and API: https://wmglab-duke.github.io/pyfibers/
- FEM/Ansys field mapping: https://wmglab-duke.github.io/pyfibers/extracellular_potentials.html
- McIntyre, Richardson & Grill (2002), original MRG model and recovery cycle:
  https://pubmed.ncbi.nlm.nih.gov/11826063/
- Vuckovic, Tosato & Struijk (2008), experiments comparing slowly rising pulses,
  depolarizing prepulses, and anodal block:
  https://pubmed.ncbi.nlm.nih.gov/18566504/
