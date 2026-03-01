import pyfibers

# Enable INFO level logging to see simulation progress
pyfibers.enable_logging()

# from pyfibers import build_fiber, FiberModel

# # create fiber model
# n_nodes = 25
# fiber = build_fiber(FiberModel.MRG_INTERPOLATION, diameter=10, n_nodes=n_nodes)


# from scipy.interpolate import interp1d

# time_step = 0.001  # milliseconds
# time_stop = 20  # milliseconds
# start, on, off = 0, 0.1, 0.2  # milliseconds
# waveform = interp1d(
#     [start, on, off, time_stop], [0, 1, 0, 0], kind="previous"
# )  # monophasic rectangular pulse

# import matplotlib.pyplot as plt
# import numpy as np

# time_steps = np.arange(0, time_stop, time_step)
# plt.plot(time_steps, waveform(time_steps))
# plt.xlim(0, 1)
# plt.title('Stimulation waveform')
# plt.xlabel('Time (ms)')
# plt.ylabel('Amplitude')
# plt.show()

# fiber.potentials = fiber.point_source_potentials(0, 250, fiber.length / 2, 1, 10)

# plt.plot(fiber.longitudinal_coordinates, fiber.potentials)
# plt.xlabel('Distance along fiber (µm)')
# plt.ylabel('Electrical potential (mV)')
# plt.title('Extracellular potentials')
# plt.show()

# from pyfibers import ScaledStim

# # Create instance of ScaledStim class
# stimulation = ScaledStim(waveform=waveform, dt=time_step, tstop=time_stop)
# print(stimulation)

# stimamp = -1.5  # technically unitless, but scales the unit (1 mA) stimulus to 1.5 mA
# ap, time = stimulation.run_sim(stimamp, fiber)
# print(f'Number of action potentials detected: {ap}')
# print(f'Time of last action potential detection: {time} ms')


import numpy as np
from pyfibers import build_fiber, FiberModel, ScaledStim
from scipy.interpolate import interp1d

# create fiber model
n_nodes = 25
fiber = build_fiber(FiberModel.MRG_INTERPOLATION, diameter=10, n_nodes=n_nodes)

# Setup for simulation
time_step = 0.001
time_stop = 10
start, on, off = 0, 0.1, 0.2
waveform = interp1d(
    [start, on, off, time_stop], [0, 1, 0, 0], kind="previous"
)  # biphasic rectangular pulse

fiber.potentials = fiber.point_source_potentials(0, 250, fiber.length / 2, 1, 10)

# Create stimulation object
stimulation = ScaledStim(waveform=waveform, dt=time_step, tstop=time_stop)

fiber.record_vm()
ap, time = stimulation.run_sim(-1.5, fiber)
print(f'Number of action potentials detected: {ap}')
print(f'Time of last action potential detection: {time} ms')

import matplotlib.pyplot as plt
import seaborn as sns

sns.set(font_scale=1.5, style='whitegrid', palette='colorblind')

plt.figure()
plt.plot(
    np.array(stimulation.time),
    list(fiber.vm[fiber.loc_index(0.9)]),
    label='end node',
    color='royalblue',
    linewidth=2,
)
plt.xlim(0, 2)
plt.legend()
plt.xlabel('Time (ms)')
plt.ylabel('$V_m$ $(mV)$')
plt.show()