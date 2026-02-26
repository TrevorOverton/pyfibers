import pyfibers

# Enable INFO level logging to see simulation progress
pyfibers.enable_logging()


from pyfibers import build_fiber, FiberModel

# create fiber model
n_nodes = 25
fiber = build_fiber(FiberModel.MRG_INTERPOLATION, diameter=10, n_nodes=n_nodes)


from scipy.interpolate import interp1d

time_step = 0.001  # milliseconds
time_stop = 20  # milliseconds
start, on, off = 0, 0.1, 0.2  # milliseconds
waveform = interp1d(
    [start, on, off, time_stop], [0, 1, 0, 0], kind="previous"
)  # monophasic rectangular pulse

import matplotlib.pyplot as plt
import numpy as np

time_steps = np.arange(0, time_stop, time_step)
plt.plot(time_steps, waveform(time_steps))
plt.xlim(0, 1)
plt.title('Stimulation waveform')
plt.xlabel('Time (ms)')
plt.ylabel('Amplitude')
plt.show()

fiber.potentials = fiber.point_source_potentials(0, 250, fiber.length / 2, 1, 10)

plt.plot(fiber.longitudinal_coordinates, fiber.potentials)
plt.xlabel('Distance along fiber (µm)')
plt.ylabel('Electrical potential (mV)')
plt.title('Extracellular potentials')
plt.show()