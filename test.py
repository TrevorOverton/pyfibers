import pyfibers

# Enable INFO level logging to see simulation progress
pyfibers.enable_logging()

from pyfibers import FiberModel
from pyfibers import build_fiber

n_nodes = 25
fiber = build_fiber(
    fiber_model=FiberModel.MRG_DISCRETE,
    diameter=10,  # um
    n_nodes=25,  # um
    temperature=37,  # C
)
print(fiber)
