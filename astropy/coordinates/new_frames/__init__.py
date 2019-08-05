"""
Alternative implementation of coordinate frames.

Objectives:

* Abolish frames with data
  Frames only contain the parameters for the transformations
* No metaprogramming or global variables to define new frames

"""
# TODO: Short-circuit transforms (for example GCRS<->ITRS or GCRS<->HCRS)
# TODO: Evaluate need for transformation graph
# TODO: Factories for custom frames or specific conventions

from .earth import CIRS, GCRS
from .enums import NutationModels
from .utils import transform
from .ecliptic import BarycentricEcliptic, HeliocentricEcliptic
from .equatorial import HCRS, ICRS

__all__ = [
    "CIRS",
    "GCRS",
    "BarycentricEcliptic",
    "HeliocentricEcliptic",
    "HCRS",
    "ICRS",
    "NutationModels",
]
