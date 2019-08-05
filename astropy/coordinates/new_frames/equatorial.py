from astropy.coordinates.builtin_frames.utils import DEFAULT_OBSTIME

from .base import HasObstime, BaseCoordinateFrame
from .helpers import apply_affine


class ICRS(BaseCoordinateFrame):
    def to_icrs(self, coords_self):
        return coords_self

    def from_icrs(self, coords_icrs):
        return coords_icrs


class HCRS(BaseCoordinateFrame, HasObstime):
    def __init__(self, obstime=DEFAULT_OBSTIME):
        self._obstime = obstime

    def to_icrs(self, coords_self):
        if coords_self.differentials:
            raise NotImplementedError
        else:
            # TODO: This should be cached! Most of the times it will be DEFAULT_OBSTIME
            from astropy.coordinates.solar_system import get_body_barycentric

            bary_sun_pos = get_body_barycentric("sun", self.obstime)

        return apply_affine(coords_self, None, bary_sun_pos)

    def from_icrs(self, coords_icrs):
        if coords_icrs.differentials:
            raise NotImplementedError
        else:
            # TODO: This should be cached! Most of the times it will be DEFAULT_OBSTIME
            from astropy.coordinates.solar_system import get_body_barycentric

            bary_sun_pos = -get_body_barycentric("sun", self.obstime)

        return apply_affine(coords_icrs, None, bary_sun_pos)
