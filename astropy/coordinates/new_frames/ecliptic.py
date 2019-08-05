from astropy import _erfa as erfa
from astropy import units as u
from astropy.coordinates.matrix_utilities import matrix_product, rotation_matrix
from astropy.coordinates.builtin_frames.utils import (
    EQUINOX_J2000,
    DEFAULT_OBSTIME,
    get_jd12,
)

from .base import HasEquinox, HasObstime, HasNutation, BaseCoordinateFrame
from .enums import NutationModels
from .helpers import apply_affine
from .equatorial import HCRS


def _mean_ecliptic_rotation_matrix(equinox):
    jd1, jd2 = get_jd12(equinox, "tt")
    rmat = erfa.ecm06(jd1, jd2)
    return rmat


def _true_ecliptic_rotation_matrix(equinox):
    # This code calls pnm06a from ERFA, which retrieves the precession
    # matrix (including frame bias) according to the IAU 2006 model, and
    # including the nutation. This family of systems is less popular
    # than the "mean" ecliptic ones, obtained using ecm06
    # (see https://github.com/astropy/astropy/pull/6508).
    jd1, jd2 = get_jd12(equinox, "tt")
    rnpb = erfa.pnm06a(jd1, jd2)
    obl = erfa.obl06(jd1, jd2) * u.radian
    return matrix_product(rotation_matrix(obl, "x"), rnpb)


class BarycentricEcliptic(BaseCoordinateFrame, HasEquinox, HasNutation):
    def __init__(self, equinox=EQUINOX_J2000, nutation=NutationModels.MEAN):
        self._equinox = equinox
        self._nutation = NutationModels(nutation)

    def to_icrs(self, coords_self):
        if self.nutation is NutationModels.MEAN:
            rotation = _mean_ecliptic_rotation_matrix(self.equinox).T
        elif self.nutation is NutationModels.TRUE:
            rotation = _true_ecliptic_rotation_matrix(self.equinox).T

        return apply_affine(coords_self, rotation)

    def from_icrs(self, coords_icrs):
        if self.nutation is NutationModels.MEAN:
            rotation = _mean_ecliptic_rotation_matrix(self.equinox)
        elif self.nutation is NutationModels.TRUE:
            rotation = _true_ecliptic_rotation_matrix(self.equinox)

        return apply_affine(coords_icrs, rotation)


class BarycentricMeanEcliptic(BarycentricEcliptic):
    def __init__(self, equinox=EQUINOX_J2000):
        super().__init__(equinox=equinox, nutation=NutationModels.MEAN)


class BarycentricTrueEcliptic(BarycentricEcliptic):
    def __init__(self, equinox=EQUINOX_J2000):
        super().__init__(equinox=equinox, nutation=NutationModels.TRUE)


class HeliocentricEcliptic(BaseCoordinateFrame, HasObstime, HasEquinox, HasNutation):
    def __init__(
        self,
        obstime=DEFAULT_OBSTIME,
        equinox=EQUINOX_J2000,
        nutation=NutationModels.MEAN,
    ):
        self._obstime = obstime
        self._equinox = equinox
        self._nutation = NutationModels(nutation)

    def to_icrs(self, coords_self):
        if self.nutation is NutationModels.MEAN:
            coords_heliocentric_equatorial = BarycentricMeanEcliptic(
                self.equinox
            ).to_icrs(coords_self)
        elif self.nutation is NutationModels.TRUE:
            coords_heliocentric_equatorial = BarycentricTrueEcliptic(
                self.equinox
            ).to_icrs(coords_self)

        return HCRS(self.obstime).to_icrs(coords_heliocentric_equatorial)

    def from_icrs(self, coords_icrs):
        coords_heliocentric_equatorial = HCRS(self.obstime).from_icrs(coords_icrs)
        if self.nutation is NutationModels.MEAN:
            coords_self = BarycentricMeanEcliptic(self.equinox).from_icrs(
                coords_heliocentric_equatorial
            )
        elif self.nutation is NutationModels.TRUE:
            coords_self = BarycentricTrueEcliptic(self.equinox).from_icrs(
                coords_heliocentric_equatorial
            )

        return coords_self


class HeliocentricMeanEcliptic(HeliocentricEcliptic):
    def __init__(self, obstime=DEFAULT_OBSTIME, equinox=EQUINOX_J2000):
        super().__init__(obstime=obstime, equinox=equinox, nutation=NutationModels.MEAN)


class HeliocentricTrueEcliptic(HeliocentricEcliptic):
    def __init__(self, obstime=DEFAULT_OBSTIME, equinox=EQUINOX_J2000):
        super().__init__(obstime=obstime, equinox=equinox, nutation=NutationModels.TRUE)
