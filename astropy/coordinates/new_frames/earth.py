from astropy import _erfa as erfa
from astropy import units as u
from astropy.coordinates import CartesianRepresentation, SphericalRepresentation
from astropy.coordinates.builtin_frames.utils import (
    DEFAULT_OBSTIME,
    get_cip,
    get_jd12,
    prepare_earth_position_vel,
)

from .base import HasObstime, HasObserver, BaseCoordinateFrame


class CIRS(BaseCoordinateFrame, HasObstime):
    def __init__(self, obstime=DEFAULT_OBSTIME):
        self._obstime = obstime

    def to_icrs(self, coords_self):
        srepr = coords_self.represent_as(SphericalRepresentation)
        cirs_ra = srepr.lon.to_value(u.radian)
        cirs_dec = srepr.lat.to_value(u.radian)

        # set up the astrometry context for ICRS<->cirs and then convert to
        # astrometric coordinate direction
        jd1, jd2 = get_jd12(self.obstime, "tt")
        x, y, s = get_cip(jd1, jd2)
        earth_pv, earth_heliocentric = prepare_earth_position_vel(self.obstime)
        astrom = erfa.apci(jd1, jd2, earth_pv, earth_heliocentric, x, y, s)
        i_ra, i_dec = erfa.aticq(cirs_ra, cirs_dec, astrom)

        # When there is a distance, apply the parallax/offset to the SSB as the
        # last step - ensures round-tripping with the icrs_to_cirs transform

        # the distance in intermedrep is *not* a real distance as it does not
        # include the offset back to the SSB
        intermedrep = SphericalRepresentation(
            lat=u.Quantity(i_dec, u.radian, copy=False),
            lon=u.Quantity(i_ra, u.radian, copy=False),
            distance=srepr.distance,
            copy=False,
        )

        astrom_eb = CartesianRepresentation(
            astrom["eb"], unit=u.au, xyz_axis=-1, copy=False
        )

        return intermedrep + astrom_eb

    def from_icrs(self, coords_icrs):
        jd1, jd2 = get_jd12(self.obstime, "tt")
        x, y, s = get_cip(jd1, jd2)
        earth_pv, earth_heliocentric = prepare_earth_position_vel(self.obstime)
        astrom = erfa.apci(jd1, jd2, earth_pv, earth_heliocentric, x, y, s)

        # When there is a distance, we first offset for parallax to get the
        # astrometric coordinate direction and *then* run the ERFA transform for
        # no parallax/PM. This ensures reversibility and is more sensible for
        # inside solar system objects
        astrom_eb = CartesianRepresentation(
            astrom["eb"], unit=u.au, xyz_axis=-1, copy=False
        )
        newcart = coords_icrs.represent_as(CartesianRepresentation) - astrom_eb

        srepr = newcart.represent_as(SphericalRepresentation)
        i_ra = srepr.lon.to_value(u.radian)
        i_dec = srepr.lat.to_value(u.radian)
        cirs_ra, cirs_dec = erfa.atciqz(i_ra, i_dec, astrom)

        newrep = SphericalRepresentation(
            lat=u.Quantity(cirs_dec, u.radian, copy=False),
            lon=u.Quantity(cirs_ra, u.radian, copy=False),
            distance=srepr.distance,
            copy=False,
        )
        return newrep


class GCRS(BaseCoordinateFrame, HasObstime, HasObserver):
    def __init__(
        self,
        obstime=DEFAULT_OBSTIME,
        obsgeoloc=CartesianRepresentation([0, 0, 0], unit=u.m),
        obsgeovel=CartesianRepresentation([0, 0, 0], unit=u.m / u.s),
    ):
        self._obstime = obstime
        self._obsgeoloc = obsgeoloc
        self._obsgeovel = obsgeovel

    def to_icrs(self, coords_self):
        srepr = coords_self.represent_as(SphericalRepresentation)
        gcrs_ra = srepr.lon.to_value(u.radian)
        gcrs_dec = srepr.lat.to_value(u.radian)

        # set up the astrometry context for ICRS<->GCRS and then convert to BCRS
        # coordinate direction
        obs_pv = erfa.pav2pv(
            self.obsgeoloc.get_xyz(xyz_axis=-1).to_value(u.m),
            self.obsgeovel.get_xyz(xyz_axis=-1).to_value(u.m / u.s),
        )

        jd1, jd2 = get_jd12(self.obstime, "tt")
        earth_pv, earth_heliocentric = prepare_earth_position_vel(self.obstime)
        astrom = erfa.apcs(jd1, jd2, obs_pv, earth_pv, earth_heliocentric)
        i_ra, i_dec = erfa.aticq(gcrs_ra, gcrs_dec, astrom)

        # When there is a distance, apply the parallax/offset to the SSB as the
        # last step - ensures round-tripping with the icrs_to_gcrs transform

        # the distance in intermedrep is *not* a real distance as it does not
        # include the offset back to the SSB
        intermedrep = SphericalRepresentation(
            lat=u.Quantity(i_dec, u.radian, copy=False),
            lon=u.Quantity(i_ra, u.radian, copy=False),
            distance=srepr.distance,
            copy=False,
        )

        astrom_eb = CartesianRepresentation(
            astrom["eb"], unit=u.au, xyz_axis=-1, copy=False
        )
        newrep = intermedrep + astrom_eb

        return newrep

    def from_icrs(self, coords_icrs):
        # first set up the astrometry context for ICRS<->GCRS. There are a few steps...
        # get the position and velocity arrays for the observatory.  Need to
        # have xyz in last dimension, and pos/vel in one-but-last.
        # (Note could use np.stack once our minimum numpy version is >=1.10.)
        obs_pv = erfa.pav2pv(
            self.obsgeoloc.get_xyz(xyz_axis=-1).to_value(u.m),
            self.obsgeovel.get_xyz(xyz_axis=-1).to_value(u.m / u.s),
        )

        # find the position and velocity of earth
        jd1, jd2 = get_jd12(self.obstime, "tt")
        earth_pv, earth_heliocentric = prepare_earth_position_vel(self.obstime)
        astrom = erfa.apcs(jd1, jd2, obs_pv, earth_pv, earth_heliocentric)

        # When there is a distance,  we first offset for parallax to get the
        # BCRS coordinate direction and *then* run the ERFA transform for no
        # parallax/PM. This ensures reversibility and is more sensible for
        # inside solar system objects
        astrom_eb = CartesianRepresentation(
            astrom["eb"], unit=u.au, xyz_axis=-1, copy=False
        )
        newcart = coords_icrs.represent_as(CartesianRepresentation) - astrom_eb

        srepr = newcart.represent_as(SphericalRepresentation)
        i_ra = srepr.lon.to_value(u.radian)
        i_dec = srepr.lat.to_value(u.radian)
        gcrs_ra, gcrs_dec = erfa.atciqz(i_ra, i_dec, astrom)

        newrep = SphericalRepresentation(
            lat=u.Quantity(gcrs_dec, u.radian, copy=False),
            lon=u.Quantity(gcrs_ra, u.radian, copy=False),
            distance=srepr.distance,
            copy=False,
        )
        return newrep
