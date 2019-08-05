class BaseCoordinateFrame:
    def to_icrs(self, coords_self):
        raise NotImplementedError

    def from_icrs(self, coords_icrs):
        raise NotImplementedError


class HasObstime:
    @property
    def obstime(self):
        return self._obstime


class HasEquinox:
    @property
    def equinox(self):
        return self._equinox


class HasNutation:
    @property
    def nutation(self):
        return self._nutation


class HasObserver:
    @property
    def obsgeoloc(self):
        return self._obsgeoloc

    @property
    def obsgeovel(self):
        return self._obsgeovel
