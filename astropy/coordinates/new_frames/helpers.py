from astropy.coordinates import CartesianRepresentation


def apply_affine(coords, left_rotation=None, offset=None, right_rotation=None):
    if left_rotation is not None:
        coords = coords.represent_as(CartesianRepresentation).transform(left_rotation)

    if offset is not None:
        coords = coords + offset

    if right_rotation is not None:
        coords = coords.represent_as(CartesianRepresentation).transform(right_rotation)

    return coords
