def transform(coords, orig_frame, dest_frame):
    coords_icrs = orig_frame.to_icrs(coords)
    coords_dest = dest_frame.from_icrs(coords_icrs)
    return coords_dest
