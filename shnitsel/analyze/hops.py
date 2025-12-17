import numpy as np


def get_hops(frames, hop_types: list[tuple[int, int]] | None = None):
    is_hop = frames.astate.st.mdiff() != 0

    res = frames.isel(frame=is_hop)
    res = res.assign_coords(
        hop_from=(frames.astate.shift({'frame': 1}, -1).isel(frame=is_hop)),
        hop_to=res.astate,
    )
    if hop_types is not None:
        acc = np.full(res.sizes['frame'], False)
        for hop_from, hop_to in hop_types:
            acc |= (res.hop_from == hop_from) & (res.hop_to == hop_to)
        res = res.isel(frame=acc)
    return res.drop_dims(['trajid_'], errors='ignore')
