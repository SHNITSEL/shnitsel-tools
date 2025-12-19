

def _positions(atXYZ:xr.DataArray, atom_idxs:Sequence):
    # FIXME: Unclear, whether this expects indices (-> .isel) or the numbers of the atoms (-> sel).
    # FIXME: Why does it call zip?
    return [
        atXYZ.sel(atom=list(idxs))
        .drop(['atNames', 'atNums'], errors='ignore')
        .rename(atom='descriptor')
        for idxs in zip(*atom_idxs)
    ]