import numpy as np
from xarray.core.indexes import (
    Index,
    PandasIndex,
    PandasMultiIndex,
    IndexSelResult,
    normalize_label,
)
from xarray.core.indexing import merge_sel_results


class VectorizableMultiIndex(PandasMultiIndex):
    from collections.abc import Iterable

    def sel(self, labels, method=None, tolerance=None) -> IndexSelResult:
        if self.index.name in labels and not set(self.index.names).isdisjoint(labels):
            # xarray.indexes.PandasMultiIndex does this check in a more complicated way.
            # For now, this should suffice:
            raise ValueError(
                f"cannot provide labels for both coordinate {self.index.name!r} (multi-index array) "
                f"and one or more coordinates among {self.index.names!r} (multi-index levels)"
            )
        elif self.index.name in labels:
            return super().sel(labels, method=method, tolerance=tolerance)

        label_arrays = {
            k: normalize_label(v, dtype=self.level_coords_dtype[k])
            for k, v in labels.items()
        }
        if any(np.atleast_1d(v).shape[0] > 1 for k, v in label_arrays.items()):
            indexer = self.index.get_locs(
                [label_arrays.get(k, slice(None)) for k in self.index.names]
            )
            # print(indexer)
            new_index = type(self)(self.index[indexer], self.dim)
            variables = new_index.create_variables()
            return IndexSelResult(
                dim_indexers={self.dim: indexer},
                indexes={d: new_index for d in [self.dim] + self.index.names},
                variables=variables,
                # drop_indexes=list(scalar_coord_values),
                # drop_coords=drop_coords,
                # rename_dims=dims_dict,
            )
        else:
            tmp = super().sel(labels, method, tolerance)
            # print(tmp)
            # tmp.drop_indexes = []
            # tmp.indexes = []
            # tmp.variables = []
            return tmp
            # return super().sel(labels, method, tolerance)

    def __repr__(self):
        return f"VectorizableMultiIndex({self.index!r})"

    def reindex_like(self, other, method=None, tolerance=None):
        print('called reindex_like')
        print(f'{self=}')
        print(f'{other=}')
        res = super().reindex_like(other)
        print(res)
        return res


class MajMinIndex(Index):
    """Warning: Experimental. Like the major and minor ticks on a graph,
    this custom index governs a MultiIndex and an associated
    simple Index.

    Each level of the MultiIndex typically contains
    many repetitions of the same value, which are disambiguated
    by the other levels. If there is some information which is
    already uniquely defined by the value of one level on its own, it
    might be convenient to put this information along a separate
    dimension. For example, if you have a dataset with multiple compounds,


    """

    def __init__(self, full_midx, sparse_idxs):
        # I think the call signature of __init__
        # is not constrained by xarray
        print('called init')
        self._full_midx = full_midx
        self._sparse_idxs = sparse_idxs

    @classmethod
    def from_variables(cls, variables, **kws):
        # The call signature of this class method
        # is constrained by xarray expectations
        print('called from_variables')
        full = {}
        sparse = {}

        for k, v in variables.items():
            if k.endswith('_'):
                sparse[k] = v
            else:
                full[k] = v

        # full_midx = PandasMultiIndex.from_variables(full, **kws)
        full_midx = VectorizableMultiIndex.from_variables(full, **kws)

        sparse_idxs = {
            k: PandasIndex.from_variables({k: v}, **kws) for k, v in sparse.items()
        }

        # print(kws)
        # print(full_midx)
        # print(sparse_idxs)

        return cls(full_midx, sparse_idxs)

    def create_variables(self, variables):
        # print(f'called create_variables({variables})')
        print('called create_variables')
        idx_variables = {}

        for index in self._sparse_idxs.values():
            idx_variables.update(index.create_variables(variables))

        return idx_variables

    def sel(self, labels):
        print('called sel')
        print(labels)
        results = []

        # transform boolean masks into label-based indexers
        # using the appropriate sparse index
        for k, indexer in labels.items():
            arr = np.asarray(indexer)
            if np.issubdtype(arr.dtype, np.bool_):
                print("Found bool index")
                print(self._sparse_idxs[k + '_'].sel({k: indexer}))
                labels[k] = self._sparse_idxs[k + '_'].to_pandas_index().values[arr]

        # res_sparse_idxs = {}
        new_sparse_idxs = {}
        for k, index in self._sparse_idxs.items():
            if k[:-1] in labels:
                this_res = index.sel({k: labels[k[:-1]]})
                results.append(this_res)

                # actually construct a suitable index for the selection result to inherit
                # we need to do this because the selection results for a PandasIndex
                # don't usually contain indexes; the recipient function somehow makes do without
                new_sparse_idxs[k] = self._sparse_idxs[k][this_res.dim_indexers[k]]

        res_full_midx = self._full_midx.sel(labels)
        new_full_midx = res_full_midx.indexes[self._full_midx.index.name]
        results.append(res_full_midx)

        res = merge_sel_results(results)

        new_index = type(
            self
        )(
            full_midx=new_full_midx,
            sparse_idxs=new_sparse_idxs,  # {k: v.indexes[k] for k, v in res_sparse_idxs.items()}
        )
        indexes = {
            d: new_index
            for d in [self._full_midx.index.name]
            + self._full_midx.index.names
            + list(self._sparse_idxs)
        }  # TODO: tidy
        res.indexes = indexes

        # res.drop_indexes = []
        # res.drop_coords = []

        # print(res.dim_indexers)
        # print(res)
        return res