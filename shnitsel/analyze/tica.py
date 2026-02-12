from typing import Generic, Hashable, TypeVar

import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import Pipeline
import xarray as xr

from shnitsel.analyze.generic import norm

OriginType = TypeVar('OriginType')
ResultType = TypeVar('ResultType')


class TICAResult(Generic[OriginType, ResultType]):
    def __init__(
        self,
        inputs: OriginType,
        reduced_dimension: Hashable,
        scale_func,
        reducer_object,
        projected_inputs: ResultType,
    ):
        if isinstance(inputs, xr.DataArray):
            assert isinstance(projected_inputs, xr.DataArray), (
                "If inputs are provided as a single data array, the results must also be a single data array"
            )
            coord_initial = [
                projected_inputs.coords['PC'],  # FIXME: Misnamed for compatibility
                inputs.coords[reduced_dimension],
            ]
            n_components = projected_inputs.sizes['PC']
            components = xr.DataArray(
                reducer_object.model.instantaneous_coefficients.T[:n_components, :],
                coords=coord_initial,
            ).assign_coords(
                TICAResult.get_extra_coords_for_loadings(inputs, reduced_dimension)
            )
        else:
            raise NotImplementedError
            # TODO (thevro): Implement

        # FIXME (thevro): Currently using "pca_" names for compatibility with biplot etc.
        self.inputs = inputs
        self.pca_mapped_dimension = reduced_dimension
        self.principal_components = components
        self.projected_inputs = projected_inputs
        self.pca_object = reducer_object

        # Incompatible with PCAResult class
        self.pca_pipeline = NotImplemented
        self.scale_func = scale_func

        # Aliases
        self.loadings = self.principal_components
        self.results = self.projected_inputs

    def get_most_significant_loadings(
        self, top_n_per: int = 5, top_n_total: int = 5
    ) -> tuple[Mapping[Hashable, xr.DataArray], xr.DataArray]:
        """Function to retrieve the most significant loadings in the
        PCA result for each individual component and in total.

        You can configure the amount of

        Parameters
        ----------
        top_n_per : int, optional
            Number of top (most significant absolute loading) n loadings per component, by default 5
        top_n_total : int, optional
            Number of overall top (i.e. most significant by 2-norm of their loadings across all PC) n features across all components, by default 5

        Returns
        -------
        tuple[Mapping[Hashable, xr.DataArray], xr.DataArray]
            First the mapping of each PC to the array holding the data of all their most significant loadings.
            Second the overall most significant loadings across all components.
        """
        loadings = self.loadings

        per_pc_results = {}
        for pc in loadings.PC.values:
            component = loadings.sel(PC=pc)

            top_n_per_local = min(component.sizes[self.pca_mapped_dimension], top_n_per)

            abs_loading = np.abs(component)
            top_arg_indices = np.argpartition(abs_loading, -top_n_per_local)[
                -top_n_per_local:
            ]
            top_arg_coords = component.coords[self.pca_mapped_dimension].values[
                top_arg_indices
            ]

            per_pc_results[pc] = component.sel(
                {self.pca_mapped_dimension: top_arg_coords}
            )

        top_n_total_local = min(loadings.sizes[self.pca_mapped_dimension], top_n_total)
        total_abs_loadings = norm(loadings, dim='PC')

        top_arg_indices = np.argpartition(total_abs_loadings, -top_n_total_local)[
            -top_n_total_local:
        ]
        top_arg_coords = loadings.coords[self.pca_mapped_dimension].values[
            top_arg_indices
        ]

        total_pc_results = loadings.sel({self.pca_mapped_dimension: top_arg_coords})
        return per_pc_results, total_pc_results

    def explain_loadings(self, top_n_per: int = 5, top_n_total: int = 5) -> str:
        """Generate a textual explanation of the top influential loadings in the PCA result.

        Tries to put the results of `get_most_significant_loadings()` into a textual form.

        Parameters
        ----------
        top_n_per : int, optional
            Number of top (most significant absolute loading) n loadings per component, by default 5
        top_n_total : int, optional
            Number of overall top (i.e. most significant by 2-norm of their loadings across all PC) n features across all components, by default 5

        Returns
        -------
        str
            A text describing the results of the principal components analysis.
        """
        per_top, total_top = self.get_most_significant_loadings(
            top_n_per=top_n_per, top_n_total=top_n_total
        )

        explanation: str = ""

        total_expl = f"Maximum contributing features overall:\n"
        for feature, indices, coeff in zip(
            total_top.descriptor.values,
            total_top.feature_indices.values,
            norm(total_top, dim='PC').values,
        ):
            total_expl += f" {feature} (weight: {coeff}) (Idxs: {indices}) \n"
        explanation += total_expl + "\n\n"

        for pc in per_top:
            loadings = per_top[pc]

            pc_expl = f"Maximum contributing features to component {pc} :\n"
            for feature, indices, coeff in zip(
                loadings.descriptor.values,
                loadings.feature_indices.values,
                loadings.values,
            ):
                pc_expl += f" {feature}  (weight: {coeff}) (Idxs: {indices}) \n"
            explanation += pc_expl + "\n"
        return explanation

    def project_array(self, other_da: xr.DataArray) -> xr.DataArray:
        scaled = self.scale_func(other_da)
        return xr.apply_ufunc(
            reducer_object.transform,
            scaled,
            input_core_dims=[[dim]],
            output_core_dims=[['PC']],
        )

    @staticmethod
    def get_extra_coords_for_loadings(
        data: xr.DataArray, dim: Hashable
    ) -> Mapping[Hashable, xr.DataArray]:
        # TODO (thevro): this static method is copied wholesale from PCAResult
        # because I didn't want to mess around with inheritance yet
        coords = {
            key: coord
            for key, coord in data.coords.items()
            if dim in coord.dims and key != dim
        }
        return coords

    def __str__(self) -> str:
        return (
            type(self).__name__
            + f" on {type(self.inputs).__name__} with {self.principal_components.sizes['PC']} components"
        )

    def __repr__(self) -> str:
        return self.__str__() + "\n" + self.explain_loadings()


def tica_direct(
    data: xr.DataArray, dim: Hashable, lagtime: int, n_components: int = 2
) -> PCAResult:
    """Wrapper function to directly apply tICA decomposition to the values in a dataarray.

    Parameters
    ----------
    data : xr.DataArray
        A DataArray with at least a dimension with a name matching `dim`
    dim : Hashable
        The name of the array-dimension to reduce (i.e. the axis along which different
        features lie)
    n_components : int, optional
        The number of independent components to return, by default 2

    Returns
    -------
    TICAResult
        The full information obtained by the fitting of the result.
        Contains the inputs for the PCA result, the principal components,
        the mapped values for the inputs, the full pipeline to apply the PCA
        transformation again to other data.

        The mapped inputs are a DataArray with the same dimensions as ``da``, except for the dimension
        indicated by `dim`, which is replaced by a dimension ``component`` of size ``n_components``.

    .. !todo
        Examples:
        ---------
        >>> tica_results1 = tica(data1, 'features')
        >>> tica_results1.projected_inputs  # See the loadings
        >>> tica_results2 = tica_results1.project_array(data2)
    """
    try:
        from deeptime.util.data import TrajectoriesDataset
        from deeptime.decomposition import TICA
    except ModuleNotFoundError as err:
        print("Please install the deeptime package")
        raise err

    if not isinstance(data, xr.DataArray):
        raise NotImplementedError()  # TODO

    data = data.transpose('frame', ...)  # Required for TrajectoriesDataset
    scaler = MinMaxScaler().fit(data)

    scale_func = lambda da: xr.apply_ufunc(
        scaler.transform,
        da,
        input_core_dims=[[dim]],
        output_core_dims=[[dim]],
    )

    scaled = scale_func(data)
    deeptime_trajs = TrajectoriesDataset.from_numpy(
        lagtime=lagtime,
        data=[scaler.transform(x) for _, x in data.groupby('atrajectory')],
    )

    # NOTE: in order to treat trajectories separately,
    # we can't use the pipeline

    # pca_object = sk_PCA(n_components=n_components)
    reducer_object = TICA(lagtime=lagtime, dim=n_components)
    reducer_object.fit_from_timeseries(deeptime_trajs)

    # pipeline = Pipeline([('scaler', scaler), ('reducer', reducer_object)])

    # pca_res: xr.DataArray = xr.apply_ufunc(
    #     pipeline.fit_transform,
    #     data,
    #     input_core_dims=[[dim]],
    #     output_core_dims=[['PC']],
    # )
    tica_res: xr.DataArray = xr.apply_ufunc(
        reducer_object.transform,
        scaled,
        input_core_dims=[[dim]],
        output_core_dims=[['PC']],
    )

    # pca_result_wrapper = PCAResult(
    #     pca_inputs=data,
    #     pca_object=pipeline[-1],
    #     pca_dimension=dim,
    #     pca_projected_inputs=pca_res,
    #     pca_pipeline=pipeline,
    # )
    tica_result_wrapper = TICAResult(
        inputs=data,
        reducer_object=reducer_object,
        reduced_dimension=dim,
        projected_inputs=tica_res,
        # pca_pipeline=pipeline,
        scale_func=scale_func,
    )

    return tica_result_wrapper
