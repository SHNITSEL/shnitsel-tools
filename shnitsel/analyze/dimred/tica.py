import logging
from typing import TYPE_CHECKING, Any, Generic, Hashable, TypeVar, overload

import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import Pipeline
import xarray as xr

from .dim_red_result import DimRedResult
from shnitsel.analyze.generic import norm
from shnitsel.core.typedefs import DimName
from shnitsel.data.tree.node import TreeNode

OriginType = TypeVar('OriginType')
ResultType = TypeVar('ResultType')

if TYPE_CHECKING:
    from deeptime.decomposition import TICA as dt_TICA


class TICAResult(
    Generic[OriginType, ResultType],
    DimRedResult[OriginType, ResultType],
):
    _scaler_object = MinMaxScaler
    _tica_object: "dt_TICA"

    def __init__(
        self,
        inputs: OriginType,
        reduced_dimension: DimName,
        scaler: MinMaxScaler,
        reducer_object: "dt_TICA",
        projected_inputs: ResultType,
        pipeline: Pipeline,
    ):
        if isinstance(inputs, xr.DataArray):
            assert isinstance(projected_inputs, xr.DataArray), (
                "If inputs are provided as a single data array, the results must also be a single data array"
            )
            coord_initial = [
                projected_inputs.coords['TIC'],
                inputs.coords[reduced_dimension],
            ]
            n_components = projected_inputs.sizes['TIC']
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
        super().__init__(
            dimred_label='TICA',
            inputs=inputs,
            mapped_dimension=reduced_dimension,
            pipeline=pipeline,
            loadings=components,
            component_dimension="TIC",
            projected_inputs=projected_inputs,
        )

        # Incompatible with PCAResult class
        # self.pca_pipeline = NotImplemented

        self._scaler_object = scaler
        self._tica_object = reducer_object

        # Aliases
        # NOTE: Removed due to new DimRed interface
        # self.loadings = self.principal_components
        # self.results = self.projected_inputs

    @overload
    def project_array(
        self, other_da: TreeNode[Any, xr.DataArray]
    ) -> TreeNode[Any, xr.DataArray]: ...
    @overload
    def project_array(self, other_da: xr.DataArray) -> xr.DataArray: ...

    def project_array(
        self, other_da: xr.DataArray | TreeNode[Any, xr.DataArray]
    ) -> xr.DataArray | TreeNode[Any, xr.DataArray]:
        """Apply the transformation encoded by this dimensionality reduction
        to the provided (set of) DataArray(s).

        Parameters
        ----------
        other_da : xr.DataArray | TreeNode[Any, xr.DataArray]
            The data to apply the transformation to.

        Returns
        -------
        xr.DataArray | TreeNode[Any, xr.DataArray]
            The transformed data
        """
        if isinstance(other_da, TreeNode):
            return other_da.map_data(self.project_array, dtype=xr.DataArray)

        if isinstance(other_da, xr.DataArray):
            if other_da.size == 0:
                return xr.zeros_like(other_da, dtype=float)

            return self.pipeline.transform(other_da)
        else:
            # numpy array
            if len(other_da.shape) < 2 or other_da.shape[-2] < 1:
                return xr.zeros_like(other_da, dtype=float)

            return self.pipeline.transform(other_da)

    @property
    def fitted_tica_object(self) -> "dt_TICA":
        return self._tica_object


def tica_direct(
    data: xr.DataArray, dim: Hashable, lagtime: int, n_components: int = 2
) -> TICAResult:
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

    pipeline, transformer = _build_tica_pipeline(
        n_components=n_components, lagtime=lagtime, reduce_dim=dim
    )

    # data = data.transpose('frame', ...)  # Required for TrajectoriesDataset
    # scaler = MinMaxScaler().fit(data)

    # scale_func = lambda da: xr.apply_ufunc(
    #     scaler.transform,
    #     da,
    #     input_core_dims=[[dim]],
    #     output_core_dims=[[dim]],
    # )

    # scaled = scale_func(data)
    # deeptime_trajs = TrajectoriesDataset.from_numpy(
    #     lagtime=lagtime,
    #     data=[scaler.transform(x) for _, x in data.groupby('atrajectory')],
    # )

    # # NOTE: in order to treat trajectories separately,
    # # we can't use the pipeline

    # # pca_object = sk_PCA(n_components=n_components)
    # reducer_object = TICA(lagtime=lagtime, dim=n_components)
    # reducer_object.fit_from_timeseries(deeptime_trajs)

    # pipeline = Pipeline([('scaler', scaler), ('reducer', reducer_object)])

    # pca_res: xr.DataArray = xr.apply_ufunc(
    #     pipeline.fit_transform,
    #     data,
    #     input_core_dims=[[dim]],
    #     output_core_dims=[['PC']],
    # )

    tica_res: xr.DataArray = pipeline.fit_transform(data)

    # pca_result_wrapper = PCAResult(
    #     pca_inputs=data,
    #     pca_object=pipeline[-1],
    #     pca_dimension=dim,
    #     pca_projected_inputs=pca_res,
    #     pca_pipeline=pipeline,
    # )
    tica_result_wrapper = TICAResult(
        inputs=data,
        reducer_object=transformer.reducer_object,
        reduced_dimension=dim,
        projected_inputs=tica_res,
        pipeline=pipeline,
        scaler=transformer.scaler,
    )

    return tica_result_wrapper


from sklearn.base import BaseEstimator, TransformerMixin


class TICATransformer(BaseEstimator, TransformerMixin):
    scaler: MinMaxScaler
    reducer_object: "dt_TICA"
    reduce_dim: DimName
    # n_components:int
    lagtime: int

    def __init__(self, n_components: int, lagtime: int, reduce_dim: DimName) -> None:
        try:
            from deeptime.decomposition import TICA as dt_TICA
        except ModuleNotFoundError as err:
            print("Please install the deeptime package")
            raise err
        
        super().__init__()
        self.scaler = MinMaxScaler()
        self.lagtime = lagtime
        self.reducer_object = dt_TICA(lagtime=lagtime, dim=n_components)
        self.reduce_dim = reduce_dim

    def scale_func(self, da: xr.DataArray | np.ndarray) -> xr.DataArray | np.ndarray:
        """Apply the initial scaling operation to the data array depending on its type"""
        if isinstance(da, xr.DataArray):
            return xr.apply_ufunc(
                self.scaler.transform,
                da,
                input_core_dims=[[self.reduce_dim]],
                output_core_dims=[[self.reduce_dim]],
            )
        else:
            return self.scaler.transform(da)

    def fit(self, X: xr.DataArray | np.ndarray, y=None):
        try:
            from deeptime.util.data import TrajectoriesDataset
            from deeptime.decomposition import TICA
        except ModuleNotFoundError:
            print(
                "Please install the deeptime package. It is required for TICA operations."
            )
            raise

        assert y is None, (
            "TICA transformer does not allow for target values to be fitted. It performs unsupervised transformation."
        )

        if isinstance(X, xr.DataArray):
            data = X.transpose('frame', ...)  # Required for TrajectoriesDataset
            self.scaler.fit(data)

            deeptime_trajs = TrajectoriesDataset.from_numpy(
                lagtime=self.lagtime,
                data=[self.scaler.transform(x) for _, x in data.groupby('atrajectory')],
            )
            self.reducer_object.fit_from_timeseries(deeptime_trajs)
        else:
            logging.warning("TICA analysis on unlabeled numpy arrays is not reliable.")
            data = X
            self.scaler.fit(data)
            deeptime_trajs = TrajectoriesDataset.from_numpy(
                lagtime=self.lagtime,
                data=self.scaler.transform(data).tolist(),
            )
            self.reducer_object.fit_from_timeseries(deeptime_trajs)

        # pipeline = Pipeline([('scaler', scaler), ('reducer', reducer_object)])

        # pca_res: xr.DataArray = xr.apply_ufunc(
        #     pipeline.fit_transform,
        #     data,
        #     input_core_dims=[[dim]],
        #     output_core_dims=[['PC']],
        # )
        return self  # The fit method typically does nothing for transformers

    def transform(self, X):
        # Your transformation logic goes here

        scaled = self.scale_func(X)
        tica_res: xr.DataArray | np.ndarray
        if isinstance(scaled, xr.DataArray):
            tica_res = xr.apply_ufunc(
                self.reducer_object.transform,
                scaled,
                input_core_dims=[[self.reduce_dim]],
                output_core_dims=[['TIC']],
            )
        else:
            tica_res = self.reducer_object.transform(scaled)

        return tica_res


def _build_tica_pipeline(
    n_components: int, lagtime: int, reduce_dim: DimName
) -> tuple[Pipeline, TICATransformer]:
    transformer = TICATransformer(n_components, lagtime, reduce_dim)
    return Pipeline([('TICA', transformer)]), transformer


tica = tica_direct
TICA = tica

time_lagged_independent_component_analysis = tica
