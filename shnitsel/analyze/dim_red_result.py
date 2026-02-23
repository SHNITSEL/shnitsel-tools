import xarray as xr

from shnitsel.data.tree.node import TreeNode
from .generic import norm
import abc
from typing import Any, Generic, Hashable, Mapping, TypeVar, overload

import numpy as np
from sklearn.pipeline import Pipeline

from shnitsel.core.typedefs import DimName

OriginType = TypeVar('OriginType')
ResultType = TypeVar('ResultType')
PredictionType = TypeVar('PredictionType')


class DimRedResult(
    Generic[OriginType, ResultType],
    abc.ABC,
):
    """Base class to hold the results of a dimensionality reduction like PCA, LDA or PLS analysis.

    Also retains input data as well as corresponding results of the decomposition.
    Input and output types are parametrized to allow for tree structures to be accurately represented.

    Provides accessors for all result meta data as well as the method `project_array(data_array)` to
    project another array of appropriate shape with dimension `mapped_dimension` to the respective representation

    Parameters
    ----------
    OriginType
        The type of the original intput data. Should either be xr.DataArray for simple types, meaning we were provided a feature array
        or a flat DataGroup with xr.DataArrays in its leaves for tree types.
    ResultType
        Matching structure to `OriginType` but with the projected and decomposed input data as data within it.
        Either an xr.DataArray or a DataGroup same as for `OriginType`.
    """

    _dimred_label: str
    _inputs: OriginType
    _pipeline: Pipeline
    _mapped_dimension: DimName
    _loadings: xr.DataArray
    _component_dimension: DimName
    _projected_inputs: ResultType

    def __init__(
        self,
        dimred_label: str,
        inputs: OriginType,
        mapped_dimension: DimName,
        pipeline: Pipeline,
        loadings: xr.DataArray,
        component_dimension: DimName,
        projected_inputs: ResultType,
    ):
        self._dimred_label = dimred_label
        self._inputs = inputs
        self._mapped_dimension = mapped_dimension
        self._pipeline = pipeline
        self._loadings = loadings
        self._component_dimension = component_dimension
        self._projected_inputs = projected_inputs

    @property
    def label(self) -> str:
        return self._dimred_label

    @property
    def inputs(self) -> OriginType:
        return self._inputs

    @property
    def mapped_dimension(self) -> DimName:
        return self._mapped_dimension

    @property
    def pipeline(self) -> Pipeline:
        return self._pipeline

    @property
    def components(self) -> xr.DataArray:
        return self.loadings

    @property
    def num_components(self) -> int:
        return self.loadings.sizes[self._component_dimension]

    @property
    def loadings(self) -> xr.DataArray:
        return self._loadings

    @property
    def component_dimension(self) -> DimName:
        return self._component_dimension

    @property
    def projected_inputs(self) -> ResultType:
        return self._projected_inputs

    @property
    def results(self) -> ResultType:
        return self.projected_inputs

    def get_most_significant_loadings(
        self, top_n_per: int = 5, top_n_total: int = 5
    ) -> tuple[Mapping[Hashable, xr.DataArray], xr.DataArray]:
        """Function to retrieve the most significant loadings in the
        dimensionality reduction result for each individual component and in total.

        You can configure the amount of

        Parameters
        ----------
        top_n_per : int, optional
            Number of top (most significant absolute loading) n loadings per component, by default 5
        top_n_total : int, optional
            Number of overall top (i.e. most significant by 2-norm of their loadings across all components)
            n features across all components, by default 5

        Returns
        -------
        tuple[Mapping[Hashable, xr.DataArray], xr.DataArray]
            First the mapping of each PC to the array holding the data of all their most significant loadings.
            Second the overall most significant loadings across all components.
        """
        loadings = self.loadings

        per_component_results = {}
        for component_label in loadings[self.component_dimension].values:
            component = loadings.sel({self.component_dimension: component_label})
            # print(component)
            # print(component.values)

            top_n_per_local = min(component.sizes[self.mapped_dimension], top_n_per)

            abs_loading = np.abs(component)
            top_arg_indices = np.argpartition(abs_loading, -top_n_per_local)[
                -top_n_per_local:
            ]
            top_arg_coords = component.coords[self.mapped_dimension].values[
                top_arg_indices
            ]

            # print(top_arg_indices)
            # print(top_arg_coords)

            per_component_results[component_label] = component.sel(
                {self.mapped_dimension: top_arg_coords}
            )

        top_n_total_local = min(loadings.sizes[self.mapped_dimension], top_n_total)
        total_abs_loadings = norm(loadings, dim=self.component_dimension)

        top_arg_indices = np.argpartition(total_abs_loadings, -top_n_total_local)[
            -top_n_total_local:
        ]
        top_arg_coords = loadings.coords[self.mapped_dimension].values[top_arg_indices]

        # print(top_arg_indices)
        # print(top_arg_coords)
        total_pc_results = loadings.sel({self.mapped_dimension: top_arg_coords})
        # print(component.feature_indices)
        return per_component_results, total_pc_results

    def explain_loadings(self, top_n_per: int = 5, top_n_total: int = 5) -> str:
        """Generate a textual explanation of the top influential loadings in the dimensionality reduction result.

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
            norm(total_top, dim=self.component_dimension).values,
        ):
            total_expl += f" {feature} (weight: {coeff}) (Idxs: {indices}) \n"
        explanation += total_expl + "\n\n"

        for component_label in per_top:
            loadings = per_top[component_label]

            component_expl = (
                f"Maximum contributing features to component {component_label} :\n"
            )
            for feature, indices, coeff in zip(
                loadings.descriptor.values,
                loadings.feature_indices.values,
                loadings.values,
            ):
                component_expl += f" {feature}  (weight: {coeff}) (Idxs: {indices}) \n"
            explanation += component_expl + "\n"
        return explanation

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

        if other_da.sizes[self.mapped_dimension] == 0:
            return xr.zeros_like(other_da, dtype=float)

        return xr.apply_ufunc(
            self.pipeline.transform,
            other_da,
            input_core_dims=[[self.mapped_dimension]],
            output_core_dims=[[self.component_dimension]],
        )

    @staticmethod
    def get_extra_coords_for_loadings(
        data: xr.DataArray, dim: Hashable
    ) -> Mapping[Hashable, xr.DataArray]:
        """Helper function to assign dropped coordinates back to
        the mapped data.

        Parameters
        ----------
        data : xr.DataArray
            The input data for the PCA from which to extract the additional coordinates
        dim : Hashable
            The dimension along which we want to gather all coordinates that may have been dropped in the pipeline

        Returns
        -------
        Mapping[Hashable, xr.DataArray]
            Mapping from coordinate name to the data array holding the coordinate values/labels
        """
        coords = {
            key: coord
            for key, coord in data.coords.items()
            if dim in coord.dims and key != dim
        }
        return coords

    def __str__(self) -> str:
        return (
            type(self).__name__
            + f" on {type(self.inputs).__name__} with {self.components.sizes[self.component_dimension]} components"
        )

    def __repr__(self) -> str:
        return self.__str__() + "\n" + self.explain_loadings()


class PredictorDimRedResult(
    Generic[OriginType, ResultType, PredictionType],
    DimRedResult[OriginType, ResultType],
):
    """Specialization of the general dimensionality reduction result
    encapsulating a reduction that provides predictive behavior.

    """

    @overload
    def predict(
        self, other_da: TreeNode[Any, xr.DataArray]
    ) -> TreeNode[Any, PredictionType]: ...
    @overload
    def predict(self, other_da: xr.DataArray) -> PredictionType: ...

    def predict(
        self, other_da: xr.DataArray | TreeNode[Any, xr.DataArray]
    ) -> PredictionType | TreeNode[Any, PredictionType]:
        """Apply the transformation encoded by this dimensionality reduction
        to the provided (set of) DataArray(s) and return the prediction of the
        underlying predictor

        Parameters
        ----------
        other_da : xr.DataArray | TreeNode[Any, xr.DataArray]
            The data to apply the transformation to and predict classification for.

        Returns
        -------
        PredictionType | TreeNode[Any, PredictionType]
            The transformed data
        """
        if isinstance(other_da, TreeNode):
            return other_da.map_data(self._predict_single)

        return self._predict_single(other_da=other_da)

    def _predict_single(self, other_da: xr.DataArray) -> PredictionType:
        if other_da.sizes[self.mapped_dimension] == 0:
            return self._empty_prediction(other_da)

        return xr.apply_ufunc(
            self.pipeline.predict,
            other_da,
            input_core_dims=[[self.mapped_dimension]],
            output_core_dims=[[self.component_dimension]],
        )

    def _empty_prediction(self, other_da: xr.DataArray) -> PredictionType:
        """Helper method to yield an empty prediction result for empty arrays

        Parameters
        ----------
        other_da : xr.DataArray
            The empty array for which a prediction was requested

        Returns
        -------
        PredictionType
            The empty prediction
        """
        return xr.zeros_like(other_da, dtype=int)
