from typing import Any, Collection, TypeAlias
from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure
from shnitsel.data.tree.data_leaf import DataLeaf
from shnitsel.data.tree.node import TreeNode
from shnitsel.vis.support.visualizeable import Visualizable

IndividualPlotsources: TypeAlias = Figure | SubFigure | Axes | Visualizable


class MultiPlot:
    _individual_plots: (
        Collection[IndividualPlotsources] | TreeNode[Any, IndividualPlotsources]
    )

    def __init__(
        self,
        individual_plots: Collection[IndividualPlotsources]
        | TreeNode[Any, IndividualPlotsources],
    ):
        self._individual_plots = individual_plots

    def plot(self):
        if isinstance(self._individual_plots, TreeNode):
            return self._individual_plots.map_filtered_nodes(
                lambda x: x.is_leaf,
                lambda x: DataLeaf(
                    data=(
                        x.path,
                        MultiPlot.plot_individual_source(x.data)
                        if x is not None
                        else None,
                    )
                ),
            ).collect_data()
        else:
            return [MultiPlot.plot_individual_source(x) for x in self._individual_plots]

    @staticmethod
    def plot_individual_source(source: IndividualPlotsources):
        if isinstance(source, Visualizable):
            return source.plot()
        else:
            return source
