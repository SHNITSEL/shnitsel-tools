import abc
from typing import Self

from matplotlib.axes import Axes
from matplotlib.figure import Figure


class Visualizable(abc.ABC):
    @abc.abstractmethod
    def plot(self, *args, **kwargs) -> Figure | Axes | Self:
        raise NotImplementedError(".plot() is not implemented for type {type(self)}")
