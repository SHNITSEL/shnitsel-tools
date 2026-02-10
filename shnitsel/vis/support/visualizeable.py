import abc
from typing import Self, Collection

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from IPython.lib.pretty import PrettyPrinter
from matplotlib.pylab import Any


class Visualizable(abc.ABC):
    @abc.abstractmethod
    def plot(self, *args, **kwargs) -> Figure | Axes | Self:
        raise NotImplementedError(".plot() is not implemented for type {type(self)}")

    def _repr_pretty_(self, p: PrettyPrinter, cycle: bool):
        """Pretty-print an object using the pretty-printer
        instance `p`.

        The Flag `cycle` indicates whether a cycle has been detected.
        If it is true, the object should avoid recursing further.

        Parameters
        ----------
        p : PrettyPrinter
            The pretty printer to use to output the pretty-printed representation.
        cycle : bool
            Flag to indicate whether a cycle in references has been detected.
        """
        if cycle:
            p.text('MyObject(...)')
        else:
            p.text('MyObject[...]')

    def _repr_svg_(self) -> str | tuple[str, dict[str, Any]]:
        """Function to obtain an svg string representing the object

        Returns
        -------
        str
            The string of the svg.
        """
        return ""

    def _repr_png_(self) -> bytes | tuple[bytes, dict[str, Any]]:
        """Function to obtain a byte stream of a png string representing the object.

        Returns
        -------
        bytes
            The bytes of the associated png
        tuple[bytes, dict[str, Any]]
            First the bytes of the png and then a metadata dict for the displaying frontend
        """
        return bytes(0)

    def _repr_jpeg_(self) -> bytes | tuple[bytes, dict[str, Any]]:
        return bytes(0)

    def _repr_html_(self) -> str | tuple[str, dict[str, Any]]:
        return ""

    def _repr_javascript_(self) -> str | tuple[str, dict[str, Any]]:
        return ""

    def _repr_markdown_(self) -> str | tuple[str, dict[str, Any]]:
        return ""

    def _repr_latex_(self) -> str | tuple[str, dict[str, Any]]:
        return ""

    def _repr_mimebundle_(
        self,
        include: Collection[str] | None = None,
        exclude: Collection[str] | None = None,
    ) -> dict[str, Any] | tuple[dict[str, Any], dict[str, Any]]:
        """Offer all available options to the callling environment in a `mime-type -> data` dict result.
        Alternatively may return one data and one metadata dict as a tuple

        Parameters
        ----------
        include : Collection[str], optional
            Requested mime types to be included in result, by default None. May be ignored.
        exclude : Collection[str], optional
            Mime types that should not be included in the result, by default None. May be ignored.

        Returns
        -------
        dict[str, Any]
            Dict holding the `mime-type -> data` mapping
        tuple[dict[str, Any], dict[str, Any]]
            The `mime-type -> data` dict with a separate visualization metadata dict (like `height` or `width` for images)
        """
        return None

    # def _ipython_display_():
    #     """Directly display an object. If this is present, all other display functions will be ignored"""
    #     ...
