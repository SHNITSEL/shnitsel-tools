from bokeh.layouts import column
from bokeh.models import ColumnDataSource  # , CustomJS
from bokeh.plotting import figure
from bokeh.io import show, output_notebook
from bokeh.settings import settings


class FrameSelector:
    selected_frame_indices: list[int] = []

    def __init__(self, df, xname, yname, title, allowed_ws_origin=None, webgl=True):
        output_notebook()
        if allowed_ws_origin is not None:
            if isinstance(allowed_ws_origin, str):
                allowed_ws_origin = [allowed_ws_origin]
            settings.allowed_ws_origin.set_value(allowed_ws_origin)

        def bkapp(doc):
            nonlocal df
            df = df.copy()
            source = ColumnDataSource(data=df)

            plot = figure(
                tools='lasso_select',  # type: ignore
                title=title,
                output_backend='webgl' if webgl else 'canvas',
            )
            plot.scatter(xname, yname, source=source, selection_color='red')

            def callback(attr, old, new):
                nonlocal self
                self.selected_frame_indices = new
                self.df_selection = df.iloc[new, :]

            source.selected.on_change('indices', callback)

            doc.add_root(column(plot))

        show(bkapp)

    # @property
    # def allowed_ws_origin(self, value):
