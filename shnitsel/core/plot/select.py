import os
from logging import warning

from bokeh.layouts import column
from bokeh.models import ColumnDataSource, CustomJS, Div
from bokeh.plotting import figure
from bokeh.io import show, output_notebook, curdoc
from bokeh.settings import settings

import pandas as pd
import xarray as xr

class FrameSelector:
    selected_frame_indices: list[int] = []

    def __init__(self, df_or_da, xname=None, yname=None, title="", allowed_ws_origin=None, webgl=True):
        output_notebook()

        if isinstance(df_or_da, pd.DataFrame):
            df = df_or_da
            da = None
        elif isinstance(df_or_da, xr.DataArray):
            da = df_or_da
            if len(da.dims) == 2:
                df = da.to_pandas()
            else:
                raise ValueError(
                    "When the first argument to FrameSelector is an "
                    "xarray.DataArray, it should have 2 dimensions, "
                    f"rather than {len(da.dims)} dimensions (namely {da.dims})."
                )
        else:
            raise TypeError(
                "The first argument to FrameSelector should be a "
                "pandas.DataFrame or a 2-dimensional xarray.DataArray"
            )

        # Column names must be strings
        for col in df.columns:
             if not isinstance(col, str):
                 df = df.rename(columns={col: str(col)})

        if xname is None:
            xname = df.columns[0]
        if yname is None:
            yname = df.columns[1]

        # Names must be strings
        xname = str(xname)
        yname = str(yname)

        self.df = df
        self.da = da
        self.xname = xname
        self.yname = yname
        self.title = title
        self.webgl = webgl

        if allowed_ws_origin is not None:
            if isinstance(allowed_ws_origin, str):
                allowed_ws_origin = [allowed_ws_origin]
            settings.allowed_ws_origin.set_value(allowed_ws_origin)
        elif 'VSCODE_PID' in os.environ:
            warning(
                "We appear to be running in VS Code and allowed_ws_origin "
                "was not provided, so setting allowed_ws_origin='*'"
            )
            settings.allowed_ws_origin.set_value('*')

        show(self.bkapp) 

    def _selected_on_change(self, attr, old, new):
        self.selected_frame_indices = new


    def bkapp(self, doc):
        source = ColumnDataSource(data=self.df)
        plot = figure(
            tools='lasso_select',  # type: ignore
            title=self.title,
            output_backend='webgl' if self.webgl else 'canvas',
        )
        plot.scatter(self.xname, self.yname, source=source, selection_color='red')

        source.selected.on_change('indices', self._selected_on_change)

        doc.add_root(column(plot))
        
    
    @property
    def df_selection(self):
        return self.df.iloc[self.selected_frame_indices, :]

    @property
    def selection(self):
        if self.da is None:
            return None
        else:
            return self.da[self.selected_frame_indices, :]


class TrajSelector(FrameSelector):
    def _selected_on_change(self, attr, old, new):
        self.selected_frame_indices = new

    def bkapp(self, doc):
        source = ColumnDataSource(data=self.df[[self.xname, self.yname]])
        assert 'trajid_time' in source.data
        plot = figure(
            tools='lasso_select,tap',  # type: ignore
            title=self.title,
            output_backend='webgl' if self.webgl else 'canvas',
        )
        source.data['_color'] = np.zeros(len(self.df))
        cmap = [
            (85, 85, 85, 1),
            (68, 119, 17, 1),
            (102, 204, 136, 1),
            (85, 85, 85, 0.02)
        ]
        scatter = plot.scatter(
            self.xname,
            self.yname,
            source=source,
            color=linear_cmap('_color', cmap, low=0, high=3),
        )

        div = Div(width=plot.width, height=10, height_policy="fixed")

        js_callback = CustomJS(
            args=dict(source=source, div=div),
            code="""
            //debugger;
            let trajids = source.selected.indices.map((i) => source.data.trajid_time[i][0]);        
            const set_unique_trajids = new Set(trajids);
            const unique_trajids = [...set_unique_trajids];

            div.text = "<span><b>trajids:</b> " + unique_trajids.join(", ") + "</span>";

            source.data['_color'] = new Uint8Array(source.data['0'].length).fill(3);

            for (let i = 0; i < source.data['0'].length; i++) {
                if (set_unique_trajids.has(source.data.trajid_time[i][0])) {
                    source.data['_color'][i] = 1;
                };
            };
            for (let i of source.selected.indices) {
                source.data['_color'][i] = 2;
            };
            source.change.emit();
        """
        )

        source.selected.js_on_change('indices', js_callback)
        source.selected.on_change('indices', self._selected_on_change)

        doc.add_root(column(plot, div))