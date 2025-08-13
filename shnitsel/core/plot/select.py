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
                self.selection = da[[], :]
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
        
        self.selected_frame_indices = []
        self.df_selection = pd.DataFrame

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

        bkapp = self._bkapp()
        show(bkapp)
    
    def _bkapp(self):
        source = ColumnDataSource(data=self.df)

        def callback(attr, old, new):
            nonlocal self
            self.selected_frame_indices = new
            self.df_selection = self.df.iloc[new, :]
            if self.da is not None:
                self.selection = self.da[new, :]

        def bkapp(doc):
            nonlocal self, source, callback
            plot = figure(
                tools='lasso_select',  # type: ignore
                title=self.title,
                output_backend='webgl' if self.webgl else 'canvas',
            )
            plot.scatter(self.xname, self.yname, source=source, selection_color='red')

            source.selected.on_change('indices', callback)

            doc.add_root(column(plot))
        
        return bkapp

class TrajSelector(FrameSelector):
    def _bkapp(self):
        source = ColumnDataSource(data=self.df)

        def callback(attr, old, new):
            nonlocal self
            self.selected_frame_indices = new
            self.df_selection = self.df.iloc[new, :]
            if self.da is not None:
                self.selection = self.da[new, :]

        def bkapp(doc):
            nonlocal self
            plot = figure(
                tools='lasso_select,tap',  # type: ignore
                title=self.title,
                output_backend='webgl' if self.webgl else 'canvas',
            )
            source2 = ColumnDataSource({k: [] for k in source.column_names})
            scatter2 = plot.scatter(self.xname, self.yname, source=source2, color='limegreen', nonselection_alpha=1)
            scatter = plot.scatter(self.xname, self.yname, source=source, selection_color='red')

            div = Div(width=plot.width, height=10, height_policy="fixed")

            js_callback = CustomJS(args=dict(source=source, plot=plot, source2=source2, div=div), code="""
                let trajids = source.selected.indices.map((i) => source.data.trajid_time[i][0]);        
                const unique_trajids = [];


                trajids.forEach((x) => {
                    if (!unique_trajids.includes(x)) {
                        unique_trajids.push(x);
                    };
                });

                div.text = "<span><b>trajids:</b> " + unique_trajids + "</span>";

                let new_indices = [];
                for (let i = 0; i < source.data['0'].length; i++) {
                    if (unique_trajids.includes(source.data.trajid_time[i][0])) {
                        new_indices.push(i);
                    }
                };
                const new_data = {};
                for (const col in source.data) {
                    new_data[col] = new_indices.map(i => source.data[col][i]);
                };
                source2.data = new_data;
                source.change.emit();
                source2.change.emit();
            """
            )

            source.selected.js_on_change('indices', js_callback)
            source.selected.on_change('indices', callback)

            doc.add_root(column(plot, div))
        
        return bkapp