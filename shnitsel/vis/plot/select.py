import os
import logging

try:
    from bokeh.io import show, output_notebook
    from bokeh.layouts import column
    from bokeh.models import ColumnDataSource, CustomJS, Div
    from bokeh.plotting import figure
    from bokeh.settings import settings
    from bokeh.transform import linear_cmap

    _bokeh_package_available = True
except ImportError:
    _bokeh_package_available = False
    logging.warning("Bokeh package missing. Plotting capabilities limited")

import numpy as np
import pandas as pd
import xarray as xr

class FrameSelector:
    selected_frame_indices: list[int] = []

    def __init__(
        self,
        data,
        data_var=None,
        *,
        dim=None,
        xname=None,
        yname=None,
        title="",
        allowed_ws_origin=None,
        webgl=True,
    ):
        if not _bokeh_package_available:
            logging.error(f"ERROR: Package <bokeh> was not found in the current environment. Please install the missing package to use the {self.__class__.__name__} class.")

        output_notebook()

        if isinstance(data, pd.DataFrame):
            df = data
            xr_obj = None
            if dim is not None:
                raise TypeError("`dim` should not be set when `data` is a DataFrame")

        elif isinstance(data, xr.DataArray):
            xr_obj = data
            if dim is None:
                dim = xr_obj.dims[0]

            if len(xr_obj.dims) != 2:
                raise ValueError(
                    "When the first argument to FrameSelector is an "
                    "xarray.DataArray, it should have 2 dimensions, "
                    f"rather than {len(xr_obj.dims)} dimensions (namely {xr_obj.dims})."
                )
            df = xr_obj.transpose(dim, ...).to_pandas()

        elif isinstance(data, xr.Dataset):
            if data_var is None:
                raise TypeError(
                    "If 'data' is an xarray.Dataset, please indicate which data_var to use"
                )
            xr_obj = data
            if dim is None:
                dim = xr_obj[data_var].dims[0]
            df = xr_obj[data_var].to_pandas()

        else:
            raise TypeError(
                "The first argument to FrameSelector should be a "
                "pandas.DataFrame, a 2-dimensional xarray.DataArray "
                "or an xarray.Dataset containing such a DataArray."
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
        self.xr_obj = xr_obj
        self.dim = dim

        self.xname = xname
        self.yname = yname
        self.title = title
        self.webgl = webgl

        if allowed_ws_origin is not None:
            if isinstance(allowed_ws_origin, str):
                allowed_ws_origin = [allowed_ws_origin]
            settings.allowed_ws_origin.set_value(allowed_ws_origin)
        elif 'VSCODE_PID' in os.environ:
            logging.warning(
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
        if self.xr_obj is None:
            return None
        else:
            return self.xr_obj[{self.dim: self.selected_frame_indices}]


class TrajSelector(FrameSelector):
    def _selected_on_change(self, attr, old, new):
        self.selected_frame_indices = new

    def bkapp(self, doc):
        source = ColumnDataSource(data=self.df[[self.xname, self.yname]])
        assert 'atrajectory_time' in source.data
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
            (85, 85, 85, 0.1)
        ]
        plot.scatter(
            self.xname,
            self.yname,
            source=source,
            color=linear_cmap('_color', cmap, low=0, high=3),
            nonselection_alpha=1,
        )

        div = Div(width=plot.width, height=10, height_policy="fixed")

        js_callback = CustomJS(
            args=dict(source=source, div=div),
            code="""
            //debugger;
            let trajids = source.selected.indices.map((i) => source.data.atrajectory_time[i][0]);        
            const set_unique_trajids = new Set(trajids);
            const unique_trajids = [...set_unique_trajids];

            div.text = "<span><b>Trajectory IDs:</b> " + unique_trajids.join(", ") + "</span>";

            source.data['_color'] = new Uint8Array(source.data['0'].length).fill(3);

            for (let i = 0; i < source.data['0'].length; i++) {
                if (set_unique_trajids.has(source.data.atrajectory_time[i][0])) {
                    source.data['_color'][i] = 1;
                };
            };
            for (let i of source.selected.indices) {
                source.data['_color'][i] = 2;
            };
            source.change.emit();
        """,
        )

        source.selected.js_on_change('indices', js_callback)
        source.selected.on_change('indices', self._selected_on_change)

        doc.add_root(column(plot, div))