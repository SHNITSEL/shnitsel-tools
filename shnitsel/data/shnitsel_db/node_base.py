import xarray as xr


class ShnitselNode(xr.DataTree):
    def __init__(self, dataset, name, level):
        # TODO FIXME Did this break in a merge somehow?
        self.super().__init__(
            dataset,
        )
