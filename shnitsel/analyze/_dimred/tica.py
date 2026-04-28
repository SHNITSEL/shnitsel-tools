def tica(da, lagtime):
    """Time-lagged independent component analysis (TICA).
    Work in progress"""
    from deeptime.util.data import TrajectoriesDataset
    from deeptime.decomposition import TICA

    scaled = scale(da).transpose('frame', ...)

    # TODO: What effect does the lagtime on the TrajectoryDataset have?
    deeptime_trajs = TrajectoriesDataset.from_numpy(
        lagtime=lagtime, data=[x for _, x in scaled.groupby('atrajectory')]
    )

    tica_reducer = TICA(lagtime=lagtime, dim=2)
    tica_reducer.fit_from_timeseries(deeptime_trajs)

    return xr.apply_ufunc(
        tica_reducer.transform,
        scaled,
        input_core_dims=[['descriptor']],
        output_core_dims=[['component']],
    ).assign_coords(component=['comp1', 'comp2'])


TICA = tica

time_lagged_independent_component_analysis = tica
