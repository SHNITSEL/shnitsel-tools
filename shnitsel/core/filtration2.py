from numbers import Number

from shnitsel.core.midx import sel_trajs


def is_dataset_stacked(obj):
    if 'frame' in obj.dims and {'trajid', 'time'}.issubset(obj.coords):
        return True
    elif {'trajid', 'time'}.issubset(mask.dims):
        return False
    else:
        raise ValueError(
            "The mask argument should be trajectories, either stacked or unstacked"
        )


def mask_from_filtranda(): ...


def mask_from_cutoffs(): ...


def mask_from_dataset():
    # done, called same name

    ########################
    # Formats we will need #
    ########################
    # two versions currently, for stacked and for unstacked
    # use `filtranda` and `thresholds` if available
    # otherwise use `good_upto`
    ...


def cutoffs_from_mask(): ...


def cutoffs_from_filtranda(): ...


def cutoffs_from_dataset(): ...


####################
# Action functions #
####################

# All the action functions take a dataset
# They can use the functions above to get the info they need


def omit(ds):
    cutoffs = cutoffs_from_dataset(ds)
    good_throughout = cutoffs['good_throughout']
    selection = good_throughout.all('criterion')
    return sel_trajs(ds, selection)


def truncate(ds): ...
def transect(ds, cutoff): ...


#
# Filtranda derivat
#


#########################
# Convenience functions #
#########################


def energy_filtranda(): ...


def sanity_check(
    frames,
    cut=False,
    *,
    units='eV',  # TODO: FIXME: Actually implement!
    etot_drift=0.2,
    etot_step=0.1,
    epot_step=0.7,
    ekin_step=0.7,
    hop_epot_step=1.0,
):
    settings = {
        k: locals()[k]
        for k in ['etot_drift', 'etot_step', 'epot_step', 'ekin_step', 'hop_epot_step']
    }
    frames = frames.assign(filtranda=energy_filtranda(frames, **settings))
    frames = get_cutoffs(frames)
    if not cut:
        return frames
    elif cut == 'truncate':
        return truncate(frames)
    elif cut == 'omit':
        return omit(frames)
    else:
        assert isinstance(cut, Number)
        return transect(frames, cut)


def filter_cleavages(
    atom_pair_identifiers: list,
    atom_pair_thresholds: dict,
): ...


###########################################
# Formats directly prerequisite for plots #
###########################################


# For check_thresholds aka plot_thresholds
def cum_max_quantiles(ds_or_da, quantiles=None):
    ...
    # Current implementation expects filtranda
    # or dataset containing filtranda


# For validity_populations
def validity_populations(ds, intersections=False): ...