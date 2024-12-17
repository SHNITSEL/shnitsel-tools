from ... import postprocess as P
from .spectra import get_spectrum, calc_spectra, get_sgroups

__all__ = ['get_spectrum', 'calc_spectra', 'get_sgroups']

# Calculate all at once
# TODO: Consider making a dataclass with all the different
# Datasets and methods to process and plot them.


def calc_all(frames):
    """
    Calculate all properties required for a datasheet and return them in a dictionary.
    Note: this function may yet be refactored as a method of a Datasheet class.
    """
    per_state = P.get_per_state(frames)
    inter_state = P.get_inter_state(frames)
    pops = P.calc_pops(frames)

    delta_E = P.time_grouped_ci(inter_state['energies'])
    noodle, hops = P.pca_and_hops(frames)

    if 'dip_trans' in frames:
        inter_state = P.assign_fosc(inter_state)
        fosc_time = P.time_grouped_ci(inter_state['fosc'])
        spectra = calc_spectra(inter_state)
        sgroups = get_sgroups(spectra)
    else:
        fosc_time = None
        spectra = None
        sgroups = None

    return {
        'per_state': per_state,
        'inter_state': inter_state,
        'pops': pops,
        'sgroups': sgroups,
        'noodle': noodle,
        'hops': hops,
        'delta_E': delta_E,
        'fosc_time': fosc_time,
        'atXYZ': frames['atXYZ'].isel(frame=0),
    }