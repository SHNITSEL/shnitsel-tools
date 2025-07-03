import numpy as np
import pandas as pd
import xarray as xr
import logging
import os
import re
import math
from itertools import product, combinations
from glob import glob
from typing import NamedTuple, Any
from tqdm.auto import tqdm
from .common import (
    get_dipoles_per_xyz,
    dip_sep,
    __atnum2symbol__,
    ConsistentValue,
    get_triangular,
)
from ..postprocess import convert_length

_re_grads = re.compile('[(](?P<nstates>[0-9]+)x(?P<natoms>[0-9]+)x3')
_re_nacs = re.compile('[(](?P<nstates>[0-9]+)x[0-9]+x(?P<natoms>[0-9]+)x3')


class IcondPath(NamedTuple):
    idx: int
    path: str
    # prefix: str | None


def nans(*dims):
    return np.full(dims, np.nan)


def list_iconds(iconds_path='./iconds/', glob_expr='**/ICOND_*'):
    dirs = glob(glob_expr, recursive=True, root_dir=iconds_path)
    if len(dirs) == 0:
        raise FileNotFoundError(
            f"The search '{glob_expr}' didn't match any directories "
            f"under {iconds_path=} "
            f"relative to working directory '{os.getcwd()}'"
        )
    names = sorted(
        [
            name
            for name in dirs
            if 'QM.out' in os.listdir(os.path.join(iconds_path, name))
        ]
    )

    return [IcondPath(int(name[6:]), os.path.join(iconds_path, name)) for name in names]


def dims_from_QM_out(f):
    # faced with redundancy, use it to ensure consistency
    nstates = ConsistentValue('nstates', weak=True)
    natoms = ConsistentValue('natoms', weak=True)

    for index, line in enumerate(f):
        if line.startswith("! 1 Hamiltonian Matrix"):
            nstates.v = int(next(f).split(' ')[0])
        elif line.startswith('! 2 Dipole Moment Matrices'):
            dim = re.split(' +', next(f).strip())
            nstates.v = int(dim[0])
        elif line.startswith('! 3 Gradient Vectors'):
            info = _re_grads.search(line)
            assert info is not None
            nstates.v, natoms.v = map(int, info.group('nstates', 'natoms'))
        elif line.startswith('! 5 Non-adiabatic couplings'):
            info = _re_nacs.search(line)
            assert info is not None
            nstates.v, natoms.v = map(int, info.group('nstates', 'natoms'))

    return nstates.v, natoms.v


def dims_from_QM_log(log):
    nstates = ConsistentValue('nstates', weak=True)
    natoms = ConsistentValue('natoms', weak=True)
    for line in log:
        if line.startswith('States:'):
            linecont = line.strip().split()
            if 'Singlet' in linecont and 'Triplet' not in linecont:
                nsinglets = int(linecont[2])
                ntriplets = 0
            elif 'Singlet' in linecont and 'Triplet' in linecont:
                nsinglets = int(linecont[1])
                ntriplets = int(linecont[3])
            elif 'Triplet' in linecont and 'Singlet' not in linecont:
                ntriplets = int(linecont[2])
                nsinglets = 0

            # calculate total number of states
            nstates.v = nsinglets + (3 * ntriplets)

        elif line.startswith('Found Geo!'):
            linecont = re.split(' ', line.strip())
            natoms.v = int(linecont[-1][0:-1])

    return nstates.v, natoms.v


def check_dims(pathlist):
    if len(pathlist) == 0:
        raise ValueError("pathlist is empty")
    nstates = ConsistentValue('nstates', ignore_none=True)
    natoms = ConsistentValue('natoms', ignore_none=True)
    for _, path in pathlist:
        try:
            with open(os.path.join(path, 'QM.out')) as f:
                nstates.v, natoms.v = dims_from_QM_out(f)
            with open(os.path.join(path, 'QM.log')) as f:
                nstates.v, natoms.v = dims_from_QM_log(f)
        except FileNotFoundError:
            pass
    return nstates.v, natoms.v


def dir_of_iconds(path='./iconds/', *, subset: (set | None) = None):
    pathlist = list_iconds(path)
    if subset is not None:
        pathlist = [icond for icond in pathlist if icond.name in subset]

    return read_iconds(pathlist)


def dirs_of_iconds(path='./iconds/', *, levels=1, subset: (set | None) = None):
    pathlist = list_iconds(path)
    if subset is not None:
        pathlist = [icond for icond in pathlist if icond.name in subset]

    return read_iconds(pathlist)


def init_iconds(indices, nstates, natoms, **res):
    template = {
        'energy': ['state'],
        'dip_all': ['state', 'state2', 'direction'],
        'dip_perm': ['state', 'direction'],
        'dip_trans': ['statecomb', 'direction'],
        'forces': ['state', 'atom', 'direction'],
        # 'has_forces': ['placeholder'],
        # 'has_forces': [],
        'phases': ['state'],
        'nacs': ['statecomb', 'atom', 'direction'],
        'atXYZ': ['atom', 'direction'],
    }

    if natoms == 0:
        # This probably means that check_dims() couldn't find natoms,
        # so we don't expect properties with an atom dimension.

        del template['forces']
        del template['nacs']

        # On the other hand, we don't worry about not knowing nstates,
        # because energy is always written.

    if indices is not None:
        niconds = len(indices)
        for varname, dims in template.items():
            dims.insert(0, 'icond')
    else:
        niconds = 0

    template['atNames'] = ['atom']

    lens = {
        # 'placeholder': 1,
        'icond': niconds,
        'state': nstates,
        'state2': nstates,
        'atom': natoms,
        'direction': 3,
        'statecomb': math.comb(nstates, 2),
    }

    coords: dict | xr.Dataset = {
        'state': (states := np.arange(nstates)),
        'state2': states,
        'atom': np.arange(natoms),
        'direction': ['x', 'y', 'z'],
        'has_forces': (
            'icond',
            x if (x := res.get('has_forces')) is not None else nans(niconds),
        ),
    }

    if indices is not None:
        coords['icond'] = indices

    coords = xr.Coordinates.from_pandas_multiindex(
        pd.MultiIndex.from_tuples(combinations(states, 2), names=['from', 'to']),
        dim='statecomb',
    ).merge(coords)

    attrs = {
        'energy': {'units': 'hartree', 'unitdim': 'Energy'},
        'e_kin': {'units': 'hartree', 'unitdim': 'Energy'},
        'dip_perm': {'long_name': "permanent dipoles", 'units': 'au'},
        'dip_trans': {'long_name': "transition dipoles", 'units': 'au'},
        'sdiag': {'long_name': 'active state (diag)'},
        'astate': {'long_name': 'active state (MCH)'},
        'forces': {'units': 'hartree/bohr', 'unitdim': 'Force'},
        'nacs': {'long_name': "nonadiabatic couplings", 'units': 'au'},
    }

    datavars = {
        varname: (
            dims,
            (
                x
                if (x := res.get(varname)) is not None
                else nans(*[lens[d] for d in dims])
            ),
            attrs[varname] if varname in attrs else {},
        )
        for varname, dims in template.items()
    }

    datavars['atNames'] = (['atom'], np.full((natoms), ''), {})

    return xr.Dataset(datavars, coords)


def read_iconds(pathlist, index=None):
    logging.info("Ensuring consistency of ICONDs dimensions")
    nstates, natoms = check_dims(pathlist)
    logging.info("Allocating Dataset for ICONDs")
    if index is None:
        index = [p.idx for p in pathlist]
    iconds = init_iconds(index, nstates, natoms)
    logging.info("Reading ICONDs data into Dataset...")

    for icond, path in tqdm(pathlist):
        with open(os.path.join(path, 'QM.out')) as f:
            parse_QM_out(f, out=iconds.sel(icond=icond))

    for icond, path in tqdm(pathlist):
        try:
            with open(os.path.join(path, 'QM.log')) as f:
                parse_QM_log_geom(f, out=iconds.sel(icond=icond))
        except FileNotFoundError:
            logging.warning(
                f"""no QM.log file found in {path}.
                This is currently used to determine geometry.
                Eventually, user-inputs will be accepted as an alternative.
                See https://github.com/SHNITSEL/db-workflow/issues/3"""
            )

    iconds.atXYZ.attrs['units'] = 'bohr'
    iconds['atXYZ'] = convert_length(iconds.atXYZ, to='angstrom')
    return iconds


def parse_QM_log(log):
    info: dict[str, Any] = {}
    for line in log:
        if line.startswith('States:'):
            linecont = re.split(' +|\t', line.strip())
            if 'Singlet' in linecont and 'Triplet' not in linecont:
                nsinglets = int(linecont[2])
                ntriplets = 0
            elif 'Singlet' in linecont and 'Triplet' in linecont:
                nsinglets = int(linecont[2])
                ntriplets = int(linecont[5])
            elif 'Triplet' in linecont and 'Singlet' not in linecont:
                ntriplets = int(linecont[2])
                nsinglets = 0

            # calculate total number of states
            nstates = nsinglets + (3 * ntriplets)

            info['nStates'] = nstates
            info['nSinglets'] = nsinglets
            info['nTriplets'] = ntriplets
            nnacs = int(nsinglets * (nsinglets - 1) / 2) + int(
                ntriplets * (ntriplets - 1) / 2
            )
            info['nNACS'] = nnacs
            info['nDipoles'] = int(nsinglets + ntriplets + nnacs)

        elif line.startswith('Method:'):
            linecont = re.split(' +|\t', line.strip())
            method = linecont[2]

            info['method'] = method

        elif line.startswith('Found Geo!'):
            linecont = re.split(' ', line.strip())
            natom = int(linecont[-1][0:-1])

            info['nAtoms'] = natom

        elif line.startswith('Geometry in Bohrs:'):
            # NB. Geometry is indeed in bohrs!
            atnames = []
            atxyz = np.zeros((natom, 3))
            for i in range(natom):
                geometry_line = re.split(' +', next(log).strip())
                atnames.append(geometry_line[0])
                atxyz[i] = [geometry_line[j] for j in range(1, 4)]

            info['atNames'] = atnames
            r = {j: i for i, j in __atnum2symbol__.items()}
            info['atNums'] = [r[i] for i in atnames]
            info['atXYZ'] = atxyz

    return info


def parse_QM_log_geom(f, out):
    # NB. Geometry is indeed in bohrs!
    while not next(f).startswith('Geometry in Bohrs:'):
        pass

    for i in range(out.sizes['atom']):
        geometry_line = next(f).strip().split()
        out['atNames'][i] = geometry_line[0]
        out['atXYZ'][i] = geometry_line[1:4]


def parse_QM_out(f, out: (xr.Dataset | None) = None):
    res: xr.Dataset | dict[str, np.ndarray]
    if out is not None:
        # write data directly into dataset
        res = out
    else:
        # write data as ndarrays into dict, then make dataset after parsing
        res = {}

    res['has_forces'] = np.array([0])
    nstates = ConsistentValue('nstates')
    natoms = ConsistentValue('natoms')

    for index, line in enumerate(f):
        if line.startswith("! 1 Hamiltonian Matrix"):
            # get number of states from dimensions of Hamiltonian
            nstates.v = int(next(f).split(' ')[0])
            if out is None:
                res['energy'] = nans(nstates.v)

            for istate in range(nstates.v):
                energyline = re.split(' +', next(f).strip())
                res['energy'][istate] = float(energyline[2 * istate])

        elif line.startswith('! 2 Dipole Moment Matrices'):
            dim = re.split(' +', next(f).strip())
            n = int(dim[0])
            m = int(dim[1])

            if out is None:
                res['dip_all'] = nans(n, m, 3)
                res['dip_perm'] = nans(n, 3)
                res['dip_trans'] = nans(math.comb(n, 2), 3)

            res['dip_all'][:, :, 0] = get_dipoles_per_xyz(f, n, m)
            next(f)
            res['dip_all'][:, :, 1] = get_dipoles_per_xyz(f, n, m)
            next(f)
            res['dip_all'][:, :, 2] = get_dipoles_per_xyz(f, n, m)

            res['dip_perm'][:], res['dip_trans'][:] = dip_sep(np.array(res['dip_all']))

        elif line.startswith('! 3 Gradient Vectors'):
            res['has_forces'] = np.array([1])

            search_res = _re_grads.search(line)
            assert search_res is not None
            get_dim = search_res.group
            nstates.v = int(get_dim('nstates'))
            natoms.v = int(get_dim('natoms'))

            if out is None:
                res['forces'] = nans(nstates.v, natoms.v, 3)

            for istate in range(nstates.v):
                next(f)
                for atom in range(natoms.v):
                    res['forces'][istate][atom] = [
                        float(entry) for entry in next(f).strip().split()
                    ]

        elif line.startswith('! 5 Non-adiabatic couplings'):
            search_res = _re_nacs.search(line)
            assert search_res is not None
            get_dim = search_res.group
            nstates.v = int(get_dim('nstates'))
            natoms.v = int(get_dim('natoms'))

            if out is None:
                res['nacs'] = nans(math.comb(nstates.v, 2), natoms.v, 3)

            nacs_all = nans(nstates.v, nstates.v, natoms.v, 3)

            for bra, ket in product(range(nstates.v), range(nstates.v)):
                # TODO info currently unused, but keep the `next(f)` no matter what!
                nac_multi = int(re.split(' +', next(f).strip())[-1])  # noqa: F841

                for atom in range(natoms.v):
                    nacs_line = re.split(' +', next(f).strip())
                    nacs_all[bra, ket, atom] = [float(n) for n in nacs_line]

            # all nacs, i.e., nacs of all singlet and triplet states
            # all diagonal elements are zero (self-coupling, e.g. S1 and S1)
            # off-diagonal elements conatin couplings of different states (e.g. S0 and S1)
            # in principle one has here the full matrix for the nacs between all singlet and triplet states
            # in the following we extract only the upper triangular elements of the matrix

            res['nacs'][:] = get_triangular(nacs_all)

        elif line.startswith('! 6 Overlap matrix'):
            nlines = int(re.split(' +', next(f).strip())[0])
            assert nlines == nstates.v

            found_overlap = False
            phasevector = np.ones((nlines))

            wvoverlap = np.zeros((nlines, nlines))
            for j in range(nlines):
                linecont = [float(n) for n in re.split(' +', next(f).strip())]
                vec = [n for n in linecont[::2]]
                assert len(vec) == nlines
                wvoverlap[j] = vec

            for istate in range(nlines):
                if np.abs(wvoverlap[istate, istate]) >= 0.5:
                    found_overlap = True
                    if wvoverlap[istate, istate] >= 0.5:
                        res['phases'][istate] = +1
                    else:
                        res['phases'][istate] = -1

            if found_overlap:
                res['phases'][:] = phasevector
                pass

        elif line.startswith('! 8 Runtime'):
            next(f)

    if out is None:
        if not res['has_forces']:
            res['forces'] = nans(natoms.v, 3)

        assert isinstance(res, dict)
        return init_iconds(indices=None, nstates=nstates.v, natoms=natoms.v, **res)

        # xr.Dataset(
        #     {
        #         'energy': (['state'], energy),
        #         'dip_all': (['state', 'state2', 'direction'], dip_all),
        #         'dip_perm': (['state', 'direction'], dip_perm),
        #         'dip_trans': (['statecomb', 'direction'], dip_trans),
        #         'forces': (['atom', 'direction'], forces),
        #         'has_forces': ([], has_forces),
        #         'phases': (['state'], phases),
        #         'nacs': (['statecomb', 'atom', 'direction'], nacs)
        #     },
        #     coords={
        #         'state': np.arange(1, nstates.v+1),
        #         'state2': np.arange(1, nstates.v+1),
        #         'atom': np.arange(natoms.v),
        #         'statecomb': np.arange(math.comb(nstates.v, 2)),
        #         'direction': ['x', 'y', 'z']
        #     }
        # )
    else:
        # all the data has already been written to `out`
        # no need to return anything
        return None

def iconds_to_frames(iconds: xr.Dataset):
    for name, var in iconds.data_vars.items():
        shape = var.data.shape
        if 0 in shape:
            raise ValueError(
                f"Variable '{name}' has shape {shape} which contains 0. "
                "Please remove this variable before converting to frames. "
                "Note: An empty variable could indicate a problem with parsing."
            )

    if 'atNames' in iconds.data_vars and 'atNames' not in iconds.coords:
        iconds = iconds.assign_coords(atNames=iconds.atNames)

    return (
        iconds.rename_dims(icond='trajid')
        .rename_vars(icond='trajid')
        .expand_dims('time')
        .assign_coords(time=('time', [0.0]))
        .stack(frame=['trajid', 'time'])
    )