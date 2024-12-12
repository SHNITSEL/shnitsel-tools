import logging, os, re, itertools
import numpy as np
import numpy.typing as npt
#import scipy.stats as st
from shnitsel.dynamic import trajparse, icondparse
from shnitsel.dynamic import postprocess
#from dboverview import plotting
##from dboverview import parsecommon
#import trajparse, icondparse, postprocess, plotting
from dataclasses import dataclass
from typing import Optional, List

NDArrayFloat = npt.NDArray[np.float64]

@dataclass
class Traj():
    nsteps: int
    natoms: int
    nsinglets: int
    ndoublets: int
    ntriplets: int
    nstates: int
    energies: NDArrayFloat
    dip_all: NDArrayFloat
    dip_perm: NDArrayFloat
    dip_trans: NDArrayFloat
    forces: NDArrayFloat
    has_forces: npt.NDArray[bool]
    phases: NDArrayFloat
    nacs: NDArrayFloat
    
    atNames: Optional[List[str]] = None
    atXYZ: Optional[NDArrayFloat] = None
    statelabels = []

    def phase_correction(self, ):
        ...

    def write(self, filename):
        ...

    def __getitem__(_, ts):
        if ts >= _.nsteps:
            raise IndexError(
                f"timestep {ts} is out of bounds; last step is {_.nsteps - 1}")
        if ts < 0:
            raise IndexError(f"timestep {ts} is out of bounds; first step is 0")

        return Traj(nsteps=1,
                    natoms=_.natoms,
                    nsinglets=_.nsinglets,
                    ndoublets=_.ndoublets,
                    ntriplets=_.ntriplets,
                    nstates=_.nstates,
                    energies=_.energies[ts],
                    dip_all=_.dip_all[ts],
                    dip_perm=_.dip_perm[ts],
                    dip_trans=_.dip_trans[ts],
                    forces=_.forces[ts],
                    has_forces=_.has_forces[ts],
                    phases=_.phases[ts],
                    nacs=_.nacs[ts]
                    )


def read():
    ...
    #return SHNITSELDB(data)
