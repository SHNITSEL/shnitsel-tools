import numpy as np
import re
import logging

__atnum2symbol__ = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 
                    6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 
                    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',                              
                    16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
                    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn',
                    26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
                    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 
                    36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
                    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh',
                    46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
                    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs',
                    56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
                    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb',
                    66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
                    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re',
                    76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
                    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
                    86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 
                    91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 
                    96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 
                    101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db',
                    106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
                    111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc',
                    116: 'Lv', 117: 'Ts', 118: 'Og' }

def get_triangular(original_array):
    """
    get_triangular - get the upper triangle of a (nstat1 x nstat2 x natoms x 3) matrix

    This function takes in a 4-dimensional numpy array (original_array) and returns a 3-dimensional numpy array (upper_tril)
    which is the upper triangle of the input matrix, obtained by excluding the diagonal elements.
    The number of steps (k) to move the diagonal above the leading diagonal is 1.
    The returned matrix has shape (len(cols), natoms, 3)

    Inputs:
    original_array: 4D numpy array of shape (nstat1, nstat2, natoms, 3) representing the input matrix

    Returns:
    upper_tril: 3D numpy array of shape (len(cols), natoms, 3) representing the upper triangle of the input matrix
    """
    # Get the indices of the upper triangle
    nstat1, nstat2, natoms, xyz = original_array.shape

    if nstat1 != nstat2:
        raise ValueError("expected square input matrix")

    rows, cols = np.triu_indices(nstat2, k=1)
    upper_tril = np.zeros((len(cols), natoms, 3))

    for i in range(len(cols)):
        me = original_array[rows[i], cols[i]]
        upper_tril[i] = me

    return upper_tril

def get_dipoles_per_xyz(file, n, m):

    dip = np.zeros((n, m))
    for istate in range(n):
        linecont = [float(i) for i in re.split(' +', next(file).strip())]
        tmp = []
        for element in linecont:
            if element != 0.0:
                tmp.append(element)
        if len(tmp) == m:
            dip[istate] = tmp

    return dip

def dip_sep(dipoles):
    """
    Separates a complete matrix of dipoles into permanent
    and transitional dipoles, removing redundancy in the process.

    Inputs:
    dipoles: 3D numpy array of shape (nstates, nstates, 3) where
        the first axis represents state before transition,
        the second axis represents state after transition and
        the third axis contains x, y and z coordinates.

    Outputs:
    dip_perm: 2D numpy array of shape (nstates, 3)
    dip_trans: 2D numpy array of shape (math.comb(nstates, 2), 3)
        in the order e.g. (for nstates = 4)
        0->1, 0->2, 0->3, 1->2, 1->3, 2->3
        where 0->1 is the transitional dipole between
        state 0 and state 1.
    """
    assert(dipoles.ndim == 3)
    nstates, check, three = dipoles.shape
    assert(nstates == check)
    assert(three == 3)
    dip_perm = np.diagonal(dipoles).T
    dip_trans = dipoles[np.triu_indices(nstates, k=1)]
    logging.debug("permanent dipoles\n" + str(dip_perm))
    logging.debug("transitional dipoles\n" +  str(dip_trans))
    return dip_perm, dip_trans

class ConsistentValue():
    defined = False
    _val = None
    name = "ConsistentValue"

    def __init__(self, name, weak=False, ignore_none=False):
        self.name = name
        self._weak = weak
        self._ignore_none = ignore_none

    @property
    def v(self):
        if self.defined:
            return self._val
        elif self._weak:
            return None
        raise AttributeError(f"{self.name}.v accessed before assignment")

    @v.setter
    def v(self, new_val):
        if self._ignore_none and new_val is None:
            return

        if self.defined and new_val != self._val:
            raise ValueError(
f"""inconsistent assignment to {self.name}:
    current value: {type(self._val).__name__} = {repr(self._val)}
    new value:  {type(new_val).__name__} = {repr(new_val)}
"""
            )

        self.defined = True
        self._val = new_val