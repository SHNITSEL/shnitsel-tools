import numpy as np

def parse_xyz(f):
    natoms = int(next(f).strip())
    # ts = 0

    atXYZ = [] #np.full((nsteps, natoms, 3), np.nan)
    atNames = [] #np.full((natoms), '')
    next(f)
    thisXYZ = np.full((natoms, 3), np.nan)
    for iatom in range(natoms):
        geometry_line = next(f).strip().split()
        # atNames[iatom] = geometry_line[0]
        # atXYZ[ts, iatom] = [float(n) for n in geometry_line[1:]]
        atNames.append(geometry_line[0])
        thisXYZ[iatom] = [float(n) for n in geometry_line[1:]]
    atXYZ.append(thisXYZ)

    for line in f:
        assert line.startswith(' '), f'line content: {line!r}'
        # ts += 1
        line = next(f)
        assert line.startswith(' '), f'line content: {line!r}'

        thisXYZ = np.full((natoms, 3), np.nan)
        for iatom, atName in enumerate(atNames):
            geometry_line = next(f).strip().split()
            assert geometry_line[0] == atName,\
                "Inconsistent atom order"
            # atXYZ[ts, iatom] = [float(n) for n in geometry_line[1:]]
            thisXYZ[iatom] = [float(n) for n in geometry_line[1:]]
        atXYZ.append(thisXYZ)

    return (atNames, np.stack(atXYZ, axis=0))