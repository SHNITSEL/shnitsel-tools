# pyrigth:off
from dataclasses import dataclass
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

import rdkit # type: ignore
from rdkit.Chem import AllChem # type: ignore
import PIL
import io

from . import (
  postprocess as P)

def figax(ax=None):
    """
    Create figure and axes-object if an axes-object is not supplied.
    """
    if ax is None:
        return plt.subplots(1,1)
    else:
        return ax.figure, ax


# NOTE difference: takes ax, original took (axs: dict) such that `'pca' in axs`
# NOTE also moved ax to go after data arguments to facilitate figax approach
def plot_noodleplot(noodle, hops=None, ax=None,
    c=None, colorbar_label=None,
    cmap=None, cnorm=None, cscale=None,
    hops_kws=None):
    fig, ax = figax(ax)

    if c is None:
        c = noodle['time']
        colorbar_label = colorbar_label or '$t$ / fs'

    cmap = cmap or mpl.colormaps['cividis_r']
    if isinstance(cmap, str):
        cmap = mpl.colormaps[cmap]
    cnorm = cnorm \
      or mpl.colors.Normalize(c.min(), c.max()) #type: ignore
    cscale = cscale or mpl.cm.ScalarMappable( #type: ignore
        norm=cnorm,
        cmap=cmap)

    for trajid, traj in noodle.groupby('trajid'):
        ctraj = c.sel(trajid=trajid)
        ax.scatter(
            traj.isel(PC=0), traj.isel(PC=1),
            c=cmap(cnorm(ctraj)),
            label=trajid, alpha=0.5, s=0.2)

    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    if hops is not None:
        hops_kws = dict(s=0.5, c='limegreen') | (hops_kws or {})
        ax.scatter(hops.isel(PC=0), hops.isel(PC=1), **hops_kws)

    # TODO facilitate custom colorbar
    fig.colorbar(cscale, ax=ax, label=colorbar_label, pad=0.02)

    return ax

# TODO: finish later!
def plot_noodleplot_lines(noodle, #hops,
  ax=None,
  cmap=None, cnorm=None, cscale=None):
    fig, ax = figax(ax)
    # One traj per line
    for trajid, traj in noodle.groupby('trajid'):

        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = mpl.collections.LineCollection()

    return ax

def get_loadings(frames):
    atXYZ = P.subtract_combinations(frames['atXYZ'], 'atom', labels=True)
    atXYZ = P.norm(atXYZ)
    _, pca_obj = P.pca(atXYZ, 'atomcomb', return_pca_object=True)
    return xr.DataArray(
      data=pca_obj.components_,
      dims=['PC', 'atomcomb'],
      coords=dict(
        atomcomb=atXYZ.atomcomb
      ),
      attrs={'natoms': frames.sizes['atom']}
    )

def plot_loadings(ax, loadings):
    for _, pcs in loadings.groupby('atomcomb'):
        pc1, pc2 = pcs.values
        ax.arrow(0, 0, pc1, pc2)
        a1, a2 = int(pcs['from']), int(pcs['to'])
        ax.text(pc1, pc2, f"{a1},{a2}")

def xyz_to_mol(atXYZ, charge=0):
    mol = rdkit.Chem.rdmolfiles.MolFromXYZBlock(P.to_xyz(atXYZ))
    rdkit.Chem.rdDetermineBonds.DetermineBonds(mol, charge)
    AllChem.Compute2DCoords(mol)
    return mol

def show_atom_numbers(atXYZ):
    mol = xyz_to_mol(atXYZ)

    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("atomLabel", str(atom.GetIdx()))

    return mol

def cluster_general(decider, n):
    clustered = np.full((n,), False)
    clusters = []
    # for each item, if it has not been clustered,
    # put those later item which have not yet been clustered
    # in a cluster with it
    for i in range(n):
        if clustered[i]:
            continue

        cluster = [i]
        for j in range(i+1, n):
            if clustered[j]:
                continue
            if decider(i, j):
                cluster.append(j)
                clustered[j] = True

        clusters.append(cluster)

    return clusters

def cluster_loadings(loadings: xr.DataArray, cutoff=0.05):
    def dist(i, j, l):
        pc1, pc2 = l.isel(atomcomb=j).values - l.isel(atomcomb=i).values
        return (pc1**2 + pc2**2)**0.5

    def decider(i, j):
        nonlocal loadings, cutoff, dist
        return dist(i, j, loadings) <= cutoff

    n = loadings.sizes['atomcomb']
    return cluster_general(decider, n)

def plot_clusters(loadings, clusters, ax=None, labels=None):
    fig, ax = figax(ax)
    for i, cluster in enumerate(clusters):
        acs = loadings.isel(atomcomb=cluster)
        x, y = acs.mean(dim='atomcomb')
        s = (labels[i] if labels is not None else
          ' '.join([f'({a1},{a2})' for a1, a2 in acs.atomcomb.values]))
        ax.arrow(0, 0, x, y)
        ax.text(x, y, s)

def highlight_pairs(mol, pairs):
    d = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo(320, 240)
    # colors = iter(mpl.colormaps['tab10'](range(10)))
    colors = iter(mpl.colormaps['rainbow'](np.linspace(0, 1, len(pairs))))

    acolors = {}
    bonds = {}
    for a1, a2 in pairs:
        if (bond := mol.GetBondBetweenAtoms(a1, a2)) is not None:
            bonds[bond.GetIdx()] = [(1, 0.5, 0.5)]
        else:
            c = tuple(next(colors))
            for a in [a1, a2]:
                if a not in acolors:
                    acolors[a] = []
                acolors[a].append(c)
            
    # d.drawOptions().fillHighlights = False
    d.drawOptions().setBackgroundColour((0.8, 0.8, 0.8, 0.5))
    d.drawOptions().padding = 0

    d.DrawMoleculeWithHighlights(
        mol,
        '',
        acolors,
        bonds,
        {}, {}, -1
    )
    d.FinishDrawing()
    return d.GetDrawingText()

def get_clusters_coords(loadings, atomcomb_clusters):
    return np.array([
      loadings.isel(atomcomb=c).mean(dim='atomcomb').values
      for c in atomcomb_clusters])

def calc_dist(point):
    x, y = point
    return (x**2 + y**2)**0.5

def calc_angle(point):
    x, y = point
    return float(np.degrees(np.arctan2(x, y)))

def separate_angles(points, min_angle=10):

    
    angles = [calc_angle(point) for point in points]

    def decider(i, j):
        nonlocal angles
        return abs(angles[i] - angles[j]) <= min_angle # degrees. Edge case: -179 and 179
    
    angle_clusters = cluster_general(decider, len(angles))
    scalefactors = {}
    for angle_cluster in angle_clusters:
        if len(angle_cluster) < 2:
            continue
        dists = np.array(
            [(idx, calc_dist(points[idx])) for idx in angle_cluster],
            dtype=[('idx', int), ('dist', float)])
        dists.sort(order='dist')
        factor = 1
        for idx, dist in dists[::-1]:
            scalefactors[idx] = factor # less extrusion for the smaller radius
            factor *= 0.8
    return scalefactors

def filter_cluster_coords(coords, n):
    radii = [(x**2 + y**2)**0.5 for x, y in coords]  
    angles = [np.degrees(np.arctan2(x, y)) for x, y in coords]
    res = set(np.argsort(radii)[-(n-2):])
    avg = np.mean(angles)
    splay = [abs(avg-angle) for angle in angles]
    return res.union(np.argsort(splay)[-2:])

def extrude(x, y, xmin, xmax, ymin, ymax):
    # for extrusion, flip negative rays into quadrant 1
    if x < 0:
        xlim = -xmin # positive
        xsgn = -1
    else:
        xlim = xmax
        xsgn = 1
    if y < 0:
        ylim = -ymin # positive
        ysgn = -1
    else:
        ylim = ymax
        ysgn = 1
    # now extrude
    x2 = abs(ylim*x/y)  # try extruding till we reach the top
    if x2 <= xlim: # have we dropped off the right?
        y2 = ylim  # if not, go with this
    else:          # but if we would have dropped off the right
        x2 = xlim  # just go as far right as possible instead   
        y2 = abs(xlim*y/x)
    return x2*xsgn, y2*ysgn

def mpl_imshow_png(ax, png, **imshow_kws):
    buffer = io.BytesIO()
    buffer.write(png)
    buffer.seek(0)
    img_array = np.array(PIL.Image.open(buffer))
    ax.axis('off')
    return ax.imshow(img_array, rasterized=True, **imshow_kws)

def plot_clusters2(ax, loadings, clusters, mol, min_angle=10, inset_scale=1, show_at_most=None):
    points = get_clusters_coords(loadings, clusters)
    if show_at_most is not None:
        indices = filter_cluster_coords(points, show_at_most)
        # clusters = [c for i, c in enumerate(clusters) if i in indices]
        # points = [p for i, p in enumerate(points) if i in indices]
    else:
        indices = range(len(clusters))
    scalefactors = separate_angles(points, min_angle)

    # else:
        # indices = range(len(clusters))
    for i, cluster in enumerate(clusters):
        acs = loadings.isel(atomcomb=cluster)
        x, y = acs.mean(dim='atomcomb')
        s = ' '.join([f'({a1},{a2})' for a1, a2 in acs.atomcomb.values])
        arrow_color = 'k' if i in indices else (0, 0, 0, 0.5)
        ax.arrow(0, 0, x, y, head_width=0.01, length_includes_head=True, color=arrow_color)

        scale = scalefactors.get(i, 1)
        
        x2, y2 = extrude(x, y, *ax.get_xlim(), *ax.get_ylim())
        x2 *= 0.8*scale
        y2 *= 0.8*scale     

        if i not in indices:
            continue

        ax.plot([x,x2], [y,y2], '--', c='darkgray', lw=0.5)

        ymin, ymax = ax.get_ylim()
        inset_size = inset_scale*np.array([7, 10])*(ymax-ymin)/65
        iax = ax.inset_axes([x2, y2, *inset_size], transform=ax.transData)
        iax.set_anchor('SW') # keep bottom-left corner of image at arrow tip!

        png = highlight_pairs(mol, acs.atomcomb.values)
        # display(Image(png))  # DEBUG
        mpl_imshow_png(iax, png)

def plot_clusters3(loadings, clusters, ax=None, labels=None, axs=None, mol=None):
    fig, ax = figax(ax)
    points = get_clusters_coords(loadings, clusters)
    if labels is None:
        labels = list('abcdefghijklmnopqrstuvwxyz')

    for i, cluster in enumerate(clusters):
        acs = loadings.isel(atomcomb=cluster)
        x, y = acs.mean(dim='atomcomb')
        s = labels[i]
        ax.arrow(0, 0, x, y, head_width=0.01, length_includes_head=True)
        
        x2, y2 = extrude(x, y, *ax.get_xlim(), *ax.get_ylim())

        ax.plot([x,x2], [y,y2], '--', c='k', lw=0.5)
        ax.text(x2, y2, s)

        if axs is not None and mol is not None:
            png = highlight_pairs(mol, acs.atomcomb.values)
            mpl_imshow_png(axs[s], png)
            axs[s].set_title(s)


##################
# New stuff for angle binning!

def get_mask(angles, theta1, theta2, seam=180):
    if not theta1 <= theta2: theta1, theta2 = theta2, theta1
    if theta1 < seam and theta2 < seam:
        mask = (angles > theta1) & (angles < theta2)
    elif theta1 < seam and theta2 > seam:
        mask = (angles > theta1) | (angles < theta2-360)
    elif theta1 > seam and theta2 > seam:
        mask = (angles > theta1-360) & (angles < theta2-360)
    else:
        print(theta1, theta2)
        mask = []
        raise ValueError()
    return mask

def circbins(angles, nbins=4, center=0):
    @dataclass
    class Sweeper:
        inner: int = 0
        outer: int = 0
    
    ppbin = len(angles)/nbins
    
    sweepers = [Sweeper(), Sweeper()] # anticlockwise; clockwise
    bins = []
    edges = []
    # first bin spreads evenly clockw and anticw
    while sweepers[0].outer + sweepers[1].outer < 180:
        sweepers[0].outer += 5
        sweepers[1].outer -= 5
        idxs = np.nonzero(
          get_mask(angles, center-sweepers[0].outer, center-sweepers[1].outer))[0]
        if len(idxs) >= ppbin:
            sweepers[0].inner = sweepers[0].outer
            sweepers[1].inner = sweepers[1].outer
            edges += [center+sweepers[0].inner, center+sweepers[1].inner]
            bins.append(idxs)
            break
    
    # intermediate bins
    for ibin in range(1, nbins-1):
        way = ibin % 2
        cur = sweepers[way] # sweepers[0] increases, sweepers[1] decreases
        sgn = [1, -1][way] # bins 1,3,5 clockw; bins 2,4,6 anti
        print(f"sweeping {["anti", "clock"][way]} from {cur.inner}°")
        while sweepers[0].outer - sweepers[1].outer < 360:
            cur.outer += sgn*10
            idxs = np.nonzero(
              get_mask(angles, center+cur.outer, center+cur.inner))[0]
            if len(idxs) >= ppbin:
                print(f"swept to {cur.outer}°")
                cur.inner = cur.outer
                bins.append(idxs)
                edges.append(center + cur.inner)
                break
    
    # last bin
    all_indices = range(len(angles))
    already_binned = np.concatenate(bins)
    bins.append(np.array([i for i in all_indices if i not in already_binned]))

    return bins, edges

def plot_bin_edges(angles, radii, bins, edges, picks, ax, labels):
    rangles = np.radians(angles)

    for e in np.radians(edges):
        ax.plot([e, e], [0, 0.4], c='gray', ls='--', lw='1')

    for a, r, s in zip(rangles[picks], radii[picks], labels[:len(picks)]):
        ax.text(a,r,s, ha='left', va='bottom', fontsize=6)

    for b, c in zip(bins, list('rgbm')):
        # ax.plot(x, y)
        # colors = ['r' if x else 'b' for x in mask]
        ax.scatter(rangles[b], radii[b], c='gray', s=5)

    ax.scatter(rangles[picks], radii[picks], c='k', s=5)

    ax.set_rlabel_position(200)