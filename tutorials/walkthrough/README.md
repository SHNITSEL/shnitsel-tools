# Walkthorugh Examples

## Step 1: Filter by geometry

Exclude trajectories where bond lengths exceed a threshold of x (C-C), xx (C-H), xx (C-N), or xx (C-H).
- [`01_filter_geoms.ipynb`]()

## Step 2: PCA across Compounds

Perform PCA across compounds, by concatenating xr.Datasets of I01, A01, A02 and A03, and use the respective pairwise distances as features for PCA.
Visualize projection on first two PCs for all compounds and highlight the points of a certain substance.

- [`02_pca_across_compounds.ipynb`]()

## Step 3: Clustering of Trajectories

In case of A03, isomerization is observed in some of the trajectories.
They can be identified by the CC=CC dihedral angles at the end-points of the trajectories.
Here we use the last 20 frames and perform PCA on this subset, to project the datapoints on two axis color highlight by dihedral.
Subsequently, we perform k-means clustering, to label them as `Z-to-Z` and `Z-to-E` (productive).

- [`03_clustering_a03.ipynb`]()

## Step 4: Kinetics across Compounds

Generate population data for the photophysical trajectories of A01, A02 and A03 (`Z-to-Z` labelled data), i.e. not the photochemical reaction of butene (`Z-to-E`labelled data, see Step 3).
Compare population curves in one plot.

- [`04_kinetics_across_compounds.ipynb`]()
  


