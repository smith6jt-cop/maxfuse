"""
Utility functions for dealing with spatial data
"""

import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN


def bind_spatial(features, nbhd, wt_on_features=0.7):
    """
    Return a new array of form [wt_on_features * features / feature_norm, (1-wt_on_features) * nbhd / nbhd_norm]

    Parameters
    ----------
    features: np.ndarray of shape (n_samples, n_features)
        Feature matrix
    nbhd: np.ndarray of shape (n_samples, n_clusters)
        Cell neighborhood composition matrix
    wt_on_features: float, default=0.7
        Weight to put on the feature matrix.

    Returns
    -------
    res: np.ndarray of shape (n_samples, n_features+n_clusters)

    """
    # normalize two kinds of info for easier tuning of weight
    feature_norm = np.linalg.norm(features)
    nbhd_norm = np.linalg.norm(nbhd)
    res = np.concatenate((
        wt_on_features * features / feature_norm,
        (1-wt_on_features) * nbhd / nbhd_norm
    ), axis=1)
    return res


def get_spatial_knn_indices(locations, n_neighbors=15, method='kd_tree'):
    """
    Compute k-nearest neighbors of locations.

    Parameters
    ----------
    locations: np.ndarray of shape (n_samples, 2)
        Data matrix
    n_neighbors: int
        Number of nearest neighbors
    method: str, default='kd_tree'
        Method to use when computing the nearest neighbors, one of ['ball_tree', 'kd_tree', 'brute']

    Returns
    -------
    knn_indices: np.ndarray of shape (n_samples, n_neighbors)
        Each row represents the knn of that sample
    """
    locations = np.array(locations)
    assert n_neighbors <= locations.shape[0]
    # k-NN indices, may be asymmetric
    _, knn_indices = NearestNeighbors(
        n_neighbors=n_neighbors, algorithm=method
    ).fit(locations).kneighbors(locations)
    return knn_indices


def get_neighborhood_composition(knn_indices, labels, log1p=False):
    """
    Compute the composition of neighbors for each sample.

    Parameters
    ----------
    knn_indices: np.ndarray of shape (n_samples, n_neighbors)
        Each row represents the knn of that sample
    labels: np.ndarray of shape (n_samples, )
        Cluster labels
    log1p: bool, default=False
        Whether to apply log1p transformation

    Returns
    -------
    comp: np.ndarray of shape (n_samples, n_neighbors)
        The composition (in proportion) of neighbors for each sample.
    """
    labels = list(labels)
    n, k = knn_indices.shape
    unique_clusters = np.unique(labels)
    n_clusters = len(unique_clusters)
    label_to_clust_idx = {label: i for i, label in enumerate(unique_clusters)}

    comp = np.zeros((n, n_clusters))
    for i, neighbors in enumerate(knn_indices):
        good_neighbors = [nb for nb in neighbors if nb != -1]
        for nb in good_neighbors:
            comp[i, label_to_clust_idx[labels[nb]]] += 1

    if log1p:
        comp = np.log1p(comp)
    return comp


def detect_tissue_regions(locations, marker_expression, marker_names,
                          marker_to_region, n_neighbors=30, min_cluster_size=10,
                          eps_quantile=0.1):
    """
    Auto-detect tissue regions using a combined approach:
    1. Classify cells by dominant marker expression
    2. Spatially cluster cells of each type using DBSCAN
    3. Assign region labels based on marker identity + spatial coherence

    Parameters
    ----------
    locations : np.ndarray of shape (n_cells, 2)
        Spatial coordinates (X, Y centroids)
    marker_expression : np.ndarray of shape (n_cells, n_markers)
        Protein marker expression values
    marker_names : list of str
        Names of markers corresponding to columns of marker_expression
    marker_to_region : dict
        Maps marker names to region names.
        e.g., {'CD20': 'B_follicle', 'CD3e': 'T_zone', 'CD68': 'Red_pulp'}
    n_neighbors : int, default=30
        Number of neighbors for spatial density estimation
    min_cluster_size : int, default=10
        Minimum number of cells to form a spatial cluster
    eps_quantile : float, default=0.1
        Quantile of k-NN distances to use as DBSCAN eps parameter

    Returns
    -------
    regions : np.ndarray of shape (n_cells,)
        Region assignment for each cell. Cells not in coherent regions
        are labeled as 'mixed'.
    region_info : dict
        Additional information about detected regions including:
        - 'n_regions': number of distinct regions found
        - 'region_counts': dict of region name to cell count
    """
    n_cells = len(locations)
    regions = np.full(n_cells, 'mixed', dtype=object)

    # Build marker name to index mapping
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    # Get available markers that map to regions
    available_markers = [m for m in marker_to_region.keys() if m in marker_to_idx]
    if not available_markers:
        return regions, {'n_regions': 0, 'region_counts': {'mixed': n_cells}}

    # Step 1: Classify cells by dominant marker
    marker_indices = [marker_to_idx[m] for m in available_markers]
    relevant_expression = marker_expression[:, marker_indices]

    # Normalize expression per marker (z-score)
    expr_mean = relevant_expression.mean(axis=0)
    expr_std = relevant_expression.std(axis=0)
    expr_std[expr_std < 1e-10] = 1.0
    normalized_expr = (relevant_expression - expr_mean) / expr_std

    # Assign each cell to marker with highest normalized expression
    # Only assign if expression is above threshold (z > 0.5)
    dominant_marker_idx = normalized_expr.argmax(axis=1)
    max_expr = normalized_expr.max(axis=1)
    high_expr_mask = max_expr > 0.5

    cell_marker_type = np.full(n_cells, 'none', dtype=object)
    for i in range(n_cells):
        if high_expr_mask[i]:
            cell_marker_type[i] = available_markers[dominant_marker_idx[i]]

    # Step 2: Compute adaptive eps for DBSCAN based on local density
    knn_dists = NearestNeighbors(n_neighbors=n_neighbors).fit(locations).kneighbors(locations)[0]
    eps = np.quantile(knn_dists[:, -1], eps_quantile)

    # Step 3: Spatially cluster cells of each marker type
    region_counts = {'mixed': 0}

    for marker in available_markers:
        region_name = marker_to_region[marker]
        marker_mask = cell_marker_type == marker

        if marker_mask.sum() < min_cluster_size:
            continue

        # Get spatial coordinates of cells with this marker
        marker_locations = locations[marker_mask]
        marker_indices_orig = np.where(marker_mask)[0]

        # DBSCAN clustering
        clustering = DBSCAN(eps=eps, min_samples=min_cluster_size // 2).fit(marker_locations)

        # Assign region to cells in valid clusters (label >= 0)
        for i, label in enumerate(clustering.labels_):
            if label >= 0:  # Not noise
                orig_idx = marker_indices_orig[i]
                regions[orig_idx] = region_name

        if region_name not in region_counts:
            region_counts[region_name] = 0
        region_counts[region_name] += (clustering.labels_ >= 0).sum()

    # Count mixed cells
    region_counts['mixed'] = (regions == 'mixed').sum()

    region_info = {
        'n_regions': len([r for r in region_counts if r != 'mixed' and region_counts[r] > 0]),
        'region_counts': region_counts
    }

    return regions, region_info


def compute_region_celltype_prior(celltype_to_region_weights, rna_labels,
                                   spatial_regions, default_weight=1.0):
    """
    Build prior distance matrix for interpolation with embedding distances.

    The prior encodes biological expectations about which cell types should
    match to which tissue regions. Lower values indicate more compatible
    matches (preferred), higher values indicate less compatible matches
    (discouraged).

    Parameters
    ----------
    celltype_to_region_weights : dict
        Nested dict mapping {celltype: {region: weight}}.
        Lower weight = more compatible (e.g., 0.1 for expected matches).
        Higher weight = less compatible (e.g., 5.0 for unlikely matches).
        Example:
            {
                'B_cell': {'B_follicle': 0.1, 'T_zone': 2.0, 'Red_pulp': 5.0},
                'T_cell': {'B_follicle': 2.0, 'T_zone': 0.1, 'Red_pulp': 3.0},
            }
    rna_labels : np.ndarray of shape (n_rna_cells,)
        Cell type labels for RNA cells
    spatial_regions : np.ndarray of shape (n_spatial_cells,)
        Region assignments for spatial (CODEX) cells
    default_weight : float, default=1.0
        Weight for unspecified celltype-region pairs (neutral)

    Returns
    -------
    prior_dist : np.ndarray of shape (n_rna_cells, n_spatial_cells)
        Prior distance matrix where prior_dist[i,j] reflects compatibility
        between RNA cell i's type and spatial cell j's region.
        Used in: final_dist = (1-w)*embed_dist + w*prior_dist
    """
    n_rna = len(rna_labels)
    n_spatial = len(spatial_regions)

    # Initialize with default weight
    prior_dist = np.full((n_rna, n_spatial), default_weight, dtype=np.float64)

    # Get unique values for efficient lookup
    unique_types = np.unique(rna_labels)
    unique_regions = np.unique(spatial_regions)

    # Build index mappings for vectorized assignment
    for celltype in unique_types:
        if celltype not in celltype_to_region_weights:
            continue

        type_mask = rna_labels == celltype
        type_indices = np.where(type_mask)[0]

        region_weights = celltype_to_region_weights[celltype]

        for region in unique_regions:
            weight = region_weights.get(region, default_weight)
            region_mask = spatial_regions == region
            region_indices = np.where(region_mask)[0]

            # Assign weight to all (celltype, region) pairs
            for rna_idx in type_indices:
                prior_dist[rna_idx, region_indices] = weight

    # Normalize to have similar scale as correlation distances (0-2 range)
    # Map weights so that: 0.1 -> ~0.1, 1.0 -> ~1.0, 5.0 -> ~1.5
    # Using log transformation for smoother scaling
    prior_dist = np.log1p(prior_dist)

    return prior_dist


def compute_neighborhood_augmented_features(features, locations, labels,
                                             n_neighbors=15, wt_on_features=0.7,
                                             log1p=False):
    """
    Augment feature matrix with spatial neighborhood composition.

    This is a convenience wrapper that combines get_spatial_knn_indices(),
    get_neighborhood_composition(), and bind_spatial() into a single call.

    Parameters
    ----------
    features : np.ndarray of shape (n_samples, n_features)
        Original feature matrix (e.g., protein expression)
    locations : np.ndarray of shape (n_samples, 2)
        Spatial coordinates (X, Y)
    labels : np.ndarray of shape (n_samples,)
        Cell type or cluster labels for computing neighborhood composition
    n_neighbors : int, default=15
        Number of spatial neighbors to consider
    wt_on_features : float, default=0.7
        Weight on original features vs neighborhood composition.
        Higher values weight features more heavily.
    log1p : bool, default=False
        Whether to log-transform neighborhood counts before combining

    Returns
    -------
    augmented : np.ndarray of shape (n_samples, n_features + n_unique_labels)
        Features concatenated with weighted neighborhood composition
    """
    # Get spatial k-NN indices
    knn_indices = get_spatial_knn_indices(locations, n_neighbors=n_neighbors)

    # Compute neighborhood composition
    nbhd_composition = get_neighborhood_composition(knn_indices, labels, log1p=log1p)

    # Combine features with neighborhood
    augmented = bind_spatial(features, nbhd_composition, wt_on_features=wt_on_features)

    return augmented
