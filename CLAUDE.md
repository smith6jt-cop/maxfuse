# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MaxFuse is a Python package for integrating single-cell datasets from different modalities with no overlapping features and/or under low signal-to-noise ratio regimes. It's designed for cross-modal data integration in "weak linkage" scenarios where linked features are few or uninformative, such as integrating spatial proteomic data with single-cell sequencing data.

The package includes MARIO (Matching with Approximate Regional Importance with Overlap) for statistical matchability testing and bipartite matching.

## Repository Structure

```
maxfuse/
├── src/maxfuse/           # Main package (src layout)
│   ├── core/              # MaxFuse core functionality
│   │   ├── model.py       # Fusor class - main entry point
│   │   ├── match_utils.py # Matching algorithms
│   │   ├── graph.py       # Graph construction & clustering
│   │   ├── utils.py       # Numerical utilities
│   │   ├── spatial_utils.py # Spatial data utilities
│   │   └── metrics.py     # Evaluation metrics
│   └── mario/             # MARIO matching module
│       ├── match.py       # Mario class
│       ├── match_utils.py # MARIO matching utilities
│       ├── cluster.py     # Clustering algorithms
│       ├── embed.py       # Embedding utilities
│       └── utils.py       # MARIO utilities
├── examples/              # Example notebooks
├── docs/                  # Sphinx documentation
├── pyproject.toml         # Package configuration
└── README.md
```

## Installation & Development

```bash
# Create environment and install
conda create -n maxfuse python=3.8
conda activate maxfuse
pip install -e .  # For development installation

# Or install from PyPI
pip install maxfuse
```

**Key dependencies:** igraph, leidenalg, numpy, pandas, scanpy, scipy, scikit-learn, matplotlib, requests

**Note:** Example notebooks in `examples/` serve as integration tests.

## Usage

```python
# Import main classes
from maxfuse import Fusor, Mario

# Or import from submodules
from maxfuse.core import spatial_utils
from maxfuse.mario import match_utils
```

## Architecture

### Core Module (`src/maxfuse/core/`)

- **`model.py`**: Contains the `Fusor` class - the main entry point for the MaxFuse pipeline. Orchestrates the entire integration workflow including batching, graph construction, matching, and embedding generation. Supports spatial awareness via region priors and neighborhood features.

- **`match_utils.py`**: Matching utilities using linear sum assignment (Hungarian algorithm). Handles initial matching via correlation distance and refined matching through iterative CCA refinement.

- **`graph.py`**: Graph construction and clustering using k-NN graphs with Leiden algorithm. Uses scanpy for graph construction and igraph/leidenalg for community detection.

- **`utils.py`**: Core numerical utilities including SVD operations, CCA embedding, centroid shrinkage, graph smoothing, and data preprocessing via scanpy.

- **`spatial_utils.py`**: Spatial data utilities for handling cell neighborhood composition, spatial k-NN indices, tissue region detection, and region-aware matching priors.

- **`metrics.py`**: Evaluation metrics including matching accuracy, FOSCTTM (fraction of samples closer than true match), and k-NN alignment scores.

### MARIO Module (`src/maxfuse/mario/`)

- **`match.py`**: Contains the `Mario` class for statistical matchability testing and bipartite matching.

- **`match_utils.py`**: MARIO-specific matching utilities.

- **`cluster.py`**: Spectral clustering and k-means implementations.

- **`embed.py`**: Embedding utilities for MARIO.

### Pipeline Flow

1. **Initialize `Fusor`** with shared and active feature arrays for both modalities
2. **`split_into_batches()`**: Partition data for scalability
3. **`construct_graphs()`**: Build k-NN graphs, optionally create metacells, cluster using Leiden
4. **`find_initial_pivots()`**: Initial matching via fuzzy smoothed embedding on shared features
5. **`refine_pivots()`**: Iterative CCA refinement on active features
6. **`filter_bad_matches()`**: Quality filtering and pivot selection
7. **`propagate()`**: Extend matching to non-pivot cells via nearest neighbor search
8. **`get_matching()`** / **`get_embedding()`**: Extract final results

### Smoothing Methods

The `Fusor` class supports two smoothing approaches (set via `method` parameter):
- **`centroid_shrinkage`**: Shrinks cells toward cluster centroids
- **`graph_smoothing`**: Shrinks cells toward neighborhood averages

## Documentation

Documentation is built with Sphinx and hosted on ReadTheDocs. Build locally:

```bash
cd docs
pip install -r requirements.txt
make html
```

## Key Data Structures

- **Matching format**: Lists of `[rows, cols, vals]` where `(rows[i], cols[i])` is a matched pair with distance/score `vals[i]`
- **Graph edges**: Lists of `[rows, cols, weights]` representing weighted k-NN graph edges
- Arrays are numpy ndarrays with shape `(n_samples, n_features)`

## Data Files

Large data files are not tracked in git but are expected in `docs/` for running example notebooks:
- `docs/1904_CC2B_cell.h5ad` - CODEX protein data
- `docs/raw_feature_bc_matrix/` - 10x scRNAseq data
- `docs/1904_gene/` - Cell Ranger output

## Archive

Manuscript analysis code is preserved in the `archive-manuscript` orphan branch. This code uses a development version of MaxFuse with slightly different API grammar and should not be used as a reference for the current API.
