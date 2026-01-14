# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MaxFuse is a Python package for integrating single-cell datasets from different modalities with no overlapping features and/or under low signal-to-noise ratio regimes. It's designed for cross-modal data integration in "weak linkage" scenarios where linked features are few or uninformative, such as integrating spatial proteomic data with single-cell sequencing data.

The package includes MARIO (Matching with Approximate Regional Importance with Overlap) for statistical matchability testing and bipartite matching.

## Repository Structure

```
repo/
├── src/
│   ├── maxfuse/               # Self-contained package (portable)
│   │   ├── __init__.py        # from maxfuse import Fusor, Mario
│   │   ├── pyproject.toml     # Package configuration
│   │   ├── README.md          # Package readme
│   │   ├── core/              # MaxFuse core functionality
│   │   │   ├── model.py       # Fusor class - main entry point
│   │   │   ├── match_utils.py # Matching algorithms
│   │   │   ├── graph.py       # Graph construction & clustering
│   │   │   ├── utils.py       # Numerical utilities
│   │   │   ├── spatial_utils.py # Spatial data utilities
│   │   │   └── metrics.py     # Evaluation metrics
│   │   └── mario/             # MARIO matching module
│   │       ├── match.py       # Mario class
│   │       ├── match_utils.py # MARIO matching utilities
│   │       ├── cluster.py     # Clustering algorithms
│   │       ├── embed.py       # Embedding utilities
│   │       └── utils.py       # MARIO utilities
│   └── mario-py/              # Original MARIO repo (reference)
├── notebooks/                 # Analysis notebooks
├── data/                      # Data files (not tracked in git)
├── results/                   # Output files (not tracked in git)
├── docs/                      # Sphinx documentation
├── media/                     # Images and media files
├── environment.yml            # Conda environment
└── CLAUDE.md                  # This file
```

## Installation & Development

```bash
# Create environment and install
conda env create -f environment.yml
conda activate maxfuse
pip install -e src/maxfuse  # For development installation

# Or install from PyPI (when published)
pip install maxfuse
```

**Key dependencies:** igraph, leidenalg, numpy, pandas, scanpy, scipy, scikit-learn, matplotlib, requests

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

## Key Data Structures

- **Matching format**: Lists of `[rows, cols, vals]` where `(rows[i], cols[i])` is a matched pair with distance/score `vals[i]`
- **Graph edges**: Lists of `[rows, cols, weights]` representing weighted k-NN graph edges
- Arrays are numpy ndarrays with shape `(n_samples, n_features)`

## Notebook Pipeline

The analysis notebooks form a connected pipeline via checkpoint files:

```
notebooks/
├── 1_preprocessing.ipynb   → results/1_preprocessing/
│   Load: data/*.tsv, data/raw_feature_bc_matrix/
│   Save: protein_adata.h5ad, rna_adata.h5ad, rna_adata_lognorm.h5ad
│
├── 2_integration.ipynb     → results/2_integration/
│   Load: results/1_preprocessing/*
│   Save: maxfuse_matching.pkl, normalized arrays, correspondence.csv
│
├── 3_visualization.ipynb   → results/3_visualization/
│   Load: results/1_preprocessing/*, results/2_integration/*
│   Save: figures, spatial clusters, aligned indices, cell types
│
└── 5_analysis.ipynb        → results/5_analysis/
    Load: results/1_preprocessing/*, results/2_integration/*, results/3_visualization/*
    Save: validation results, cross-validation, cell type expression profiles
```

### Running the Pipeline

1. Run `1_preprocessing.ipynb` first (loads raw data, QC filtering)
2. Run `2_integration.ipynb` (MaxFuse/MARIO integration)
3. Run `3_visualization.ipynb` (visualizations, spatial mapping)
4. Run `5_analysis.ipynb` (statistical validation, gene expression inference)

Each notebook checks for required input files and provides clear error messages if prerequisites are missing.

### Checkpoint File Formats

| Data Type | Format | Location |
|-----------|--------|----------|
| AnnData objects | `.h5ad` | `results/*/` |
| Matching results | `.pkl`, `.csv` | `results/2_integration/` |
| NumPy arrays | `.npy` | `results/2_integration/` |
| Parameters | `.json` | `results/*/` |

## Data Files

Data files are stored in `data/` (not tracked in git):
- CODEX protein data (.h5ad files, .tsv from QuPath)
- scRNAseq data (10x format: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
- Conversion tables (.csv)

Results are stored in `results/` (not tracked in git), organized by notebook.

## Integration Validation Best Practices

### Cross-Validation Direction
When validating RNA-protein integration accuracy, the correct question is:
- **Correct**: "Can RNA expression predict protein levels?" (RNA → Protein)
- **Wrong**: "Can protein predict RNA?" (small panel can't explain large panel)

The rich RNA panel (18k genes) should predict the sparse protein panel (26 markers), not vice versa.

### Key Metrics
1. **Spearman correlation** (shared features): Validates that matched cells have correlated expression
2. **R² from cross-validation**: Measures how well RNA predicts protein (expect 5-35% for good markers)
3. **Permutation p-value**: Tests if matching is better than random (this is the critical metric)

### Interpreting Results
- Low R² (4-7%) is normal for RNA→protein prediction due to post-transcriptional regulation
- Permutation p-value < 0.01 confirms the integration captured real biological signal
- Best performing markers: B cell (MS4A1/CD20), macrophage (CD68), T cell (CD3E)

## Archive

Manuscript analysis code is preserved in the `archive-manuscript` orphan branch.

## Skills Registry

This project uses a skills registry for Claude memory persistence at `.skills_registry/`.

### Commands
- `/advise` - Search skills before starting work (finds related experiments, what worked/failed)
- `/retrospective` - Save learnings from current session as a new skill

### Available Skills
- `plugins/maxfuse/repo-reorganization/` - Python package with pyproject.toml inside package directory
- `plugins/maxfuse/region-aware-matching/` - Spatial region-aware cell matching for CODEX/scRNAseq
- `plugins/scientific/notebook-checkpoint-pattern/` - Connecting notebooks via checkpoint files
- `plugins/scientific/plotly-figurewidget-interactive/` - Fix for Plotly FigureWidget interactive dashboard bugs
- `plugins/scientific/project-data-separation/` - Separating repository code from user data
- `plugins/scientific/sparse-expression-visualization/` - Bar plots vs boxplots for sparse scRNA-seq data

### Adding New Skills
1. Copy `templates/experiment-skill-template/` to `plugins/maxfuse/your-skill-name/`
2. Update `plugin.json` with trigger conditions
3. Fill in `SKILL.md` with goal, what worked, what failed, exact parameters
4. Push to the Skills_Registry repo
