# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## ⚠️ CRITICAL: Never Use Cell Numbers

**NEVER say "cell 38", "Cell 12", or any cell index.** Cell numbers are:
- NOT visible to users in Jupyter
- Change when cells are added/deleted
- Impossible for users to find

**Always refer to cells by:**
- First line of code: "the cell starting with `# PRE-FILTER:`"
- Section header: "the cell under '## Load Data'"
- Purpose: "the pre-filter cell", "the normalization cell"
- Variable: "the cell that creates `protein_shared`"

This applies to ALL communication - tool calls can use indices internally, but messages to users must use descriptive references.

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

## Notebook Cell References

See **⚠️ CRITICAL: Never Use Cell Numbers** at the top of this file. This rule is absolute - there are no exceptions.

## Editing Large Notebooks (Token Limit Workarounds)

Large notebooks with many cells or outputs can exceed Claude's token limit, causing incomplete reads or silent failures. Use these strategies:

### Strategy 1: Direct Cell Editing (Preferred)

**Don't read the whole notebook.** Use the `NotebookEdit` tool directly with cell index:

```
User: "Edit cell 39 in notebooks/2_integration.ipynb - change cca_components from 25 to 7"
Claude: [Uses NotebookEdit with cell index 39, no Read needed]
```

To find cell indices, run:
```bash
python scripts/list_cells.py notebooks/2_integration.ipynb
```

### Strategy 2: User Provides Cell Content

Instead of asking Claude to read the notebook:
```
User: "Here's the cell I want to modify:
[paste cell code]
Change cca_components from 25 to 7"
```

### Strategy 3: Read Specific Line Ranges

For partial notebook inspection:
```
Read notebooks/2_integration.ipynb lines 500-600
```

### Strategy 4: Describe the Cell Location

Reference cells by surrounding context:
```
User: "Edit the cell under '## Step 7: MaxFuse Integration' that starts with 'fusor = Fusor('"
```

### Helper Script: `scripts/list_cells.py`

Lists all cells with indices and first lines (without full content):
```bash
# List cells in a notebook
python scripts/list_cells.py notebooks/2_integration.ipynb

# Add cell IDs to old notebooks (upgrades to nbformat 4.5)
python scripts/list_cells.py notebooks/*.ipynb --add-ids
```

Output format:
```
Idx  | Cell ID      | Type     | First Line
-----|--------------|----------|----------------------------------
0    | -            | markdown | # Integration: MaxFuse and MARIO
1    | -            | code     | import numpy as np
39   | -            | code     | fusor = Fusor(
```

### MCP Notebook Tools

The `notebook` MCP server (`mcp-jupyter`) provides:
- `list_cells` - List cells with indices (no full content)
- `get_cell_source` - Read specific cell by index
- `edit_cell_source` - Modify cell by index
- `insert_cell` / `delete_cell` - Add/remove cells

These tools allow targeted edits without loading the entire notebook into context.

### When Claude Must Read a Large Notebook

If full context is needed:
1. **Split into smaller notebooks** - Separate preprocessing, integration, visualization
2. **Extract large code to `.py` modules** - Keep cells as thin wrappers
3. **Clear outputs** before reading - Run `jupyter nbconvert --clear-output --inplace notebook.ipynb`

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

## MaxFuse Parameter Tuning

### Critical: Scale Alignment

**BOTH modalities must be z-scored to the same scale before integration:**
- RNA: normalize_total → log1p → z-score (mean≈0, std≈1, range [-5, 5])
- Protein: z-score using StandardScaler (mean≈0, std≈1, range [-5, 5])

If protein is used "as-is" without z-scoring, RNA will dominate matching due to variance mismatch.

### Parameter Scaling by Panel Size

When protein panel size changes, parameters must be adjusted. **Critical rule: All SVD components must be < n_features.**

| Parameter | 26 Markers | 59 Markers | 19 Markers | Rule |
|-----------|-----------|-----------|------------|------|
| CCA Components | 25 | 7 | **10** | `min(n_shared - 1, 10)` for small panels |
| SVD for CCA (RNA) | 40 | 50 | 50 | Increase with more protein info |
| SVD for CCA (Protein) | None | 35 | **15** | `min(15, n_prot - 1)` |
| SVD for Graphs (Protein) | 15 | 30 | **15** | `min(15, n_prot - 1)` |
| SVD for Initial Pivots | 25/20 | 20/18 | **15/15** | `min(15, n_shared - 1)` |
| Metacell Size | 2 | 2 | **1** | Disable for <25 features |

**CCA Component Rule**: `cca_components = min(n_shared - 1, 10)` for small panels (<25 features)

**Metacell Rule**: With <25 features, metacell_size=2 provides minimal noise reduction. Use `metacell_size=1` (disabled) or increase to 5+ if noise is genuinely problematic.

Using too many CCA/SVD components with small panels causes overfitting (trivially perfect correlations).

### Diagnostic Checks

1. **Scale alignment**: Both modalities should have mean≈0, std≈0.8-1.0
2. **Canonical correlations**: Should decay smoothly, not all >0.95 (overfitting sign)
3. **UMAP alignment**: Matched pairs should be close, not crossing clusters

## Pre-filtering for CD45+ RNA Data

When RNA data is CD45+ sorted (immune cells only), protein data may contain non-immune cells that shouldn't match.

### Two-Stage Pre-filtering

```python
# Stage 1: Marker-specific z-score filtering
# Remove cells where specific markers have extreme low values (e.g., detection failures)
MARKER_ZSCORE_THRESHOLDS = {
    'MPO': -2.5,  # Filter cells where MPO z-score < -2.5
    # Add more markers as needed
}

for marker, threshold in MARKER_ZSCORE_THRESHOLDS.items():
    marker_idx = marker_names.index(marker)
    keep_mask &= protein_shared[:, marker_idx] >= threshold

# Stage 2: Immune cell filtering (for Pancreas tissue)
# Keep only protein cells with immune signature above threshold
immune_markers = ['CD3E', 'MS4A1', 'CD68', 'PTPRC']
immune_score = protein_shared[:, immune_marker_idx].sum(axis=1)
threshold = np.percentile(immune_score[is_pln], 25)  # Use pLN as reference
keep_mask &= (is_pln) | (immune_score > threshold)
```

### Why Pre-filter?
- Non-immune cells create "bad mode" in bimodal matching scores
- Cross-tissue matches (pLN RNA → Pancreas epithelial) are biologically meaningless
- Filtering BEFORE Fusor creation is critical - arrays must be updated before integration

### Cell Execution Order Matters
After pre-filtering, **re-run all subsequent cells** to use filtered arrays. The Fusor and tissue priors must see the reduced data.

## Diagnosing Bimodal Matching Scores

If matching scores show two distinct modes (e.g., peaks at 0.2 and 0.8):

### Diagnostic Checklist
1. **Tissue mismatch**: Check cross-tabulation of RNA tissue vs protein tissue in matches
2. **Expression levels**: Compare total protein expression between score modes
3. **Shared feature correlation**: Calculate per-pair correlation across shared features
4. **CCA effect**: Compare initial vs refined score distributions (kurtosis)

### Common Causes
| Symptom | Cause | Fix |
|---------|-------|-----|
| Cross-tissue matches in bad mode | Tissue imbalance (e.g., 476k Pancreas vs 21k pLN) | Increase region prior weight to 0.7+ |
| Bad mode has low/negative expression | Non-immune cells in protein data | Pre-filter by immune score |
| CCA increases bimodality | CCA correctly separates matchable/unmatchable | Use GMM threshold to filter bad matches |
| Specific marker at -3.0 | Detection failure or normalization artifact | Filter cells by marker z-score threshold |

## Integration Validation Best Practices

### Cross-Validation Direction
When validating RNA-protein integration accuracy, the correct question is:
- **Correct**: "Can RNA expression predict protein levels?" (RNA → Protein)
- **Wrong**: "Can protein predict RNA?" (small panel can't explain large panel)

The rich RNA panel (18k genes) should predict the protein panel (59 markers), not vice versa.

### Key Metrics
1. **Spearman correlation** (shared features): Validates that matched cells have correlated expression
2. **R² from cross-validation**: Measures how well RNA predicts protein (expect 5-35% for good markers)
3. **Permutation p-value**: Tests if matching is better than random (this is the critical metric)

### Interpreting Results
- Low R² (4-7%) is normal for RNA→protein prediction due to post-transcriptional regulation
- Permutation p-value < 0.01 confirms the integration captured real biological signal
- Best performing markers: B cell (MS4A1/CD20), macrophage (CD68), T cell (CD3E)

## Critical: Notebook Data Consistency

When modifying data sources in notebooks, **check the ENTIRE notebook** for dependent variables:

### Common Pitfall: Mixed Data Sources
If switching from one data source to another (e.g., raw protein data → gated protein data), ALL arrays derived from that source must be updated:
- `protein_shared` and `protein_active` must come from the SAME source
- Shape mismatches (e.g., 1.9M vs 1.2M cells) indicate mixed data sources
- Search for ALL references to the old variable name before making changes

### Gated Protein Data
- Location: `data/6551_leiden_umap.h5ad` (1.2M cells, 59 markers)
- Created by: `6551_Analysis.ipynb` using scimap gating with marker-specific thresholds
- Transformation: raw → scimap gate → log1p → sc.pp.scale()
- **Do NOT guess detection thresholds** - they were set per-marker in napari/scimap

### Never Guess Parameters
When unclear about thresholds, scaling, or transformations:
1. Find the notebook/script that created the data
2. Read the actual processing steps
3. Use those exact values
4. If source is unclear, ASK - don't try multiple guesses

### Before Changing Any Data Source
1. Grep for all references to the variable being replaced
2. Trace which arrays are derived from it
3. Update ALL downstream uses
4. Verify shapes match across modalities

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
- `plugins/maxfuse/parameter-scaling/` - MaxFuse parameter tuning for different protein panel sizes (19, 26, 59+ markers)
- `plugins/maxfuse/cross-modal-normalization/` - Scale alignment for RNA-protein integration (BOTH must be z-scored)
- `plugins/maxfuse/bimodal-score-diagnosis/` - Diagnosing bimodal matching scores (tissue mismatch, non-immune cells, marker issues)
- `plugins/scientific/notebook-checkpoint-pattern/` - Connecting notebooks via checkpoint files
- `plugins/scientific/notebook-cell-references/` - How to refer to notebook cells (NEVER use cell numbers)
- `plugins/scientific/notebook-data-consistency/` - Check ALL dependent arrays when changing data sources
- `plugins/scientific/plotly-figurewidget-interactive/` - Fix for Plotly FigureWidget interactive dashboard bugs
- `plugins/scientific/project-data-separation/` - Separating repository code from user data
- `plugins/scientific/sparse-expression-visualization/` - Bar plots vs boxplots for sparse scRNA-seq data

### Adding New Skills
1. Copy `templates/experiment-skill-template/` to `plugins/maxfuse/your-skill-name/`
2. Update `plugin.json` with trigger conditions
3. Fill in `SKILL.md` with goal, what worked, what failed, exact parameters
4. Push to the Skills_Registry repo
