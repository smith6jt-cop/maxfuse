# Notebooks CLAUDE.md

Instructions specific to the analysis notebooks in this folder.

## Notebook Pipeline

```
1_preprocessing.ipynb  → results/1_preprocessing/
2_integration.ipynb    → results/2_integration/
3_visualization.ipynb  → results/3_visualization/
4_analysis.ipynb       → results/4_analysis/
```

Run notebooks in order. Each saves checkpoints that subsequent notebooks load.

## Key Variables in 2_integration.ipynb

### Normalized Arrays (for MaxFuse input)
- `rna_shared` / `protein_shared` - Normalized shared features
- `rna_active` / `protein_active` - Normalized active features (HVGs for RNA, all non-excluded for protein)

### Pre-normalization Arrays (for visualization)
- `rna_shared_raw` - Raw counts
- `rna_after_log` - After log1p transformation
- `protein_shared_raw` - Log-transformed fluorescence (from `protein_adata.layers['log']`)

### Detection Masks
- `rna_detection_mask` - Boolean, True where count > 0
- `protein_detection_mask` - Boolean, True where gate-positive

### Constants
- `ZERO_VALUE = -1.0` - Fixed value for undetected RNA features (tuned to reduce bimodality)
- `EXCLUDED_MARKERS` - List of markers to exclude (non-immune markers for CD45+ RNA data)

## Normalization Strategy

### RNA (Sparse Data)
Detection-aware normalization:
1. `normalize_total` (library size normalization)
2. `log1p` (variance stabilization)
3. Per-feature z-score on **detected cells only**
4. Undetected cells → `ZERO_VALUE`

### Protein (Continuous Data)
Standard z-score on LOG-TRANSFORMED fluorescence:
1. Extract log layer: `protein_adata.layers['log']` (NOT gated 0-1 values)
2. StandardScaler z-score (mean=0, std=1)
3. Gate-based detection mask for comparisons only

**Key insight**: Protein data has no true zeros (it's fluorescence intensity), so detection-aware normalization is wrong for protein. Use standard z-score.

**CRITICAL**: Use the LOG layer, NOT gated 0-1 values. The gated values are right-skewed (56% < 0 after z-score, skewness=0.62), causing bimodal matching. Log values have better symmetry (53% < 0, skewness=0.30).

## Large Cell Ratio Matching (>100:1 protein:RNA)

When protein cells >> RNA cells, use these parameter adjustments:

```python
# Sqrt-scaled matching ratio (not linear)
ratio = n_protein / n_rna
matching_ratio = min(100, max(10, int(np.sqrt(ratio)) + 5))

# GMM-guided pivot filtering (not fixed 20%)
pivot_filter_prop = min(0.2, bad_mode_fraction + 0.02)

# Reduced propagation filter
propagate_filter_prop = 0.05  # Not 0.1
```

See skill: `.skills_registry/plugins/maxfuse/large-cell-ratio-matching/`

## Visualization Best Practices

### Scatter Plot Labels
Use `adjustText` to prevent label overlap:
```python
from adjustText import adjust_text

texts = []
for i, name in enumerate(names):
    texts.append(ax.text(x[i], y[i], name, fontsize=7))
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5))
```

### Per-Feature Expression Plots
**Problem**: Z-scored means are ~0 by design, making bar charts useless.

**Solution**: Use pre-normalization values with rank-scaling:
```python
from scipy.stats import rankdata

# Compute pre-normalization means
rna_log_means = [rna_after_log[detected, i].mean() for i in range(n_features)]
prot_raw_means = [protein_shared_raw[:, i].mean() for i in range(n_features)]

# Rank-scale to 0-1 for visualization
rna_ranks = rankdata(rna_log_means, method='average')
rna_scaled = (rna_ranks - 1) / (n_features - 1)
```

### Detection Rate Scatter
**Problem**: Low detection rates cluster in bottom-left corner.

**Solution**: Sqrt-scale both axes:
```python
rna_det_sqrt = np.sqrt(rna_det_rates)
prot_det_sqrt = np.sqrt(prot_det_rates)
ax.scatter(rna_det_sqrt, prot_det_sqrt, ...)
ax.set_xlabel('RNA detection (√%)')
```

### Bar Charts with Many Features
Always rotate x-axis labels:
```python
ax.set_xticklabels([f[:6] for f in feature_names], rotation=45, ha='right', fontsize=7)
```

### Histogram Auto-scaling
Don't hardcode bin ranges - use actual data:
```python
bins = np.linspace(data.min(), data.max(), 50)
ax.hist(data, bins=bins, ...)
```

## Marker Exclusions (CD45+ RNA Data)

When RNA data is CD45+ only (immune cells), exclude non-immune markers:
```python
EXCLUDED_MARKERS = [
    'Vimentin', 'Podoplanin', 'aSMA', 'CD31', 'Cytokeratin',
    'E-Cadherin', 'Pan-Cytokeratin', 'EpCAM', 'Collagen-IV',
    'Synaptophysin', 'Chromogranin-A', 'CD35', 'B2M', 'Na-K-ATPase',
    'PCNA', 'Ki67', 'p53', 'Caspase-3', 'DNA1', 'DNA2', 'Background'
]
```

Apply consistently in BOTH:
1. Shared feature matching cell
2. Protein active features cell

## Pre-filtering Cells

When RNA is CD45+ sorted, add pre-filtering BEFORE Fusor creation:

### Marker-Specific Z-Score Filtering
Filter cells where specific markers have extreme values without excluding the marker:
```python
MARKER_ZSCORE_THRESHOLDS = {
    'MPO': -2.5,  # Remove cells where MPO < -2.5 (keeps MPO as a feature)
}
for marker, threshold in MARKER_ZSCORE_THRESHOLDS.items():
    marker_idx = marker_names.index(marker)
    keep_mask &= protein_shared[:, marker_idx] >= threshold
```

### Immune Cell Filtering
Remove non-immune Pancreas cells:
```python
immune_markers = ['CD3E', 'MS4A1', 'CD68', 'PTPRC']
immune_score = protein_shared[:, marker_idx].sum(axis=1)
threshold = np.percentile(immune_score[is_pln], 25)
keep_mask &= (is_pln) | (immune_score > threshold)
```

### Critical: Cell Execution Order
After pre-filtering updates `protein_shared`, `protein_active`, `protein_adata`:
1. **Re-run ALL subsequent cells** - Fusor, batching, priors, etc.
2. Add verification prints to confirm filtered data is used:
```python
print(f"protein_shared: {protein_shared.shape}")
if protein_shared.shape[0] > 400000:
    print("WARNING: Filter may not have run!")
```

## Common Issues

### "Protein distribution shows bar at -3"
Old cached values from detection-aware normalization. **Fix**: Restart kernel, re-run all cells.

### "Only 77% RNA coverage after integration"
`matching_ratio` too high for large cell ratios. **Fix**: Use sqrt-scaling (see above).

### "Per-Feature Mean plot shows only flat lines"
Plotting z-scored means which are ~0 by design. **Fix**: Use pre-normalization values with rank-scaling.

### "Labels overlapping in scatter plot"
**Fix**: Use `adjustText` library (see above).

### "Bimodal matching scores (peaks at 0.2 and 0.8)"
Usually tissue mismatch or non-immune cells. **Diagnose**: Check tissue cross-tabulation, expression levels by mode. **Fix**: Pre-filter + increase region prior weight.

### "Filtered data not being used"
Cell execution order issue. **Fix**: After pre-filter cell, use "Run All Below" or manually run each subsequent cell in order.

### "Gap in z-scored values between -3.0 and -2.8"
Specific marker has detection issues. **Diagnose**: Check which marker has most values at -3.0. **Fix**: Add marker-specific z-score threshold to pre-filter.

### "Protein distribution peaked below zero / right-skewed"
Using gated 0-1 values instead of log-transformed fluorescence. **Diagnose**: Check `(protein_shared < 0).mean()` - if >55%, you're using gated values. **Fix**: Use `protein_adata.layers['log']` for normalization instead of `.X` (gated values). The log layer has lower skewness (0.30 vs 0.62) and better centering.

### "Bimodal initial pivot scores"
Often caused by scale mismatch between sparse RNA and dense protein. **Diagnose**: Check RNA sparsity - if >80% of values are at ZERO_VALUE, that's the issue. **Fix**:
1. Use log layer for protein normalization (reduces skewness)
2. Tune ZERO_VALUE: -1.0 works better than -3.0 for sparse RNA (mode separation 0.131 vs 0.176)
3. The smaller gap between undetected RNA and protein improves correlation-based matching

### "Conflicting match statistics between cells"
**CRITICAL**: MaxFuse has TWO matching stages with different statistics:

| Stage | What it measures | Typical count | Typical quality |
|-------|------------------|---------------|-----------------|
| **Pivot** | High-confidence anchor matches | ~1,500 RNA | ~80% "good" (score≥0.5) |
| **Propagated** | Extended matches via kNN | ~5,500 RNA | ~65% "good" |

**Example confusion:**
- Cell A (pivots): "1159/1419 RNA have good matches" (81.7%)
- Cell B (propagated): "60.5% RNA coverage" (5588/9236)

These are BOTH correct but measure DIFFERENT stages! The notebook now includes a comparison cell after propagation that shows both stages side-by-side.

**When comparing stats, always check:**
1. Is this pivot-level or propagated-level?
2. What's the denominator (pivot count vs total RNA)?
3. What threshold defines "good" (0.5 is standard)?
