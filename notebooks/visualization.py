# Visualize DETECTION-AWARE normalization
# Skip if raw CODEX data not available

if not CODEX_RAW_AVAILABLE:
    print("Skipping: Raw CODEX data not available (run from preprocessing to analyze)")
else:
    # Visualize DETECTION-AWARE normalization
    # Updated to use Cell Median-based detection for protein data
    
    fig, axes = plt.subplots(4, 4, figsize=(20, 16))
    
    feature_names = list(rna_protein_correspondence[:, 0])
    protein_names = list(rna_protein_correspondence[:, 1])
    n_features = rna_shared_raw.shape[1]
    
    # Define x_pos and width for bar charts
    x_pos = np.arange(n_features)
    width = 0.35
    
    # Pick a feature with moderate expression
    detection_rates_rna = [(rna_shared_raw[:, i] > 0).mean() for i in range(n_features)]
    good_features = [i for i, d in enumerate(detection_rates_rna) if 0.2 < d < 0.8]
    if good_features:
        best_feat_idx = good_features[len(good_features)//2]
    else:
        best_feat_idx = np.argmax(detection_rates_rna)
    feat_name = feature_names[best_feat_idx]
    prot_name = protein_names[best_feat_idx]
    
    # For protein, use Cell Median-based detection (from cell 26)
    # Build protein detection rates from Cell Median
    protein_to_median_col = {}
    for prot in protein_names:
        for col in codex_df.columns:
            if prot in col and 'Cell:' in col and 'Median' in col:
                protein_to_median_col[prot] = col
                break
    
    # Calculate detection rates: RNA uses >0, Protein uses Cell Median > 0.5
    rna_det_rates = [(rna_shared_raw[:, i] > 0).mean() * 100 for i in range(n_features)]
    prot_det_rates = []
    for i, pname in enumerate(protein_names):
        if pname in protein_to_median_col:
            median_vals = codex_df[protein_to_median_col[pname]].values
            prot_det_rates.append((median_vals > 0.5).mean() * 100)
        else:
            # Fallback to raw >0
            prot_det_rates.append((protein_shared_raw[:, i] > 0).mean() * 100)
    
    rna_det = rna_det_rates[best_feat_idx]
    prot_det = prot_det_rates[best_feat_idx]
    
    print(f"Example feature: {feat_name} / {prot_name}")
    print(f"  RNA detection: {rna_det:.1f}% ({int(rna_det/100 * rna_shared_raw.shape[0]):,} cells)")
    print(f"  Protein detection (Cell Median > 0): {prot_det:.1f}% ({int(prot_det/100 * protein_shared_raw.shape[0]):,} cells)")
    
    ZERO_VALUE = -2.5  # Must match value used in normalization
    
    # ============================================================
    # Row 1: RNA transformation pipeline
    # ============================================================
    ax = axes[0, 0]
    raw_vals = rna_shared_raw[:, best_feat_idx]
    ax.hist(raw_vals, bins=50, alpha=0.7, color='steelblue', edgecolor='white')
    ax.set_title(f'RNA: Raw Counts\n(zeros: {(raw_vals==0).mean()*100:.0f}%)', fontsize=10)
    ax.set_xlabel('Count')
    ax.set_ylabel('Cells')
    
    ax = axes[0, 1]
    log_vals = rna_shared_after_log[:, best_feat_idx]
    zeros = log_vals == 0
    nonzeros = ~zeros
    ax.hist(log_vals[nonzeros], bins=30, alpha=0.7, color='coral', edgecolor='white', label=f'Non-zero ({nonzeros.sum():,})')
    if zeros.sum() > 0:
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=2, label=f'Zeros ({zeros.sum():,})')
    ax.set_title('RNA: After log1p\n(zeros at 0)', fontsize=10)
    ax.set_xlabel('log1p(count)')
    ax.legend(fontsize=8)
    
    ax = axes[0, 2]
    norm_vals = rna_shared_after_scale[:, best_feat_idx]
    nonzero_mask = norm_vals > ZERO_VALUE + 0.1
    ax.hist(norm_vals[nonzero_mask], bins=30, alpha=0.7, color='forestgreen', edgecolor='white', 
            label=f'Detected ({nonzero_mask.sum():,})')
    ax.axvline(x=ZERO_VALUE, color='red', linestyle='--', linewidth=2, 
               label=f'Not detected ({(~nonzero_mask).sum():,})')
    ax.set_title('RNA: Detection-Aware Normalized\n(zeros → fixed value)', fontsize=10)
    ax.set_xlabel('Normalized value')
    ax.legend(fontsize=8)
    
    ax = axes[0, 3]
    ax.axis('off')
    ax.text(0.1, 0.85, 'RNA Pipeline:', fontsize=12, fontweight='bold', transform=ax.transAxes)
    ax.text(0.1, 0.70, '1. normalize_total (library size)', fontsize=10, transform=ax.transAxes)
    ax.text(0.1, 0.58, '2. log1p (variance stabilization)', fontsize=10, transform=ax.transAxes)
    ax.text(0.1, 0.46, '3. Detection-aware normalization:', fontsize=10, transform=ax.transAxes)
    ax.text(0.15, 0.34, '• Non-zeros → rank → normal quantiles', fontsize=9, transform=ax.transAxes)
    ax.text(0.15, 0.22, f'• Zeros → fixed value ({ZERO_VALUE})', fontsize=9, transform=ax.transAxes)
    
    # ============================================================
    # Row 2: Protein transformation pipeline
    # ============================================================
    ax = axes[1, 0]
    raw_vals = protein_shared_raw[:, best_feat_idx]
    # Show Cell Median-based detection
    if prot_name in protein_to_median_col:
        median_vals = codex_df[protein_to_median_col[prot_name]].values
        is_bg = median_vals <= 0.5
        pct_bg = is_bg.mean() * 100
    else:
        pct_bg = (raw_vals == 0).mean() * 100
    ax.hist(raw_vals, bins=50, alpha=0.7, color='darkorange', edgecolor='white')
    ax.set_title(f'Protein: Raw MFI\n(background: {pct_bg:.0f}% via Cell Median)', fontsize=10)
    ax.set_xlabel('Mean Fluorescence Intensity')
    ax.set_ylabel('Cells')
    
    ax = axes[1, 1]
    arcsinh_vals = protein_shared_after_arcsinh[:, best_feat_idx]
    # Use Cell Median for zero classification
    if prot_name in protein_to_median_col:
        median_vals = codex_df[protein_to_median_col[prot_name]].values
        zeros = median_vals <= 0.5
    else:
        zeros = protein_shared_raw[:, best_feat_idx] == 0
    nonzeros = ~zeros
    ax.hist(arcsinh_vals[nonzeros], bins=30, alpha=0.7, color='purple', edgecolor='white', label=f'Detected ({nonzeros.sum():,})')
    if zeros.sum() > 0:
        # Show background cells as a separate histogram
        ax.hist(arcsinh_vals[zeros], bins=30, alpha=0.5, color='gray', edgecolor='white', label=f'Background ({zeros.sum():,})')
    ax.set_title('Protein: After arcsinh(x/5)\n(Cell Median-based detection)', fontsize=10)
    ax.set_xlabel('arcsinh(MFI/5)')
    ax.legend(fontsize=8)
    
    ax = axes[1, 2]
    norm_vals = protein_shared_after[:, best_feat_idx]
    nonzero_mask = norm_vals > ZERO_VALUE + 0.1
    ax.hist(norm_vals[nonzero_mask], bins=30, alpha=0.7, color='forestgreen', edgecolor='white',
            label=f'Detected ({nonzero_mask.sum():,})')
    ax.axvline(x=ZERO_VALUE, color='red', linestyle='--', linewidth=2,
               label=f'Not detected ({(~nonzero_mask).sum():,})')
    ax.set_title('Protein: Detection-Aware Normalized\n(background → fixed value)', fontsize=10)
    ax.set_xlabel('Normalized value')
    ax.legend(fontsize=8)
    
    ax = axes[1, 3]
    ax.axis('off')
    ax.text(0.1, 0.85, 'Protein Pipeline:', fontsize=12, fontweight='bold', transform=ax.transAxes)
    ax.text(0.1, 0.70, '1. arcsinh(x/5) (standard for cytometry)', fontsize=10, transform=ax.transAxes)
    ax.text(0.15, 0.58, '• Linear near zero, log-like at high', fontsize=9, transform=ax.transAxes)
    ax.text(0.1, 0.46, '2. Detection via Cell Median:', fontsize=10, transform=ax.transAxes)
    ax.text(0.15, 0.34, '• Median > 0 → detected', fontsize=9, transform=ax.transAxes)
    ax.text(0.15, 0.22, '• Median = 0 → background', fontsize=9, transform=ax.transAxes)
    
    # ============================================================
    # Row 3: Distribution comparison - NON-ZERO VALUES ONLY
    # ============================================================
    
    # Get non-zero values only
    rna_nonzero = rna_shared_after_scale[rna_shared_after_scale > ZERO_VALUE + 0.1]
    prot_nonzero = protein_shared_after[protein_shared_after > ZERO_VALUE + 0.1]
    
    ax = axes[2, 0]
    bins = np.linspace(-3, 3, 50)
    ax.hist(rna_nonzero, bins=bins, alpha=0.6, density=True, label=f'RNA ({len(rna_nonzero):,})', color='steelblue')
    ax.hist(prot_nonzero, bins=bins, alpha=0.6, density=True, label=f'Protein ({len(prot_nonzero):,})', color='darkorange')
    ax.set_title('Non-Zero Values Only\n(should overlap well)', fontsize=10)
    ax.set_xlabel('Normalized value')
    ax.set_ylabel('Density')
    ax.legend(fontsize=8)
    
    # Box plot of non-zero values
    ax = axes[2, 1]
    rna_sample = rna_nonzero[::max(1, len(rna_nonzero)//5000)]
    prot_sample = prot_nonzero[::max(1, len(prot_nonzero)//5000)]
    bp = ax.boxplot([rna_sample, prot_sample], labels=['RNA', 'Protein'], patch_artist=True)
    bp['boxes'][0].set_facecolor('steelblue')
    bp['boxes'][1].set_facecolor('darkorange')
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
    ax.set_ylabel('Normalized value (non-zero only)')
    ax.set_title('Overall Distribution (Non-Zero Values)', fontsize=10)
    
    # Per-feature mean of NON-ZERO values (should be ~0)
    ax = axes[2, 2]
    rna_means_nz = []
    prot_means_nz = []
    for i in range(n_features):
        rna_nz = rna_shared_after_scale[:, i][rna_shared_after_scale[:, i] > ZERO_VALUE + 0.1]
        prot_nz = protein_shared_after[:, i][protein_shared_after[:, i] > ZERO_VALUE + 0.1]
        rna_means_nz.append(rna_nz.mean() if len(rna_nz) > 0 else 0)
        prot_means_nz.append(prot_nz.mean() if len(prot_nz) > 0 else 0)
    
    ax.bar(x_pos - width/2, rna_means_nz, width, label='RNA', alpha=0.8, color='steelblue')
    ax.bar(x_pos + width/2, prot_means_nz, width, label='Protein', alpha=0.8, color='darkorange')
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Feature index')
    ax.set_ylabel('Mean (non-zero only)')
    ax.set_title('Per-Feature Mean (Non-Zero Values)', fontsize=10)
    ax.legend(fontsize=8)
    ax.set_xticks(x_pos)
    
    # Per-feature std of NON-ZERO values (should be ~1)
    ax = axes[2, 3]
    rna_stds_nz = []
    prot_stds_nz = []
    for i in range(n_features):
        rna_nz = rna_shared_after_scale[:, i][rna_shared_after_scale[:, i] > ZERO_VALUE + 0.1]
        prot_nz = protein_shared_after[:, i][protein_shared_after[:, i] > ZERO_VALUE + 0.1]
        rna_stds_nz.append(rna_nz.std() if len(rna_nz) > 1 else 0)
        prot_stds_nz.append(prot_nz.std() if len(prot_nz) > 1 else 0)
    
    ax.bar(x_pos - width/2, rna_stds_nz, width, label='RNA', alpha=0.8, color='steelblue')
    ax.bar(x_pos + width/2, prot_stds_nz, width, label='Protein', alpha=0.8, color='darkorange')
    ax.axhline(y=1, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Feature index')
    ax.set_ylabel('Std Dev (non-zero only)')
    ax.set_title('Per-Feature Std Dev (Non-Zero Values)', fontsize=10)
    ax.legend(fontsize=8)
    ax.set_xticks(x_pos)
    
    # ============================================================
    # Row 4: Detection rates and summary
    # ============================================================
    
    # Detection rate comparison
    ax = axes[3, 0]
    ax.bar(x_pos - width/2, rna_det_rates, width, label='RNA', alpha=0.8, color='steelblue')
    ax.bar(x_pos + width/2, prot_det_rates, width, label='Protein', alpha=0.8, color='darkorange')
    ax.set_xlabel('Feature index')
    ax.set_ylabel('Detection rate (%)')
    ax.set_title('Per-Feature Detection Rate', fontsize=10)
    ax.legend(fontsize=8)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f[:6] for f in feature_names], rotation=45, ha='right', fontsize=8)
    
    # Detection rate scatter
    ax = axes[3, 1]
    ax.scatter(rna_det_rates, prot_det_rates, s=50, alpha=0.7)
    for i, fname in enumerate(feature_names):
        ax.annotate(fname[:6], (rna_det_rates[i], prot_det_rates[i]), fontsize=7)
    ax.plot([0, 100], [0, 100], 'r--', alpha=0.5)
    ax.set_xlabel('RNA detection rate (%)')
    ax.set_ylabel('Protein detection rate (%)')
    ax.set_title('Detection Rate Comparison', fontsize=10)
    
    # Empty placeholder
    ax = axes[3, 2]
    ax.axis('off')
    
    # Summary text
    ax = axes[3, 3]
    ax.axis('off')
    summary = f"""DETECTION-AWARE NORMALIZATION SUMMARY
    {"="*45}
    
    Cell counts:
      RNA:     {rna_shared_raw.shape[0]:>10,} cells
      Protein: {protein_shared_raw.shape[0]:>10,} cells
    
    Detection method:
      RNA: count > 0
      Protein: Cell Median > 0 (histogram matching)
    
    Non-zero values after normalization:
      RNA:     mean={np.mean(rna_means_nz):.3f}, std={np.mean(rna_stds_nz):.3f}
      Protein: mean={np.mean(prot_means_nz):.3f}, std={np.mean(prot_stds_nz):.3f}
    
    Zero handling:
      All "not detected" set to: {ZERO_VALUE}
      (Below all detected values)
    """
    ax.text(0.0, 0.95, summary, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plt.suptitle('Detection-Aware Normalization Results', fontsize=14, y=1.01)
    plt.show()
    
    print("\n" + "="*70)
    print("KEY INSIGHT:")
    print("="*70)
    print("RNA detection: count > 0")
    print("Protein detection: Cell Median > 0 (accounts for histogram matching)")
    print("This ensures 'not detected' is biologically meaningful in both modalities.")