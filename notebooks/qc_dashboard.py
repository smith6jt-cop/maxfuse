"""
Interactive QC Dashboard for scRNA-seq preprocessing.

This module provides a linked slider-based QC dashboard using Plotly
with proper log-scale support for all axes.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import ipywidgets as widgets
from IPython.display import display, clear_output


def get_renderer():
    """Auto-detect environment and set appropriate Plotly renderer."""
    try:
        from IPython import get_ipython
        shell = get_ipython().__class__.__name__
        if 'ZMQInteractiveShell' in shell:
            import os
            if 'VSCODE' in os.environ.get('TERM_PROGRAM', ''):
                return 'vscode'
            return 'notebook'
    except:
        pass
    return 'browser'


def log_histogram(data, n_bins=100, min_val=None):
    """
    Compute histogram with log-spaced bins for log-scale display.

    Parameters
    ----------
    data : array-like
        Data values (positive values only for log scale)
    n_bins : int
        Number of bins
    min_val : float, optional
        Minimum value for binning. If None, uses data minimum.
        Values below this are excluded.

    Returns
    -------
    bin_centers : ndarray
        Center of each bin (geometric mean of edges)
    counts : ndarray
        Count in each bin
    bin_edges : ndarray
        Bin edges (log-spaced)
    """
    data = np.asarray(data)

    # Filter to positive values for log scale
    if min_val is None:
        min_val = data[data > 0].min() if (data > 0).any() else 1

    data_valid = data[data >= min_val]
    if len(data_valid) == 0:
        return np.array([]), np.array([]), np.array([])

    max_val = data_valid.max()

    # Log-spaced bin edges
    bin_edges = np.logspace(np.log10(min_val), np.log10(max_val), n_bins + 1)
    counts, _ = np.histogram(data_valid, bins=bin_edges)

    # Geometric mean for bin centers (appropriate for log scale)
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])

    return bin_centers, counts, bin_edges


class QCDashboard:
    """
    Interactive QC dashboard with linked threshold sliders.

    Parameters
    ----------
    rna_adata : AnnData
        AnnData object with QC metrics in .obs:
        - 'total_counts': UMI counts per cell
        - 'n_genes_by_counts': genes detected per cell
        - 'pct_counts_mt': mitochondrial percentage
        - 'tissue' (optional): tissue type for coloring
    min_umi_threshold : int, default 5
        Minimum UMI count for a valid cell
    max_points : int, default 50000
        Maximum points to display in scatter plots (subsampled if exceeded)
    tissue_col : str, default 'tissue'
        Column name in .obs containing tissue type for coloring

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.read_h5ad('my_data.h5ad')
    >>> sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    >>> dashboard = QCDashboard(adata)
    >>> dashboard.display()
    """

    # Default color palette for tissue types
    TISSUE_COLORS = [
        '#1f77b4',  # blue
        '#ff7f0e',  # orange
        '#2ca02c',  # green
        '#d62728',  # red
        '#9467bd',  # purple
        '#8c564b',  # brown
        '#e377c2',  # pink
        '#7f7f7f',  # gray
        '#bcbd22',  # olive
        '#17becf',  # cyan
    ]

    def __init__(self, rna_adata, min_umi_threshold=5, max_points=50000, tissue_col='tissue'):
        self.obs = rna_adata.obs
        self.min_umi_threshold = min_umi_threshold
        self.max_points = max_points
        self.tissue_col = tissue_col

        # Extract data arrays
        self.total_counts = self.obs['total_counts'].values
        self.n_genes = self.obs['n_genes_by_counts'].values
        self.pct_mt = self.obs['pct_counts_mt'].values

        # Extract tissue info if available
        if tissue_col in self.obs.columns:
            self.tissue = self.obs[tissue_col].values
            self.tissue_types = sorted(self.obs[tissue_col].unique())
            self.tissue_color_map = {
                t: self.TISSUE_COLORS[i % len(self.TISSUE_COLORS)]
                for i, t in enumerate(self.tissue_types)
            }
            self.has_tissue = True
        else:
            self.tissue = None
            self.tissue_types = []
            self.tissue_color_map = {}
            self.has_tissue = False

        # Create valid data mask
        self.valid_mask = (
            ~np.isnan(self.total_counts) &
            ~np.isnan(self.n_genes) &
            ~np.isnan(self.pct_mt) &
            (self.total_counts >= min_umi_threshold)
        )

        # Calculate data ranges from valid data
        valid_counts = self.total_counts[self.valid_mask]
        valid_genes = self.n_genes[self.valid_mask]
        valid_mt = self.pct_mt[self.valid_mask]

        self.umi_data_min = int(np.floor(valid_counts.min()))
        self.umi_data_max = int(np.ceil(valid_counts.max()))
        self.genes_data_min = int(np.floor(valid_genes.min()))
        self.genes_data_max = int(np.ceil(valid_genes.max()))
        self.mt_data_min = 0.01  # Small positive value for log scale
        self.mt_data_max = float(np.ceil(valid_mt.max()))

        # Subsample for scatter plots
        self._create_scatter_mask()

        # Pre-compute histograms with log-spaced bins
        self._compute_histograms()

        # Initialize sliders and outputs
        self._create_widgets()

        # Set Plotly renderer
        pio.renderers.default = get_renderer()

    def _create_scatter_mask(self):
        """Create subsampling mask for scatter plots."""
        n_valid = self.valid_mask.sum()
        if n_valid > self.max_points:
            valid_indices = np.where(self.valid_mask)[0]
            subsample_idx = np.random.choice(valid_indices, self.max_points, replace=False)
            self.scatter_mask = np.zeros(len(self.obs), dtype=bool)
            self.scatter_mask[subsample_idx] = True
        else:
            self.scatter_mask = self.valid_mask.copy()

        # Pre-extract scatter data
        self.scatter_umi = self.total_counts[self.scatter_mask]
        self.scatter_mt = self.pct_mt[self.scatter_mask]
        self.scatter_genes = self.n_genes[self.scatter_mask]

        # Extract tissue for scatter points
        if self.has_tissue:
            self.scatter_tissue = self.tissue[self.scatter_mask]
        else:
            self.scatter_tissue = None

    def _compute_histograms(self):
        """Pre-compute histograms with log-spaced bins."""
        valid_umi = self.total_counts[self.valid_mask]
        valid_genes = self.n_genes[self.valid_mask]
        valid_mt = self.pct_mt[self.valid_mask]

        # MT% histogram - log-spaced bins, handle zeros
        mt_positive = valid_mt[valid_mt > 0]
        self.mt_hist_centers, self.mt_hist_counts, _ = log_histogram(
            mt_positive, n_bins=80, min_val=0.01
        )

        # UMI histogram - log-spaced bins
        self.umi_hist_centers, self.umi_hist_counts, _ = log_histogram(
            valid_umi, n_bins=80, min_val=max(1, self.umi_data_min)
        )

        # Genes histogram - log-spaced bins
        self.genes_hist_centers, self.genes_hist_counts, _ = log_histogram(
            valid_genes, n_bins=80, min_val=max(1, self.genes_data_min)
        )

    def _create_widgets(self):
        """Create slider widgets."""
        style = {'description_width': '80px'}
        layout = widgets.Layout(width='350px')

        self.mt_slider = widgets.FloatSlider(
            value=self.mt_data_max,
            min=0.1,
            max=self.mt_data_max,
            step=0.5,
            description='Max MT%:',
            readout_format='.1f',
            style=style,
            layout=layout
        )

        self.min_umi_slider = widgets.IntSlider(
            value=self.umi_data_min,
            min=self.umi_data_min,
            max=self.umi_data_max // 2,
            step=10,
            description='Min UMI:',
            style=style,
            layout=layout
        )

        self.max_umi_slider = widgets.IntSlider(
            value=self.umi_data_max,
            min=50,
            max=self.umi_data_max,
            step=50,
            description='Max UMI:',
            style=style,
            layout=layout
        )

        self.min_genes_slider = widgets.IntSlider(
            value=self.genes_data_min,
            min=self.genes_data_min,
            max=self.genes_data_max // 2,
            step=25,
            description='Min Genes:',
            style=style,
            layout=layout
        )

        self.max_genes_slider = widgets.IntSlider(
            value=self.genes_data_max,
            min=25,
            max=self.genes_data_max,
            step=25,
            description='Max Genes:',
            style=style,
            layout=layout
        )

        # Output widgets
        self.fig_output = widgets.Output()
        self.cell_count_output = widgets.HTML(value='')

        # Connect observers
        for slider in [self.mt_slider, self.min_umi_slider, self.max_umi_slider,
                       self.min_genes_slider, self.max_genes_slider]:
            slider.observe(self._update_dashboard, names='value')

    def _create_figure(self, mt_max, min_umi, max_umi, min_genes, max_genes):
        """Create the QC dashboard figure with current threshold values."""
        # Dynamic subplot titles based on whether tissue info is available
        genes_title = 'Genes vs Counts (by tissue)' if self.has_tissue else 'Genes vs Counts (colored by MT%)'
        mt_umi_title = 'MT% vs UMI (by tissue)' if self.has_tissue else 'MT% vs UMI Counts'

        fig = make_subplots(
            rows=2, cols=3,
            subplot_titles=(
                'MT% Distribution',
                mt_umi_title,
                genes_title,
                'UMI Distribution',
                'Genes Distribution',
                'Cells Passing Filters'
            ),
            specs=[[{}, {}, {}], [{}, {}, {"type": "domain"}]],
            vertical_spacing=0.15,
            horizontal_spacing=0.08
        )

        # 1. MT% histogram (log-log scale)
        fig.add_trace(go.Scatter(
            x=self.mt_hist_centers,
            y=self.mt_hist_counts,
            mode='lines',
            fill='tozeroy',
            fillcolor='rgba(70, 130, 180, 0.4)',
            line=dict(color='steelblue', width=1),
            showlegend=False
        ), row=1, col=1)

        # 2. MT% vs UMI scatter (log-log) - colored by tissue
        # Filter scatter points for positive MT% (required for log scale)
        mt_positive_mask = self.scatter_mt > 0

        if self.has_tissue:
            # Add separate trace per tissue type for legend
            for tissue_type in self.tissue_types:
                tissue_mask = mt_positive_mask & (self.scatter_tissue == tissue_type)
                fig.add_trace(go.Scattergl(
                    x=self.scatter_umi[tissue_mask],
                    y=self.scatter_mt[tissue_mask],
                    mode='markers',
                    marker=dict(size=2, opacity=0.4, color=self.tissue_color_map[tissue_type]),
                    name=tissue_type,
                    legendgroup=tissue_type,
                    showlegend=True,
                    hovertemplate=f'{tissue_type}<br>UMI: %{{x}}<br>MT%: %{{y:.1f}}<extra></extra>'
                ), row=1, col=2)
        else:
            fig.add_trace(go.Scattergl(
                x=self.scatter_umi[mt_positive_mask],
                y=self.scatter_mt[mt_positive_mask],
                mode='markers',
                marker=dict(size=2, opacity=0.3, color='steelblue'),
                showlegend=False,
                hovertemplate='UMI: %{x}<br>MT%%: %{y:.1f}<extra></extra>'
            ), row=1, col=2)

        # 3. Genes vs Counts (log-log) - colored by tissue
        if self.has_tissue:
            for tissue_type in self.tissue_types:
                tissue_mask = self.scatter_tissue == tissue_type
                fig.add_trace(go.Scattergl(
                    x=self.scatter_umi[tissue_mask],
                    y=self.scatter_genes[tissue_mask],
                    mode='markers',
                    marker=dict(size=2, opacity=0.4, color=self.tissue_color_map[tissue_type]),
                    name=tissue_type,
                    legendgroup=tissue_type,
                    showlegend=False,  # Legend already shown from plot 2
                    hovertemplate=f'{tissue_type}<br>UMI: %{{x}}<br>Genes: %{{y}}<extra></extra>'
                ), row=1, col=3)
        else:
            fig.add_trace(go.Scattergl(
                x=self.scatter_umi,
                y=self.scatter_genes,
                mode='markers',
                marker=dict(
                    size=2, opacity=0.5,
                    color=self.scatter_mt,
                    colorscale='RdYlBu_r',
                    cmin=0, cmax=min(80, self.mt_data_max),
                    colorbar=dict(title='MT%', x=1.02, len=0.4, y=0.8)
                ),
                showlegend=False,
                hovertemplate='UMI: %{x}<br>Genes: %{y}<br>MT%%: %{marker.color:.1f}<extra></extra>'
            ), row=1, col=3)

        # 4. UMI histogram (log-log scale)
        fig.add_trace(go.Scatter(
            x=self.umi_hist_centers,
            y=self.umi_hist_counts,
            mode='lines',
            fill='tozeroy',
            fillcolor='rgba(70, 130, 180, 0.4)',
            line=dict(color='steelblue', width=1),
            showlegend=False
        ), row=2, col=1)

        # 5. Genes histogram (log-log scale)
        fig.add_trace(go.Scatter(
            x=self.genes_hist_centers,
            y=self.genes_hist_counts,
            mode='lines',
            fill='tozeroy',
            fillcolor='rgba(70, 130, 180, 0.4)',
            line=dict(color='steelblue', width=1),
            showlegend=False
        ), row=2, col=2)

        # 6. Calculate passing cells and add pie chart
        mask = (
            self.valid_mask &
            (self.total_counts >= min_umi) &
            (self.total_counts <= max_umi) &
            (self.n_genes >= min_genes) &
            (self.n_genes <= max_genes) &
            (self.pct_mt <= mt_max)
        )
        n_pass = mask.sum()
        n_fail = self.valid_mask.sum() - n_pass

        fig.add_trace(go.Pie(
            values=[n_pass, n_fail],
            labels=['Pass', 'Fail'],
            marker_colors=['#90EE90', '#FFB6C1'],
            hole=0.4,
            showlegend=False,
            textinfo='percent+label'
        ), row=2, col=3)

        # Add threshold lines using annotations (works reliably with log scale)
        shapes = self._create_threshold_shapes(mt_max, min_umi, max_umi, min_genes, max_genes)
        fig.update_layout(shapes=shapes)

        # Update axes to log scale
        # Row 1, Col 1: MT% histogram
        fig.update_xaxes(type='log', title_text='MT%', row=1, col=1)
        fig.update_yaxes(type='log', title_text='Count', row=1, col=1)

        # Row 1, Col 2: MT% vs UMI
        fig.update_xaxes(type='log', title_text='UMI Counts', row=1, col=2)
        fig.update_yaxes(type='log', title_text='MT%', row=1, col=2)

        # Row 1, Col 3: Genes vs Counts
        fig.update_xaxes(type='log', title_text='UMI Counts', row=1, col=3)
        fig.update_yaxes(type='log', title_text='Genes', row=1, col=3)

        # Row 2, Col 1: UMI histogram
        fig.update_xaxes(type='log', title_text='UMI Counts', row=2, col=1)
        fig.update_yaxes(type='log', title_text='Count', row=2, col=1)

        # Row 2, Col 2: Genes histogram
        fig.update_xaxes(type='log', title_text='Genes', row=2, col=2)
        fig.update_yaxes(type='log', title_text='Count', row=2, col=2)

        layout_kwargs = dict(
            height=600,
            width=1400,
            title_text='Interactive QC Dashboard - Drag sliders to adjust thresholds',
            title_x=0.5
        )

        # Add legend configuration if tissue coloring is used
        if self.has_tissue:
            layout_kwargs['legend'] = dict(
                orientation='h',
                yanchor='bottom',
                y=1.02,
                xanchor='center',
                x=0.5,
                title='Tissue'
            )

        fig.update_layout(**layout_kwargs)

        return fig, n_pass

    def _create_threshold_shapes(self, mt_max, min_umi, max_umi, min_genes, max_genes):
        """
        Create threshold line shapes that work correctly with log scale.

        For log-scale axes, we use data coordinates (xref='x', yref='y')
        with appropriate axis ranges.
        """
        shapes = []

        # Get axis ranges for domain-based positioning
        # These are approximate but work for line positioning
        mt_y_max = self.mt_hist_counts.max() * 2 if len(self.mt_hist_counts) > 0 else 1000
        umi_y_max = self.umi_hist_counts.max() * 2 if len(self.umi_hist_counts) > 0 else 1000
        genes_y_max = self.genes_hist_counts.max() * 2 if len(self.genes_hist_counts) > 0 else 1000

        # MT% threshold lines
        # Row 1, Col 1: MT% histogram - vertical line at mt_max
        shapes.append(dict(
            type='line',
            x0=mt_max, x1=mt_max,
            y0=0.5, y1=mt_y_max,
            xref='x', yref='y',
            line=dict(color='red', dash='dash', width=2)
        ))

        # Row 1, Col 2: MT% vs UMI - horizontal line at mt_max
        shapes.append(dict(
            type='line',
            x0=self.umi_data_min, x1=self.umi_data_max,
            y0=mt_max, y1=mt_max,
            xref='x2', yref='y2',
            line=dict(color='red', dash='dash', width=2)
        ))

        # Min UMI lines (green)
        # Row 1, Col 2: MT% vs UMI - vertical
        shapes.append(dict(
            type='line',
            x0=min_umi, x1=min_umi,
            y0=0.01, y1=self.mt_data_max,
            xref='x2', yref='y2',
            line=dict(color='green', dash='dash', width=2)
        ))
        # Row 1, Col 3: Genes vs UMI - vertical
        shapes.append(dict(
            type='line',
            x0=min_umi, x1=min_umi,
            y0=self.genes_data_min, y1=self.genes_data_max,
            xref='x3', yref='y3',
            line=dict(color='green', dash='dash', width=2)
        ))
        # Row 2, Col 1: UMI histogram - vertical
        shapes.append(dict(
            type='line',
            x0=min_umi, x1=min_umi,
            y0=0.5, y1=umi_y_max,
            xref='x4', yref='y4',
            line=dict(color='green', dash='dash', width=2)
        ))

        # Max UMI lines (orange)
        # Row 1, Col 2: MT% vs UMI - vertical
        shapes.append(dict(
            type='line',
            x0=max_umi, x1=max_umi,
            y0=0.01, y1=self.mt_data_max,
            xref='x2', yref='y2',
            line=dict(color='orange', dash='dash', width=2)
        ))
        # Row 1, Col 3: Genes vs UMI - vertical
        shapes.append(dict(
            type='line',
            x0=max_umi, x1=max_umi,
            y0=self.genes_data_min, y1=self.genes_data_max,
            xref='x3', yref='y3',
            line=dict(color='orange', dash='dash', width=2)
        ))
        # Row 2, Col 1: UMI histogram - vertical
        shapes.append(dict(
            type='line',
            x0=max_umi, x1=max_umi,
            y0=0.5, y1=umi_y_max,
            xref='x4', yref='y4',
            line=dict(color='orange', dash='dash', width=2)
        ))

        # Min genes lines (purple)
        # Row 1, Col 3: Genes vs UMI - horizontal
        shapes.append(dict(
            type='line',
            x0=self.umi_data_min, x1=self.umi_data_max,
            y0=min_genes, y1=min_genes,
            xref='x3', yref='y3',
            line=dict(color='purple', dash='dash', width=2)
        ))
        # Row 2, Col 2: Genes histogram - vertical
        shapes.append(dict(
            type='line',
            x0=min_genes, x1=min_genes,
            y0=0.5, y1=genes_y_max,
            xref='x5', yref='y5',
            line=dict(color='purple', dash='dash', width=2)
        ))

        # Max genes lines (brown)
        # Row 1, Col 3: Genes vs UMI - horizontal
        shapes.append(dict(
            type='line',
            x0=self.umi_data_min, x1=self.umi_data_max,
            y0=max_genes, y1=max_genes,
            xref='x3', yref='y3',
            line=dict(color='saddlebrown', dash='dash', width=2)
        ))
        # Row 2, Col 2: Genes histogram - vertical
        shapes.append(dict(
            type='line',
            x0=max_genes, x1=max_genes,
            y0=0.5, y1=genes_y_max,
            xref='x5', yref='y5',
            line=dict(color='saddlebrown', dash='dash', width=2)
        ))

        return shapes

    def _update_dashboard(self, change=None):
        """Update the dashboard when sliders change."""
        mt_max = self.mt_slider.value
        min_umi = self.min_umi_slider.value
        max_umi = self.max_umi_slider.value
        min_genes = self.min_genes_slider.value
        max_genes = self.max_genes_slider.value

        fig, n_pass = self._create_figure(mt_max, min_umi, max_umi, min_genes, max_genes)
        pct = 100 * n_pass / self.valid_mask.sum()

        # Update cell count display
        self.cell_count_output.value = f"""
        <div style="font-size: 14px; padding: 10px; background: linear-gradient(to right, #90EE90 {pct}%, #FFB6C1 {pct}%);
                    border-radius: 5px; text-align: center;">
            <b>Cells passing:</b> {n_pass:,} / {self.valid_mask.sum():,} ({pct:.1f}%)
        </div>
        """

        # Update figure
        with self.fig_output:
            clear_output(wait=True)
            fig.show()

    def display(self):
        """Display the interactive dashboard."""
        print(f"Total barcodes: {len(self.obs):,}")
        print(f"Valid cells (no NaN, >= {self.min_umi_threshold} UMI): {self.valid_mask.sum():,}")
        print(f"Excluded (NaN or low UMI): {(~self.valid_mask).sum():,}")
        print(f"\nData ranges:")
        print(f"  UMI: {self.umi_data_min:,} - {self.umi_data_max:,}")
        print(f"  Genes: {self.genes_data_min:,} - {self.genes_data_max:,}")
        print(f"  MT%: {self.mt_data_min:.2f} - {self.mt_data_max:.1f}")

        if self.has_tissue:
            tissue_counts = {t: (self.tissue == t).sum() for t in self.tissue_types}
            print(f"\nTissue types ({len(self.tissue_types)}):")
            for t, count in tissue_counts.items():
                print(f"  {t}: {count:,} cells")
        else:
            print(f"\nNo '{self.tissue_col}' column found - scatters colored uniformly/by MT%")

        # Create layout
        slider_box = widgets.VBox([
            widgets.HTML('<h4 style="margin:5px 0">Filter Thresholds</h4>'),
            widgets.HBox([self.min_umi_slider, self.max_umi_slider]),
            widgets.HBox([self.min_genes_slider, self.max_genes_slider]),
            widgets.HBox([self.mt_slider]),
            self.cell_count_output
        ], layout=widgets.Layout(padding='10px', border='1px solid #ddd', border_radius='5px'))

        display(slider_box)
        display(self.fig_output)

        # Initial render
        self._update_dashboard()

        print("""
Threshold Legend:
  - Red (dashed): Max MT%
  - Green (dashed): Min UMI
  - Orange (dashed): Max UMI
  - Purple (dashed): Min Genes
  - Brown (dashed): Max Genes
""")

    def get_thresholds(self):
        """
        Return current threshold values as a dictionary.

        Useful for programmatically applying the selected filters.

        Returns
        -------
        dict
            Dictionary with keys: 'min_umi', 'max_umi', 'min_genes',
            'max_genes', 'max_mt'
        """
        return {
            'min_umi': self.min_umi_slider.value,
            'max_umi': self.max_umi_slider.value,
            'min_genes': self.min_genes_slider.value,
            'max_genes': self.max_genes_slider.value,
            'max_mt': self.mt_slider.value
        }

    def get_passing_mask(self):
        """
        Return boolean mask for cells passing current filters.

        Returns
        -------
        ndarray
            Boolean array, True for cells passing all current threshold filters.
        """
        t = self.get_thresholds()
        return (
            self.valid_mask &
            (self.total_counts >= t['min_umi']) &
            (self.total_counts <= t['max_umi']) &
            (self.n_genes >= t['min_genes']) &
            (self.n_genes <= t['max_genes']) &
            (self.pct_mt <= t['max_mt'])
        )
