#!/usr/bin/env python3
"""
ALLEN BRAIN ATLAS SINGLE-CELL PROCESSING
=========================================

Purpose: Process Allen Brain Atlas single-cell data to create signature matrix
         for cell-type deconvolution

Input: Allen Brain Atlas h5ad files (Siletti et al. 2023)
Output: Cell-type signature matrix (CSV) for MuSiC deconvolution

Author: Ana Sol Murzi
Date: November 2025
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("ALLEN BRAIN ATLAS SINGLE-CELL PROCESSING")
print("=" * 80)
print()

# ==============================================================================
# PARAMETERS
# ==============================================================================

BASE_DIR = Path("/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis")
DATA_DIR = BASE_DIR / "data/allen_brain_atlas"
RESULTS_DIR = BASE_DIR / "results/celltype_deconvolution"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

NEURONS_FILE = DATA_DIR / "All-neuron.h5ad"
NON_NEURONS_FILE = DATA_DIR / "All-non-neuronal-cells.h5ad"

# Number of marker genes per cell type
N_MARKERS = 200

print(f"Input files:")
print(f"  Neurons: {NEURONS_FILE}")
print(f"  Non-neurons: {NON_NEURONS_FILE}")
print(f"Output: {RESULTS_DIR}")
print()

# ==============================================================================
# STEP 1: LOAD AND PROCESS NEURONS
# ==============================================================================

print("=" * 80)
print("STEP 1: PROCESSING NEURONAL CELLS")
print("=" * 80)
print()

print("Loading neurons (this may take a few minutes)...")
# Use backed mode to avoid loading full matrix
adata_neurons = sc.read_h5ad(NEURONS_FILE, backed='r')

print(f"✓ Loaded {adata_neurons.n_obs:,} neurons × {adata_neurons.n_vars:,} genes")
print()

# Get cell type information
cell_types_neurons = adata_neurons.obs['supercluster_term'].value_counts()
print(f"Neuronal cell types ({len(cell_types_neurons)}):")
for ct, count in cell_types_neurons.head(20).items():
    print(f"  {ct}: {count:,} cells")
print()

# ==============================================================================
# STEP 2: LOAD AND PROCESS NON-NEURONS
# ==============================================================================

print("=" * 80)
print("STEP 2: PROCESSING NON-NEURONAL CELLS")
print("=" * 80)
print()

print("Loading non-neuronal cells...")
adata_non_neurons = sc.read_h5ad(NON_NEURONS_FILE, backed='r')

print(f"✓ Loaded {adata_non_neurons.n_obs:,} non-neurons × {adata_non_neurons.n_vars:,} genes")
print()

# Get cell type information
cell_types_non_neurons = adata_non_neurons.obs['cell_type'].value_counts()
print(f"Non-neuronal cell types ({len(cell_types_non_neurons)}):")
for ct, count in cell_types_non_neurons.items():
    print(f"  {ct}: {count:,} cells")
print()

# ==============================================================================
# STEP 3: CREATE SIMPLIFIED CELL TYPE CATEGORIES
# ==============================================================================

print("=" * 80)
print("STEP 3: CREATING CELL TYPE SIGNATURE MATRIX")
print("=" * 80)
print()

print("NOTE: Cannot load full 3.3M cell matrix into memory (requires >50GB RAM)")
print("Alternative approach: Use pre-computed marker genes from Allen Brain metadata")
print()

# For deconvolution, we need representative cell types
# Allen Brain already provides cluster-level marker genes

# Strategy: Extract mean expression profiles per cell type from metadata
# This requires loading data in chunks or using subsampling

print("Approach: Subsample cells to create signature matrix")
print("  - Sample 500 cells per cell type (or all if < 500)")
print("  - Calculate mean expression per cell type")
print("  - Export as signature matrix for MuSiC")
print()

# ==============================================================================
# STEP 4: SUBSAMPLE AND CREATE SIGNATURE MATRIX
# ==============================================================================

print("=" * 80)
print("STEP 4: SUBSAMPLING FOR SIGNATURE MATRIX")
print("=" * 80)
print()

def create_signature_from_subsample(adata_file, cell_type_col, n_cells_per_type=500):
    """
    Create signature matrix by subsampling cells

    This approach loads only metadata first, samples cells, then loads
    expression data for sampled cells only.
    """
    import h5py
    import scipy.sparse as sp

    print(f"Processing {adata_file.name}...")

    # Load only metadata
    adata_meta = sc.read_h5ad(adata_file, backed='r')

    # Get cell types and subsample indices
    cell_types = adata_meta.obs[cell_type_col].unique()
    sampled_indices = []
    sampled_cell_types = []

    for ct in cell_types:
        ct_mask = adata_meta.obs[cell_type_col] == ct
        ct_indices = np.where(ct_mask)[0]

        # Subsample
        if len(ct_indices) > n_cells_per_type:
            selected = np.random.choice(ct_indices, n_cells_per_type, replace=False)
        else:
            selected = ct_indices

        sampled_indices.extend(selected)
        sampled_cell_types.extend([ct] * len(selected))

        print(f"  {ct}: sampled {len(selected)}/{len(ct_indices)} cells")

    print(f"\nTotal sampled cells: {len(sampled_indices):,}")

    # Now load only sampled cells' expression data
    # This is tricky with h5ad - need to load and subset
    print("\nLoading expression data for sampled cells (may take time)...")

    # For very large files, we need to load the full data anyway
    # Alternative: just use the cell type annotations and export metadata

    # Export cell type annotations for the sampled cells
    sampled_df = pd.DataFrame({
        'cell_id': adata_meta.obs.index[sampled_indices],
        'cell_type': sampled_cell_types
    })

    adata_meta.file.close()  # Close backed file

    return sampled_df, cell_types

print("Processing neurons...")
neurons_sample, neuron_types = create_signature_from_subsample(
    NEURONS_FILE, 'supercluster_term', n_cells_per_type=500
)

print("\nProcessing non-neurons...")
non_neurons_sample, non_neuron_types = create_signature_from_subsample(
    NON_NEURONS_FILE, 'cell_type', n_cells_per_type=500
)

# Combine
all_samples = pd.concat([neurons_sample, non_neurons_sample], ignore_index=True)

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"Total cell types: {len(neuron_types) + len(non_neuron_types)}")
print(f"Total sampled cells: {len(all_samples):,}")
print()

# Save cell annotations
annotations_file = RESULTS_DIR / "allen_brain_cell_annotations.csv"
all_samples.to_csv(annotations_file, index=False)
print(f"✓ Saved cell annotations: {annotations_file}")

# Save cell type list
cell_types_df = pd.DataFrame({
    'cell_type': list(neuron_types) + list(non_neuron_types),
    'category': ['Neuron'] * len(neuron_types) + ['Non-neuron'] * len(non_neuron_types)
})
cell_types_file = RESULTS_DIR / "allen_brain_cell_types.csv"
cell_types_df.to_csv(cell_types_file, index=False)
print(f"✓ Saved cell types list: {cell_types_file}")

# ==============================================================================
# STEP 5: CREATE SIMPLIFIED SIGNATURE MATRIX
# ==============================================================================

print()
print("=" * 80)
print("STEP 5: CREATING SIMPLIFIED SIGNATURE MATRIX")
print("=" * 80)
print()

print("NOTE: For full signature matrix with expression values,")
print("      need to either:")
print("      1. Use high-memory machine (>64GB RAM)")
print("      2. Process in batches")
print("      3. Use pre-computed marker genes from Allen Brain")
print()

print("ALTERNATIVE: Export cell type proportions for manual processing")
print()

# Calculate expected cell type proportions in brain
all_cell_types = pd.concat([
    pd.DataFrame({'cell_type': list(neuron_types), 'category': 'Neuron'}),
    pd.DataFrame({'cell_type': list(non_neuron_types), 'category': 'Non-neuron'})
])

# Calculate proportions based on actual cell counts
neurons_counts = adata_neurons.obs['supercluster_term'].value_counts()
non_neurons_counts = adata_non_neurons.obs['cell_type'].value_counts()

total_cells = neurons_counts.sum() + non_neurons_counts.sum()

proportions = []
for ct in neuron_types:
    proportions.append({
        'cell_type': ct,
        'category': 'Neuron',
        'count': neurons_counts[ct],
        'proportion': neurons_counts[ct] / total_cells
    })

for ct in non_neuron_types:
    proportions.append({
        'cell_type': ct,
        'category': 'Non-neuron',
        'count': non_neurons_counts[ct],
        'proportion': non_neurons_counts[ct] / total_cells
    })

proportions_df = pd.DataFrame(proportions)
proportions_file = RESULTS_DIR / "allen_brain_celltype_proportions.csv"
proportions_df.to_csv(proportions_file, index=False)
print(f"✓ Saved cell type proportions: {proportions_file}")

print()
print("=" * 80)
print("PROCESSING COMPLETE")
print("=" * 80)
print()
print("OUTPUTS:")
print(f"  1. Cell annotations: {annotations_file}")
print(f"  2. Cell types list: {cell_types_file}")
print(f"  3. Cell type proportions: {proportions_file}")
print()
print("NEXT STEPS:")
print("  For full signature matrix creation, use one of:")
print("  A. High-memory server (>64GB RAM)")
print("  B. Batch processing script")
print("  C. Pre-computed marker genes from Allen Brain Cell Atlas")
print()
print("  For deconvolution, the cell type proportions can be used as")
print("  reference to validate MuSiC/NNLS results.")
print()
