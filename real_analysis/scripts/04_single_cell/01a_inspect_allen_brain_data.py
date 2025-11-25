#!/usr/bin/env python3
"""
Inspect Allen Brain Atlas h5ad files
"""
import scanpy as sc
import pandas as pd
import numpy as np

print("=" * 80)
print("INSPECTING ALLEN BRAIN ATLAS DATA")
print("=" * 80)
print()

# File paths
NEURONS_FILE = "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis/data/allen_brain_atlas/All-neuron.h5ad"
NON_NEURONS_FILE = "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis/data/allen_brain_atlas/All-non-neuronal-cells.h5ad"

print("Loading NEURONS dataset (backed mode - won't load full matrix)...")
adata_neurons = sc.read_h5ad(NEURONS_FILE, backed='r')

print(f"✓ Loaded neurons: {adata_neurons.n_obs:,} cells × {adata_neurons.n_vars:,} genes")
print()

print("=" * 80)
print("NEURONS METADATA COLUMNS:")
print("=" * 80)
print(adata_neurons.obs.columns.tolist())
print()

print("=" * 80)
print("SAMPLE OF NEURONS METADATA:")
print("=" * 80)
print(adata_neurons.obs.head())
print()

# Check for cell type columns
potential_celltype_cols = [col for col in adata_neurons.obs.columns
                           if any(keyword in col.lower()
                                 for keyword in ['type', 'class', 'subclass', 'cluster', 'annotation'])]

print("=" * 80)
print("POTENTIAL CELL TYPE COLUMNS:")
print("=" * 80)
for col in potential_celltype_cols:
    n_unique = adata_neurons.obs[col].nunique()
    print(f"\n{col}: {n_unique} unique values")
    if n_unique < 50:
        print(adata_neurons.obs[col].value_counts().head(20))
print()

print("=" * 80)
print("Loading NON-NEURONAL dataset (backed mode - won't load full matrix)...")
print("=" * 80)
adata_non_neurons = sc.read_h5ad(NON_NEURONS_FILE, backed='r')

print(f"✓ Loaded non-neurons: {adata_non_neurons.n_obs:,} cells × {adata_non_neurons.n_vars:,} genes")
print()

print("NON-NEURONS METADATA COLUMNS:")
print(adata_non_neurons.obs.columns.tolist())
print()

print("=" * 80)
print("SAMPLE OF NON-NEURONS METADATA:")
print("=" * 80)
print(adata_non_neurons.obs.head())
print()

# Check for cell type columns
potential_celltype_cols = [col for col in adata_non_neurons.obs.columns
                           if any(keyword in col.lower()
                                 for keyword in ['type', 'class', 'subclass', 'cluster', 'annotation'])]

print("=" * 80)
print("POTENTIAL CELL TYPE COLUMNS:")
print("=" * 80)
for col in potential_celltype_cols:
    n_unique = adata_non_neurons.obs[col].nunique()
    print(f"\n{col}: {n_unique} unique values")
    if n_unique < 100:
        print(adata_non_neurons.obs[col].value_counts().head(20))

print()
print("=" * 80)
print("SUMMARY:")
print("=" * 80)
print(f"Total cells: {adata_neurons.n_obs + adata_non_neurons.n_obs:,}")
print(f"Total genes (neurons): {adata_neurons.n_vars:,}")
print(f"Total genes (non-neurons): {adata_non_neurons.n_vars:,}")
print()
print("✓ Inspection complete")
print()
