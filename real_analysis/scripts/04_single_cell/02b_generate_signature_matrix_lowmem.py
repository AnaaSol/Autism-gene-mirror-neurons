#!/usr/bin/env python3
"""
GENERATE SIGNATURE MATRIX - LOW MEMORY/DISK MODE
================================================
Purpose: Calculate mean gene expression per cell type by processing
         h5ad files in chunks, avoiding memory/disk overflow.

Input: Siletti et al. 2023 .h5ad files (Neurons & Non-neurons)
Output: gtex_singlecell_signature_matrix.csv (Ready for MuSiC)

Author: Ana Sol Murzi (Optimized for Low Resource Environment)
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import scipy.sparse as sp
import warnings
import gc

warnings.filterwarnings('ignore')

print("=" * 80)
print("GENERATING SIGNATURE MATRIX (LOW MEMORY MODE)")
print("=" * 80)

# PATHS
BASE_DIR = Path("/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis")
DATA_DIR = BASE_DIR / "data/allen_brain_atlas"
RESULTS_DIR = BASE_DIR / "results/celltype_deconvolution"
NEURONS_FILE = DATA_DIR / "All-neuron.h5ad"
NON_NEURONS_FILE = DATA_DIR / "All-non-neuronal-cells.h5ad"

# CONFIG
CHUNK_SIZE = 10000  # Number of cells to load at once

def compute_mean_expression(file_path, cell_type_col, prefix=""):
    print(f"\nProcessing: {file_path.name}")
    print(f"Target column: {cell_type_col}")
    
    # 1. Load Metadata Only
    adata_backed = sc.read_h5ad(file_path, backed='r')
    total_cells = adata_backed.n_obs
    var_names = adata_backed.var_names.to_numpy().astype(str)
    n_vars = len(var_names)
    
    # Get unique cell types
    cell_types = adata_backed.obs[cell_type_col].unique().tolist()
    print(f"Found {len(cell_types)} cell types.")
    
    # 2. Initialize accumulators (Stores Sum of expression)
    # Shape: (n_cell_types, n_genes)
    type_sums = pd.DataFrame(0.0, index=cell_types, columns=var_names)
    type_counts = pd.Series(0, index=cell_types)
    
    # 3. Iterate in Chunks
    print(f"Processing {total_cells:,} cells in chunks of {CHUNK_SIZE}...")
    
    start = 0
    chunk_num = 0
    
    while start < total_cells:
        end = min(start + CHUNK_SIZE, total_cells)
        
        # Load ONLY this chunk into memory
        chunk = adata_backed[start:end]
        
        # Get X (expression) and obs (metadata) for chunk
        # Ideally raw counts or normalized. Usually Siletti data is raw counts in X.
        # Note: If X is sparse, we densify only when adding to sum to save RAM
        X_chunk = chunk.X
        
        # Get cell types for this chunk
        chunk_types = chunk.obs[cell_type_col].values
        
        # Aggregate per cell type in this chunk
        for ct in cell_types:
            # Boolean mask for this cell type in current chunk
            mask = (chunk_types == ct)
            count = np.sum(mask)
            
            if count > 0:
                # Sum expression for these cells
                if sp.issparse(X_chunk):
                    # Efficient sparse sum
                    sum_expr = X_chunk[mask].sum(axis=0).A1 
                else:
                    sum_expr = X_chunk[mask].sum(axis=0)
                
                # Add to global accumulators
                type_sums.loc[ct] += sum_expr
                type_counts[ct] += count
        
        # Clean up memory
        del chunk, X_chunk
        gc.collect()
        
        start = end
        chunk_num += 1
        if chunk_num % 10 == 0:
            print(f"  Processed {end:,} cells...")

    adata_backed.file.close()
    
    # 4. Calculate Means
    print("Calculating means...")
    # Avoid division by zero
    type_counts = type_counts.replace(0, 1) 
    signature_df = type_sums.div(type_counts, axis=0)
    
    return signature_df.T  # Transpose to (Genes x CellTypes)

# ==============================================================================
# EXECUTION
# ==============================================================================

# 1. Process Neurons (Using 'supercluster_term')
# Using supercluster_term gives us "MGE interneuron", "L2/3 IT", etc.
neuron_sig = compute_mean_expression(NEURONS_FILE, 'supercluster_term')

# 2. Process Non-Neurons (Using 'cell_type')
# Using cell_type gives us "Astrocyte", "Oligodendrocyte", etc.
non_neuron_sig = compute_mean_expression(NON_NEURONS_FILE, 'cell_type')

print("\nMerging matrices...")

# 3. Merge and Handle Missing Genes
# Align indexes (genes) - keep only common genes or union
common_genes = neuron_sig.index.intersection(non_neuron_sig.index)
print(f"Common genes: {len(common_genes):,}")

final_sig = pd.concat([
    neuron_sig.loc[common_genes], 
    non_neuron_sig.loc[common_genes]
], axis=1)

print(f"Final Matrix Shape: {final_sig.shape}")

# 4. Save
output_file = RESULTS_DIR / "gtex_singlecell_signature_matrix.csv" # Keeping name for compatibility
print(f"Saving to {output_file}...")
final_sig.to_csv(output_file)

print("\n" + "="*80)
print("✓ SUCCESS! Signature Matrix Created.")
print("="*80)
print(f"Disk space used: Only ~{final_sig.memory_usage().sum()/1e6:.2f} MB")
