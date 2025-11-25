#!/usr/bin/env python3
"""
EMPIRICAL REGULATOR SELECTION (STRICT MODE - 18 MIRROR GENES)
Select ONLY the Top 10 regulators for each target gene to build a focused network.
"""

import pandas as pd
import numpy as np
import gzip
import sys
from scipy import stats
from pathlib import Path
import heapq

# --- CONFIGURACIÓN ---
BASE_DIR = Path("/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis")
GCT_FILE = BASE_DIR / "data/gtex_bulk/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
ATTR_FILE = BASE_DIR / "data/gtex_bulk/GTEx_Samples.txt"
OUTPUT_FILE = BASE_DIR / "results/causal_networks/empirical_genes_list.txt"

# TU LISTA MAESTRA DE GENES (15 Genes Espejo)
TARGET_GENES = [
    "ALDOA", "CHD8", "CACNA1C", "CNTNAP2", "FOXP2", 
    "MAPK3", "OXTR", "TAOK2", "TBX6", "YPEL3", 
    "KCTD13", "MECP2", "PPP4C", "SHANK3", "THEMIS"
]

# Vamos a ser un poco más estrictos para mantener la red manejable
# 18 genes x 5 reguladores = 90 nodos (Ideal para visualización)
# Si ponemos 10, serían ~180 nodos (aún manejable, pero denso)
TOP_N = 3

# Lista amplia de potenciales reguladores para filtrar el genoma
KNOWN_TFS_SET = set([
    'CTCF', 'SP1', 'REST', 'GATA3', 'FOXA1', 'NFKB1', 'RELA', 'JUN', 'FOS', 
    'MYC', 'STAT3', 'ESR1', 'NR3C1', 'POU5F1', 'SOX2', 'TP53', 'E2F1', 
    'YY1', 'CREB1', 'SMAD3', 'HIF1A', 'SREBF1', 'PAX6', 'NEUROD1', 'OLIG2',
    'TBR1', 'FOXP1', 'FOXP2', 'CHD8', 'MECP2', 'RBFOX1', 'ELF4', 'MEF2C',
    'NPAS4', 'BDNF', 'FMR1', 'ADNP', 'ASH1L', 'BCL11A', 'CTNNB1', 'DLX1',
    'DLX2', 'EGR1', 'EMX1', 'FOXG1', 'GABPA', 'HDAC1', 'HDAC2', 'HES1',
    'LHX2', 'LHX6', 'MEF2A', 'NANOG', 'NEUROG2', 'NKX2-1', 'NR4A2', 'OTX2',
    'PAX5', 'PBX1', 'POU3F2', 'RUNX1', 'SATB2', 'SIX3', 'SOX9', 'STAT5A',
    'TAL1', 'TCF4', 'TFE3', 'USF1', 'ZEB1', 'ZEB2', 'ZNF', 'POGZ',
    'EP300', 'CREBBP', 'KMT2A', 'EZH2', 'SUZ12', 'STAG2', 'RAD21', 'SMC3', 'FOXP1'
]) 

print("="*80)
print(f"BÚSQUEDA EMPÍRICA DE REGULADORES ({len(TARGET_GENES)} GENES)")
print("="*80)

# 1. OBTENER MUESTRAS DE CEREBRO
print("1. Identificando muestras de cerebro...")
samples_df = pd.read_csv(ATTR_FILE, sep="\t", low_memory=False)
brain_samples_ids = set(samples_df[samples_df['SMTSD'].str.contains("Brain", na=False)]['SAMPID'])
print(f"   Total muestras cerebro: {len(brain_samples_ids)}")

# 2. LEER TARGET GENES
print("2. Extrayendo expresión de genes objetivo...")
target_expression = {}
sample_indices = []

with gzip.open(GCT_FILE, 'rt') as f:
    f.readline()
    f.readline()
    header = f.readline().strip().split('\t')
    
    for i, sid in enumerate(header[2:]):
        if sid in brain_samples_ids:
            sample_indices.append(i)
            
    print(f"   Indices mapeados: {len(sample_indices)}")
    
    for line in f:
        parts = line.strip().split('\t')
        gene_name = parts[1]
        
        if gene_name in TARGET_GENES:
            values = [float(parts[2+i]) for i in sample_indices]
            target_expression[gene_name] = np.array(values)

found = list(target_expression.keys())
print(f"   Genes encontrados ({len(found)}/{len(TARGET_GENES)}): {found}")

# 3. BARRIDO GENÓMICO
print("3. Buscando Top reguladores por gen...")
top_candidates = {t: [] for t in found}

with gzip.open(GCT_FILE, 'rt') as f:
    for _ in range(3): f.readline()
    
    for i, line in enumerate(f):
        parts = line.strip().split('\t')
        gene_sym = parts[1]
        
        # Filtro: TF conocido o familia relevante
        is_potential_tf = (gene_sym in KNOWN_TFS_SET) or \
                          any(x in gene_sym for x in ['ZNF', 'FOX', 'SOX', 'HOX', 'STAT', 'HDAC', 'KMT'])
        
        if not is_potential_tf or gene_sym in TARGET_GENES:
            continue
            
        values = np.array([float(parts[2+idx]) for idx in sample_indices])
        if np.std(values) == 0: continue

        for target in found:
            target_vals = target_expression[target]
            if np.std(target_vals) == 0: continue
            
            corr, _ = stats.pearsonr(np.log1p(values), np.log1p(target_vals))
            abs_corr = abs(corr)
            
            if abs_corr > 0.3:
                heap = top_candidates[target]
                if len(heap) < TOP_N:
                    heapq.heappush(heap, (abs_corr, gene_sym))
                else:
                    if abs_corr > heap[0][0]:
                        heapq.heappushpop(heap, (abs_corr, gene_sym))

        if i % 5000 == 0:
            print(f"   Analizados {i} genes...", end='\r')

# 4. CONSOLIDAR Y GUARDAR
print("\n4. Consolidando lista final...")
final_regulators = set()

for target, candidates in top_candidates.items():
    sorted_candidates = sorted(candidates, key=lambda x: x[0], reverse=True)
    # Descomentar para ver detalles
    # print(f"\n--- {target} ---")
    for corr, reg in sorted_candidates:
        # print(f"  {reg}: {corr:.3f}")
        final_regulators.add(reg)

final_list = list(final_regulators) + found
print(f"\nTotal nodos para la red: {len(final_list)}")

with open(OUTPUT_FILE, 'w') as f:
    f.write("\n".join(final_list))

print(f"✓ Lista guardada en: {OUTPUT_FILE}")
