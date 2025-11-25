#!/usr/bin/env python3
import gzip
import os
import sys
from pathlib import Path

# CONFIGURACIÓN
# Ajusta la ruta a tu archivo GTEx descargado
INPUT_FILE = "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis/data/gtex_bulk/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
OUTPUT_DIR = "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis/data/gtex_bulk/splits"
BATCH_SIZE = 500  # Muestras por archivo (Seguro para 4GB RAM)

print("="*80)
print(f"DIVIDIENDO ARCHIVO GTEx GIGANTE EN LOTES DE {BATCH_SIZE} MUESTRAS")
print("="*80)

# Crear directorio de salida
os.makedirs(OUTPUT_DIR, exist_ok=True)

if not os.path.exists(INPUT_FILE):
    print(f"Error: No encuentro el archivo {INPUT_FILE}")
    sys.exit(1)

print(f"Leyendo: {os.path.basename(INPUT_FILE)}...")

# Abrimos el archivo gigante
with gzip.open(INPUT_FILE, 'rt') as f_in:
    # 1. Saltar las primeras 2 líneas de metadatos del formato GCT
    version = f_in.readline()
    dims = f_in.readline()
    
    # 2. Leer encabezado (Samples)
    header_line = f_in.readline().strip()
    header_parts = header_line.split('\t')
    
    # Estructura: [0]=Name, [1]=Description, [2:]=Samples
    gene_col_name = header_parts[0]  # "Name"
    sample_names = header_parts[2:]
    total_samples = len(sample_names)
    
    print(f"Total muestras detectadas: {total_samples:,}")
    print(f"Total lotes estimados: {(total_samples // BATCH_SIZE) + 1}")
    
    # 3. Preparar gestores de archivos para los lotes
    batch_files = []
    
    # Calculamos cuántos archivos necesitamos
    for i in range(0, total_samples, BATCH_SIZE):
        batch_idx = i // BATCH_SIZE
        fname = os.path.join(OUTPUT_DIR, f"batch_{batch_idx:02d}.csv.gz")
        
        # Abrimos archivo para escribir comprimido
        f_out = gzip.open(fname, 'wt')
        
        # Escribimos el encabezado: "Name,Sample1,Sample2..."
        # Tomamos el subconjunto de muestras para este lote
        current_samples = sample_names[i : i + BATCH_SIZE]
        f_out.write(gene_col_name + "," + ",".join(current_samples) + "\n")
        
        batch_files.append(f_out)
        
    print(f"✓ Archivos temporales creados. Procesando genes...")

    # 4. Procesar genes línea por línea
    gene_count = 0
    for line in f_in:
        parts = line.strip().split('\t')
        gene_id = parts[0]
        # Saltamos parts[1] (Description) para ahorrar espacio
        expression_values = parts[2:]
        
        # Repartir los valores en los archivos correspondientes
        for batch_idx, f_out in enumerate(batch_files):
            start = batch_idx * BATCH_SIZE
            end = start + BATCH_SIZE
            
            # Extraer valores para este lote
            batch_values = expression_values[start:end]
            
            # Escribir línea: "ENSG...,0.1,0.5,..."
            f_out.write(gene_id + "," + ",".join(batch_values) + "\n")
            
        gene_count += 1
        if gene_count % 5000 == 0:
            print(f"  Procesados {gene_count:,} genes...", end='\r')

    # 5. Cerrar todo
    for f in batch_files:
        f.close()

print(f"\n✓ ¡LISTO! Procesados {gene_count:,} genes.")
print(f"Los lotes están en: {OUTPUT_DIR}")
