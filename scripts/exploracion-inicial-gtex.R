# ==============================================================================
# GTEx - ANÃLISIS FINAL LIMPIO Y ORDENADO
# Script completo para anÃ¡lisis por tejidos cerebrales
# ==============================================================================

library(data.table)
library(dplyr)

# Configurar data.table para usar menos memoria
setDTthreads(2)

cat("=== ANÃLISIS GTEx OPTIMIZADO PARA MEMORIA ===\n\n")

# ==============================================================================
# FUNCIÃ“N 1: CARGAR SOLO METADATOS
# ==============================================================================

load_metadata_only <- function() {
  cat("Cargando solo metadatos (ligero)...\n")
  
  sample_metadata <- fread("~/Desktop/autism_mirror_neurons/data/raw/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
  
  # Definir tejidos de interÃ©s
  brain_tissues <- c(
    "Brain - Frontal Cortex (BA9)",
    "Brain - Cortex",
    "Brain - Anterior cingulate cortex (BA24)",
    "Brain - Cerebellar Hemisphere"
  )
  
  # Filtrar solo muestras cerebrales
  brain_samples <- sample_metadata[SMTSD %in% brain_tissues]
  
  cat(paste("Muestras cerebrales identificadas:", nrow(brain_samples), "\n"))
  
  return(list(
    brain_samples = brain_samples,
    all_metadata = sample_metadata,
    brain_tissues = brain_tissues
  ))
}

# ==============================================================================
# FUNCIÃ“N 2: ANÃLISIS POR TEJIDO (CORREGIDO)
# ==============================================================================

analyze_by_tissue_light <- function(brain_expr_data, metadata_info) {
  cat("\n=== ANÃLISIS POR TEJIDO (CORREGIDO) ===\n")
  
  # Preparar info de muestras
  brain_samples <- metadata_info$brain_samples
  sample_columns <- names(brain_expr_data)[-(1:2)]
  
  # Crear mapeo muestra -> tejido (solo para muestras que tenemos en expresiÃ³n)
  sample_tissue_map <- brain_samples[SAMPID %in% sample_columns, .(SAMPID, SMTSD)]
  setnames(sample_tissue_map, c("sample_id", "tissue"))
  
  cat(paste("Total muestras con expresiÃ³n:", length(sample_columns), "\n"))
  cat(paste("Muestras mapeadas a tejidos:", nrow(sample_tissue_map), "\n"))
  
  # Verificar distribuciÃ³n por tejido
  tissue_counts <- sample_tissue_map[, .N, by = tissue]
  cat("DistribuciÃ³n por tejido:\n")
  print(tissue_counts)
  
  # AnÃ¡lisis por tejido
  results_by_tissue <- list()
  
  for(current_tissue in unique(sample_tissue_map$tissue)) {  # <- VARIABLE CORRECTA
    cat(paste("\n--- Procesando:", current_tissue, "---\n"))
    
    # Obtener muestras de este tejido ESPECÃFICO
    tissue_samples <- sample_tissue_map[tissue == current_tissue]$sample_id  # <- COMPARACIÃ“N CORRECTA
    cat(paste("Muestras identificadas para este tejido:", length(tissue_samples), "\n"))
    
    # Verificar que las muestras existen en los datos
    available_tissue_samples <- intersect(tissue_samples, sample_columns)
    cat(paste("Muestras disponibles en datos:", length(available_tissue_samples), "\n"))
    
    if(length(available_tissue_samples) == 0) {
      cat("Â¡Sin muestras disponibles! Saltando tejido.\n")
      next
    }
    
    # Extraer datos de este tejido ESPECÃFICO
    tissue_cols <- c("Name", "Description", available_tissue_samples)
    tissue_data <- brain_expr_data[, ..tissue_cols]
    
    # Aplicar QC especÃ­fico por tejido
    cat("Aplicando QC especÃ­fico por tejido...\n")
    expr_matrix <- as.matrix(tissue_data[, -(1:2)])
    
    # Filtrar genes con expresiÃ³n en al menos 50% de muestras de este tejido
    min_samples_tissue <- ceiling(ncol(expr_matrix) * 0.5)
    expressed_per_gene <- rowSums(expr_matrix > 1)  # TPM > 1
    keep_genes_tissue <- expressed_per_gene >= min_samples_tissue
    
    # Aplicar filtro
    tissue_data_filtered <- tissue_data[keep_genes_tissue]
    expr_matrix_filtered <- as.matrix(tissue_data_filtered[, -(1:2)])
    
    cat(paste("Genes antes de QC tisular:", nrow(tissue_data), "\n"))
    cat(paste("Genes despuÃ©s de QC tisular:", nrow(tissue_data_filtered), "\n"))
    
    # Normalizar (log2)
    expr_log2 <- log2(expr_matrix_filtered + 1)
    
    # Crear resultado final
    tissue_normalized <- data.table(
      gene_id = tissue_data_filtered$Name,
      expr_log2
    )
    
    results_by_tissue[[current_tissue]] <- tissue_normalized
    
    # EstadÃ­sticas finales
    cat(paste("- Genes finales:", nrow(tissue_normalized), "\n"))
    cat(paste("- Muestras finales:", ncol(expr_log2), "\n"))
    cat(paste("- ExpresiÃ³n media:", round(mean(expr_log2), 3), "\n"))
    cat(paste("- ExpresiÃ³n mediana:", round(median(expr_log2), 3), "\n"))
  }
  
  return(results_by_tissue)
}

# ==============================================================================
# FUNCIÃ“N 3: GUARDAR RESULTADOS
# ==============================================================================

save_optimized_results <- function(results_by_tissue, metadata_info) {
  cat("\n=== GUARDANDO RESULTADOS ===\n")
  
  # Crear directorios
  if(!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
  if(!dir.exists("data/processed/by_tissue_optimized")) dir.create("data/processed/by_tissue_optimized", recursive = TRUE)
  
  # Guardar por tejido
  summary_stats <- data.table()
  
  for(tissue_name in names(results_by_tissue)) {
    # Limpiar nombre
    tissue_clean <- gsub("[^A-Za-z0-9]", "_", tissue_name)
    
    # Guardar datos
    fwrite(results_by_tissue[[tissue_name]], 
           paste0("data/processed/by_tissue_optimized/", tissue_clean, "_normalized.txt"),
           sep = "\t", quote = FALSE)
    
    # EstadÃ­sticas
    tissue_data <- results_by_tissue[[tissue_name]]
    expr_matrix <- as.matrix(tissue_data[, -1])
    
    summary_stats <- rbind(summary_stats, data.table(
      Tejido = tissue_name,
      Genes = nrow(tissue_data),
      Muestras = ncol(expr_matrix),
      Expr_Media = round(mean(expr_matrix), 3),
      Expr_Mediana = round(median(expr_matrix), 3),
      Expr_Min = round(min(expr_matrix), 3),
      Expr_Max = round(max(expr_matrix), 3)
    ))
    
    cat(paste("Guardado:", tissue_clean, "\n"))
  }
  
  # Guardar resumen
  fwrite(summary_stats, "data/processed/summary_optimized.txt", sep = "\t", quote = FALSE)
  
  return(summary_stats)
}

# ==============================================================================
# EJECUCIÃ“N PRINCIPAL EN ORDEN CORRECTO
# ==============================================================================

cat("Iniciando anÃ¡lisis optimizado...\n")

# PASO 1: Cargar metadatos
metadata_info <- load_metadata_only()

# PASO 2: Cargar y procesar datos de expresiÃ³n
cat("Usando estrategia de matriz completa...\n")
cat("Cargando datos de expresiÃ³n...\n")

expr_data <- fread("~/Desktop/autism_mirror_neurons/data/raw/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", 
                   skip = 2, header = TRUE)

# Obtener muestras cerebrales
brain_sample_ids <- metadata_info$brain_samples$SAMPID
expr_sample_cols <- names(expr_data)[-(1:2)]
available_brain_samples <- intersect(brain_sample_ids, expr_sample_cols)

# Filtrar solo columnas necesarias
brain_expr_cols <- c("Name", "Description", available_brain_samples)
brain_expr_data <- expr_data[, ..brain_expr_cols]
rm(expr_data); gc()  # Limpiar memoria

# PASO 3: Aplicar QC general
cat("Aplicando control de calidad...\n")
gene_info <- brain_expr_data[, .(Name, Description)]
expr_matrix <- as.matrix(brain_expr_data[, -(1:2)])

# QC: genes expresados en >20% muestras con TPM > 0.1
expressed_samples <- rowSums(expr_matrix > 0.1)
min_samples_threshold <- ceiling(ncol(expr_matrix) * 0.2)
keep_genes <- expressed_samples >= min_samples_threshold

brain_expr_qc <- brain_expr_data[keep_genes]
cat(paste("Genes despuÃ©s de QC:", nrow(brain_expr_qc), "\n"))

rm(brain_expr_data, expr_matrix); gc()

# PASO 4: Analizar por tejido
results_by_tissue <- analyze_by_tissue_light(brain_expr_qc, metadata_info)

# PASO 5: Guardar resultados
if(exists("results_by_tissue") && length(results_by_tissue) > 0) {
  summary_stats <- save_optimized_results(results_by_tissue, metadata_info)
  
  # Resumen final
  cat("\n=== RESUMEN FINAL CORRECTO ===\n")
  print(summary_stats)
  
  cat("\nArchivos generados en data/processed/by_tissue_optimized/\n")
  cat("AnÃ¡lisis completado exitosamente.\n")
} else {
  cat("ERROR: No se pudo crear results_by_tissue\n")
}
