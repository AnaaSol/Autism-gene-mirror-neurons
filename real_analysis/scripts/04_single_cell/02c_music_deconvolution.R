# 02_music_deconvolution.R
# MODO: NNLS (Non-Negative Least Squares) por Lotes
# PROPÓSITO: Deconvolución usando Matriz de Firmas pre-calculada
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(Matrix)
    library(tidyverse)
    library(data.table)
    library(nnls) # Usamos el motor NNLS directo
})

cat("================================================================\n")
cat("DECONVOLUCIÓN NNLS POR LOTES (REAL DATA)\n")
cat("================================================================\n\n")

# --- RUTAS ---
BASE_DIR <- "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis"
SIG_FILE <- file.path(BASE_DIR, "results/celltype_deconvolution/gtex_singlecell_signature_matrix.csv")
SPLIT_DIR <- file.path(BASE_DIR, "data/gtex_bulk/splits")
OUTPUT_FILE <- file.path(BASE_DIR, "results/celltype_deconvolution/gtex_celltype_proportions.csv")

# 1. CARGAR FIRMAS SINGLE-CELL
cat("1. Cargando firmas (Allen Brain Atlas)...\n")
if(!file.exists(SIG_FILE)) stop("Falta matriz de firmas")
sig_matrix <- read.csv(SIG_FILE, row.names = 1, check.names = FALSE)
cat("   Dimensiones:", dim(sig_matrix)[1], "genes x", dim(sig_matrix)[2], "tipos celulares\n")

# 2. LISTAR LOTES
batch_files <- list.files(SPLIT_DIR, pattern = "batch_.*.csv.gz", full.names = TRUE)
if(length(batch_files) == 0) stop("No hay lotes en: ", SPLIT_DIR)
cat("2. Encontrados", length(batch_files), "lotes para procesar.\n\n")

# --- FUNCIÓN AUXILIAR NNLS ---
run_nnls <- function(bulk_mat, sig_mat) {
    # bulk_mat: Genes x Muestras
    # sig_mat:  Genes x TiposCelulares
    
    # Inicializar matriz de resultados
    n_samples <- ncol(bulk_mat)
    n_types <- ncol(sig_mat)
    results <- matrix(0, nrow = n_samples, ncol = n_types)
    colnames(results) <- colnames(sig_mat)
    rownames(results) <- colnames(bulk_mat)
    
    # Iterar por cada muestra (columna)
    for(j in 1:n_samples) {
        y <- bulk_mat[, j]
        # Resolver: minimizar ||y - sig * x|| sujeto a x >= 0
        fit <- nnls::nnls(as.matrix(sig_mat), y)
        coefs <- coef(fit)
        
        # Normalizar para que sumen 1 (Proporciones)
        if(sum(coefs) > 0) {
            coefs <- coefs / sum(coefs)
        }
        results[j, ] <- coefs
    }
    return(as.data.frame(results))
}

# 3. PROCESAR CADA LOTE
all_results <- list()

for(i in seq_along(batch_files)) {
    bf <- batch_files[i]
    cat(sprintf("--- Procesando Lote %d/%d: %s ---\n", i, length(batch_files), basename(bf)))
    
    # Leer lote (Forzando CSV)
    bulk_data <- data.table::fread(bf, data.table = FALSE, sep = ",")
    rownames(bulk_data) <- bulk_data[[1]] 
    bulk_matrix <- as.matrix(bulk_data[, -1]) 
    
    # Limpiar IDs
    rownames(bulk_matrix) <- sub("\\..*", "", rownames(bulk_matrix))
    
    # Cruzar genes
    common <- intersect(rownames(bulk_matrix), rownames(sig_matrix))
    
    if(length(common) > 0) {
        bulk_subset <- bulk_matrix[common, , drop=FALSE]
        sig_subset <- as.matrix(sig_matrix[common, ])
        
        cat("   Genes comunes:", length(common), "| Muestras:", ncol(bulk_subset), "\n")
        
        tryCatch({
            # EJECUTAR NNLS
            res <- run_nnls(bulk_subset, sig_subset)
            res$SampleID <- rownames(res)
            all_results[[i]] <- res
            cat("   ✓ Éxito\n")
            
        }, error = function(e) {
            cat("   ⚠️ Falló este lote:", e$message, "\n")
        })
    } else {
        cat("   ⚠️ Sin genes en común.\n")
    }
    
    rm(bulk_data, bulk_matrix, bulk_subset, sig_subset); gc()
}

# 4. UNIR Y GUARDAR
if(length(all_results) > 0) {
    cat("\n3. Uniendo resultados...\n")
    final_df <- do.call(rbind, all_results)
    final_df <- final_df %>% select(SampleID, everything())
    
    write.csv(final_df, OUTPUT_FILE, row.names = FALSE)
    cat("✓ DECONVOLUCIÓN COMPLETA. Guardado en:", OUTPUT_FILE, "\n")
} else {
    cat("\n❌ No se generaron resultados.\n")
}
