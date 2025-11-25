# ==============================================================================
# COLOCALIZACIÓN ESTÁNDAR EN BULK (CORREGIDA: SIN DUPLICADOS)
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(coloc)
})

cat("================================================================\n")
cat("COLOCALIZACIÓN BULK TISSUE (REAL DATA - CLEANED)\n")
cat("================================================================\n\n")

# --- RUTAS ---
BASE_DIR <- "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis"
RESULTS_DIR <- file.path(BASE_DIR, "results/colocalization_standard")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

EQTL_FILE <- file.path(BASE_DIR, "results/eqtl_extracted/real_gtex_eqtls.csv")
GWAS_FILE <- file.path(BASE_DIR, "data/gwas/iPSYCH-PGC_ASD_Nov2017")
OUTPUT_FILE <- file.path(RESULTS_DIR, "bulk_coloc_all_variants.csv")

# Tamaños de muestra
tissue_sample_sizes <- list(
    "Amygdala" = 129, "Anterior cingulate cortex (BA24)" = 147,
    "Caudate (basal ganglia)" = 194, "Cerebellar Hemisphere" = 175,
    "Cerebellum" = 209, "Cortex" = 205, "Frontal Cortex (BA9)" = 175,
    "Hippocampus" = 165, "Hypothalamus" = 170,
    "Nucleus accumbens (basal ganglia)" = 202, "Putamen (basal ganglia)" = 170,
    "Spinal cord (cervical c-1)" = 126, "Substantia nigra" = 114
)

# 1. CARGAR DATOS
cat("1. Cargando datos...\n")
if(!file.exists(EQTL_FILE)) stop("Falta real_gtex_eqtls.csv")
eqtls_all <- read.csv(EQTL_FILE)

# Parsear posiciones y LIMPIEZA INICIAL DE DUPLICADOS
# Si un SNP aparece dos veces para el mismo gen/tejido, nos quedamos con el menor P-val
eqtls_all <- eqtls_all %>%
    mutate(
        chr_clean = sub("chr", "", sub("_.*", "", variant_id)),
        pos_hg38 = as.numeric(sub(".*_([0-9]+)_.*_.*_.*", "\\1", variant_id))
    ) %>%
    group_by(gene_symbol, tissue, variant_id) %>%
    slice_min(order_by = pval_nominal, n = 1, with_ties = FALSE) %>%
    ungroup()

if(!file.exists(GWAS_FILE)) stop("Falta archivo GWAS")
gwas_raw <- fread(GWAS_FILE, select = c("CHR", "BP", "P", "OR", "SE", "SNP"))
# Limpieza de GWAS por si acaso
gwas_raw <- gwas_raw[!duplicated(gwas_raw$SNP), ]

cat("   Datos cargados y limpios.\n")

results_list <- list()
unique_pairs <- eqtls_all %>% select(gene_symbol, tissue) %>% distinct()

cat(sprintf("   Se analizarán %d pares Gen-Tejido.\n", nrow(unique_pairs)))

# 2. BUCLE DE ANÁLISIS
for(i in 1:nrow(unique_pairs)) {
    gene <- unique_pairs$gene_symbol[i]
    tissue <- unique_pairs$tissue[i]
    
    cat(sprintf("\n--- Analizando %s en %s ---\n", gene, tissue))
    
    pair_eqtls <- eqtls_all %>% filter(gene_symbol == gene, tissue == tissue)
    if(nrow(pair_eqtls) < 5) { cat("   Pocas variantes (<5). Saltando.\n"); next }
    
    # Ventana GWAS
    chr_target <- pair_eqtls$chr_clean[1]
    if(chr_target == "X") chr_gwas <- 23 else chr_gwas <- chr_target
    min_pos <- min(pair_eqtls$pos_hg38) - 500000 
    max_pos <- max(pair_eqtls$pos_hg38) + 500000
    
    gwas_region <- gwas_raw[CHR == chr_gwas & BP >= min_pos & BP <= max_pos]
    if(nrow(gwas_region) == 0) { cat("   Sin datos GWAS.\n"); next }
    
    # Match y Filtrado de Duplicados Final
    top_eqtl <- pair_eqtls[which.min(pair_eqtls$pval_nominal), ]
    top_gwas_local <- gwas_region[which.min(abs(gwas_region$BP - top_eqtl$pos_hg38)), ]
    offset <- top_eqtl$pos_hg38 - top_gwas_local$BP
    
    beta1 <- c(); varbeta1 <- c(); p1 <- c(); maf1 <- c(); snp_ids <- c()
    beta2 <- c(); varbeta2 <- c(); p2 <- c()
    
    for(k in 1:nrow(pair_eqtls)) {
        eqtl_row <- pair_eqtls[k,]
        target_pos_hg19 <- eqtl_row$pos_hg38 - offset
        match_idx <- which.min(abs(gwas_region$BP - target_pos_hg19))
        dist <- abs(gwas_region$BP[match_idx] - target_pos_hg19)
        
        if(dist < 1000) {
            gwas_row <- gwas_region[match_idx, ]
            # Check duplicado en vectores (por si acaso)
            if(!(eqtl_row$variant_id %in% snp_ids)) {
                beta1 <- c(beta1, eqtl_row$slope)
                varbeta1 <- c(varbeta1, eqtl_row$slope_se^2)
                p1 <- c(p1, eqtl_row$pval_nominal)
                maf1 <- c(maf1, eqtl_row$maf)
                snp_ids <- c(snp_ids, eqtl_row$variant_id)
                
                beta2 <- c(beta2, log(gwas_row$OR))
                varbeta2 <- c(varbeta2, gwas_row$SE^2)
                p2 <- c(p2, gwas_row$P)
            }
        }
    }
    
    cat(sprintf("   Variantes únicas encontradas: %d\n", length(snp_ids)))
    
    if(length(snp_ids) < 5) { cat("   Insuficientes variantes comunes.\n"); next }
    
    # Ejecutar COLOC
    N_gtex <- tissue_sample_sizes[[tissue]]
    if(is.null(N_gtex)) N_gtex <- 150
    
    maf_clean <- maf1
    maf_clean[maf_clean == 0] <- 0.001
    sdY_vec <- sqrt(varbeta1) * sqrt(2 * N_gtex * maf_clean * (1-maf_clean))
    
    ds1 <- list(
        beta = beta1, varbeta = varbeta1, pvalues = p1,
        N = N_gtex, type = "quant", 
        snp = snp_ids, sdY = sdY_vec
    )
    
    ds2 <- list(
        beta = beta2, varbeta = varbeta2, pvalues = p2,
        N = 46350, type = "cc", s = 0.4,
        snp = snp_ids
    )
    
    tryCatch({
        res <- coloc.abf(ds1, ds2)
        pp4 <- res$summary["PP.H4.abf"]
        pp3 <- res$summary["PP.H3.abf"]
        cat(sprintf("   PP.H4: %.4f | PP.H3: %.4f\n", pp4, pp3))
        
        results_list[[length(results_list)+1]] <- data.frame(
            Gene = gene, Tissue = tissue,
            N_SNPs = length(snp_ids),
            PP_H3 = pp3, PP_H4 = pp4,
            Result = ifelse(pp4 > 0.75, "COLOC", "NO_COLOC")
        )
    }, error = function(e) { cat("   Error coloc:", e$message, "\n") })
}

# Guardar
if(length(results_list) > 0) {
    final_df <- do.call(rbind, results_list)
    write.csv(final_df, OUTPUT_FILE, row.names = FALSE)
    cat("\n✓ Guardado en:", OUTPUT_FILE, "\n")
}
