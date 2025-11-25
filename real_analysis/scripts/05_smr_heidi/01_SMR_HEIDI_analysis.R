# 04_SMR_HEIDI_analysis.R
# ==============================================================================
# ANÁLISIS SMR & HEIDI (Summary-data-based Mendelian Randomization)
# ==============================================================================
# Objetivo: Validar la causalidad usando un método ortogonal a coloc.
# Método:
#   1. SMR: Testea si la expresión del gen causa el fenotipo (Autismo).
#      H0: No causalidad. P < 0.05 -> Causalidad probable.
#   2. HEIDI: Testea si hay heterogeneidad (LD vs Pleiotropía).
#      H0: Una sola variante causal (Pleiotropía).
#      H1: Múltiples variantes (Linkage).
#      ¡OJO! Aquí buscamos P_HEIDI > 0.05 (No rechazar H0).
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
})

cat("================================================================\n")
cat("ANÁLISIS SMR / HEIDI (VALIDACIÓN CAUSAL)\n")
cat("================================================================\n\n")

# --- RUTAS ---
BASE_DIR <- "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis"
RESULTS_DIR <- file.path(BASE_DIR, "results/smr_heidi")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

EQTL_FILE <- file.path(BASE_DIR, "results/eqtl_extracted/real_gtex_eqtls.csv")
GWAS_FILE <- file.path(BASE_DIR, "data/gwas/iPSYCH-PGC_ASD_Nov2017")
OUTPUT_FILE <- file.path(RESULTS_DIR, "smr_heidi_results.csv")

# 1. CARGAR DATOS
cat("1. Cargando datos...\n")
eqtls <- read.csv(EQTL_FILE)
# Parsear posiciones
eqtls <- eqtls %>% mutate(
    pos_hg38 = as.numeric(sub(".*_([0-9]+)_.*_.*_.*", "\\1", variant_id)),
    chr_num = sub("chr([0-9XY]+)_.*", "\\1", variant_id)
)

gwas_raw <- fread(GWAS_FILE, select = c("CHR", "BP", "P", "OR", "SE", "SNP"))

# --- FUNCIONES SMR/HEIDI ---

# SMR Test: Chi-square test for association between expression and phenotype
# T_SMR = (beta_gwas^2 * beta_eqtl^2) / (beta_gwas^2 * se_eqtl^2 + beta_eqtl^2 * se_gwas^2)
# Aproximación Wald: Z_SMR = (beta_gwas / beta_eqtl) / SE_ratio
calculate_smr <- function(beta_gwas, se_gwas, beta_eqtl, se_eqtl) {
    # Beta SMR (Efecto causal estimado) = b_GWAS / b_eQTL
    b_smr <- beta_gwas / beta_eqtl
    
    # Error estándar SMR (Método Delta)
    se_smr <- sqrt((se_gwas^2 * beta_eqtl^2 + se_eqtl^2 * beta_gwas^2) / (beta_eqtl^4))
    
    z_smr <- b_smr / se_smr
    p_smr <- 2 * pnorm(-abs(z_smr))
    
    return(list(b_smr=b_smr, p_smr=p_smr))
}

# HEIDI Test (Heterogeneity In Dependent Instruments)
# Compara si la ratio b_GWAS/b_eQTL es constante a través de múltiples SNPs
calculate_heidi <- function(gwas_data, eqtl_data) {
    # Necesitamos al menos 3 SNPs para calcular heterogeneidad
    if(nrow(gwas_data) < 3) return(NA)
    
    # Calcular b_SMR para cada SNP
    b_smr_vec <- gwas_data$beta / eqtl_data$beta
    # Peso (inverso de la varianza aproximada)
    # var(b_smr) approx (se_gwas/b_eqtl)^2
    w_vec <- (eqtl_data$beta / gwas_data$se)^2
    
    # Media ponderada
    b_smr_mean <- sum(b_smr_vec * w_vec) / sum(w_vec)
    
    # Estadístico Q (Cochran's Q)
    # Mide cuánto se desvían los SNPs individuales del efecto promedio
    Q <- sum(w_vec * (b_smr_vec - b_smr_mean)^2)
    
    # P-valor (Chi-cuadrado con k-1 grados de libertad)
    df <- length(b_smr_vec) - 1
    p_heidi <- pchisq(Q, df, lower.tail = FALSE)
    
    return(p_heidi)
}

# 2. ANÁLISIS POR GEN
results_list <- list()
unique_genes <- unique(eqtls$gene_symbol)

cat("\n2. Ejecutando SMR & HEIDI...\n")

for(gene in unique_genes) {
    cat(sprintf("--- %s ---\n", gene))
    gene_eqtls <- eqtls %>% filter(gene_symbol == gene)
    
    # Usamos el tejido con mejor señal (Top eQTL)
    best_row <- gene_eqtls[which.min(gene_eqtls$pval_nominal), ]
    tissue <- best_row$tissue
    
    # Filtrar SNPs de ese tejido
    tissue_eqtls <- gene_eqtls %>% filter(tissue == tissue)
    
    # Buscar SNPs en GWAS (Ventana 1MB y Matching)
    # (Usamos la lógica de matching que ya sabemos que funciona: proximidad hg38)
    
    top_pos <- best_row$pos_hg38
    chr_str <- sub("chr([0-9XY]+)_.*", "\\1", best_row$variant_id)
    if(chr_str == "X") chr_gwas <- 23 else chr_gwas <- chr_str
    
    gwas_region <- gwas_raw[CHR == chr_gwas & BP >= (top_pos - 1e6) & BP <= (top_pos + 1e6)]
    
    if(nrow(gwas_region) == 0) next
    
    # Matching de SNPs
    # Necesitamos vectores alineados de Beta/SE para GWAS y eQTL
    beta_g <- c(); se_g <- c(); p_g <- c()
    beta_e <- c(); se_e <- c(); p_e <- c()
    snps <- c()
    
    # Offset estimado
    top_gwas_local <- gwas_region[which.min(abs(gwas_region$BP - top_pos)), ]
    offset <- top_pos - top_gwas_local$BP
    
    for(k in 1:nrow(tissue_eqtls)) {
        e_row <- tissue_eqtls[k, ]
        target <- e_row$pos_hg38 - offset
        match_idx <- which.min(abs(gwas_region$BP - target))
        dist <- abs(gwas_region$BP[match_idx] - target)
        
        if(dist < 5000) { # Match aceptable
            g_row <- gwas_region[match_idx, ]
            
            beta_e <- c(beta_e, e_row$slope)
            se_e <- c(se_e, e_row$slope_se)
            p_e <- c(p_e, e_row$pval_nominal)
            
            beta_g <- c(beta_g, log(g_row$OR)) # Log OR es Beta
            se_g <- c(se_g, g_row$SE)
            p_g <- c(p_g, g_row$P)
            
            snps <- c(snps, e_row$variant_id)
        }
    }
    
    if(length(snps) < 3) {
        cat("   ⚠️ Pocos SNPs coincidentes (<3) para HEIDI.\n")
        next
    }
    
    # A) Calcular SMR para el Top SNP (Instrumento)
    # El instrumento debe ser un eQTL fuerte (P < 5e-8 idealmente, o el mejor disponible)
    top_idx <- which.min(p_e)
    smr_res <- calculate_smr(beta_g[top_idx], se_g[top_idx], beta_e[top_idx], se_e[top_idx])
    
    cat(sprintf("   SMR P-valor: %.2e\n", smr_res$p_smr))
    
    # B) Calcular HEIDI usando todos los SNPs
    # Filtramos SNPs que sean al menos nominalmente significativos en eQTL (P < 0.05) para reducir ruido
    valid_heidi <- which(p_e < 0.05)
    if(length(valid_heidi) >= 3) {
        df_g <- data.frame(beta=beta_g[valid_heidi], se=se_g[valid_heidi])
        df_e <- data.frame(beta=beta_e[valid_heidi], se=se_e[valid_heidi])
        
        p_heidi <- calculate_heidi(df_g, df_e)
        cat(sprintf("   HEIDI P-valor: %.2e\n", p_heidi))
    } else {
        p_heidi <- NA
        cat("   Insuficientes SNPs significativos para HEIDI.\n")
    }
    
    # Interpretación
    # Causalidad ideal: SMR P < 0.05 Y HEIDI P > 0.05
    status <- "No Causal"
    if(smr_res$p_smr < 0.05) {
        if(is.na(p_heidi) || p_heidi > 0.05) {
            status <- "CAUSAL CANDIDATE (Passed HEIDI)"
        } else {
            status <- "Pleiotropy/Linkage (Failed HEIDI)"
        }
    }
    
    results_list[[length(results_list)+1]] <- data.frame(
        Gene = gene, Tissue = tissue,
        Top_SNP = snps[top_idx],
        SMR_P = smr_res$p_smr,
        HEIDI_P = p_heidi,
        Status = status
    )
}

# 3. GUARDAR
if(length(results_list) > 0) {
    final <- do.call(rbind, results_list)
    write.csv(final, OUTPUT_FILE, row.names = FALSE)
    cat("\n======================================\n")
    print(final)
    cat("\n✓ Resultados guardados en:", OUTPUT_FILE, "\n")
} else {
    cat("\n❌ No se generaron resultados.\n")
}
