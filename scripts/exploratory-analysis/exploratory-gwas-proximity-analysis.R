# ===============================================================================
# ANÃLISIS DE COLOCALIZACIÃ“N - GENES DE NEURONAS ESPEJO
# ===============================================================================

suppressMessages({
  library(data.table)
  library(dplyr)
})

cat("=== ANÃLISIS DE COLOCALIZACIÃ“N ===\n")
cat("Genes de neuronas espejo vs GWAS autismo Grove 2019\n\n")

# 1. CARGAR DATOS
cat("1. Cargando datos...\n")

# Cargar GWAS Grove 2019
gwas <- fread("~/Desktop/autism_mirror_neurons/data/raw/Grove2019_ASD_GWAS.txt")
cat("  GWAS:", nrow(gwas), "SNPs cargados\n")

# Definir genes de neuronas espejo con coordenadas genÃ³micas precisas (hg19/GRCh37)
mirror_genes <- data.table(
  gene_symbol = c("FOXP2", "FOXP1", "CNTNAP2", "CHD8", "CACNA1C", "THEMIS"),
  chromosome = c(7, 3, 7, 14, 12, 6),
  start_pos = c(114086327, 71011832, 145818786, 21385455, 2162089, 42857000),
  end_pos = c(114693772, 71585540, 147891751, 21543618, 2798436, 42950000),
  evidence_level = c("High", "Medium", "Medium", "Low", "Low", "Medium"),
  function_category = c("Transcription", "Transcription", "Synaptic", "Chromatin", "Ion_channel", "Development")
)

cat("  Genes candidatos:", nrow(mirror_genes), "\n\n")

# 2. FUNCIÃ“N DE COLOCALIZACIÃ“N
cat("2. Definiendo ventanas cis y analizando colocalizaciÃ³n...\n")

analyze_colocalization <- function(gene_info, gwas_data, cis_window = 1e6) {
  
  gene_name <- gene_info$gene_symbol
  chr <- gene_info$chromosome
  gene_start <- gene_info$start_pos
  gene_end <- gene_info$end_pos
  
  cat("  Analizando", gene_name, "- Chr", chr, "\n")
  
  # Definir ventana cis (Â±1Mb del gen)
  window_start <- max(1, gene_start - cis_window)
  window_end <- gene_end + cis_window
  
  # Filtrar SNPs GWAS en la ventana cis
  cis_snps <- gwas_data[
    CHR == chr & 
    BP >= window_start & 
    BP <= window_end &
    !is.na(P) & P > 0
  ]
  
  if (nrow(cis_snps) == 0) {
    cat("    âš ï¸ Sin SNPs en ventana cis\n")
    return(NULL)
  }
  
  cat("    ðŸ“Š", nrow(cis_snps), "SNPs en ventana cis\n")
  
  # Calcular distancia al centro del gen
  gene_center <- (gene_start + gene_end) / 2
  cis_snps[, distance_to_gene := abs(BP - gene_center)]
  
  # Calcular evidencia de asociaciÃ³n
  cis_snps[, neg_log_p := -log10(P)]
  cis_snps[, z_score := sign(log(OR)) * sqrt(qchisq(1-P, 1))]
  
  # Score de colocalizaciÃ³n basado en proximidad y significancia
  cis_snps[, coloc_score := neg_log_p * exp(-distance_to_gene / 100000)]
  
  # Identificar top SNPs
  top_snps <- cis_snps[order(-coloc_score)][1:min(10, nrow(cis_snps))]
  
  # EstadÃ­sticas del gen
  gene_stats <- list(
    gene = gene_name,
    chromosome = chr,
    n_cis_snps = nrow(cis_snps),
    min_p = min(cis_snps$P),
    max_coloc_score = max(cis_snps$coloc_score),
    mean_distance = mean(cis_snps$distance_to_gene),
    top_snp = top_snps$SNP[1],
    top_snp_p = top_snps$P[1],
    top_snp_distance = top_snps$distance_to_gene[1]
  )
  
  # Categorizar evidencia
  if (gene_stats$min_p < 1e-5) {
    evidence_level <- "Strong"
  } else if (gene_stats$min_p < 1e-3) {
    evidence_level <- "Moderate" 
  } else if (gene_stats$min_p < 0.05) {
    evidence_level <- "Weak"
  } else {
    evidence_level <- "Minimal"
  }
  
  gene_stats$evidence_level <- evidence_level
  
  cat("    ðŸŽ¯ Min P-value:", sprintf("%.2e", gene_stats$min_p), 
      "- Evidencia:", evidence_level, "\n")
  
  return(list(
    stats = gene_stats,
    top_snps = top_snps,
    all_cis_snps = cis_snps
  ))
}

# 3. EJECUTAR ANÃLISIS PARA TODOS LOS GENES
cat("\n3. Ejecutando anÃ¡lisis de colocalizaciÃ³n...\n")
cat("==========================================\n")

colocalization_results <- list()
gene_summary <- data.table()

for (i in 1:nrow(mirror_genes)) {
  gene_info <- mirror_genes[i]
  
  result <- analyze_colocalization(gene_info, gwas)
  
  if (!is.null(result)) {
    colocalization_results[[gene_info$gene_symbol]] <- result
    
    # Agregar a resumen
    gene_summary <- rbind(gene_summary, as.data.table(result$stats))
  }
  cat("\n")
}

# 4. GENERAR RESUMEN CONSOLIDADO
cat("4. RESUMEN DE RESULTADOS\n")
cat("========================\n")

if (nrow(gene_summary) > 0) {
  # Ordenar por evidencia
  gene_summary_ordered <- gene_summary[order(min_p)]
  
  cat("Genes de neuronas espejo - Evidencia de asociaciÃ³n con autismo:\n\n")
  
  for (i in 1:nrow(gene_summary_ordered)) {
    gene_data <- gene_summary_ordered[i]
    cat(sprintf("ðŸ§¬ %-8s (Chr%d): P=%-10s, Evidencia=%-8s, SNPs_cis=%-4d\n",
                gene_data$gene, gene_data$chromosome, 
                sprintf("%.2e", gene_data$min_p),
                gene_data$evidence_level, gene_data$n_cis_snps))
  }
  
  cat("\nDistribuciÃ³n de evidencia:\n")
  evidence_counts <- table(gene_summary$evidence_level)
  for (level in names(evidence_counts)) {
    cat("  ", level, ":", evidence_counts[level], "genes\n")
  }
  
  # Top 3 genes con mejor evidencia
  cat("\nðŸ† TOP 3 GENES CON MEJOR EVIDENCIA:\n")
  top_3 <- gene_summary_ordered[1:min(3, nrow(gene_summary_ordered))]
  
  for (i in 1:nrow(top_3)) {
    gene <- top_3[i]
    cat(sprintf("%d. %s - P=%.2e (Chr%d, %d SNPs cis)\n",
                i, gene$gene, gene$min_p, gene$chromosome, gene$n_cis_snps))
    
    # Mostrar top SNP
    if (gene$gene %in% names(colocalization_results)) {
      top_snp_info <- colocalization_results[[gene$gene]]$top_snps[1]
      cat(sprintf("   Top SNP: %s (P=%.2e, dist=%.0f kb)\n",
                  top_snp_info$SNP, top_snp_info$P, 
                  top_snp_info$distance_to_gene/1000))
    }
  }
  
} else {
  cat("âš ï¸ No se encontraron resultados de colocalizaciÃ³n\n")
}

# 5. GUARDAR RESULTADOS
cat("\n5. Guardando resultados...\n")

# Guardar resumen de genes
fwrite(gene_summary, "colocalization_gene_summary.csv")
cat("  âœ… Resumen guardado: colocalization_gene_summary.csv\n")

# Guardar detalles de top SNPs para cada gen
for (gene_name in names(colocalization_results)) {
  top_snps_file <- paste0("colocalization_", gene_name, "_top_snps.csv")
  fwrite(colocalization_results[[gene_name]]$top_snps, top_snps_file)
}

cat("  âœ… Detalles de SNPs guardados para cada gen\n")

# 6. ANÃLISIS FUNCIONAL BÃSICO
cat("\n6. ANÃLISIS FUNCIONAL\n")
cat("====================\n")

if (nrow(gene_summary) > 0) {
  # Genes con evidencia fuerte/moderada
  significant_genes <- gene_summary[evidence_level %in% c("Strong", "Moderate")]
  
  if (nrow(significant_genes) > 0) {
    cat("Genes con evidencia significativa para estudios funcionales:\n")
    
    for (i in 1:nrow(significant_genes)) {
      gene_name <- significant_genes[i]$gene
      original_info <- mirror_genes[gene_symbol == gene_name]
      
      cat(sprintf("ðŸ”¬ %s:\n", gene_name))
      cat(sprintf("   FunciÃ³n: %s\n", original_info$function_category))
      cat(sprintf("   Evidencia previa: %s\n", original_info$evidence_level))
      cat(sprintf("   P-value mÃ­nimo: %.2e\n", significant_genes[i]$min_p))
      cat("   RecomendaciÃ³n: Priorizar para validaciÃ³n experimental\n\n")
    }
  } else {
    cat("Sin genes con evidencia fuerte/moderada identificados\n")
    cat("RecomendaciÃ³n: Expandir ventana cis o usar threshold menos conservador\n")
  }
}

cat("=== ANÃLISIS DE COLOCALIZACIÃ“N COMPLETADO ===\n")
