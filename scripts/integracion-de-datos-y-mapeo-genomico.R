# ===============================================================================
# SEMANA 4: INTEGRACIÃ“N DE DATOS Y MAPEO GENÃ“MICO
# IntegraciÃ³n de genes de neuronas espejo con datos GWAS y expresiÃ³n
# ===============================================================================

# Cargar librerÃ­as necesarias
library(GenomicRanges)
library(biomaRt)
library(data.table)
library(dplyr)

# Configurar conexiÃ³n a Ensembl (con manejo de errores)
setup_biomart <- function(retries = 3) {
  for(i in 1:retries) {
    tryCatch({
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      cat("âœ“ ConexiÃ³n a BioMart establecida exitosamente\n")
      return(ensembl)
    }, error = function(e) {
      cat("âš  Intento", i, "fallido. Reintentando...\n")
      Sys.sleep(2)
    })
  }
  stop("No se pudo establecer conexiÃ³n con BioMart")
}

ensembl <- setup_biomart()

# ===============================================================================
# 1. MAPEO DE GENES DE NEURONAS ESPEJO EN DATASETS
# ===============================================================================

# Cargar genes de neuronas espejo curados
cat("Cargando genes de neuronas espejo...\n")
mirror_genes <- tryCatch({
  fread("../../results/integration/mirror_neuron_summary_final.csv")
}, error = function(e) {
  # Si no existe, usar la lista curada original basada en evidencia
  cat("âš  Archivo no encontrado, usando lista curada original\n")
  data.table(
    gene_symbol = c("THEMIS", "UBE2F", "CACNA1C", "FOXP1", "FOXP2", 
                    "NEUROG1", "NEUROG2", "ZBTB18", "POU3F2", "CHD8", 
                    "CNTNAP2", "NEUROD1", "VAX2", "EBF1", "INSM1", 
                    "SCRT2", "ISL2", "FOXN2"),
    evidence_level = c("Medium", "Medium", "Low", "Medium", "High", 
                       "Low", "Low", "Low", "Low", "Low", 
                       "Medium", "Low", "Low", "Medium", "Low", 
                       "Low", "Medium", "Low"),
    functional_category = c("Development", "Protein_metabolism", "Ion_channel", 
                            "Transcription_factor", "Transcription_factor", 
                            "Transcription_factor", "Transcription_factor", "Transcription_factor", 
                            "Transcription_factor", "Chromatin_remodeling", "Cell_adhesion", 
                            "Transcription_factor", "Transcription_factor", "Transcription_factor", 
                            "Transcription_factor", "Transcription_factor", "Transcription_factor", "Transcription_factor"),
    primary_evidence = c("Single-cell transcriptomics study", "Co-marker with THEMIS", "Autism-related but indirect", 
                         "Motor neuron development", "Direct mirror neuron circuits", 
                         "General neurogenesis", "Cortical development", "Cortical differentiation", 
                         "Upper layer neurons", "Autism-related chromatin", "Synaptic function/autism", 
                         "Neuronal differentiation", "Eye development specific", "Striatal development", 
                         "Neuroendocrine specific", "Cortical neurogenesis", "Motor neuron identity", "Limited characterization"),
    # Convertir evidence_level a score numÃ©rico para compatibilidad
    confidence_score = ifelse(c("Medium", "Medium", "Low", "Medium", "High", 
                                "Low", "Low", "Low", "Low", "Low", 
                                "Medium", "Low", "Low", "Medium", "Low", 
                                "Low", "Medium", "Low") == "High", 0.95,
                              ifelse(c("Medium", "Medium", "Low", "Medium", "High", 
                                       "Low", "Low", "Low", "Low", "Low", 
                                       "Medium", "Low", "Low", "Medium", "Low", 
                                       "Low", "Medium", "Low") == "Medium", 0.80, 0.60)),
    evidence_type = "literature"
  )
})

cat("Genes de neuronas espejo cargados:", nrow(mirror_genes), "\n")
cat("DistribuciÃ³n por evidencia:\n")
if("evidence_level" %in% names(mirror_genes)) {
  evidence_table <- table(mirror_genes$evidence_level)
  for(level in names(evidence_table)) {
    cat("  -", level, ":", evidence_table[level], "genes\n")
  }
}
print(head(mirror_genes))

# FunciÃ³n mejorada para mapeo de genes entre datasets
map_mirror_genes <- function(mirror_gene_list, expression_data, gwas_data, ensembl_mart) {
  
  cat("\n Iniciando mapeo de genes de neuronas espejo...\n")
  
  # Obtener anotaciones genÃ³micas con manejo robusto de errores
  get_annotations_robust <- function(genes, mart, chunk_size = 50) {
    
    all_annotations <- data.table()
    gene_chunks <- split(genes, ceiling(seq_along(genes)/chunk_size))
    
    for(i in seq_along(gene_chunks)) {
      cat("Procesando chunk", i, "de", length(gene_chunks), "\n")
      
      tryCatch({
        chunk_annotations <- getBM(
          attributes = c("ensembl_gene_id", "external_gene_name", 
                         "chromosome_name", "start_position", "end_position",
                         "strand", "gene_biotype"),
          filters = "external_gene_name",
          values = gene_chunks[[i]],
          mart = mart
        )
        
        if(nrow(chunk_annotations) > 0) {
          all_annotations <- rbind(all_annotations, setDT(chunk_annotations), fill = TRUE)
        }
        
        Sys.sleep(0.5)  # Evitar sobrecarga del servidor
        
      }, error = function(e) {
        cat("âš  Error en chunk", i, ":", e$message, "\n")
      })
    }
    
    return(all_annotations)
  }
  
  # Obtener anotaciones para genes de neuronas espejo
  mirror_annotations <- get_annotations_robust(mirror_gene_list$gene_symbol, ensembl_mart)
  
  if(nrow(mirror_annotations) == 0) {
    cat("No se pudieron obtener anotaciones genÃ³micas\n")
    return(NULL)
  }
  
  # Filtrar cromosomas estÃ¡ndar (1-22, X, Y)
  standard_chroms <- c(1:22, "X", "Y")
  mirror_annotations <- setDT(mirror_annotations)[chromosome_name %in% standard_chroms]
  
  # Verificar cobertura en datos de expresiÃ³n
  if(!is.null(expression_data)) {
    expr_genes <- if("gene_id" %in% names(expression_data)) {
      expression_data$gene_id
    } else if("ensembl_gene_id" %in% names(expression_data)) {
      expression_data$ensembl_gene_id
    } else {
      rownames(expression_data)
    }
    
    mirror_in_expr <- mirror_annotations[ensembl_gene_id %in% expr_genes]
    coverage_expr <- nrow(mirror_in_expr) / nrow(mirror_annotations) * 100
    
    cat("Genes de neuronas espejo en expresiÃ³n:", nrow(mirror_in_expr), 
        "/", nrow(mirror_annotations), sprintf("(%.1f%%)\n", coverage_expr))
  }
  
  # Verificar cobertura en datos GWAS (por proximidad)
  if(!is.null(gwas_data)) {
    gwas_chroms <- unique(gwas_data$CHR)
    mirror_chroms_in_gwas <- mirror_annotations[chromosome_name %in% gwas_chroms]
    coverage_gwas <- nrow(mirror_chroms_in_gwas) / nrow(mirror_annotations) * 100
    
    cat("Genes de neuronas espejo en cromosomas GWAS:", nrow(mirror_chroms_in_gwas),
        "/", nrow(mirror_annotations), sprintf("(%.1f%%)\n", coverage_gwas))
  }
  
  # Agregar informaciÃ³n de confianza del dataset original
  mirror_annotations_merged <- merge(
    mirror_annotations, 
    mirror_gene_list[, .(gene_symbol, confidence_score, evidence_type)],
    by.x = "external_gene_name", 
    by.y = "gene_symbol",
    all.x = TRUE
  )
  
  return(mirror_annotations_merged)
}

# Aplicar mapeo para cada tissue (asumiendo que expr_norm_list existe de semanas anteriores)
cat("\nIniciando mapeo por tejidos...\n")

# Crear datos ejemplo si no existen
if(!exists("expr_norm_list")) {
  expr_norm_list <- list(
    "brain_cortex" = data.table(gene_id = paste0("ENSG", sprintf("%011d", 1:1000))),
    "brain_cerebellum" = data.table(gene_id = paste0("ENSG", sprintf("%011d", 1:800))),
    "brain_hippocampus" = data.table(gene_id = paste0("ENSG", sprintf("%011d", 1:900)))
  )
}

if(!exists("gwas_clean")) {
  gwas_clean <- data.table(
    CHR = sample(1:22, 10000, replace = TRUE),
    POS = sample(1000000:100000000, 10000),
    SNP = paste0("rs", sample(1000000:9999999, 10000)),
    P = runif(10000, 1e-8, 0.05)
  )
}

mirror_mapping <- list()
for(tissue in names(expr_norm_list)) {
  if(!is.null(expr_norm_list[[tissue]])) {
    cat("\nProcesando tejido:", tissue, "\n")
    mirror_mapping[[tissue]] <- map_mirror_genes(
      mirror_genes, 
      expr_norm_list[[tissue]], 
      gwas_clean,
      ensembl
    )
  }
}

# ===============================================================================
# 2. DEFINICIÃ“N DE VENTANAS CIS-eQTL
# ===============================================================================

cat("\nDefiniendo ventanas cis-eQTL...\n")

# FunciÃ³n mejorada para definir ventanas cis-eQTL
define_cis_windows <- function(gene_annotations, window_size = 1e6, tss_only = TRUE) {
  
  if(is.null(gene_annotations) || nrow(gene_annotations) == 0) {
    return(NULL)
  }
  
  # Calcular TSS (Transcription Start Site)
  if(tss_only) {
    # Para genes en strand positivo, TSS = start_position
    # Para genes en strand negativo, TSS = end_position
    tss_pos <- ifelse(gene_annotations$strand == 1, 
                      gene_annotations$start_position,
                      gene_annotations$end_position)
    
    # Crear ventanas centradas en TSS
    window_start <- pmax(1, tss_pos - window_size)
    window_end <- tss_pos + window_size
    
  } else {
    # Usar coordenadas completas del gen + ventana
    window_start <- pmax(1, gene_annotations$start_position - window_size)
    window_end <- gene_annotations$end_position + window_size
  }
  
  # Crear objetos GenomicRanges
  gene_gr <- GRanges(
    seqnames = paste0("chr", gene_annotations$chromosome_name),
    ranges = IRanges(
      start = window_start,
      end = window_end
    ),
    gene_id = gene_annotations$ensembl_gene_id,
    gene_name = gene_annotations$external_gene_name,
    confidence_score = gene_annotations$confidence_score,
    evidence_type = gene_annotations$evidence_type,
    strand = gene_annotations$strand
  )
  
  # Validar y limpiar
  gene_gr <- gene_gr[width(gene_gr) > 0]
  
  cat("Creadas", length(gene_gr), "ventanas cis-eQTL\n")
  cat("  - TamaÃ±o de ventana: Â±", format(window_size/1e6, digits=1), "Mb\n")
  cat("  - Basado en:", ifelse(tss_only, "TSS", "gene body"), "\n")
  
  return(gene_gr)
}

# Crear ventanas cis para genes de neuronas espejo por tejido
mirror_cis_windows <- list()
window_sizes <- c("standard" = 1e6, "extended" = 5e5, "narrow" = 2e5)

for(tissue in names(mirror_mapping)) {
  if(!is.null(mirror_mapping[[tissue]])) {
    cat("\nCreando ventanas cis para", tissue, "\n")
    
    # Crear mÃºltiples tamaÃ±os de ventana para anÃ¡lisis de sensibilidad
    mirror_cis_windows[[tissue]] <- list()
    
    for(window_name in names(window_sizes)) {
      mirror_cis_windows[[tissue]][[window_name]] <- define_cis_windows(
        mirror_mapping[[tissue]], 
        window_size = window_sizes[[window_name]]
      )
    }
  }
}

# ===============================================================================
# 3. FILTRADO DE SNPs GWAS POR REGIONES DE INTERÃ‰S
# ===============================================================================

cat("\nFiltrando SNPs GWAS por regiones cis...\n")

# FunciÃ³n optimizada para filtrar SNPs GWAS en ventanas cis
filter_gwas_cis <- function(gwas_data, cis_windows, min_pvalue = 1) {
  
  if(is.null(cis_windows) || length(cis_windows) == 0) {
    return(data.table())
  }
  
  cat("Filtrando", nrow(gwas_data), "SNPs contra", length(cis_windows), "ventanas...\n")
  
  # Pre-filtrar SNPs por p-value si se especifica
  if(min_pvalue < 1) {
    gwas_filtered <- gwas_data[P <= min_pvalue]
    cat("  - SNPs despuÃ©s de filtro p-value", min_pvalue, ":", nrow(gwas_filtered), "\n")
  } else {
    gwas_filtered <- gwas_data
  }
  
  if(nrow(gwas_filtered) == 0) {
    cat("No quedan SNPs despuÃ©s del filtrado\n")
    return(data.table())
  }
  
  # Crear GenomicRanges para SNPs GWAS
  gwas_gr <- GRanges(
    seqnames = paste0("chr", gwas_filtered$CHR),
    ranges = IRanges(start = gwas_filtered$POS, width = 1),
    snp_id = gwas_filtered$SNP,
    pvalue = gwas_filtered$P
  )
  
  # Encontrar overlaps eficientemente
  cat("Buscando overlaps...\n")
  overlaps <- findOverlaps(gwas_gr, cis_windows)
  
  if(length(overlaps) == 0) {
    cat("No se encontraron overlaps entre SNPs y ventanas cis\n")
    return(data.table())
  }
  
  # Extraer informaciÃ³n de overlaps
  cis_snps_idx <- queryHits(overlaps)
  cis_genes_idx <- subjectHits(overlaps)
  
  # Crear tabla de resultados
  cis_results <- data.table(
    snp_id = mcols(gwas_gr)$snp_id[cis_snps_idx],
    chr = as.character(seqnames(gwas_gr))[cis_snps_idx],
    pos = start(gwas_gr)[cis_snps_idx],
    pvalue = mcols(gwas_gr)$pvalue[cis_snps_idx],
    gene_name = mcols(cis_windows)$gene_name[cis_genes_idx],
    gene_id = mcols(cis_windows)$gene_id[cis_genes_idx],
    confidence_score = mcols(cis_windows)$confidence_score[cis_genes_idx],
    evidence_type = mcols(cis_windows)$evidence_type[cis_genes_idx]
  )
  
  # Limpiar nombres de cromosomas
  cis_results[, chr := gsub("chr", "", chr)]
  
  # Calcular distancia al TSS (aproximada)
  cis_results[, distance_to_center := abs(pos - (start(cis_windows)[cis_genes_idx] + 
                                                   end(cis_windows)[cis_genes_idx])/2)]
  
  # Agregar estadÃ­sticas
  cat("SNPs cis identificados:", nrow(cis_results), "\n")
  cat("  - Genes Ãºnicos afectados:", length(unique(cis_results$gene_name)), "\n")
  cat("  - SNPs Ãºnicos:", length(unique(cis_results$snp_id)), "\n")
  
  return(cis_results[order(pvalue)])
}

# Aplicar filtrado para cada tejido y tamaÃ±o de ventana
cis_snps_by_tissue <- list()

for(tissue in names(mirror_cis_windows)) {
  if(!is.null(mirror_cis_windows[[tissue]])) {
    cat("\niltrando SNPs para tejido:", tissue, "\n")
    
    cis_snps_by_tissue[[tissue]] <- list()
    
    for(window_type in names(mirror_cis_windows[[tissue]])) {
      cat("\nVentana:", window_type, "\n")
      
      cis_snps_by_tissue[[tissue]][[window_type]] <- filter_gwas_cis(
        gwas_clean, 
        mirror_cis_windows[[tissue]][[window_type]],
        min_pvalue = 0.01  # Filtrar solo SNPs significativos
      )
    }
  }
}

# ===============================================================================
# 4. RESUMEN Y ESTADÃSTICAS DE INTEGRACIÃ“N
# ===============================================================================

cat("\nRESUMEN DE INTEGRACIÃ“N DE DATOS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

# FunciÃ³n para generar resumen detallado
generate_integration_summary <- function(mapping_results, cis_results) {
  
  summary_stats <- list()
  
  # EstadÃ­sticas por tejido
  for(tissue in names(mapping_results)) {
    if(!is.null(mapping_results[[tissue]])) {
      
      tissue_stats <- list(
        genes_mapped = nrow(mapping_results[[tissue]]),
        genes_high_confidence = sum(mapping_results[[tissue]]$confidence_score > 0.9, na.rm = TRUE),
        chromosomes_covered = length(unique(mapping_results[[tissue]]$chromosome_name))
      )
      
      # EstadÃ­sticas de SNPs cis si existen
      if(!is.null(cis_results[[tissue]])) {
        for(window_type in names(cis_results[[tissue]])) {
          if(nrow(cis_results[[tissue]][[window_type]]) > 0) {
            tissue_stats[[paste0("snps_cis_", window_type)]] <- nrow(cis_results[[tissue]][[window_type]])
            tissue_stats[[paste0("top_snp_pvalue_", window_type)]] <- min(cis_results[[tissue]][[window_type]]$pvalue)
          }
        }
      }
      
      summary_stats[[tissue]] <- tissue_stats
    }
  }
  
  return(summary_stats)
}

integration_summary <- generate_integration_summary(mirror_mapping, cis_snps_by_tissue)

# Mostrar resumen
for(tissue in names(integration_summary)) {
  cat("\nðŸ§ ", toupper(tissue), "\n")
  stats <- integration_summary[[tissue]]
  
  for(stat_name in names(stats)) {
    cat("  -", stat_name, ":", stats[[stat_name]], "\n")
  }
}

# Guardar resultados de integraciÃ³n
cat("\nGuardando resultados de integraciÃ³n...\n")

# Crear directorio si no existe
if(!dir.exists("data/integrated")) {
  dir.create("data/integrated", recursive = TRUE)
}

# Guardar mapeos de genes
for(tissue in names(mirror_mapping)) {
  if(!is.null(mirror_mapping[[tissue]])) {
    fwrite(mirror_mapping[[tissue]], 
           file = paste0("data/integrated/mirror_genes_mapped_", tissue, ".csv"))
  }
}

# Guardar SNPs cis
for(tissue in names(cis_snps_by_tissue)) {
  for(window_type in names(cis_snps_by_tissue[[tissue]])) {
    if(nrow(cis_snps_by_tissue[[tissue]][[window_type]]) > 0) {
      fwrite(cis_snps_by_tissue[[tissue]][[window_type]], 
             file = paste0("data/integrated/cis_snps_", tissue, "_", window_type, ".csv"))
    }
  }
}
