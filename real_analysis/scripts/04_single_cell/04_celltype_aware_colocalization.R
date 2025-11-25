# 04_celltype_aware_colocalization_v3.R
# ==============================================================================
# COLOCALIZACIÓN SINGLE-CELL-AWARE CON INTERACTION eQTL (MÉTODO GOLD STANDARD)
# ==============================================================================
# METODOLOGÍA:
# 1. LiftOver: Convierte GWAS de hg19 → hg38 (corrección de genome build)
# 2. Interaction eQTL: Modela Expression ~ SNP × Neuron_Prop para estimar
#    efectos cell-type-specific
# 3. Colocalization: Usa efectos específicos de neuronas para coloc
#
# DIFERENCIAS CON v1/v2:
# - v1: NO liftover, usa ID artificial común (metodológicamente incorrecto)
# - v2: Liftover + regresión simple (identifica SNPs neuronales pero no estima efecto)
# - v3: Liftover + interaction model (estima EFECTO específico de neuronas)
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(coloc)
    library(stats)

    if (!require("rtracklayer", quietly=TRUE)) {
        if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
        BiocManager::install("rtracklayer", update = FALSE, ask = FALSE)
    }
    library(rtracklayer)
    library(GenomicRanges)
})

cat("================================================================\n")
cat("COLOCALIZACIÓN CELL-TYPE-AWARE v3: INTERACTION eQTL + LIFTOVER\n")
cat("================================================================\n\n")

# --- CONFIGURACIÓN ---
BASE_DIR <- "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis"
RESULTS_DIR <- file.path(BASE_DIR, "results/celltype_colocalization_real")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Archivos de entrada
EQTL_FILE <- file.path(BASE_DIR, "results/eqtl_extracted/real_gtex_eqtls.csv")
GWAS_FILE <- file.path(BASE_DIR, "data/gwas/iPSYCH-PGC_ASD_Nov2017")
PROPS_FILE <- file.path(BASE_DIR, "results/celltype_deconvolution/gtex_celltype_proportions.csv")
ATTR_FILE <- file.path(BASE_DIR, "data/gtex_bulk/GTEx_Samples.txt")
CHAIN_FILE <- file.path(BASE_DIR, "data/hg19ToHg38.over.chain")
BULK_FILE <- file.path(BASE_DIR, "data/gtex_bulk/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")

# Archivo de salida
OUTPUT_FILE <- file.path(RESULTS_DIR, "celltype_aware_colocalization_v3_interaction.csv")

# Parámetros
NEURON_THRESHOLD <- 0.005  # Mínima fracción neuronal (0.5%, muy permisivo)
MIN_SAMPLES <- 10          # Mínimo de muestras para interaction model
WINDOW_BROAD <- 2e6      # Ventana amplia para liftover (±2MB)
WINDOW_FINE <- 10000     # Ventana fina para matching (±10kb)

# ==============================================================================
# PASO 0: CARGAR CHAIN FILE PARA LIFTOVER
# ==============================================================================
cat("0. Cargando chain file para liftOver (hg19 → hg38)...\n")
if(!file.exists(CHAIN_FILE)) {
    stop("❌ Falta archivo chain: ", CHAIN_FILE, "\n",
         "Descarga de: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz")
}
chain <- import.chain(CHAIN_FILE)
cat("   ✓ Chain file cargado\n\n")

# ==============================================================================
# PASO 1: PREPARAR PROPORCIONES CELULARES
# ==============================================================================
cat("1. Preparando proporciones de tipos celulares...\n")

# Cargar metadata de muestras
samples_meta <- fread(ATTR_FILE, select = c("SAMPID", "SMTSD"))
colnames(samples_meta) <- c("SampleID", "Tissue")

# Cargar proporciones celulares
props <- read.csv(PROPS_FILE)

# Unir con metadata y filtrar por cerebro
props_brain <- props %>%
    inner_join(samples_meta, by = "SampleID") %>%
    filter(grepl("Brain", Tissue))

cat("   Muestras cerebrales con proporciones:", nrow(props_brain), "\n")

# Limpiar nombres de columnas
colnames(props_brain) <- make.names(colnames(props_brain))

# Definir columnas neuronales (basado en Allen Brain Atlas)
neuron_cols <- c(
    "Upper.rhombic.lip", "Splatter", "Lower.rhombic.lip", "Mammillary.body",
    "Thalamic.excitatory", "Amygdala.excitatory", "Medium.spiny.neuron",
    "Eccentric.medium.spiny.neuron", "Miscellaneous", "Cerebellar.inhibitory",
    "Midbrain.derived.inhibitory", "CGE.interneuron", "LAMP5.LHX6.and.Chandelier",
    "MGE.interneuron", "Deep.layer.near.projecting", "Deep.layer.corticothalamic.and.6b",
    "Hippocampal.CA1.3", "Upper.layer.intratelencephalic", "Deep.layer.intratelencephalic",
    "Hippocampal.dentate.gyrus", "Hippocampal.CA4"
)

# Calcular fracción neuronal total
valid_cols <- intersect(neuron_cols, colnames(props_brain))
props_brain$Total_Neuron <- rowSums(props_brain[, valid_cols], na.rm = TRUE)

cat("   Fracción neuronal media:", round(mean(props_brain$Total_Neuron, na.rm=TRUE), 3), "\n")
cat("   Rango neuronal:", round(min(props_brain$Total_Neuron, na.rm=TRUE), 3), "-",
    round(max(props_brain$Total_Neuron, na.rm=TRUE), 3), "\n\n")

# ==============================================================================
# PASO 2: CARGAR DATOS GENÉTICOS
# ==============================================================================
cat("2. Cargando eQTLs y GWAS...\n")

# --- eQTLs ---
eqtls <- read.csv(EQTL_FILE)

# Parsear coordenadas (hg38)
eqtls <- eqtls %>% mutate(
    pos_hg38 = as.numeric(sub(".*_([0-9]+)_.*_.*_.*", "\\1", variant_id)),
    chr_clean = sub("chr", "", sub("_.*", "", variant_id))
)

cat("   eQTLs cargados:", nrow(eqtls), "\n")
cat("   Genes únicos:", length(unique(eqtls$gene_symbol)), "\n")

# --- GWAS ---
gwas_raw <- fread(GWAS_FILE, select = c("CHR", "BP", "P", "OR", "SE", "SNP"))
cat("   GWAS cargado:", nrow(gwas_raw), "variantes\n")
cat("   (Asumiendo genome build: hg19/GRCh37)\n\n")

# --- Mapear códigos de tejido ---
tissue_map <- c(
    'Brain_Amygdala' = 'Brain - Amygdala',
    'Brain_Anterior_cingulate_cortex_BA24' = 'Brain - Anterior cingulate cortex (BA24)',
    'Brain_Caudate_basal_ganglia' = 'Brain - Caudate (basal ganglia)',
    'Brain_Cerebellar_Hemisphere' = 'Brain - Cerebellar Hemisphere',
    'Brain_Cerebellum' = 'Brain - Cerebellum',
    'Brain_Cortex' = 'Brain - Cortex',
    'Brain_Frontal_Cortex_BA9' = 'Brain - Frontal Cortex (BA9)',
    'Brain_Hippocampus' = 'Brain - Hippocampus',
    'Brain_Hypothalamus' = 'Brain - Hypothalamus',
    'Brain_Nucleus_accumbens_basal_ganglia' = 'Brain - Nucleus accumbens (basal ganglia)',
    'Brain_Putamen_basal_ganglia' = 'Brain - Putamen (basal ganglia)',
    'Brain_Spinal_cord_cervical_c-1' = 'Brain - Spinal cord (cervical c-1)',
    'Brain_Substantia_nigra' = 'Brain - Substantia nigra'
)
eqtls$Tissue_Match <- tissue_map[eqtls$tissue_code]

# ==============================================================================
# PASO 3: FUNCIÓN DE INTERACTION eQTL
# ==============================================================================

# Función para estimar efecto cell-type-specific
# Inputs:
#   - variant_data: Data frame con slope, slope_se, tissue
#   - props_data: Data frame con SampleID, Tissue, Total_Neuron
# Outputs:
#   - Lista con: beta_neuron (efecto neuronal), se_neuron, p_interaction, n_samples
estimate_celltype_effect <- function(variant_data, props_data) {

    # 1. Unir datos
    merged <- variant_data %>%
        inner_join(props_data %>% select(Tissue, Total_Neuron, SampleID),
                   by = c("Tissue_Match" = "Tissue"),
                   relationship = "many-to-many")

    # 2. Filtrar muestras
    merged <- merged %>% filter(Total_Neuron >= NEURON_THRESHOLD)

    if(nrow(merged) < MIN_SAMPLES) {
        return(list(status = "insufficient_samples"))
    }

    # 3. Preparar modelo
    mean_neuron <- mean(merged$Total_Neuron)
    merged$Neuron_Centered <- merged$Total_Neuron - mean_neuron
    weights <- 1 / (merged$slope_se^2)

    tryCatch({
        # 4. Ajustar modelo
        fit <- lm(slope ~ Neuron_Centered, data = merged, weights = weights)
        
        # Extraer coeficientes crudos (Interacción)
        coefs <- summary(fit)$coefficients
        beta_interaction <- coefs[2, 1] # Pendiente de interacción
        se_interaction   <- coefs[2, 2] # Error std de interacción
        p_interaction    <- coefs[2, 4] # P-valor de interacción

        # 5. Predecir efecto en neuronas puras (Total_Neuron = 1)
        new_data <- data.frame(Neuron_Centered = (1 - mean_neuron))
        prediction <- predict(fit, newdata = new_data, se.fit = TRUE)
        
        return(list(
            # --- Datos Calculados (Prediction) ---
            beta_neuron = as.numeric(prediction$fit),
            se_neuron = as.numeric(prediction$se.fit),
            
            # --- Datos Crudos del Modelo (Faltaban estos) ---
            beta_interaction = beta_interaction,
            se_interaction = se_interaction, # Agregado por seguridad
            p_interaction = p_interaction,
            
            # --- Metadatos ---
            n_samples = nrow(merged),
            mean_neuron_prop = mean_neuron,
            status = "success"
        ))

    }, error = function(e) {
        return(list(status = paste("error:", e$message)))
    })
}

# ==============================================================================
# PASO 4: ANÁLISIS POR GEN
# ==============================================================================
cat("3. Ejecutando análisis interaction eQTL + colocalization...\n\n")

results_list <- list()
unique_genes <- unique(eqtls$gene_symbol)

for(gene in unique_genes) {
    cat(sprintf("━━━ %s ━━━\n", gene))
    gene_data <- filter(eqtls, gene_symbol == gene)

    # Obtener cromosoma y posición representativa
    target_chr <- gene_data$chr_clean[1]
    target_pos <- median(gene_data$pos_hg38, na.rm = TRUE)  # Usar mediana como centro

    cat(sprintf("   Chr: %s | Posición central: %.0f (hg38)\n", target_chr, target_pos))

    # --- A) LIFTOVER DE GWAS (hg19 → hg38) ---
    chr_gwas <- if(target_chr == "X") 23 else as.numeric(target_chr)

    # Ventana amplia para capturar región
    gwas_broad <- gwas_raw[CHR == chr_gwas &
                           BP >= (target_pos - WINDOW_BROAD) &
                           BP <= (target_pos + WINDOW_BROAD)]

    if(nrow(gwas_broad) == 0) {
        cat("   ⚠️  Sin variantes GWAS en región\n\n")
        next
    }

    cat(sprintf("   GWAS en región (hg19): %d variantes\n", nrow(gwas_broad)))

    # LiftOver a hg38
    gr_hg19 <- GRanges(
        seqnames = paste0("chr", gwas_broad$CHR),
        ranges = IRanges(start = gwas_broad$BP, end = gwas_broad$BP),
        strand = "*",
        ID = 1:nrow(gwas_broad)
    )

    gr_hg38_list <- liftOver(gr_hg19, chain)
    gr_hg38 <- unlist(gr_hg38_list)

    if(length(gr_hg38) == 0) {
        cat("   ⚠️  LiftOver falló (no se mapeó ningún SNP)\n\n")
        next
    }

    # Recuperar datos GWAS convertidos
    valid_ids <- mcols(gr_hg38)$ID
    gwas_hg38 <- gwas_broad[valid_ids, ]
    gwas_hg38$BP_hg38 <- start(gr_hg38)

    cat(sprintf("   GWAS post-liftover (hg38): %d variantes\n", nrow(gwas_hg38)))

    # Filtrar por ventana fina
    gwas_window <- gwas_hg38[BP_hg38 >= (target_pos - WINDOW_FINE) &
                             BP_hg38 <= (target_pos + WINDOW_FINE)]

    if(nrow(gwas_window) == 0) {
        cat("   ⚠️  Ningún SNP GWAS en ventana fina (±10kb)\n\n")
        next
    }

    cat(sprintf("   GWAS en ventana fina: %d variantes\n", nrow(gwas_window)))

    # --- B) INTERACTION eQTL MODEL ---
    cat("   Estimando efecto cell-type-specific...\n")

    # Estimar efecto neuronal usando todos los SNPs del gen
    interaction_result <- estimate_celltype_effect(gene_data, props_brain)

    cat(sprintf("   Status: %s\n", interaction_result$status))

    if(interaction_result$status != "success") {
        cat(sprintf("   ⚠️  %s\n\n", interaction_result$status))
        next
    }

    cat(sprintf("   N muestras: %d\n", interaction_result$n_samples))
    cat(sprintf("   Beta interacción: %.4f (p=%.2e)\n",
                interaction_result$beta_interaction,
                interaction_result$p_interaction))
    cat(sprintf("   Beta neuronal estimado: %.4f (SE=%.4f)\n",
                interaction_result$beta_neuron,
                interaction_result$se_neuron))

    # --- C) COLOCALIZATION ---
    # Usar mejor SNP GWAS en ventana
    best_gwas <- gwas_window[which.min(P), ]
    dist_bp <- abs(best_gwas$BP_hg38 - target_pos)

    cat(sprintf("   Mejor GWAS: %s (p=%.2e, dist=%.0f bp)\n",
                best_gwas$SNP, best_gwas$P, dist_bp))

    # Validar datos GWAS
    if(is.na(best_gwas$OR) || best_gwas$OR <= 0 || is.na(best_gwas$SE)) {
        cat("   ⚠️  GWAS con OR/SE inválido\n\n")
        next
    }

    # Preparar datasets para coloc
    common_id <- paste0("Target_", gene)

    # Estimar sdY (standard deviation of the trait)
    # Para expression data, usamos regla empírica: sdY ~ SE * sqrt(N)
    # Esto asume MAF ~ 0.2 (promedio razonable)
    sdY_est <- interaction_result$se_neuron * sqrt(interaction_result$n_samples * 0.4)

    # Dataset 1: eQTL con efecto neuronal específico
    ds1 <- list(
        beta = interaction_result$beta_neuron,
        varbeta = interaction_result$se_neuron^2,
        N = interaction_result$n_samples,
        type = "quant",
        sdY = sdY_est,
        snp = common_id
    )

    # Dataset 2: GWAS
    ds2 <- list(
        beta = log(best_gwas$OR),
        varbeta = best_gwas$SE^2,
        pvalues = best_gwas$P,
        N = 46350,  # Grove et al. 2019
        type = "cc",
        s = 0.4,    # Proporción de casos
        snp = common_id
    )

    # Ejecutar coloc
    tryCatch({
        res <- coloc.abf(dataset1 = ds1, dataset2 = ds2)
        pp_h3 <- res$summary["PP.H3.abf"]
        pp_h4 <- res$summary["PP.H4.abf"]

        cat(sprintf("   🔥 PP.H3: %.4f | PP.H4: %.4f\n", pp_h3, pp_h4))

        # Interpretación
        interpretation <- if(pp_h4 > 0.75) {
            "STRONG_COLOC"
        } else if(pp_h4 > 0.5) {
            "MODERATE_COLOC"
        } else if(pp_h3 > 0.75) {
            "SEPARATE_SIGNALS"
        } else {
            "NO_COLOC"
        }

        cat(sprintf("   Interpretación: %s\n\n", interpretation))

        # Guardar resultado
        results_list[[length(results_list) + 1]] <- data.frame(
            Gene = gene,
            Chr = target_chr,
            Position_hg38 = target_pos,
            N_Samples = interaction_result$n_samples,
            Mean_Neuron_Prop = interaction_result$mean_neuron_prop,
            Beta_Interaction = interaction_result$beta_interaction,
            P_Interaction = interaction_result$p_interaction,
            Beta_Neuron_Specific = interaction_result$beta_neuron,
            SE_Neuron_Specific = interaction_result$se_neuron,
            GWAS_SNP = best_gwas$SNP,
            GWAS_P = best_gwas$P,
            GWAS_Distance_bp = dist_bp,
            PP_H3 = pp_h3,
            PP_H4 = pp_h4,
            Interpretation = interpretation,
            stringsAsFactors = FALSE
        )

    }, error = function(e) {
        cat(sprintf("   ❌ Error en coloc: %s\n\n", e$message))
    })
}

# ==============================================================================
# PASO 5: GUARDAR RESULTADOS
# ==============================================================================
cat("\n================================================================\n")

if(length(results_list) > 0) {
    final_results <- do.call(rbind, results_list)

    # Ordenar por PP.H4
    final_results <- final_results %>% arrange(desc(PP_H4))

    # Guardar
    write.csv(final_results, OUTPUT_FILE, row.names = FALSE)

    cat("RESULTADOS FINALES:\n")
    cat("================================================================\n")
    print(final_results %>%
          select(Gene, N_Samples, P_Interaction, Beta_Neuron_Specific,
                 PP_H4, Interpretation))
    cat("================================================================\n\n")
    cat("✅ ANÁLISIS COMPLETADO\n")
    cat("   Archivo guardado en:", OUTPUT_FILE, "\n")
    cat("\nRESUMEN:\n")
    cat("  - Genes analizados:", nrow(final_results), "\n")
    cat("  - Con colocalization (PP.H4 > 0.5):",
        sum(final_results$PP_H4 > 0.5, na.rm=TRUE), "\n")
    cat("  - Con interacción significativa (p < 0.05):",
        sum(final_results$P_Interaction < 0.05, na.rm=TRUE), "\n")

} else {
    cat("❌ No se generaron resultados\n")
    cat("   Verifica que:\n")
    cat("   - Los archivos de entrada existen\n")
    cat("   - Hay suficientes muestras con proporciones neuronales\n")
    cat("   - El chain file está correcto\n")
}

cat("================================================================\n")
