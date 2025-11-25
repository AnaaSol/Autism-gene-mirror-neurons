# 03_celltype_eqtl_interaction.R
# MODO: VISUALIZACIÓN DIRECTA (Para genes con 1 solo tejido)
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
})

cat("================================================================\n")
cat("ANÁLISIS DE INTERACCIÓN: VISUALIZACIÓN DIRECTA\n")
cat("================================================================\n\n")

BASE_DIR <- "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis"
RESULTS_DIR <- file.path(BASE_DIR, "results/celltype_eqtls")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

PROPS_FILE <- file.path(BASE_DIR, "results/celltype_deconvolution/gtex_celltype_proportions.csv")
EQTL_FILE <- file.path(BASE_DIR, "results/eqtl_extracted/real_gtex_eqtls.csv")
PLOT_FILE <- file.path(RESULTS_DIR, "eqtl_tissue_specificity.pdf")
SUMMARY_FILE <- file.path(RESULTS_DIR, "eqtl_summary.csv")

# 1. CARGAR METADATOS Y PROPORCIONES
cat("1. Preparando datos de tejidos...\n")
ATTR_FILE <- file.path(BASE_DIR, "data/gtex_bulk/GTEx_Samples.txt")
samples_meta <- fread(ATTR_FILE, select = c("SAMPID", "SMTSD"))
colnames(samples_meta) <- c("SampleID", "Tissue")

props <- read.csv(PROPS_FILE)
props_brain <- props %>%
    inner_join(samples_meta, by = "SampleID") %>%
    filter(grepl("Brain", Tissue))

# 2. CALCULAR NEURONAS PROMEDIO
# Nombres limpios para R
colnames(props_brain) <- make.names(colnames(props_brain))

# Columnas de neuronas (Allen Brain Atlas)
neuron_cols <- c(
    "Upper.rhombic.lip", "Splatter", "Lower.rhombic.lip", "Mammillary.body", 
    "Thalamic.excitatory", "Amygdala.excitatory", "Medium.spiny.neuron", 
    "Eccentric.medium.spiny.neuron", "Miscellaneous", "Cerebellar.inhibitory", 
    "Midbrain.derived.inhibitory", "CGE.interneuron", "LAMP5.LHX6.and.Chandelier", 
    "MGE.interneuron", "Deep.layer.near.projecting", "Deep.layer.corticothalamic.and.6b", 
    "Hippocampal.CA1.3", "Upper.layer.intratelencephalic", "Deep.layer.intratelencephalic", 
    "Hippocampal.dentate.gyrus", "Hippocampal.CA4"
)
valid_cols <- intersect(neuron_cols, colnames(props_brain))
props_brain$Total_Neuron <- rowSums(props_brain[, valid_cols], na.rm = TRUE)

tissue_stats <- props_brain %>%
    group_by(Tissue) %>%
    summarise(
        Mean_Neuron = mean(Total_Neuron, na.rm=TRUE),
        N_Samples = n()
    )

# 3. CARGAR eQTLs
cat("2. Cargando eQTLs...\n")
eqtls <- read.csv(EQTL_FILE)

# Mapeo de nombres
eqtls$Tissue_Match <- NA
eqtls$Tissue_Match[eqtls$tissue_code == "Brain_Frontal_Cortex_BA9"] <- "Brain - Frontal Cortex (BA9)"
eqtls$Tissue_Match[eqtls$tissue_code == "Brain_Cerebellar_Hemisphere"] <- "Brain - Cerebellar Hemisphere"
eqtls$Tissue_Match[eqtls$tissue_code == "Brain_Cortex"] <- "Brain - Cortex"
eqtls$Tissue_Match[eqtls$tissue_code == "Brain_Anterior_cingulate_cortex_BA24"] <- "Brain - Anterior cingulate cortex (BA24)"

# Unir
plot_data <- eqtls %>%
    inner_join(tissue_stats, by = c("Tissue_Match" = "Tissue")) %>%
    group_by(gene_symbol, Tissue_Match) %>%
    slice_min(order_by = pval_nominal, n = 1) %>% # Mejor SNP por gen/tejido
    ungroup()

cat("   Puntos a graficar:", nrow(plot_data), "\n")

# 4. GRAFICAR (BARRAS / PUNTOS)
cat("3. Generando gráfico resumen...\n")

pdf(PLOT_FILE, width = 10, height = 7)

# Gráfico 1: ¿Dónde es significativo cada gen y cuántas neuronas hay ahí?
p1 <- ggplot(plot_data, aes(x = gene_symbol, y = abs(slope))) +
    geom_bar(stat = "identity", aes(fill = Mean_Neuron), color="black", width=0.7) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Fracción Neuronal\ndel Tejido") +
    geom_text(aes(label = Tissue_Match), vjust = -0.5, size = 3, angle = 0) +
    labs(
        title = "Fuerza del eQTL y Contexto Neuronal",
        subtitle = "Cada barra es el tejido donde el eQTL es significativo (P < 0.05)",
        x = "Gen de Interés",
        y = "Tamaño del Efecto (|Slope|)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(face="bold"))

print(p1)

# Gráfico 2: Scatter de TODOS los genes juntos (Contexto Global)
# Si juntamos todos los genes, ¿vemos tendencia general?
if(nrow(plot_data) > 2) {
    p2 <- ggplot(plot_data, aes(x = Mean_Neuron, y = abs(slope))) +
        geom_point(aes(color = gene_symbol), size = 5) +
        geom_text(aes(label = gene_symbol), vjust = -1, size=3) +
        geom_smooth(method = "lm", se = FALSE, color = "gray", linetype="dashed", alpha=0.5) +
        labs(
            title = "Tendencia Global: Efecto eQTL vs Abundancia Neuronal",
            subtitle = "Comparación entre distintos genes y sus tejidos preferidos",
            x = "Proporción de Neuronas en el Tejido",
            y = "Fuerza del eQTL (|Slope|)"
        ) +
        theme_bw()
    print(p2)
}

dev.off()
write.csv(plot_data, SUMMARY_FILE, row.names = FALSE)

cat("\n✓ ANÁLISIS COMPLETADO.\n")
cat("  Gráfico guardado en:", PLOT_FILE, "\n")
cat("  Datos guardados en:", SUMMARY_FILE, "\n")
cat("\nCONCLUSIÓN PRELIMINAR:\n")
print(plot_data %>% select(gene_symbol, Tissue_Match, Mean_Neuron, slope))
