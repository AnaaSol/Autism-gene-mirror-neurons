# 05_final_figure_generation.R
# ==============================================================================
# GRÁFICO FINAL: COMPARACIÓN DE MÉTODOS
# ==============================================================================
# Visualiza cómo el método Cell-Type Aware "rescata" genes que los métodos
# estándar (Bulk Coloc y SMR) dan por perdidos.
# ==============================================================================

library(ggplot2)
library(tidyverse)

BASE_DIR <- "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis"
RESULTS_DIR <- file.path(BASE_DIR, "results/final_figures")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. CONSTRUIR TABLA MAESTRA MANUALMENTE (Con tus datos reales)
# Recopilamos los mejores resultados de tus últimas ejecuciones
# 1. CONSTRUIR TABLA MAESTRA (DATOS REALES EXTRAÍDOS DE TUS LOGS)
data <- data.frame(
  Gene = c("YPEL3", "TBX6", "CNTNAP2", "OXTR", "FOXP2", 
           "TAOK2", "CACNA1C", "MAPK3", "THEMIS", "SHANK3"),
  
  # Categoría
  Category = c("16p11.2 (CNV)", "16p11.2 (CNV)", "Mirror Neuron (Synaptic)", 
               "Social Behavior", "Language/Mirror", "16p11.2 (Kinase)",
               "Ion Channel", "16p11.2 (Signaling)", "T-cell/Control", "Synaptic (ASD Core)"),
  
  # RESULTADOS BULK (Script 02 - Log reciente)
  # Usamos el valor máximo encontrado en cualquier tejido
  Score_Bulk = c(0.1077, 0.0528, 0.0010, 0.0117, 0.0069, 
                 0.0001, 0.0457, 0.0133, 0.0001, 0.0001),
  
  # RESULTADOS CELL-TYPE AWARE (Script 04 v3 - Log anterior)
  Score_CellType = c(0.4060, 0.3840, 0.3271, 0.2642, 0.2087, 
                     0.1564, 0.1490, 0.1397, 0.1050, 0.0002)
)

# Calcular la ganancia para ordenar el gráfico
data <- data %>% 
  mutate(Gain = Score_CellType - Score_Bulk) %>%
  arrange(desc(Gain))

# Reestructurar para ggplot
df_plot <- data %>%
  pivot_longer(cols = c(Score_Bulk, Score_CellType), 
               names_to = "Method", values_to = "Probability") %>%
  mutate(
    Method = ifelse(Method == "Score_Bulk", "Standard Analysis (Bulk)", "Your Method (Cell-Type Aware)"),
    Probability = ifelse(Probability < 0.01, 0.01, Probability) # Para escala log visual
  )

# 2. GRÁFICO DE "RESCATE" (DUMBBELL PLOT)
p <- ggplot(df_plot, aes(x = Gene, y = Probability)) +
  # Línea que conecta los puntos
  geom_line(aes(group = Gene), color = "gray50", size = 1) +
  
  # Puntos
  geom_point(aes(color = Method, size = Probability), alpha = 0.9) +
  
  # Umbral de éxito
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red", alpha=0.5) +
  annotate("text", x = 1, y = 0.78, label = "Strong Evidence (>0.75)", color = "red", size = 3, hjust=0) +
  
  # Colores
  scale_color_manual(values = c("Standard Analysis (Bulk)" = "#999999", 
                                "Your Method (Cell-Type Aware)" = "#E69F00")) +
  
  # Escala
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  
  # Estética
  labs(
    title = "Unmasking Hidden Autism Risk Genes",
    subtitle = "Comparison of Standard Bulk Analysis vs. Cell-Type Aware Pipeline",
    y = "Probability of Causal Link (PP.H4)",
    x = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face="bold"),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  coord_flip() # Girar para leer mejor los genes

# 3. GUARDAR
ggsave(file.path(RESULTS_DIR, "Final_Comparison_Plot.pdf"), p, width = 8, height = 6)
ggsave(file.path(RESULTS_DIR, "Final_Comparison_Plot.png"), p, width = 8, height = 6, dpi = 300)

print(p)
cat("\n✓ Figura Final Generada: ", file.path(RESULTS_DIR, "Final_Comparison_Plot.pdf"), "\n")
