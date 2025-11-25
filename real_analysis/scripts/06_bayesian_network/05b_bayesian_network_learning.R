# 05_bayesian_network_learning.R
# ==============================================================================
# INFERENCIA CAUSAL + VISUALIZACIÓN CIRCULAR (CLEAN)
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(bnlearn)
    library(data.table)
    # Librería nueva para gráficos bonitos
    if (!require("igraph")) install.packages("igraph", repos="https://cloud.r-project.org")
    library(igraph)
})

cat("================================================================\n")
cat("REDES BAYESIANAS: VISUALIZACIÓN CIRCULAR\n")
cat("================================================================\n\n")

# --- RUTAS ---
BASE_DIR <- "/home/ana/Desktop/Autism-gene-mirror-neurons/real_analysis"
RESULTS_DIR <- file.path(BASE_DIR, "results/causal_networks")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

EXPR_FILE <- file.path(BASE_DIR, "data/gtex_bulk/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
METADATA_FILE <- file.path(BASE_DIR, "data/gtex_bulk/GTEx_Samples.txt")
EMPIRICAL_LIST_FILE <- file.path(RESULTS_DIR, "empirical_genes_list.txt")

# 1. CARGAR GENES
cat("1. Cargando lista de genes...\n")
if(!file.exists(EMPIRICAL_LIST_FILE)) stop("Falta empirical_genes_list.txt")

gene_list <- readLines(EMPIRICAL_LIST_FILE)
gene_list <- gene_list[gene_list != ""]

# Definir roles para colorear después
target_genes <- c("ALDOA", "CHD8", "CACNA1C", "CNTNAP2", "FOXP2", "MAPK3", "OXTR", "TAOK2", "TBX6", "YPEL3", "KCTD13", "MECP2", "PPP4C", "SHANK3", "THEMIS")
regulator_genes <- setdiff(gene_list, target_genes)

cat("   Total nodos:", length(gene_list), "\n")

# 2. PREPARAR DATOS (FILTRADO RÁPIDO)
cat("2. Cargando expresión (Cerebro + Genes VIP)...\n")

samples_meta <- fread(METADATA_FILE, select = c("SAMPID", "SMTSD"))
brain_samples <- samples_meta$SAMPID[grepl("Brain", samples_meta$SMTSD)]

# Usamos grep para leer solo lo necesario (Ahorro de RAM)
gene_pattern <- paste(gene_list, collapse = "|")
cmd <- paste0("zcat ", EXPR_FILE, " | grep -E 'Name|Description|", gene_pattern, "'")
expr_data <- fread(cmd = cmd) 

expr_data <- expr_data %>% 
    filter(Description %in% gene_list) %>%
    select(-Name) %>% 
    distinct(Description, .keep_all = TRUE)

# Transponer
expr_mat <- t(as.matrix(expr_data[, -1]))
colnames(expr_mat) <- expr_data$Description

# Filtrar muestras
valid_samples <- intersect(rownames(expr_mat), brain_samples)
expr_clean <- as.data.frame(expr_mat[valid_samples, ])
expr_clean <- expr_clean %>% select(any_of(gene_list))

cat("   Matriz lista:", nrow(expr_clean), "muestras x", ncol(expr_clean), "genes\n")

# 3. RESTRICCIONES (BLACKLIST)
cat("3. Aplicando reglas biológicas...\n")
# Prohibido: Target -> Regulador
blacklist <- expand.grid(from = target_genes, to = regulator_genes)
blacklist <- data.frame(from = blacklist$from, to = blacklist$to)
# Filtrar solo genes que realmente existen en la matriz cargada
blacklist <- blacklist[blacklist$from %in% colnames(expr_clean) & 
                       blacklist$to %in% colnames(expr_clean), ]

# 4. APRENDIZAJE (BOOTSTRAP)
cat("4. Aprendiendo estructura (Bootstrap x 50)...\n")
# Usamos bootstrap para encontrar solo las conexiones ROBUSTAS
# Esto limpia mucho el gráfico (elimina ruido)
str_res <- boot.strength(expr_clean, R = 50, algorithm = "hc", 
                         algorithm.args = list(blacklist = blacklist))

# Nos quedamos solo con arcos que aparecen >70% de las veces (Alta confianza)
avg_network <- averaged.network(str_res, threshold = 0.70)

# 5. VISUALIZACIÓN CIRCULAR CON IGRAPH
cat("5. Generando gráfico circular...\n")

# Convertir a objeto igraph
g <- as.igraph(avg_network)

# --- ESTÉTICA DEL GRÁFICO ---

# A) Colores de Nodos
# Creamos un vector de colores
V(g)$color <- ifelse(V(g)$name %in% target_genes, 
                     "#FF7F0E", # Naranja (Autismo Targets)
                     "#ADD8E6") # Azul Claro (Reguladores)

V(g)$label.color <- "black"
V(g)$frame.color <- "gray"
V(g)$size <- 20 # Tamaño del nodo

# B) Estilo de Aristas (Flechas)
# Grosor basado en la fuerza de la conexión (si está disponible)
# Aquí usamos un grosor estándar pero curvamos las líneas
E(g)$color <- "gray50"
E(g)$width <- 1.5
E(g)$arrow.size <- 0.5
E(g)$curved <- 0.2 # Curvatura para ver la dirección

# C) Layout Circular
layout_circ <- layout_in_circle(g)

# GUARDAR PDF
PLOT_FILE <- file.path(RESULTS_DIR, "causal_network_circular.pdf")
pdf(PLOT_FILE, width = 12, height = 12)

par(mar=c(1,1,3,1)) # Márgenes
plot(g, 
     layout = layout_circ,
     main = "Red Regulatoria Causal (Cerebro Humano)\nNaranja: Genes Neuronas Espejo | Azul: Reguladores",
     vertex.label.font = 2,      # Negrita
     vertex.label.cex = 0.8,     # Tamaño letra
     vertex.label.dist = 1.5,    # Etiquetas fuera del nodo
     vertex.label.degree = -pi/2 # Ajuste de posición
)

# Añadir leyenda
legend("topleft", 
       legend = c("Genes Autismo (Targets)", "Reguladores (TFs)"), 
       col = c("#FF7F0E", "#ADD8E6"), 
       pch = 19, pt.cex = 2, bty = "n")

dev.off()

# Exportar tabla de conexiones fuertes para el paper
strong_edges <- str_res[str_res$strength > 0.7 & str_res$direction > 0.5, ]
write.csv(strong_edges, file.path(RESULTS_DIR, "robust_connections.csv"), row.names = FALSE)

cat("\n✓ GRÁFICO CIRCULAR CREADO: ", PLOT_FILE, "\n")
cat("✓ Tabla de conexiones guardada.\n")
