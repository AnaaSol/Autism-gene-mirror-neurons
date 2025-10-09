library(data.table)
library(ggplot2)

# Cargar datos GWAS
gwas_grove <- fread("~/Desktop/autism_mirror_neurons/data/raw/Grove2019_ASD_GWAS.txt")

# ExploraciÃ³n bÃ¡sica
print(paste("NÃºmero de SNPs:", nrow(gwas_grove)))
print(paste("Columnas disponibles:", paste(colnames(gwas_grove), collapse=", ")))

# Verificar estructura
head(gwas_grove)
summary(gwas_grove)

# Chequear valores faltantes
missing_summary <- gwas_grove[, lapply(.SD, function(x) sum(is.na(x)))]
print(missing_summary)

# Filtros de calidad estÃ¡ndar para GWAS
gwas_qc <- function(gwas_data) {
  
  # 1. Filtro por MAF (Minor Allele Frequency)
  if("MAF" %in% colnames(gwas_data)) {
    gwas_data <- gwas_data[MAF > 0.01]
    cat("Filtro MAF > 0.01 aplicado\n")
  }
  
  # 2. Filtro por INFO score (calidad de imputaciÃ³n)
  if("INFO" %in% colnames(gwas_data)) {
    gwas_data <- gwas_data[INFO > 0.8]
    cat("Filtro INFO > 0.8 aplicado\n")
  }
  
  # 3. Remover SNPs con HWE muy bajo
  if("HWE_P" %in% colnames(gwas_data)) {
    gwas_data <- gwas_data[HWE_P > 1e-6]
    cat("Filtro HWE_P > 1e-6 aplicado\n")
  }
  
  # 4. Remover valores faltantes en columnas crÃ­ticas
  # CAMBIO AQUÃ: "POS" por "BP"
  required_cols <- c("CHR", "BP", "A1", "A2", "P")
  gwas_data <- gwas_data[complete.cases(gwas_data[, ..required_cols])]
  
  return(gwas_data)
}

# Ejecutar el control de calidad
gwas_clean <- gwas_qc(gwas_grove)
print(paste("SNPs despuÃ©s de QC:", nrow(gwas_clean)))

# Verificar cuÃ¡ntos SNPs se filtraron
print(paste("SNPs originales:", nrow(gwas_grove)))
print(paste("SNPs despuÃ©s de filtro INFO > 0.8:", nrow(gwas_clean)))
print(paste("SNPs removidos:", nrow(gwas_grove) - nrow(gwas_clean)))

# Verificar distribuciÃ³n de valores P
summary(gwas_clean$P)
