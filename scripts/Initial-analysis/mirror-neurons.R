# Crear lista curada de genes de neuronas espejo basada en revision bibliografica
mirror_genes <- data.frame(
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
  # Agregar informaciÃ³n bibliogrÃ¡fica de soporte
  primary_evidence = c("Single-cell transcriptomics study", "Co-marker with THEMIS", "Autism-related but indirect", 
                       "Motor neuron development", "Direct mirror neuron circuits", 
                       "General neurogenesis", "Cortical development", "Cortical differentiation", 
                       "Upper layer neurons", "Autism-related chromatin", "Synaptic function/autism", 
                       "Neuronal differentiation", "Eye development specific", "Striatal development", 
                       "Neuroendocrine specific", "Cortical neurogenesis", "Motor neuron identity", "Limited characterization"),
  notes = c("Identified as molecular marker in cluster 85", "Protein metabolism regulator", "Calcium channel, indirect evidence", 
            "Motor cortex expression", "Strong speech/language connection", 
            "Glutamatergic specification", "Epigenome remodeler", "Transcriptional repressor", 
            "POU domain factor", "CHD8 haploinsufficiency", "FOXP2-regulated, frontal expression", 
            "Pioneer transcription factor", "Ventral eye specification", "Early B-cell factor", 
            "Monoaminergic specification", "Neurogenesis antagonist", "LIM-homeodomain", "Poorly characterized")
)

# Guardar archivo con metadatos adicionales
write.csv(mirror_genes, "../../results/integration/mirror_neuron_genes_curated_evidence_based.csv", row.names=FALSE)

# Generar resumen de evidencia
cat("Resumen de niveles de evidencia para genes de neuronas espejo:\n")
cat("High evidence: ", sum(mirror_genes$evidence_level == "High"), " genes\n")
cat("Medium evidence: ", sum(mirror_genes$evidence_level == "Medium"), " genes\n") 
cat("Low evidence: ", sum(mirror_genes$evidence_level == "Low"), " genes\n")

# Mostrar genes con evidencia alta y media
cat("\nGenes con evidencia ALTA:\n")
print(mirror_genes[mirror_genes$evidence_level == "High", c("gene_symbol", "functional_category", "primary_evidence")])

cat("\nGenes con evidencia MEDIA:\n")
print(mirror_genes[mirror_genes$evidence_level == "Medium", c("gene_symbol", "functional_category", "primary_evidence")])
