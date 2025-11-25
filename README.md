# Genetic Analysis: Autism and Mirror Neuron System

Integrated analysis of the genetic architecture of Autism Spectrum Disorder (ASD) focused on genes associated with the mirror neuron system, using GWAS and gene expression data.

## Overview

This project investigates the relationship between genetic variants associated with autism and the expression of genes implicated in the mirror neuron system through colocalization analysis, causal inference, and cell-type-specific analyses.

### Data Sources

- **GWAS**: iPSYCH-PGC ASD (Grove et al. 2019) - 18,381 cases, 27,969 controls
- **eQTLs**: GTEx Analysis V8 (dbGaP phs000424.v8.p2) - 13 brain tissues
- **Single-cell**: Allen Brain Atlas - Human M1 cortex (neuronal and non-neuronal)
- **LD**: 1000 Genomes Phase 3 (CEU population, GRCh38)

### Genes Analyzed

13 candidate genes related to mirror neurons and autism risk:

**Core mirror neuron system genes:**
- FOXP2 (language and communication)
- CNTNAP2 (neuronal and synaptic development)

**Autism risk genes with neuronal function:**
- CHD8, SHANK3, MECP2 (chromatin regulation/synapses)
- CACNA1C, OXTR (calcium/oxytocin signaling)
- MAPK3, TAOK2, PPP4C (signaling cascades)
- ALDOA, TBX6, YPEL3, THEMIS, KCTD13 (metabolism/development)

## Project Structure

```
real_analysis/
├── data/
│   ├── gwas/                    # iPSYCH-PGC ASD GWAS summary statistics
│   ├── gtex_eqtl/              # GTEx v8 eQTL data (brain tissues)
│   ├── gtex_bulk/              # GTEx bulk RNA-seq expression data
│   └── allen_brain_atlas/      # Single-cell RNA-seq (M1 cortex)
│
├── scripts/
│   ├── 01_data_extraction/     # GTEx eQTL extraction
│   ├── 02_ld_calculation/      # LD calculation using LDlink API
│   ├── 03_colocalization/      # Bayesian colocalization (bulk tissue)
│   ├── 04_single_cell/         # Cell-type specific analysis
│   ├── 05_smr_heidi/          # Summary-based Mendelian Randomization
│   ├── 06_bayesian_network/    # Causal network learning
│   └── final_figure_generation.R
│
└── results/
    ├── eqtl_extracted/         # Extracted eQTLs by gene/tissue
    ├── ld_analysis/            # LD matrices by gene
    ├── colocalization_standard/ # Bulk colocalization results
    ├── celltype_deconvolution/ # Cell type proportions (MuSiC)
    ├── celltype_eqtls/         # Cell-type aware eQTLs
    ├── celltype_colocalization_real/ # Cell-type specific colocalization
    ├── smr_heidi/              # Causal inference results
    ├── causal_networks/        # Learned Bayesian networks
    ├── validation/             # Validation analyses
    └── final_figures/          # Final manuscript figures
```

## Analysis Pipeline

### 01. eQTL Data Extraction

**Script:** `scripts/01_data_extraction/extract_real_gtex_eqtls.py`

Extracts significant eQTL associations from GTEx v8 for the 13 target genes in brain tissues.

**Output:**
- 2,686 total eQTL associations
- 697 unique variants
- 13 brain tissues
- P-value range: 1.52e-30 to 1.80e-04

### 02. Linkage Disequilibrium (LD) Calculation

**Script:** `scripts/02_ld_calculation/calculate_ld_ldlink.py`

Calculates LD matrices using LDlink API (NIH/NCI) for eQTL variants, utilizing 1000 Genomes Phase 3 data (CEU population).

**Metrics:**
- Variants analyzed: 702 (filtered by availability in 1000G)
- Reference population: CEU (European)
- r² range: 0.0 - 1.0
- High LD pairs (r² > 0.8) identified per gene

### 03. Bayesian Colocalization (Bulk Tissue)

**Scripts:**
- `scripts/03_colocalization/01_prepare_gwas_data.py` - GWAS data preparation
- `scripts/03_colocalization/02_bayesian_colocalization.R` - Colocalization analysis

Applies Bayesian colocalization method (Giambartolomei et al. 2014) to evaluate whether GWAS and eQTLs share causal signals.

**Hypotheses evaluated:**
- H0: No association
- H1: GWAS only
- H2: eQTL only
- H3: Both signals, different causal variants
- H4: Both signals, shared causal variant (colocalization)

**Criterion:** PP.H4 > 0.75 for positive colocalization

### 04. Cell-Type Specific Analysis

**Scripts:**
- `01a_inspect_allen_brain_data.py` - Allen Brain Atlas data inspection
- `01b_process_allen_brain_singlecell.py` - Single-cell data processing
- `02a_split_gtex.py` - GTEx data splitting for processing
- `02b_generate_signature_matrix_lowmem.py` - Signature matrix generation
- `02c_music_deconvolution.R` - Cell proportion deconvolution (MuSiC)
- `03_celltype_eqtl_interaction.R` - eQTL x cell type interaction analysis
- `04_celltype_aware_colocalization.R` - Cell-type adjusted colocalization

**Pipeline:**
1. Deconvolution of bulk tissue into cell type proportions
2. Generation of gene signature matrix by cell type
3. eQTL analysis with cell_type x genotype interaction
4. Colocalization adjusted for cellular composition

**Cell types analyzed:**
- Excitatory neurons
- Inhibitory neurons
- Astrocytes
- Oligodendrocytes
- Microglia
- Endothelial cells

### 05. Summary-based Mendelian Randomization (SMR-HEIDI)

**Script:** `scripts/05_smr_heidi/01_SMR_HEIDI_analysis.R`

Applies SMR-HEIDI (Zhu et al. 2016) for causal inference between gene expression and autism phenotype, distinguishing pleiotropy from causality.

**Tests:**
- **SMR test**: Evaluates causal effect of gene expression on ASD
- **HEIDI test**: Distinguishes causality from pleiotropy/LD
- Criterion: P_SMR < 0.05 and P_HEIDI > 0.01

### 06. Bayesian Network Learning

**Scripts:**
- `05a_empirical_gene_selection.py` - Empirical gene selection
- `05b_bayesian_network_learning.R` - Network structure learning

Learns causal networks between genes using Bayesian structure learning algorithms (bnlearn) integrating expression and genotype data.

**Methods:**
- Hill-climbing with BIC score
- Tabu search
- Bootstrap for arc stability

### 07. Final Figure Generation

**Script:** `scripts/final_figure_generation.R`

Generates integrated visualizations of all analyses for publication.

## Main Results

### eQTL Statistics by Gene

| Gene      | eQTLs | Variants | Tissues | Min P-value |
|----------|-------|-----------|---------|-------------|
| YPEL3    | 754   | 119       | 6       | 1.52e-30    |
| TBX6     | 639   | 103       | 12      | 4.08e-15    |
| MAPK3    | 593   | 143       | 6       | 1.07e-19    |
| OXTR     | 304   | 70        | 13      | 2.51e-19    |
| CACNA1C  | 271   | 146       | 4       | 2.98e-17    |
| CNTNAP2  | 43    | 43        | 3       | 5.16e-06    |
| FOXP2    | 37    | 36        | 2       | 1.70e-05    |
| CHD8     | 20    | 19        | 1       | 2.68e-09    |
| ALDOA    | 11    | 10        | 2       | 2.09e-05    |
| TAOK2    | 8     | 7         | 1       | -           |
| SHANK3   | 2     | 2         | 1       | 1.10e-05    |
| THEMIS   | 2     | 2         | 1       | 5.10e-06    |
| PPP4C    | 2     | 2         | 1       | -           |

**Note:** KCTD13 and MECP2 showed no significant eQTLs in GTEx v8 brain tissues.

### Linkage Disequilibrium

Genes with high average LD (mean r² > 0.6):
- FOXP2: r² = 0.921
- TAOK2: r² = 0.877
- YPEL3: r² = 0.644
- CHD8: r² = 0.644

Genes with low LD (complex haplotype structure):
- CNTNAP2: r² = 0.175
- OXTR: r² = 0.182

## System Requirements

### Software

**Python 3.8+:**
- pandas >= 1.3.0
- numpy >= 1.21.0
- scipy >= 1.7.0
- scanpy >= 1.8.0 (single-cell analysis)
- anndata >= 0.8.0
- requests (LDlink API)

**R 4.1+:**
- tidyverse
- data.table
- coloc
- MuSiC (deconvolution)
- bnlearn (Bayesian networks)
- ggplot2
- patchwork

### External Data Required

1. **GTEx v8 eQTL data**:
   - Brain_*.signif_variant_gene_pairs.txt.gz

2. **GTEx v8 bulk expression**:
   - Gene_TPM matrices per tissue

3. **Allen Brain Atlas single-cell**:
   - human_M1_10x_cells.h5ad

4. **iPSYCH-PGC GWAS**:
   - iPSYCH-PGC_ASD_Nov2017

## Running the Pipeline

### Recommended Sequence

```bash
# 1. eQTL extraction
cd scripts/01_data_extraction
python extract_real_gtex_eqtls.py

# 2. LD calculation
cd ../02_ld_calculation
python calculate_ld_ldlink.py

# 3. Bulk colocalization
cd ../03_colocalization
python 01_prepare_gwas_data.py
Rscript 02_bayesian_colocalization.R

# 4. Single-cell analysis
cd ../04_single_cell
python 01a_inspect_allen_brain_data.py
python 01b_process_allen_brain_singlecell.py
python 02a_split_gtex.py
python 02b_generate_signature_matrix_lowmem.py
Rscript 02c_music_deconvolution.R
Rscript 03_celltype_eqtl_interaction.R
Rscript 04_celltype_aware_colocalization.R

# 5. SMR-HEIDI
cd ../05_smr_heidi
Rscript 01_SMR_HEIDI_analysis.R

# 6. Bayesian networks
cd ../06_bayesian_network
python 05a_empirical_gene_selection.py
Rscript 05b_bayesian_network_learning.R

# 7. Final figures
cd ../
Rscript final_figure_generation.R
```


## Methodological Considerations

### Genomic Coordinates

- **GTEx v8**: hg38
- **GWAS iPSYCH-PGC**: hg19
- **Conversion**: Performed via liftOver when necessary
- All coordinates reported in final results are in **GRCh38**

### Quality Control

**eQTLs:**
- Only GTEx-significant associations (FDR < 0.05)
- Duplicate filtering (same SNP-gene-tissue)
- Exclusion of complex multiallelic variants

**GWAS:**
- Standard QC applied by Grove et al. 2019
- INFO > 0.9, MAF > 0.01
- No GWAS significance threshold applied to maximize power in colocalization

**LD:**
- Only variants available in 1000G Phase 3
- CEU population (n=99) as European reference
- Variants with call rate < 0.95 excluded

### Limitations

1. **Sample size**: Analyses are limited by GTEx brain tissue sample sizes (n=100-200 per tissue)
2. **Bulk tissue**: Cellular heterogeneity may mask cell-type specific effects (mitigated with cell-type aware analysis)
3. **Population**: Analysis limited to European ancestry (both GWAS and LD reference)
4. **Causality**: Colocalization and SMR analyses suggest but do not prove causality

## Key References

### Data

- **GTEx Consortium** (2020). The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science, 369(6509), 1318-1330.
- **Grove et al.** (2019). Identification of common genetic risk variants for autism spectrum disorder. Nature Genetics, 51(3), 431-444.

### Methods

- **Giambartolomei et al.** (2014). Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genetics, 10(5), e1004383.
- **Zhu et al.** (2016). Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nature Genetics, 48(5), 481-487.
- **Wang et al.** (2019). Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nature Communications, 10(1), 380.

## Author

Ana Sol Murzi

## Contact

For questions about the analysis or data access, contact me at anasolm26@gmail.com.

## License

This project uses controlled-access data. Scripts are open source but data requires access requests to:
- GTEx Portal (gtexportal.org)
- iPSYCH (ipsych.dk)
- Allen Brain Atlas (brain-map.org)
