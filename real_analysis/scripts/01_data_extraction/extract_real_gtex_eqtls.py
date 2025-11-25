#!/usr/bin/env python3
"""
EXTRACT REAL GTEx v8 eQTLs FOR TARGET GENES
============================================

Purpose: Extract significant eQTL associations from GTEx v8 for genes
         associated with mirror neuron function and autism risk.

Author: Ana Sol Murzi (with guidance from Claudio)
Date: November 2025

Input:
    - GTEx v8 .signif_variant_gene_pairs.txt.gz files
    - Brain tissues: Frontal Cortex BA9, Anterior Cingulate BA24, Cortex, Cerebellar Hemisphere

Output:
    - CSV file with real eQTL associations
    - Quality control report
    - Variant list for LD calculation

Data Source: GTEx Portal (https://gtexportal.org/)
             GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2)

CRITICAL: This script uses ONLY real data from GTEx. No simulations.
"""

import pandas as pd
import numpy as np
import gzip
from pathlib import Path
import logging
from collections import defaultdict
import sys

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('../../results/eqtl_extracted/extraction.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class RealGTExExtractor:
    """
    Extract real eQTL data from GTEx v8 for target genes.
    NO SIMULATIONS. ONLY REAL DATA.
    """

    def __init__(self, base_dir="../../"):
        self.base_dir = Path(base_dir).resolve()
        self.data_dir = self.base_dir / "data" / "gtex_eqtl" / "GTEx_v8"
        self.results_dir = self.base_dir / "results" / "eqtl_extracted"

        # Create results directory
        self.results_dir.mkdir(parents=True, exist_ok=True)

        # Target genes with ENSEMBL IDs from GTEx egenes files
        self.target_genes = {
            'CACNA1C': {
                'ensembl': 'ENSG00000151067',
                'chr': 'chr12',
                'description': 'Calcium voltage-gated channel subunit alpha1 C'
            },
            'CHD8': {
                'ensembl': 'ENSG00000100888',
                'chr': 'chr14',
                'description': 'Chromodomain helicase DNA binding protein 8'
            },
            'CNTNAP2': {
                'ensembl': 'ENSG00000174469',
                'chr': 'chr7',
                'description': 'Contactin associated protein 2'
            },
            'FOXP2': {
                'ensembl': 'ENSG00000128573',
                'chr': 'chr7',
                'description': 'Forkhead box P2'
            },
            'THEMIS': {
                'ensembl': 'ENSG00000172673',
                'chr': 'chr6',
                'description': 'Thymocyte selection associated'
            },
            'SHANK3': {
                'ensembl': 'ENSG00000251322',
                'chr': 'chr22',
                'description': 'SH3 and multiple ankyrin repeat domains 3'
            },
            'OXTR': {
                'ensembl': 'ENSG00000180914',
                'chr': 'chr3',
                'description': 'Oxytocin receptor'
            },
            'MECP2': {
                'ensembl': 'ENSG00000169057',
                'chr': 'chrX',
                'description': 'Methyl-CpG binding protein 2'
            },
            # --- 16p11.2 CNV REGION DRIVERS ---
            'MAPK3': {
                'ensembl': 'ENSG00000102882',
                'chr': 'chr16',
                'description': 'Mitogen-activated protein kinase 3'
            },
            'KCTD13': {
                'ensembl': 'ENSG00000174943',
                'chr': 'chr16',
                'description': 'Potassium channel tetramerization domain containing 13'
            },
            'TAOK2': {
                'ensembl': 'ENSG00000149930',
                'chr': 'chr16',
                'description': 'TAO kinase 2'
            },
            'ALDOA': {
                'ensembl': 'ENSG00000149925',
                'chr': 'chr16',
                'description': 'Aldolase, fructose-bisphosphate A'
            },
            'PPP4C': {
                'ensembl': 'ENSG00000149923',
                'chr': 'chr16',
                'description': 'Protein phosphatase 4 catalytic subunit'
            },
            'TBX6': {
                'ensembl': 'ENSG00000149922',
                'chr': 'chr16',
                'description': 'T-box transcription factor 6'
            },
            'YPEL3': {
                'ensembl': 'ENSG00000090238',
                'chr': 'chr16',
                'description': 'Yippee like 3'
            }
        }

        # Brain tissues of interest for mirror neuron analysis
        self.brain_tissues = {
            'Brain_Amygdala': 'Amygdala',
            'Brain_Anterior_cingulate_cortex_BA24': 'Anterior cingulate cortex (BA24)',
            'Brain_Caudate_basal_ganglia': 'Caudate (basal ganglia)',
            'Brain_Cerebellar_Hemisphere': 'Cerebellar Hemisphere',
            'Brain_Cerebellum': 'Cerebellum',
            'Brain_Cortex': 'Cortex',
            'Brain_Frontal_Cortex_BA9': 'Frontal Cortex (BA9)',
            'Brain_Hippocampus': 'Hippocampus',
            'Brain_Hypothalamus': 'Hypothalamus',
            'Brain_Nucleus_accumbens_basal_ganglia': 'Nucleus accumbens (basal ganglia)',
            'Brain_Putamen_basal_ganglia': 'Putamen (basal ganglia)',
            'Brain_Spinal_cord_cervical_c-1': 'Spinal cord (cervical c-1)',
            'Brain_Substantia_nigra': 'Substantia nigra'
        }

        logger.info(f"Initialized RealGTExExtractor")
        logger.info(f"Data directory: {self.data_dir}")
        logger.info(f"Target genes: {len(self.target_genes)}")
        logger.info(f"Target tissues: {len(self.brain_tissues)}")

    def verify_data_files(self):
        """Verify that GTEx data files exist"""
        logger.info("\n=== VERIFYING DATA FILES ===")

        missing_files = []
        found_files = []

        for tissue_code in self.brain_tissues.keys():
            eqtl_file = self.data_dir / f"{tissue_code}.v8.signif_variant_gene_pairs.txt.gz"
            egenes_file = self.data_dir / f"{tissue_code}.v8.egenes.txt.gz"

            if eqtl_file.exists():
                found_files.append(str(eqtl_file))
                logger.info(f"✓ Found: {eqtl_file.name}")
            else:
                missing_files.append(str(eqtl_file))
                logger.warning(f"✗ Missing: {eqtl_file.name}")

        if missing_files:
            logger.error(f"Missing {len(missing_files)} data files!")
            logger.error("Cannot proceed without GTEx data files.")
            return False

        logger.info(f"\n✓ All {len(found_files)} required files found")
        return True

    def extract_eqtls_for_gene(self, tissue_code, gene_symbol, gene_info):
        """
        Extract all significant eQTLs for a specific gene from GTEx data.

        Args:
            tissue_code: GTEx tissue code (e.g., 'Brain_Frontal_Cortex_BA9')
            gene_symbol: Gene symbol (e.g., 'CACNA1C')
            gene_info: Dictionary with gene metadata

        Returns:
            DataFrame with real eQTL associations
        """
        eqtl_file = self.data_dir / f"{tissue_code}.v8.signif_variant_gene_pairs.txt.gz"
        ensembl_id = gene_info['ensembl']

        logger.info(f"  Searching for {gene_symbol} ({ensembl_id}) in {tissue_code}...")

        # Read GTEx eQTL file and filter for target gene
        eqtls = []

        try:
            with gzip.open(eqtl_file, 'rt') as f:
                header = f.readline().strip().split('\t')

                # Find column indices
                variant_idx = header.index('variant_id')
                gene_idx = header.index('gene_id')
                pval_idx = header.index('pval_nominal')
                slope_idx = header.index('slope')
                slope_se_idx = header.index('slope_se')
                maf_idx = header.index('maf')

                # Read line by line to avoid loading entire file
                for line in f:
                    fields = line.strip().split('\t')
                    gene_id = fields[gene_idx]

                    # Check if this is our target gene (match ENSEMBL ID without version)
                    if gene_id.split('.')[0] == ensembl_id.split('.')[0]:
                        eqtls.append({
                            'gene_symbol': gene_symbol,
                            'gene_id': gene_id,
                            'variant_id': fields[variant_idx],
                            'pval_nominal': float(fields[pval_idx]),
                            'slope': float(fields[slope_idx]),
                            'slope_se': float(fields[slope_se_idx]),
                            'maf': float(fields[maf_idx]),
                            'tissue': self.brain_tissues[tissue_code],
                            'tissue_code': tissue_code
                        })

        except FileNotFoundError:
            logger.warning(f"  File not found: {eqtl_file}")
            return pd.DataFrame()
        except Exception as e:
            logger.error(f"  Error reading {eqtl_file}: {e}")
            return pd.DataFrame()

        if eqtls:
            df = pd.DataFrame(eqtls)
            logger.info(f"  Found {len(df)} significant eQTLs for {gene_symbol}")
            return df
        else:
            logger.info(f"  No eQTLs found for {gene_symbol}")
            return pd.DataFrame()

    def extract_all_eqtls(self):
        """
        Extract eQTLs for all target genes across all brain tissues.
        This is the main extraction function.
        """
        logger.info("\n=== EXTRACTING REAL eQTLs FROM GTEx v8 ===")
        logger.info("NOTE: Using ONLY real data from GTEx. NO SIMULATIONS.\n")

        all_eqtls = []
        gene_tissue_counts = defaultdict(lambda: defaultdict(int))

        for gene_symbol, gene_info in self.target_genes.items():
            logger.info(f"\nProcessing {gene_symbol}...")

            for tissue_code, tissue_name in self.brain_tissues.items():
                eqtls_df = self.extract_eqtls_for_gene(tissue_code, gene_symbol, gene_info)

                if not eqtls_df.empty:
                    all_eqtls.append(eqtls_df)
                    gene_tissue_counts[gene_symbol][tissue_code] = len(eqtls_df)

        if not all_eqtls:
            logger.error("No eQTLs extracted! Check data files and gene IDs.")
            return None

        # Combine all eQTLs
        combined_eqtls = pd.concat(all_eqtls, ignore_index=True)

        logger.info(f"\n=== EXTRACTION COMPLETE ===")
        logger.info(f"Total eQTL associations extracted: {len(combined_eqtls)}")
        logger.info(f"Unique variants: {combined_eqtls['variant_id'].nunique()}")
        logger.info(f"Genes with eQTLs: {combined_eqtls['gene_symbol'].nunique()}")

        return combined_eqtls, gene_tissue_counts

    def quality_control(self, eqtls_df):
        """
        Perform quality control checks on extracted eQTLs.
        """
        logger.info("\n=== QUALITY CONTROL ===")

        # Check 1: P-value distribution
        logger.info(f"\nP-value statistics:")
        logger.info(f"  Min p-value: {eqtls_df['pval_nominal'].min():.2e}")
        logger.info(f"  Median p-value: {eqtls_df['pval_nominal'].median():.2e}")
        logger.info(f"  Max p-value: {eqtls_df['pval_nominal'].max():.2e}")

        # Check 2: Variant ID format (should be chr_pos_ref_alt_b38)
        sample_variants = eqtls_df['variant_id'].head(5).tolist()
        logger.info(f"\nSample variant IDs (should be chr_pos_ref_alt_b38):")
        for var in sample_variants:
            logger.info(f"  {var}")

        # Check 3: Effect sizes
        logger.info(f"\nEffect size (slope) statistics:")
        logger.info(f"  Mean |slope|: {eqtls_df['slope'].abs().mean():.3f}")
        logger.info(f"  Range: [{eqtls_df['slope'].min():.3f}, {eqtls_df['slope'].max():.3f}]")

        # Check 4: MAF distribution
        logger.info(f"\nMinor Allele Frequency statistics:")
        logger.info(f"  Mean MAF: {eqtls_df['maf'].mean():.3f}")
        logger.info(f"  Rare variants (MAF < 0.05): {(eqtls_df['maf'] < 0.05).sum()}")
        logger.info(f"  Common variants (MAF ≥ 0.05): {(eqtls_df['maf'] >= 0.05).sum()}")

        # Check 5: Per-gene summary
        logger.info(f"\neQTLs per gene:")
        for gene in eqtls_df['gene_symbol'].unique():
            n_eqtls = (eqtls_df['gene_symbol'] == gene).sum()
            n_tissues = eqtls_df[eqtls_df['gene_symbol'] == gene]['tissue'].nunique()
            logger.info(f"  {gene}: {n_eqtls} eQTLs across {n_tissues} tissues")

        return True

    def save_results(self, eqtls_df, gene_tissue_counts):
        """Save extracted eQTLs and summary statistics"""
        logger.info("\n=== SAVING RESULTS ===")

        # Save main eQTL file
        output_file = self.results_dir / "real_gtex_eqtls.csv"
        eqtls_df.to_csv(output_file, index=False)
        logger.info(f"✓ Saved eQTLs to: {output_file}")

        # Save variant list (for LD calculation)
        variant_list_file = self.results_dir / "variant_list_for_ld.txt"
        unique_variants = eqtls_df['variant_id'].unique()
        with open(variant_list_file, 'w') as f:
            for variant in sorted(unique_variants):
                f.write(f"{variant}\n")
        logger.info(f"✓ Saved variant list ({len(unique_variants)} variants): {variant_list_file}")

        # Save summary report
        summary_file = self.results_dir / "extraction_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("REAL GTEx v8 eQTL EXTRACTION SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Date: {pd.Timestamp.now()}\n")
            f.write(f"Source: GTEx Analysis V8\n")
            f.write(f"Genome Build: GRCh38 (hg38/b38)\n\n")

            f.write(f"TOTAL STATISTICS:\n")
            f.write(f"  Total eQTL associations: {len(eqtls_df)}\n")
            f.write(f"  Unique variants: {eqtls_df['variant_id'].nunique()}\n")
            f.write(f"  Genes analyzed: {eqtls_df['gene_symbol'].nunique()}\n")
            f.write(f"  Tissues: {eqtls_df['tissue'].nunique()}\n\n")

            f.write(f"PER-GENE SUMMARY:\n")
            for gene in sorted(self.target_genes.keys()):
                if gene in gene_tissue_counts:
                    total = sum(gene_tissue_counts[gene].values())
                    f.write(f"\n  {gene}: {total} total eQTLs\n")
                    for tissue, count in gene_tissue_counts[gene].items():
                        f.write(f"    - {tissue}: {count}\n")
                else:
                    f.write(f"\n  {gene}: No eQTLs found\n")

            f.write(f"\nP-VALUE STATISTICS:\n")
            f.write(f"  Min: {eqtls_df['pval_nominal'].min():.2e}\n")
            f.write(f"  Median: {eqtls_df['pval_nominal'].median():.2e}\n")
            f.write(f"  Max: {eqtls_df['pval_nominal'].max():.2e}\n")

        logger.info(f"✓ Saved summary report: {summary_file}")

        return True

    def run(self):
        """Execute the complete extraction pipeline"""
        logger.info("\n" + "=" * 70)
        logger.info("REAL GTEx v8 eQTL EXTRACTION PIPELINE")
        logger.info("NO SIMULATIONS - ONLY REAL DATA")
        logger.info("=" * 70 + "\n")

        # Step 1: Verify data files exist
        if not self.verify_data_files():
            logger.error("Data verification failed. Cannot proceed.")
            return False

        # Step 2: Extract eQTLs
        result = self.extract_all_eqtls()
        if result is None:
            logger.error("Extraction failed. No eQTLs found.")
            return False

        eqtls_df, gene_tissue_counts = result

        # Step 3: Quality control
        self.quality_control(eqtls_df)

        # Step 4: Save results
        self.save_results(eqtls_df, gene_tissue_counts)

        logger.info("\n" + "=" * 70)
        logger.info("EXTRACTION PIPELINE COMPLETED SUCCESSFULLY")
        logger.info("=" * 70 + "\n")

        return True

def main():
    """Main execution function"""
    extractor = RealGTExExtractor()
    success = extractor.run()

    if success:
        print("\n✓ Real eQTL extraction completed successfully!")
        print(f"✓ Results saved to: {extractor.results_dir}")
        print("\nNext steps:")
        print("  1. Review extraction_summary.txt")
        print("  2. Proceed to LD calculation (script 02)")
        return 0
    else:
        print("\n✗ Extraction failed. Check log for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
