#!/usr/bin/env python3
"""
CALCULATE REAL LD MATRICES USING LDlink API
===========================================

Purpose: Calculate pairwise linkage disequilibrium (r²) between eQTL variants
         using real genotype data from 1000 Genomes via LDlink API.

Author: Ana Sol Murzi (with guidance from Claudio)
Date: November 2025

Input:
    - List of variants from GTEx eQTL extraction
    - Format: chr_pos_ref_alt_b38 (GRCh38/hg38)

Output:
    - LD matrices (r²) for each gene
    - Summary statistics
    - QC visualizations

Data Source: LDlink (NIH/NCI)
             https://ldlink.nci.nih.gov/
             1000 Genomes Phase 3, CEU population (European ancestry)

CRITICAL: This uses REAL genotype data, not proximity-based estimates.
"""

import pandas as pd
import numpy as np
import requests
import time
import json
from pathlib import Path
import logging
import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('../../results/ld_analysis/ld_calculation.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class LDlinkCalculator:
    """
    Calculate LD matrices using LDlink API with real 1000 Genomes data.
    NO PROXIMITY ESTIMATES. ONLY REAL LD FROM GENOTYPES.
    """

    def __init__(self, base_dir="../../"):
        self.base_dir = Path(base_dir).resolve()
        self.eqtl_dir = self.base_dir / "results" / "eqtl_extracted"
        self.ld_dir = self.base_dir / "results" / "ld_analysis"

        # Create LD results directory
        self.ld_dir.mkdir(parents=True, exist_ok=True)

        # LDlink API settings
        self.ldlink_base_url = "https://ldlink.nih.gov/LDlinkRest"
        self.api_token = "7afc6202f9a0"  # LDlink API token
        self.population = "CEU"  # European ancestry (matches Grove et al. GWAS)
        self.genome_build = "grch38_high_coverage"  # Our GTEx variants are in GRCh38

        # Rate limiting
        self.request_delay = 1.0  # seconds between requests
        self.max_retries = 3

        logger.info("Initialized LDlinkCalculator")
        logger.info(f"Population: {self.population} (European)")
        logger.info(f"Genome build: {self.genome_build}")
        logger.info(f"Results directory: {self.ld_dir}")

    def load_eqtl_variants(self):
        """Load eQTL variants from extraction results"""
        logger.info("\n=== LOADING eQTL VARIANTS ===")

        eqtl_file = self.eqtl_dir / "real_gtex_eqtls.csv"

        if not eqtl_file.exists():
            logger.error(f"eQTL file not found: {eqtl_file}")
            return None

        eqtls_df = pd.read_csv(eqtl_file)
        logger.info(f"Loaded {len(eqtls_df)} eQTL associations")
        logger.info(f"Unique variants: {eqtls_df['variant_id'].nunique()}")
        logger.info(f"Genes: {eqtls_df['gene_symbol'].nunique()}")

        return eqtls_df

    def convert_to_rsid(self, variant_id):
        """
        Convert GTEx variant ID (chr_pos_ref_alt_b38) to rsID using LDlink SNPchip.

        Args:
            variant_id: GTEx format (e.g., 'chr12_2188187_C_T_b38')

        Returns:
            rsID or original variant_id if conversion fails
        """
        try:
            # Parse variant ID
            parts = variant_id.replace('_b38', '').split('_')
            if len(parts) < 4:
                return variant_id

            chrom = parts[0].replace('chr', '')
            pos = parts[1]

            # Query LDlink SNPchip to get rsID
            url = f"{self.ldlink_base_url}/snpchip"
            params = {
                'chr': chrom,
                'pos': pos,
                'genome_build': self.genome_build
            }

            time.sleep(self.request_delay)
            response = requests.get(url, params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()
                if 'rsid' in data and data['rsid']:
                    logger.debug(f"Converted {variant_id} → {data['rsid']}")
                    return data['rsid']

            logger.debug(f"Could not convert {variant_id}, using original ID")
            return variant_id

        except Exception as e:
            logger.warning(f"Error converting {variant_id}: {e}")
            return variant_id

    def query_ldmatrix(self, variants, gene_name):
        """
        Query LDlink LDmatrix API for pairwise LD between variants.

        Args:
            variants: List of variant IDs or rsIDs
            gene_name: Gene name for logging

        Returns:
            DataFrame with LD matrix
        """
        logger.info(f"\n  Querying LDmatrix for {len(variants)} variants...")

        # LDmatrix API endpoint with token
        url = f"{self.ldlink_base_url}/ldmatrix?token={self.api_token}"

        # Prepare variant list (newline separated)
        snps_string = "\n".join(variants)

        # API parameters (JSON format for POST)
        payload = {
            'snps': snps_string,
            'pop': self.population,
            'r2_d': 'r2',  # Request r² values
            'genome_build': self.genome_build
        }

        headers = {
            'Content-Type': 'application/json'
        }

        try:
            # Make request
            time.sleep(self.request_delay)
            logger.info(f"  Sending request to LDlink API...")

            response = requests.post(url, json=payload, headers=headers, timeout=120)

            if response.status_code != 200:
                logger.error(f"  API error: Status code {response.status_code}")
                logger.error(f"  Response: {response.text}")
                return None

            # Parse response
            logger.info(f"  Parsing LD matrix...")
            lines = response.text.strip().split('\n')

            # Check for error messages
            if any('error' in line.lower() for line in lines[:5]):
                logger.error(f"  API returned error: {lines[0]}")
                return None

            # Parse tab-delimited matrix
            # First line is header with variant IDs
            header = lines[0].split('\t')
            variant_ids = header[1:]  # Skip first column (row names)

            # Parse data rows
            ld_data = []
            row_names = []

            for line in lines[1:]:
                if not line.strip():
                    continue

                fields = line.split('\t')
                row_names.append(fields[0])
                # Handle 'NA', '-', and empty values
                ld_values = []
                for x in fields[1:]:
                    if x in ['-', 'NA', '']:
                        ld_values.append(np.nan)
                    else:
                        try:
                            ld_values.append(float(x))
                        except ValueError:
                            ld_values.append(np.nan)
                ld_data.append(ld_values)

            # Create DataFrame
            ld_matrix = pd.DataFrame(
                ld_data,
                index=row_names,
                columns=variant_ids
            )

            logger.info(f"  ✓ Received LD matrix: {ld_matrix.shape}")
            logger.info(f"  ✓ Mean r²: {ld_matrix.mean().mean():.3f}")
            logger.info(f"  ✓ High LD pairs (r² > 0.8): {(ld_matrix > 0.8).sum().sum() // 2}")

            return ld_matrix

        except requests.exceptions.Timeout:
            logger.error(f"  Request timeout for {gene_name}")
            return None
        except Exception as e:
            logger.error(f"  Error querying LDmatrix: {e}")
            return None

    def calculate_ld_by_gene(self, eqtls_df):
        """
        Calculate LD matrices for each gene separately.

        Args:
            eqtls_df: DataFrame with eQTL data

        Returns:
            Dictionary of LD matrices by gene
        """
        logger.info("\n=== CALCULATING LD MATRICES BY GENE ===")
        logger.info("Using LDlink API with 1000 Genomes CEU population")
        logger.info("NOTE: This uses REAL genotypes, not proximity estimates\n")

        ld_matrices = {}

        for gene in eqtls_df['gene_symbol'].unique():
            logger.info(f"\nProcessing {gene}...")

            # Get variants for this gene
            gene_variants = eqtls_df[eqtls_df['gene_symbol'] == gene]['variant_id'].unique()
            logger.info(f"  Variants: {len(gene_variants)}")

            if len(gene_variants) < 2:
                logger.warning(f"  Skipping {gene}: need at least 2 variants for LD")
                continue

            if len(gene_variants) > 200:
                logger.warning(f"  {gene} has {len(gene_variants)} variants (max 200)")
                logger.info(f"  Using top 200 by significance...")

                # Get top 200 by p-value
                gene_eqtls = eqtls_df[eqtls_df['gene_symbol'] == gene].copy()
                gene_eqtls = gene_eqtls.sort_values('pval_nominal')
                gene_variants = gene_eqtls.head(200)['variant_id'].unique()

            # Convert to rsIDs (LDlink prefers rsIDs but can handle chr:pos)
            logger.info(f"  Preparing variants...")
            variants_for_api = []

            for var in gene_variants:
                # LDlink can accept chr:pos format
                # Convert chr12_2188187_C_T_b38 → 12:2188187
                parts = var.replace('_b38', '').split('_')
                if len(parts) >= 2:
                    chrom = parts[0].replace('chr', '')
                    pos = parts[1]
                    api_format = f"{chrom}:{pos}"
                    variants_for_api.append(api_format)
                else:
                    variants_for_api.append(var)

            # Query API
            ld_matrix = self.query_ldmatrix(variants_for_api, gene)

            if ld_matrix is not None:
                # Map back to original variant IDs
                id_mapping = dict(zip(variants_for_api, gene_variants))

                ld_matrix.index = ld_matrix.index.map(lambda x: id_mapping.get(x, x))
                ld_matrix.columns = ld_matrix.columns.map(lambda x: id_mapping.get(x, x))

                ld_matrices[gene] = ld_matrix
                logger.info(f"  ✓ Saved LD matrix for {gene}")
            else:
                logger.error(f"  ✗ Failed to get LD matrix for {gene}")

            # Be respectful of API rate limits
            time.sleep(2)

        logger.info(f"\n=== LD CALCULATION COMPLETE ===")
        logger.info(f"Successfully calculated LD for {len(ld_matrices)} genes")

        return ld_matrices

    def quality_control(self, ld_matrices):
        """
        Perform quality control on LD matrices.
        """
        logger.info("\n=== QUALITY CONTROL ===")

        for gene, ld_matrix in ld_matrices.items():
            logger.info(f"\n{gene}:")
            logger.info(f"  Matrix size: {ld_matrix.shape[0]} × {ld_matrix.shape[1]}")
            logger.info(f"  Mean r²: {ld_matrix.mean().mean():.3f}")
            logger.info(f"  Max r²: {ld_matrix.max().max():.3f}")
            logger.info(f"  High LD pairs (r² > 0.8): {(ld_matrix > 0.8).sum().sum() // 2}")
            logger.info(f"  Low LD pairs (r² < 0.2): {(ld_matrix < 0.2).sum().sum() // 2}")

            # Check for missing values
            missing = ld_matrix.isna().sum().sum()
            if missing > 0:
                logger.warning(f"  Missing values: {missing}")

        return True

    def create_ld_heatmaps(self, ld_matrices):
        """Create heatmap visualizations of LD matrices"""
        logger.info("\n=== CREATING LD HEATMAPS ===")

        for gene, ld_matrix in ld_matrices.items():
            try:
                fig, ax = plt.subplots(figsize=(12, 10))

                # Create heatmap
                sns.heatmap(
                    ld_matrix,
                    cmap='RdYlBu_r',
                    vmin=0,
                    vmax=1,
                    square=True,
                    cbar_kws={'label': 'r²'},
                    ax=ax
                )

                ax.set_title(f'LD Matrix (r²) for {gene}', fontsize=14, fontweight='bold')
                ax.set_xlabel('Variants')
                ax.set_ylabel('Variants')

                # Save
                output_file = self.ld_dir / f"ld_heatmap_{gene}.png"
                plt.tight_layout()
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                plt.close()

                logger.info(f"  ✓ Saved heatmap: {output_file.name}")

            except Exception as e:
                logger.error(f"  Error creating heatmap for {gene}: {e}")

    def save_results(self, ld_matrices):
        """Save LD matrices and summary statistics"""
        logger.info("\n=== SAVING RESULTS ===")

        # Save each gene's LD matrix as CSV
        for gene, ld_matrix in ld_matrices.items():
            output_file = self.ld_dir / f"ld_matrix_{gene}.csv"
            ld_matrix.to_csv(output_file)
            logger.info(f"  ✓ Saved: {output_file.name}")

        # Save summary statistics
        summary_file = self.ld_dir / "ld_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("REAL LD CALCULATION SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Date: {pd.Timestamp.now()}\n")
            f.write(f"Method: LDlink API (NIH/NCI)\n")
            f.write(f"Reference: 1000 Genomes Phase 3\n")
            f.write(f"Population: {self.population} (European)\n")
            f.write(f"Genome Build: {self.genome_build}\n\n")

            f.write("PER-GENE SUMMARY:\n")
            for gene, ld_matrix in ld_matrices.items():
                f.write(f"\n{gene}:\n")
                f.write(f"  Variants: {ld_matrix.shape[0]}\n")
                f.write(f"  Mean r²: {ld_matrix.mean().mean():.3f}\n")
                f.write(f"  Max r²: {ld_matrix.max().max():.3f}\n")
                f.write(f"  High LD (r² > 0.8): {(ld_matrix > 0.8).sum().sum() // 2} pairs\n")

        logger.info(f"  ✓ Saved summary: {summary_file.name}")

        return True

    def run(self):
        """Execute the complete LD calculation pipeline"""
        logger.info("\n" + "=" * 70)
        logger.info("REAL LD CALCULATION USING LDlink API")
        logger.info("1000 Genomes Phase 3 - CEU Population")
        logger.info("NO PROXIMITY ESTIMATES - ONLY REAL GENOTYPES")
        logger.info("=" * 70 + "\n")

        # Step 1: Load variants
        eqtls_df = self.load_eqtl_variants()
        if eqtls_df is None:
            return False

        # Step 2: Calculate LD matrices
        ld_matrices = self.calculate_ld_by_gene(eqtls_df)

        if not ld_matrices:
            logger.error("No LD matrices calculated!")
            return False

        # Step 3: Quality control
        self.quality_control(ld_matrices)

        # Step 4: Create visualizations
        self.create_ld_heatmaps(ld_matrices)

        # Step 5: Save results
        self.save_results(ld_matrices)

        logger.info("\n" + "=" * 70)
        logger.info("LD CALCULATION COMPLETED SUCCESSFULLY")
        logger.info("=" * 70 + "\n")

        return True

def main():
    """Main execution function"""
    calculator = LDlinkCalculator()
    success = calculator.run()

    if success:
        print("\n✓ LD calculation completed successfully!")
        print(f"✓ Results saved to: {calculator.ld_dir}")
        print("\nNext steps:")
        print("  1. Review LD matrices and heatmaps")
        print("  2. Proceed to colocalization analysis (script 03)")
        return 0
    else:
        print("\n✗ LD calculation failed. Check log for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
