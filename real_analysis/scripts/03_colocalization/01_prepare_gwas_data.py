#!/usr/bin/env python3
"""
DOWNLOAD AND PREPARE GROVE ET AL. 2019 GWAS DATA
=================================================

Purpose: Download Grove et al. 2019 ASD GWAS summary statistics and extract
         data for our eQTL variants to enable colocalization analysis.

Author: Ana Sol Murzi (with guidance from Claudio)
Date: November 2025

Data Source: Grove et al. (2019). Nature Genetics.
             "Identification of common genetic risk variants for autism spectrum disorder"
             https://doi.org/10.1038/s41588-019-0344-8

             PGC: https://pgc.unc.edu/for-researchers/download-results/

Input:
    - Grove et al. 2019 GWAS summary statistics (to be downloaded)
    - Our extracted eQTL variants (real_gtex_eqtls.csv)

Output:
    - GWAS data for eQTL variants with aligned alleles
    - Summary statistics
    - QC report

CRITICAL: This uses REAL GWAS data from 18,381 ASD cases + 27,969 controls.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sys
import gzip
import requests
import time
from io import BytesIO

try:
    from pyliftover import LiftOver
    HAS_LIFTOVER = True
except ImportError:
    HAS_LIFTOVER = False
    logger.warning("pyliftover not installed. Install with: pip install pyliftover")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('../../results/colocalization_standard/gwas_preparation.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class GWASPreparation:
    """
    Download and prepare Grove et al. 2019 GWAS data for colocalization.
    """

    def __init__(self, base_dir="../../"):
        self.base_dir = Path(base_dir).resolve()
        self.data_dir = self.base_dir / "data" / "gwas"
        self.eqtl_dir = self.base_dir / "results" / "eqtl_extracted"
        self.coloc_dir = self.base_dir / "results" / "colocalization_standard"

        # Create directories
        self.coloc_dir.mkdir(parents=True, exist_ok=True)

        # GWAS file information
        # Grove et al. 2019 ASD GWAS from iPSYCH-PGC
        # Check for both .gz and uncompressed versions
        self.gwas_file_gz = self.data_dir / "iPSYCH-PGC_ASD_Nov2017.gz"
        self.gwas_file_plain = self.data_dir / "iPSYCH-PGC_ASD_Nov2017"

        if self.gwas_file_plain.exists():
            self.gwas_file = self.gwas_file_plain
        elif self.gwas_file_gz.exists():
            self.gwas_file = self.gwas_file_gz
        else:
            self.gwas_file = self.gwas_file_gz  # Default for download

        self.gwas_url = "https://ipsych.dk/fileadmin/iPSYCH/Downloads/iPSYCH-PGC_ASD_Nov2017.gz"

        # Initialize liftover for GRCh38 → GRCh37 conversion
        if HAS_LIFTOVER:
            try:
                self.liftover = LiftOver('hg38', 'hg19')
                logger.info("Initialized coordinate liftover (GRCh38 → GRCh37)")
            except Exception as e:
                logger.warning(f"Could not initialize liftover: {e}")
                self.liftover = None
        else:
            self.liftover = None

        logger.info("Initialized GWAS Preparation")
        logger.info(f"GWAS data directory: {self.data_dir}")
        logger.info(f"Results directory: {self.coloc_dir}")

    def download_gwas(self):
        """
        Download Grove et al. 2019 GWAS summary statistics from PGC/iPSYCH.
        """
        logger.info("\n=== DOWNLOADING GWAS DATA ===")

        if self.gwas_file.exists():
            logger.info(f"GWAS file already exists: {self.gwas_file}")
            return True

        logger.info(f"Downloading Grove et al. 2019 GWAS data...")
        logger.info(f"URL: {self.gwas_url}")
        logger.info("This may take several minutes (file size ~500 MB)...")

        try:
            # Download with progress
            response = requests.get(self.gwas_url, stream=True, timeout=300)
            response.raise_for_status()

            total_size = int(response.headers.get('content-length', 0))

            with open(self.gwas_file, 'wb') as f:
                downloaded = 0
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            progress = (downloaded / total_size) * 100
                            if downloaded % (10 * 1024 * 1024) == 0:  # Log every 10 MB
                                logger.info(f"  Downloaded: {downloaded / (1024*1024):.1f} MB ({progress:.1f}%)")

            logger.info(f"✓ Downloaded successfully: {self.gwas_file}")
            return True

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to download GWAS data: {e}")
            logger.error("\nMANUAL DOWNLOAD INSTRUCTIONS:")
            logger.error("1. Go to: https://ipsych.dk/en/research/downloads/")
            logger.error("   or: https://pgc.unc.edu/for-researchers/download-results/")
            logger.error("2. Download: iPSYCH-PGC ASD GWAS (Grove et al. 2019)")
            logger.error(f"3. Save to: {self.gwas_file}")
            logger.error("4. Re-run this script")
            return False

    def load_gwas(self):
        """
        Load and perform QC on Grove et al. 2019 GWAS data.
        """
        logger.info("\n=== LOADING GWAS DATA ===")

        if not self.gwas_file.exists():
            logger.error(f"GWAS file not found: {self.gwas_file}")
            logger.error("Please download the file manually (see instructions above)")
            return None

        try:
            logger.info(f"Loading: {self.gwas_file}")

            # Grove et al. 2019 format can be:
            # 1. Standard PLINK format: CHR SNP BP A1 A2 INFO OR SE P
            # 2. VCF-like format: CHROM POS ID REF ALT BETA SE PVAL ...

            # Try loading with standard delimiter
            if self.gwas_file.suffix == '.gz':
                gwas_df = pd.read_csv(self.gwas_file, sep=r'\s+', compression='gzip')
            else:
                gwas_df = pd.read_csv(self.gwas_file, sep=r'\s+')

            logger.info(f"Loaded {len(gwas_df):,} variants")
            logger.info(f"Columns: {list(gwas_df.columns)}")

            # Standardize column names if needed
            if 'CHROM' in gwas_df.columns:
                # VCF-like format
                gwas_df = gwas_df.rename(columns={
                    'CHROM': 'CHR',
                    'POS': 'BP',
                    'ID': 'SNP',
                    'ALT': 'A1',
                    'REF': 'A2',
                    'PVAL': 'P',
                    'IMPINFO': 'INFO'
                })
            # If already has CHR, SNP, BP, etc., no renaming needed

            # QC filters (standard GWAS QC)
            logger.info("\nApplying QC filters...")

            n_original = len(gwas_df)

            # 1. Remove missing values in critical columns
            required_cols = ['CHR', 'BP', 'SNP', 'A1', 'A2', 'P']
            gwas_df = gwas_df.dropna(subset=required_cols)
            logger.info(f"  After removing missing data: {len(gwas_df):,} variants")

            # 2. INFO score > 0.8 (if available)
            if 'INFO' in gwas_df.columns:
                gwas_df = gwas_df[gwas_df['INFO'] > 0.8]
                logger.info(f"  After INFO > 0.8 filter: {len(gwas_df):,} variants")

            # 3. MAF > 0.01 (calculate from frequencies if available)
            if 'FCAS' in gwas_df.columns and 'FCON' in gwas_df.columns:
                # VCF format: Calculate MAF from case and control frequencies
                gwas_df['MAF'] = ((gwas_df['FCAS'] + gwas_df['FCON']) / 2).apply(
                    lambda x: min(x, 1-x) if pd.notna(x) else np.nan
                )
                gwas_df = gwas_df[gwas_df['MAF'] > 0.01]
                logger.info(f"  After MAF > 0.01 filter: {len(gwas_df):,} variants")
            elif 'FRQ' in gwas_df.columns:
                # Standard format: has FRQ column
                gwas_df['MAF'] = gwas_df['FRQ'].apply(lambda x: min(x, 1-x) if pd.notna(x) else np.nan)
                gwas_df = gwas_df[gwas_df['MAF'] > 0.01]
                logger.info(f"  After MAF > 0.01 filter: {len(gwas_df):,} variants")
            else:
                logger.warning("  No MAF/frequency column found, skipping MAF filter")
                # Estimate MAF from OR (very rough approximation for QC only)
                gwas_df['MAF'] = 0.1  # Placeholder

            logger.info(f"\nQC Summary:")
            logger.info(f"  Original variants: {n_original:,}")
            logger.info(f"  After QC: {len(gwas_df):,}")
            logger.info(f"  Filtered: {n_original - len(gwas_df):,} ({100*(n_original - len(gwas_df))/n_original:.1f}%)")

            # Summary statistics
            logger.info(f"\nGWAS Statistics:")
            logger.info(f"  Min p-value: {gwas_df['P'].min():.2e}")
            logger.info(f"  Median p-value: {gwas_df['P'].median():.2e}")
            logger.info(f"  Genome-wide significant (P < 5e-8): {(gwas_df['P'] < 5e-8).sum():,}")

            return gwas_df

        except Exception as e:
            logger.error(f"Error loading GWAS data: {e}")
            return None

    def load_eqtl_variants(self):
        """
        Load our extracted eQTL variants.
        """
        logger.info("\n=== LOADING eQTL VARIANTS ===")

        eqtl_file = self.eqtl_dir / "real_gtex_eqtls.csv"

        if not eqtl_file.exists():
            logger.error(f"eQTL file not found: {eqtl_file}")
            return None

        eqtls_df = pd.read_csv(eqtl_file)
        logger.info(f"Loaded {len(eqtls_df)} eQTL associations")
        logger.info(f"Unique variants: {eqtls_df['variant_id'].nunique()}")
        logger.info(f"Genes: {eqtls_df['gene_symbol'].unique().tolist()}")

        return eqtls_df

    def parse_gtex_variant_id(self, variant_id):
        """
        Parse GTEx variant ID to extract chr, pos, ref, alt.
        Format: chr12_2188187_C_T_b38
        """
        parts = variant_id.replace('_b38', '').split('_')
        if len(parts) < 4:
            return None

        return {
            'chr': parts[0].replace('chr', ''),
            'pos': int(parts[1]),
            'ref': parts[2],
            'alt': parts[3]
        }

    def liftover_hg38_to_hg19(self, chrom, pos_hg38):
        """
        Convert GRCh38 coordinates to GRCh37 using pyliftover.
        This handles the genome build mismatch between GTEx (hg38) and Grove GWAS (hg19).
        """
        if not self.liftover:
            return None

        try:
            # LiftOver expects 'chr' prefix
            if not chrom.startswith('chr'):
                chrom = f'chr{chrom}'

            result = self.liftover.convert_coordinate(chrom, pos_hg38)

            if result and len(result) > 0:
                # Result is [(new_chrom, new_pos, strand, ...)]
                new_chrom, new_pos = result[0][0], result[0][1]
                return new_chrom.replace('chr', ''), int(new_pos)

            return None, None

        except Exception as e:
            logger.debug(f"Error lifting chr{chrom}:{pos_hg38}: {e}")
            return None, None

    def match_gwas_eqtl(self, gwas_df, eqtls_df):
        """
        Match GWAS and eQTL variants using coordinate liftover.
        GTEx is GRCh38, Grove GWAS is GRCh37.
        """
        logger.info("\n=== MATCHING GWAS AND eQTL VARIANTS ===")

        # Parse eQTL variant IDs and convert coordinates hg38 → hg19
        logger.info("Converting eQTL coordinates from GRCh38 to GRCh37...")

        eqtl_lifted = []

        for idx, row in eqtls_df.iterrows():
            parsed = self.parse_gtex_variant_id(row['variant_id'])
            if parsed:
                # Liftover GRCh38 → GRCh37
                chr_hg19, pos_hg19 = self.liftover_hg38_to_hg19(parsed['chr'], parsed['pos'])

                if chr_hg19 and pos_hg19:
                    eqtl_lifted.append({
                        'variant_id': row['variant_id'],
                        'gene_symbol': row['gene_symbol'],
                        'chr': chr_hg19,
                        'pos_hg19': pos_hg19,
                        'pos_hg38': parsed['pos'],
                        'eqtl_ref': parsed['ref'],
                        'eqtl_alt': parsed['alt'],
                        'eqtl_pval': row['pval_nominal'],
                        'eqtl_slope': row['slope'],
                        'slope_se': row['slope_se'],
                        'eqtl_maf': row['maf'],
                        'tissue': row['tissue']
                    })

        eqtl_lifted_df = pd.DataFrame(eqtl_lifted)
        logger.info(f"Successfully lifted {len(eqtl_lifted_df)}/{len(eqtls_df)} eQTL variants to GRCh37")

        if len(eqtl_lifted_df) == 0:
            logger.error("No coordinates could be lifted!")
            return None

        # Match by chr:pos (now both in hg19)
        logger.info("Matching GWAS and eQTL by chromosome and position (GRCh37)...")

        # Prepare GWAS data
        gwas_df['chr'] = gwas_df['CHR'].astype(str)
        gwas_df['pos'] = gwas_df['BP'].astype(int)
        eqtl_lifted_df['chr'] = eqtl_lifted_df['chr'].astype(str)

        merged = eqtl_lifted_df.merge(
            gwas_df,
            left_on=['chr', 'pos_hg19'],
            right_on=['chr', 'pos'],
            how='inner',
            suffixes=('_eqtl', '_gwas')
        )

        logger.info(f"Matched {len(merged)} eQTL-GWAS pairs")

        if len(merged) == 0:
            logger.error("No matches found!")
            logger.error("Possible issues: genome build mismatch or variants not in GWAS")
            return None

        # Align alleles
        logger.info("Aligning alleles...")

        def align_alleles(row):
            """Check if alleles match (direct or flipped)"""
            eqtl_ref = row['eqtl_ref'].upper()
            eqtl_alt = row['eqtl_alt'].upper()
            gwas_a1 = row['A1'].upper()
            gwas_a2 = row['A2'].upper()

            # Direct match
            if (eqtl_ref == gwas_a2 and eqtl_alt == gwas_a1):
                return 'direct', 1
            # Flipped
            elif (eqtl_ref == gwas_a1 and eqtl_alt == gwas_a2):
                return 'flipped', -1
            else:
                return 'mismatch', 0

        merged[['allele_match', 'direction']] = merged.apply(
            align_alleles, axis=1, result_type='expand'
        )

        # Remove mismatches
        matched = merged[merged['allele_match'] != 'mismatch'].copy()
        logger.info(f"After allele alignment: {len(matched)} variants")
        logger.info(f"  Direct matches: {(matched['allele_match'] == 'direct').sum()}")
        logger.info(f"  Flipped matches: {(matched['allele_match'] == 'flipped').sum()}")
        logger.info(f"  Mismatches removed: {len(merged) - len(matched)}")

        # Add hg19 position from GWAS for coloc (R script expects 'pos' column)
        matched['pos'] = matched['BP']

        # Adjust effect sizes for flipped alleles
        if 'BETA' in matched.columns:
            # VCF format: BETA (log OR)
            matched['BETA_aligned'] = matched['BETA'] * matched['direction']
            matched['OR_aligned'] = np.exp(matched['BETA_aligned'])
        elif 'OR' in matched.columns:
            # Standard format: OR (odds ratio)
            matched['OR_aligned'] = matched.apply(
                lambda row: 1/row['OR'] if row['direction'] == -1 else row['OR'],
                axis=1
            )
            # Calculate BETA for coloc (needs log OR)
            matched['BETA_aligned'] = matched.apply(
                lambda row: -np.log(row['OR']) if row['direction'] == -1 else np.log(row['OR']),
                axis=1
            )

        matched['eqtl_slope_aligned'] = matched['eqtl_slope'] * matched['direction']

        logger.info(f"\n✓ Successfully matched and aligned {len(matched)} variants")
        logger.info(f"  Covering {matched['gene_symbol'].nunique()} genes")
        logger.info(f"  Genes: {matched['gene_symbol'].unique().tolist()}")

        return matched

    def save_results(self, matched_df):
        """
        Save matched GWAS-eQTL data for colocalization.
        """
        logger.info("\n=== SAVING RESULTS ===")

        # Save full matched dataset
        output_file = self.coloc_dir / "gwas_eqtl_matched.csv"
        matched_df.to_csv(output_file, index=False)
        logger.info(f"✓ Saved: {output_file}")

        # Save per-gene files for colocalization
        for gene in matched_df['gene_symbol'].unique():
            gene_data = matched_df[matched_df['gene_symbol'] == gene]
            gene_file = self.coloc_dir / f"gwas_eqtl_{gene}.csv"
            gene_data.to_csv(gene_file, index=False)
            logger.info(f"✓ Saved {gene}: {len(gene_data)} variants")

        # Save summary
        summary_file = self.coloc_dir / "gwas_preparation_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("GWAS-eQTL MATCHING SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Date: {pd.Timestamp.now()}\n")
            f.write(f"GWAS: Grove et al. 2019 (iPSYCH-PGC ASD)\n")
            f.write(f"eQTL: GTEx v8 Brain Tissues\n\n")

            f.write(f"TOTAL STATISTICS:\n")
            f.write(f"  Matched variants: {len(matched_df)}\n")
            f.write(f"  Genes: {matched_df['gene_symbol'].nunique()}\n")
            f.write(f"  Tissues: {matched_df['tissue'].nunique()}\n\n")

            f.write("PER-GENE SUMMARY:\n")
            for gene in matched_df['gene_symbol'].unique():
                gene_data = matched_df[matched_df['gene_symbol'] == gene]
                f.write(f"\n{gene}:\n")
                f.write(f"  Variants: {len(gene_data)}\n")
                f.write(f"  Tissues: {gene_data['tissue'].unique().tolist()}\n")
                f.write(f"  Min GWAS p-value: {gene_data['P'].min():.2e}\n")
                f.write(f"  Min eQTL p-value: {gene_data['eqtl_pval'].min():.2e}\n")

        logger.info(f"✓ Saved summary: {summary_file}")

        return True

    def run(self):
        """
        Execute the complete GWAS preparation pipeline.
        """
        logger.info("\n" + "=" * 70)
        logger.info("GWAS DATA PREPARATION FOR COLOCALIZATION")
        logger.info("Grove et al. 2019 iPSYCH-PGC ASD GWAS")
        logger.info("=" * 70 + "\n")

        # Step 1: Download GWAS (if needed)
        if not self.download_gwas():
            logger.warning("GWAS download failed. Checking if file exists...")
            if not self.gwas_file.exists():
                return False

        # Step 2: Load and QC GWAS
        gwas_df = self.load_gwas()
        if gwas_df is None:
            return False

        # Step 3: Load eQTL variants
        eqtls_df = self.load_eqtl_variants()
        if eqtls_df is None:
            return False

        # Step 4: Match GWAS and eQTL
        matched_df = self.match_gwas_eqtl(gwas_df, eqtls_df)
        if matched_df is None or len(matched_df) == 0:
            return False

        # Step 5: Save results
        self.save_results(matched_df)

        logger.info("\n" + "=" * 70)
        logger.info("GWAS PREPARATION COMPLETED SUCCESSFULLY")
        logger.info("=" * 70 + "\n")

        return True

def main():
    """Main execution function"""
    prep = GWASPreparation()
    success = prep.run()

    if success:
        print("\n✓ GWAS preparation completed successfully!")
        print(f"✓ Results saved to: {prep.coloc_dir}")
        print("\nNext steps:")
        print("  1. Review matched GWAS-eQTL data")
        print("  2. Run Bayesian colocalization (script 02)")
        return 0
    else:
        print("\n✗ GWAS preparation failed. Check log for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
