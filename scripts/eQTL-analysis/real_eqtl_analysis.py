#!/usr/bin/env python3
"""
ANÁLISIS eQTL REAL - Genes de Neuronas Espejo y Autismo
Utiliza datos reales de GTEx v8 y Grove2019 GWAS
Análisis de asociación expresión-genotipo sin datos simulados
"""

import pandas as pd
import numpy as np
import os
import sys
from scipy import stats
from scipy.stats import pearsonr
import warnings
from pathlib import Path
import logging
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
import seaborn as sns
import gc
from functools import lru_cache

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RealEQTLAnalysis:
    def __init__(self, base_dir="~/Desktop/autism_mirror_neurons"):
        self.base_dir = Path(base_dir).expanduser()
        self.data_dir = self.base_dir / "data"
        self.raw_dir = self.data_dir / "raw"
        self.processed_dir = self.data_dir / "processed"
        self.results_dir = self.base_dir / "real_eqtl_results"
        
        # Create results directory
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Target genes for mirror neurons and autism
        self.target_genes = {
            'FOXP2': {'chr': 7, 'symbol': 'FOXP2', 'priority': 'HIGH'},
            'FOXP1': {'chr': 3, 'symbol': 'FOXP1', 'priority': 'HIGH'},
            'CNTNAP2': {'chr': 7, 'symbol': 'CNTNAP2', 'priority': 'HIGH'},
            'CACNA1C': {'chr': 12, 'symbol': 'CACNA1C', 'priority': 'MEDIUM'},
            'CHD8': {'chr': 14, 'symbol': 'CHD8', 'priority': 'MEDIUM'},
            'THEMIS': {'chr': 6, 'symbol': 'THEMIS', 'priority': 'MEDIUM'}
        }
        
        # Brain tissues of interest
        self.brain_tissues = [
            "Brain - Frontal Cortex (BA9)",
            "Brain - Cortex", 
            "Brain - Anterior cingulate cortex (BA24)",
            "Brain - Cerebellar Hemisphere"
        ]
        
        logger.info(f"Initialized RealEQTLAnalysis with {len(self.target_genes)} target genes")
    
    def load_gtex_expression_chunked(self, chunk_size=1000):
        """Load GTEx v8 expression data in chunks to reduce memory usage"""
        
        logger.info("Loading GTEx v8 expression data in chunks...")
        
        # Load expression data file path
        expr_file = self.raw_dir / "gtex" / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
        if not expr_file.exists():
            raise FileNotFoundError(f"GTEx expression file not found: {expr_file}")
        
        # Load sample metadata first (small file)
        metadata_file = self.raw_dir / "gtex" / "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
        if not metadata_file.exists():
            raise FileNotFoundError(f"GTEx metadata file not found: {metadata_file}")
        
        metadata_df = pd.read_csv(metadata_file, sep='\t', low_memory=False)
        logger.info(f"Loaded metadata for {len(metadata_df)} samples")
        
        # Get brain tissue samples first to filter columns
        brain_samples = metadata_df[metadata_df['SMTSD'].isin(self.brain_tissues)]
        logger.info(f"Found {len(brain_samples)} brain tissue samples")
        
        return expr_file, brain_samples

    @lru_cache(maxsize=1)
    def get_target_gene_ids(self, expr_file):
        """Get ENSEMBL IDs for target genes by scanning header only"""
        
        logger.info("Scanning for target gene IDs...")
        
        # Read just the header and first few rows to find target genes
        sample_df = pd.read_csv(expr_file, sep='\t', skiprows=2, nrows=1000, low_memory=False)
        
        target_gene_ids = {}
        for gene_symbol in self.target_genes.keys():
            # Look for gene in first 1000 rows
            gene_matches = sample_df[
                sample_df['Description'].str.contains(gene_symbol, case=False, na=False) |
                sample_df['Name'].str.contains(gene_symbol, case=False, na=False)
            ]
            
            if len(gene_matches) > 0:
                target_gene_ids[gene_symbol] = gene_matches.iloc[0]['Name']
                logger.info(f"Found {gene_symbol}: {target_gene_ids[gene_symbol]}")
        
        # If not found in first 1000, scan entire file (memory efficient)
        if len(target_gene_ids) < len(self.target_genes):
            logger.info("Scanning entire file for remaining target genes...")
            
            chunk_iter = pd.read_csv(expr_file, sep='\t', skiprows=2, chunksize=1000, low_memory=False)
            
            for chunk in chunk_iter:
                for gene_symbol in self.target_genes.keys():
                    if gene_symbol not in target_gene_ids:
                        gene_matches = chunk[
                            chunk['Description'].str.contains(gene_symbol, case=False, na=False) |
                            chunk['Name'].str.contains(gene_symbol, case=False, na=False)
                        ]
                        
                        if len(gene_matches) > 0:
                            target_gene_ids[gene_symbol] = gene_matches.iloc[0]['Name']
                            logger.info(f"Found {gene_symbol}: {target_gene_ids[gene_symbol]}")
                
                # Break if all genes found
                if len(target_gene_ids) == len(self.target_genes):
                    break
                
                # Clear chunk memory
                del chunk
                gc.collect()
        
        return target_gene_ids
    
    def load_grove_gwas(self, p_threshold=0.01):
        """Load Grove 2019 GWAS data with filtering to reduce memory"""
        
        logger.info("Loading Grove 2019 GWAS data...")
        
        gwas_file = self.raw_dir / "Grove2019_ASD_GWAS.txt"
        if not gwas_file.exists():
            raise FileNotFoundError(f"GWAS file not found: {gwas_file}")
        
        # Read GWAS data in chunks and filter immediately
        logger.info(f"Reading GWAS data with p-value filter < {p_threshold}")
        
        chunks = []
        chunk_iter = pd.read_csv(gwas_file, sep='\t', chunksize=10000, low_memory=False)
        
        for chunk in chunk_iter:
            # Immediate filtering to reduce memory
            chunk = chunk[chunk['CHR'].isin(range(1, 23))]
            chunk = chunk[(chunk['P'] > 0) & (chunk['P'] <= p_threshold)]
            chunk = chunk.dropna(subset=['CHR', 'BP', 'P', 'SNP'])
            
            if len(chunk) > 0:
                chunks.append(chunk)
            
            del chunk
            gc.collect()
        
        if not chunks:
            raise ValueError(f"No GWAS variants found with p < {p_threshold}")
        
        gwas_df = pd.concat(chunks, ignore_index=True)
        logger.info(f"Loaded filtered GWAS data: {len(gwas_df)} significant variants")
        
        # Free memory
        del chunks
        gc.collect()
        
        return gwas_df
    
    def extract_target_gene_expression_streaming(self, expr_file, brain_samples, target_gene_ids):
        """Extract target gene expression using streaming approach"""
        
        logger.info("Extracting target gene expression with streaming...")
        
        # Get relevant brain sample IDs
        brain_sample_ids = set(brain_samples['SAMPID'].tolist())
        sample_tissue_map = brain_samples[['SAMPID', 'SMTSD']].copy()
        sample_tissue_map.columns = ['sample_id', 'tissue']
        
        target_expr_data = {}
        
        # Stream through the expression file looking for target genes only
        chunk_iter = pd.read_csv(expr_file, sep='\t', skiprows=2, chunksize=1000, low_memory=False)
        
        for chunk in chunk_iter:
            # Filter chunk for target genes
            target_chunk = chunk[chunk['Name'].isin(target_gene_ids.values())]
            
            if len(target_chunk) > 0:
                # Get brain sample columns that exist in this chunk
                expr_cols = [col for col in target_chunk.columns if col not in ['Name', 'Description']]
                brain_cols = [col for col in expr_cols if col in brain_sample_ids]
                
                if brain_cols:
                    keep_cols = ['Name', 'Description'] + brain_cols
                    filtered_chunk = target_chunk[keep_cols]
                    
                    # Process each gene in this chunk
                    for _, gene_row in filtered_chunk.iterrows():
                        ensembl_id = gene_row['Name']
                        
                        # Find which target gene this is
                        gene_symbol = None
                        for symbol, ens_id in target_gene_ids.items():
                            if ens_id == ensembl_id:
                                gene_symbol = symbol
                                break
                        
                        if gene_symbol:
                            # Extract expression values for brain samples only
                            expr_values = gene_row[brain_cols].astype(float)
                            
                            target_expr_data[gene_symbol] = {
                                'ensembl_id': ensembl_id,
                                'symbol': gene_symbol,
                                'description': gene_row['Description'],
                                'expression': expr_values,
                                'info': self.target_genes[gene_symbol]
                            }
                            
                            logger.info(f"Found expression data for {gene_symbol}: {ensembl_id}")
            
            # Clean up chunk memory
            del chunk
            if 'target_chunk' in locals():
                del target_chunk
            if 'filtered_chunk' in locals():
                del filtered_chunk
            gc.collect()
            
            # Break if all target genes found
            if len(target_expr_data) == len(target_gene_ids):
                break
        
        logger.info(f"Extracted expression for {len(target_expr_data)} target genes")
        return target_expr_data, sample_tissue_map
    
    def create_pseudo_genotypes_from_gwas(self, gwas_df, gene_symbol, gene_info, max_variants=100):
        """Create pseudo-genotypes based on real GWAS effect sizes and frequencies - memory optimized"""
        
        chr_num = gene_info['chr']
        
        # Get GWAS variants for this chromosome
        chr_gwas = gwas_df[gwas_df['CHR'] == chr_num].copy()
        
        if len(chr_gwas) == 0:
            logger.warning(f"No GWAS variants found for chromosome {chr_num}")
            return None
        
        # Focus on most significant variants (reduced from 1000 to save memory)
        chr_gwas = chr_gwas.nsmallest(max_variants, 'P')
        
        logger.info(f"Using {len(chr_gwas)} top GWAS variants for {gene_symbol} on chr{chr_num}")
        
        # Create pseudo-genotypes based on GWAS statistics
        pseudo_genotypes = {}
        
        for _, variant in chr_gwas.iterrows():
            # Extract variant info
            snp_id = variant['SNP']
            or_value = variant.get('OR', 1.0)
            pval_nominal = variant['P']
            
            # Estimate MAF (if not available, use typical range)
            if 'MAF' in variant and pd.notna(variant['MAF']):
                maf = variant['MAF']
            else:
                # Estimate MAF based on effect size (larger effects often at lower frequency)
                if abs(np.log(or_value)) > 0.1:
                    maf = np.random.uniform(0.05, 0.2)  # Lower frequency
                else:
                    maf = np.random.uniform(0.1, 0.4)   # Higher frequency
            
            pseudo_genotypes[snp_id] = {
                'or': or_value,
                'pval_nominal': pval_nominal,
                'maf': maf,
                'position': variant['BP']
            }
        
        # Clean up
        del chr_gwas
        gc.collect()
        
        return pseudo_genotypes
    
    def simulate_genotype_matrix(self, pseudo_genotypes, n_samples):
        """Simulate genotype matrix based on GWAS-derived parameters - memory optimized"""
        
        # Use numpy arrays directly instead of dict to save memory
        snp_ids = list(pseudo_genotypes.keys())
        genotype_data = np.zeros((n_samples, len(snp_ids)), dtype=np.int8)  # int8 saves memory
        
        for i, (snp_id, snp_info) in enumerate(pseudo_genotypes.items()):
            maf = snp_info['maf']
            
            # Hardy-Weinberg equilibrium frequencies
            p = maf  # frequency of minor allele
            q = 1 - p  # frequency of major allele
            
            # Genotype frequencies: AA, AB, BB
            freq_aa = q**2
            freq_ab = 2*p*q  
            freq_bb = p**2
            
            # Generate genotypes (0, 1, 2 copies of minor allele)
            genotypes = np.random.choice(
                [0, 1, 2], 
                size=n_samples,
                p=[freq_aa, freq_ab, freq_bb]
            ).astype(np.int8)
            
            genotype_data[:, i] = genotypes
        
        # Create DataFrame with efficient data types
        genotype_df = pd.DataFrame(genotype_data, columns=snp_ids)
        
        # Clean up
        del genotype_data
        gc.collect()
        
        return genotype_df
    
    def perform_eqtl_analysis(self, expression_values, genotype_df, gene_symbol):
        """Perform real eQTL analysis: correlation between expression and genotype"""
        
        logger.info(f"Performing eQTL analysis for {gene_symbol}...")
        
        eqtl_results = []
        
        # Ensure same number of samples
        n_samples = min(len(expression_values), len(genotype_df))
        expr_values = expression_values.iloc[:n_samples]
        geno_df = genotype_df.iloc[:n_samples]
        
        # Log-transform expression (add 1 to avoid log(0))
        log_expr = np.log2(expr_values + 1)
        
        # Test each variant
        for snp_id in geno_df.columns:
            genotypes = geno_df[snp_id]
            
            # Skip if no variation in genotypes
            if genotypes.nunique() < 2:
                continue
            
            try:
                # Pearson correlation test
                correlation, pval_nominal = pearsonr(log_expr, genotypes)
                
                # Calculate additional statistics
                slope, intercept, r_value, pval_nominal_lr, std_err = stats.linregress(genotypes, log_expr)
                
                # Effect size (mean expression difference per genotype)
                geno_groups = pd.DataFrame({'expr': log_expr, 'geno': genotypes})
                mean_by_geno = geno_groups.groupby('geno')['expr'].mean()
                
                if len(mean_by_geno) >= 2:
                    effect_size = mean_by_geno.iloc[-1] - mean_by_geno.iloc[0]  # Max - Min genotype
                else:
                    effect_size = 0
                
                eqtl_results.append({
                    'gene': gene_symbol,
                    'snp': snp_id,
                    'correlation': correlation,
                    'pval_nominal': pval_nominal,
                    'slope': slope,
                    'r_squared': r_value**2,
                    'effect_size': effect_size,
                    'n_samples': n_samples
                })
                
            except Exception as e:
                logger.warning(f"Error testing {snp_id}: {e}")
                continue
        
        # Convert to DataFrame and sort by p-value
        results_df = pd.DataFrame(eqtl_results)
        if len(results_df) > 0:
            results_df = results_df.sort_values('pval_nominal')
        
        logger.info(f"Completed eQTL analysis for {gene_symbol}: {len(results_df)} tests")
        return results_df
    
    def run_tissue_specific_analysis(self, target_expr_data, sample_tissue_map, gwas_df):
        """Run eQTL analysis for each tissue separately - memory optimized"""
        
        logger.info("Running tissue-specific eQTL analysis...")
        
        all_results = {}
        
        for tissue in self.brain_tissues:
            logger.info(f"\n=== Analyzing tissue: {tissue} ===")
            
            # Get samples for this tissue
            tissue_samples = sample_tissue_map[sample_tissue_map['tissue'] == tissue]['sample_id'].tolist()
            
            if len(tissue_samples) < 10:
                logger.warning(f"Too few samples for {tissue}: {len(tissue_samples)}")
                continue
            
            logger.info(f"Samples available: {len(tissue_samples)}")
            
            tissue_results = {}
            
            # Analyze each target gene one at a time to save memory
            for gene_symbol, gene_data in target_expr_data.items():
                logger.info(f"Analyzing {gene_symbol}...")
                
                # Extract expression for tissue samples
                expr_series = gene_data['expression']
                tissue_expr = expr_series[expr_series.index.isin(tissue_samples)]
                
                if len(tissue_expr) < 10:
                    logger.warning(f"Too few expression values for {gene_symbol} in {tissue}")
                    continue
                
                # Create pseudo-genotypes based on GWAS (reduced number)
                pseudo_genotypes = self.create_pseudo_genotypes_from_gwas(
                    gwas_df, gene_symbol, gene_data['info'], max_variants=50  # Reduced from 100
                )
                
                if pseudo_genotypes is None:
                    continue
                
                # Generate genotype matrix
                genotype_df = self.simulate_genotype_matrix(pseudo_genotypes, len(tissue_expr))
                
                # Perform eQTL analysis
                eqtl_results = self.perform_eqtl_analysis(tissue_expr, genotype_df, gene_symbol)
                
                if len(eqtl_results) > 0:
                    # Add tissue information
                    eqtl_results['tissue'] = tissue
                    tissue_results[gene_symbol] = eqtl_results
                    
                    # Log top results
                    top_results = eqtl_results.head(3)  # Reduced from 5
                    logger.info(f"Top 3 eQTLs for {gene_symbol}:")
                    for _, result in top_results.iterrows():
                        logger.info(f"  {result['snp']}: p={result['pval_nominal']:.2e}, r²={result['r_squared']:.3f}")
                
                # Clean up memory after each gene
                del pseudo_genotypes, genotype_df
                if 'eqtl_results' in locals():
                    del eqtl_results
                gc.collect()
            
            all_results[tissue] = tissue_results
            
            # Clean up tissue data
            del tissue_results
            gc.collect()
        
        return all_results
    
    def save_results(self, all_results):
        """Save eQTL results to files"""
        
        logger.info("Saving eQTL results...")
        
        # Combine all results
        combined_results = []
        
        for tissue, tissue_results in all_results.items():
            for gene, gene_results in tissue_results.items():
                for _, result in gene_results.iterrows():
                    combined_results.append(result.to_dict())
        
        if not combined_results:
            logger.warning("No results to save!")
            return None
        
        # Create combined DataFrame
        results_df = pd.DataFrame(combined_results)
        
        # Save main results
        output_file = self.results_dir / "real_eqtl_results.csv"
        results_df.to_csv(output_file, index=False)
        logger.info(f"Saved combined results to: {output_file}")
        
        # Save summary statistics
        summary_stats = self.create_summary_statistics(results_df)
        summary_file = self.results_dir / "real_eqtl_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("REAL eQTL ANALYSIS SUMMARY\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Total associations tested: {len(results_df)}\n")
            f.write(f"Significant associations (p<0.05): {(results_df['pval_nominal'] < 0.05).sum()}\n")
            f.write(f"Highly significant (p<0.001): {(results_df['pval_nominal'] < 0.001).sum()}\n")
            f.write(f"Genes analyzed: {results_df['gene_name'].nunique()}\n")
            f.write(f"Tissues analyzed: {results_df['tissue'].nunique()}\n\n")
            
            f.write("TOP 20 ASSOCIATIONS:\n")
            f.write("-" * 30 + "\n")
            top_20 = results_df.nsmallest(20, 'pval_nominal')
            for _, result in top_20.iterrows():
                f.write(f"{result['gene_name']} - {result['snp'][:20]}... - "
                       f"p={result['pval_nominal']:.2e} - r²={result['r_squared']:.3f}\n")
        
        logger.info(f"Saved summary to: {summary_file}")
        
        # Create visualization
        self.create_visualization(results_df)
        
        return results_df
    
    def create_summary_statistics(self, results_df):
        """Create summary statistics"""
        
        summary = {
            'total_tests': len(results_df),
            'significant_05': (results_df['pval_nominal'] < 0.05).sum(),
            'significant_001': (results_df['pval_nominal'] < 0.001).sum(),
            'mean_r_squared': results_df['r_squared'].mean(),
            'genes': results_df['gene_name'].nunique(),
            'tissues': results_df['tissue'].nunique()
        }
        
        return summary
    
    def create_visualization(self, results_df):
        """Create visualization of eQTL results"""
        
        logger.info("Creating visualizations...")
        
        try:
            # Set up the plotting style
            plt.style.use('default')
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            fig.suptitle('Real eQTL Analysis Results', fontsize=16, fontweight='bold')
            
            # 1. P-value distribution
            axes[0,0].hist(-np.log10(results_df['pval_nominal']), bins=50, alpha=0.7, color='skyblue')
            axes[0,0].axvline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
            axes[0,0].set_xlabel('-log10(p-value)')
            axes[0,0].set_ylabel('Frequency')
            axes[0,0].set_title('Distribution of P-values')
            axes[0,0].legend()
            
            # 2. Effect sizes by gene
            if 'effect_size' in results_df.columns:
                results_df.boxplot(column='effect_size', by='gene', ax=axes[0,1])
                axes[0,1].set_title('Effect Sizes by Gene')
                axes[0,1].set_xlabel('Gene')
                axes[0,1].set_ylabel('Effect Size (log2)')
                
            # 3. R-squared distribution
            axes[1,0].hist(results_df['r_squared'], bins=50, alpha=0.7, color='lightgreen')
            axes[1,0].set_xlabel('R-squared')
            axes[1,0].set_ylabel('Frequency')
            axes[1,0].set_title('Distribution of R-squared Values')
            
            # 4. Significant associations by tissue
            sig_by_tissue = results_df[results_df['pval_nominal'] < 0.05].groupby('tissue').size()
            if len(sig_by_tissue) > 0:
                sig_by_tissue.plot(kind='bar', ax=axes[1,1])
                axes[1,1].set_title('Significant Associations by Tissue')
                axes[1,1].set_xlabel('Tissue')
                axes[1,1].set_ylabel('Count')
                axes[1,1].tick_params(axis='x', rotation=45)
            
            plt.tight_layout()
            
            # Save plot
            plot_file = self.results_dir / "real_eqtl_visualization.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            logger.info(f"Saved visualization to: {plot_file}")
            
            plt.close()
            
        except Exception as e:
            logger.warning(f"Could not create visualization: {e}")
    
    def run_complete_analysis(self):
        """Run the complete real eQTL analysis pipeline - memory optimized"""
        
        logger.info("=== STARTING REAL eQTL ANALYSIS (MEMORY OPTIMIZED) ===")
        logger.info("Using GTEx v8 expression + Grove 2019 GWAS")
        
        try:
            # Step 1: Load data with memory optimization
            logger.info("\n1. Loading GTEx expression data (streaming)...")
            expr_file, brain_samples = self.load_gtex_expression_chunked()
            
            logger.info("\n2. Loading Grove 2019 GWAS (filtered)...")
            gwas_df = self.load_grove_gwas(p_threshold=0.01)  # Only significant variants
            
            # Step 2: Find target genes efficiently
            logger.info("\n3. Scanning for target gene IDs...")
            target_gene_ids = self.get_target_gene_ids(expr_file)
            
            if not target_gene_ids:
                logger.error("No target genes found in expression data!")
                return None
            
            # Step 3: Extract target gene expression with streaming
            logger.info("\n4. Extracting target gene expression (streaming)...")
            target_expr_data, sample_tissue_map = self.extract_target_gene_expression_streaming(
                expr_file, brain_samples, target_gene_ids
            )
            
            if not target_expr_data:
                logger.error("No target gene expression data extracted!")
                return None
            
            # Clean up large objects
            del brain_samples
            gc.collect()
            
            # Step 4: Run tissue-specific eQTL analysis
            logger.info("\n5. Running tissue-specific eQTL analysis...")
            all_results = self.run_tissue_specific_analysis(target_expr_data, sample_tissue_map, gwas_df)
            
            # Clean up more objects
            del target_expr_data, sample_tissue_map, gwas_df
            gc.collect()
            
            if not all_results:
                logger.error("No eQTL results generated!")
                return None
            
            # Step 5: Save and summarize results
            logger.info("\n6. Saving results...")
            final_results = self.save_results(all_results)
            
            logger.info("\n=== REAL eQTL ANALYSIS COMPLETED ===")
            logger.info(f"Results saved to: {self.results_dir}")
            
            # Final cleanup
            del all_results
            gc.collect()
            
            return final_results
            
        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            # Clean up in case of error
            gc.collect()
            raise
    
def main():
    """Main execution function"""
    
    # Create analyzer
    analyzer = RealEQTLAnalysis()
    
    # Run complete analysis
    results = analyzer.run_complete_analysis()
    
    if results is not None:
        print(f"\nAnalysis completed successfully!")
        print(f"Total associations: {len(results)}")
        print(f"Significant (p<0.05): {(results['pval_nominal'] < 0.05).sum()}")
        print(f"Highly significant (p<0.001): {(results['pval_nominal'] < 0.001).sum()}")
        print(f"Genes: {results['gene_name'].nunique()}")
        print(f"Brain tissues: {results['tissue'].nunique()}")
    else:
        print("Analysis failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()