#!/usr/bin/env python3
"""
Fresh Real Network Analysis - Autism Mirror Neurons
Análisis de red basado únicamente en datos reales GTEx
Version fresh con salida actualizada
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
import logging
from itertools import combinations
from scipy.stats import pearsonr
import json

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class FreshRealNetworkAnalysis:
    """
    Análisis de red basado en datos reales de eQTLs
    """
    
    def __init__(self):
        self.base_dir = Path("/home/ana/Desktop/autism_mirror_neurons")
        self.eqtl_file = self.base_dir / "results" / "real_eqtl_results_improved" / "gtex_real_eqtl_results_improved.csv"
        self.output_dir = Path("/home/ana/Desktop/autism_mirror_neurons/results/python_fresh")
        self.output_dir.mkdir(exist_ok=True)
        self.plots_dir = Path("/home/ana/Desktop/autism_mirror_neurons/plots/python_fresh")
        
        # Network for storing connections
        self.network = nx.Graph()
        self.connections = []
        
        # Gene info
        self.gene_info = {
            'FOXP2': {'color': '#E74C3C', 'priority': 'HIGH'},
            'FOXP1': {'color': '#3498DB', 'priority': 'HIGH'},
            'CNTNAP2': {'color': '#2ECC71', 'priority': 'HIGH'},
            'CACNA1C': {'color': '#F39C12', 'priority': 'MEDIUM'},
            'CHD8': {'color': '#9B59B6', 'priority': 'MEDIUM'},
            'THEMIS': {'color': '#1ABC9C', 'priority': 'MEDIUM'}
        }
    
    def load_data(self):
        """Load real eQTL data"""
        logger.info("Loading fresh real data...")
        
        if not self.eqtl_file.exists():
            raise FileNotFoundError(f"eQTL file not found: {self.eqtl_file}")
        
        self.eqtl_df = pd.read_csv(self.eqtl_file)
        logger.info(f"Loaded {len(self.eqtl_df)} real eQTL associations")
        
        # Get gene list
        self.genes = sorted(self.eqtl_df['gene_name'].unique())
        logger.info(f"Genes found: {self.genes}")
    
    def analyze_coexpression(self):
        """Analyze co-expression patterns based on effect sizes"""
        logger.info("Analyzing co-expression...")
        
        connections = 0
        
        # Create effect size matrix by gene and tissue
        pivot_df = self.eqtl_df.pivot_table(
            index=['variant_pos', 'tissue_clean'], 
            columns='gene_name', 
            values='slope', 
            aggfunc='mean'
        ).fillna(0)
        
        # Calculate correlations between genes
        for gene1, gene2 in combinations(self.genes, 2):
            if gene1 in pivot_df.columns and gene2 in pivot_df.columns:
                
                # Get non-zero values for both genes
                mask = (pivot_df[gene1] != 0) & (pivot_df[gene2] != 0)
                if mask.sum() >= 5:  # At least 5 shared data points
                    
                    corr, p_val = pearsonr(pivot_df[gene1][mask], pivot_df[gene2][mask])
                    
                    if abs(corr) > 0.7 and p_val < 0.05:
                        
                        # Find shared tissues
                        shared_tissues = []
                        for tissue in self.eqtl_df['tissue_clean'].unique():
                            if (gene1 in self.eqtl_df[self.eqtl_df['tissue_clean'] == tissue]['gene_name'].values and 
                                gene2 in self.eqtl_df[self.eqtl_df['tissue_clean'] == tissue]['gene_name'].values):
                                shared_tissues.append(tissue)
                        
                        connection = {
                            'gene1': gene1,
                            'gene2': gene2,
                            'connection_type': 'Co-expression',
                            'weight': abs(corr) * 3,  # Scale weight
                            'evidence_strength': 'Strong' if abs(corr) > 0.8 else 'Moderate',
                            'details': f"Correlation={corr:.3f} (p={p_val:.3f}) in {len(shared_tissues)} tissues",
                            'color': '#3498DB'
                        }
                        
                        self.connections.append(connection)
                        self.network.add_edge(gene1, gene2, **connection)
                        connections += 1
        
        logger.info(f"Found {connections} co-expression connections")
        return connections
    
    def analyze_statistical_similarity(self):
        """Analyze statistical similarity between genes"""
        logger.info("Analyzing statistical similarity...")
        
        connections = 0
        
        # Create statistical profiles for each gene
        gene_profiles = {}
        for gene in self.genes:
            gene_data = self.eqtl_df[self.eqtl_df['gene_name'] == gene]
            
            if len(gene_data) > 0:
                profile = {
                    'mean_effect': gene_data['slope'].mean(),
                    'mean_significance': gene_data['neg_log10_pval'].mean(),
                    'n_tissues': gene_data['tissue_clean'].nunique(),
                    'n_eqtls': len(gene_data)
                }
                gene_profiles[gene] = profile
        
        # Compare profiles between gene pairs
        for gene1, gene2 in combinations(self.genes, 2):
            if gene1 in gene_profiles and gene2 in gene_profiles:
                
                prof1 = gene_profiles[gene1]
                prof2 = gene_profiles[gene2]
                
                # Calculate similarity score
                effect_sim = 1 - abs(prof1['mean_effect'] - prof2['mean_effect']) / 2
                sig_sim = 1 - abs(prof1['mean_significance'] - prof2['mean_significance']) / 20
                tissue_sim = min(prof1['n_tissues'], prof2['n_tissues']) / max(prof1['n_tissues'], prof2['n_tissues'])
                
                overall_sim = (effect_sim + sig_sim + tissue_sim) / 3
                
                if overall_sim > 0.6:
                    connection = {
                        'gene1': gene1,
                        'gene2': gene2,
                        'connection_type': 'Statistical Similarity',
                        'weight': overall_sim * 2,
                        'evidence_strength': 'Strong' if overall_sim > 0.8 else 'Moderate',
                        'details': f"Profile similarity={overall_sim:.3f}",
                        'color': '#F39C12'
                    }
                    
                    self.connections.append(connection)
                    self.network.add_edge(gene1, gene2, **connection)
                    connections += 1
        
        logger.info(f"Found {connections} statistical similarity connections")
        return connections
    
    def analyze_coloc_evidence(self):
        """Analyze colocalización evidence based on significance patterns"""
        logger.info("Analyzing colocalización evidence...")
        
        connections = 0
        
        # Calculate coloc evidence scores for each gene
        gene_coloc_scores = {}
        for gene in self.genes:
            gene_data = self.eqtl_df[self.eqtl_df['gene_name'] == gene]
            
            if len(gene_data) > 0:
                # Score based on significance and tissue distribution
                score = (gene_data['neg_log10_pval'].sum() / len(gene_data)) * gene_data['tissue_clean'].nunique()
                gene_coloc_scores[gene] = score
        
        # Find genes with similar evidence levels
        for gene1, gene2 in combinations(self.genes, 2):
            if gene1 in gene_coloc_scores and gene2 in gene_coloc_scores:
                
                score1 = gene_coloc_scores[gene1]
                score2 = gene_coloc_scores[gene2]
                
                # Calculate similarity in evidence strength
                if min(score1, score2) > 5:  # Both have reasonable evidence
                    
                    ratio = min(score1, score2) / max(score1, score2)
                    
                    if ratio > 0.6:  # Similar evidence levels
                        connection = {
                            'gene1': gene1,
                            'gene2': gene2,
                            'connection_type': 'Coloc Evidence',
                            'weight': ratio * 4,
                            'evidence_strength': 'Strong' if ratio > 0.8 else 'Moderate',
                            'details': f"Similar coloc evidence (scores: {score1:.1f}, {score2:.1f})",
                            'color': '#2ECC71'
                        }
                        
                        self.connections.append(connection)
                        self.network.add_edge(gene1, gene2, **connection)
                        connections += 1
        
        logger.info(f"Found {connections} colocalización connections")
        return connections
    
    def save_network_results(self):
        """Save network analysis results"""
        logger.info("Saving fresh network results...")
        
        # Convert connections to DataFrame
        if self.connections:
            connections_df = pd.DataFrame(self.connections)
            
            # Remove duplicates
            connections_df = connections_df.drop_duplicates(
                subset=['gene1', 'gene2', 'connection_type']
            )
            
            # Save to CSV
            output_file = self.output_dir / "fresh_network_connections.csv"
            connections_df.to_csv(output_file, index=False)
            logger.info(f"Saved {len(connections_df)} unique connections to {output_file}")
            
            return connections_df
        else:
            logger.warning("No connections found to save")
            return pd.DataFrame()
    
    def create_network_visualization(self):
        """Create network visualization"""
        logger.info("Creating fresh network visualization...")
        
        if len(self.network.nodes()) == 0:
            logger.warning("No network to visualize")
            return
        
        plt.figure(figsize=(14, 10))
        
        # Set up layout
        pos = nx.spring_layout(self.network, seed=42, k=3, iterations=50)
        
        # Draw nodes
        node_colors = [self.gene_info.get(node, {}).get('color', 'gray') for node in self.network.nodes()]
        node_sizes = [1000 + len(self.eqtl_df[self.eqtl_df['gene_name'] == node]) * 10 
                     for node in self.network.nodes()]
        
        nx.draw_networkx_nodes(self.network, pos, 
                              node_color=node_colors, 
                              node_size=node_sizes,
                              alpha=0.8)
        
        # Draw edges by type
        edge_types = {
            'Co-expression': '#3498DB',
            'Statistical Similarity': '#F39C12', 
            'Coloc Evidence': '#2ECC71'
        }
        
        for edge_type, color in edge_types.items():
            edges = [(u, v) for u, v, d in self.network.edges(data=True) 
                    if d.get('connection_type') == edge_type]
            
            if edges:
                weights = [self.network[u][v]['weight'] for u, v in edges]
                nx.draw_networkx_edges(self.network, pos, edgelist=edges,
                                     edge_color=color, width=weights, alpha=0.7,
                                     label=edge_type)
        
        # Draw labels
        nx.draw_networkx_labels(self.network, pos, font_size=12, font_weight='bold')
        
        plt.title('Fresh Autism Mirror Neuron Gene Network\nBased on Real GTEx eQTL Data', 
                 fontsize=16, fontweight='bold', pad=20)
        plt.legend(loc='upper right')
        plt.axis('off')
        
        # Save plot
        plt.tight_layout()
        plt.savefig(self.plots_dir / "fresh_network_visualization.png", 
                   dpi=300, bbox_inches='tight')
        plt.savefig(self.plots_dir / "fresh_network_visualization.pdf", 
                   bbox_inches='tight')
        
        logger.info("Fresh network visualization saved")
    
    def generate_network_summary(self):
        """Generate network analysis summary"""
        logger.info("Generating fresh network summary...")
        
        summary = {
            'total_genes': len(self.genes),
            'total_eqtls': len(self.eqtl_df),
            'total_connections': len(self.connections),
            'genes_in_network': len(self.network.nodes()),
            'connection_types': {},
            'gene_degrees': {}
        }
        
        # Count connection types
        for conn in self.connections:
            conn_type = conn['connection_type']
            summary['connection_types'][conn_type] = summary['connection_types'].get(conn_type, 0) + 1
        
        # Calculate node degrees
        for node in self.network.nodes():
            summary['gene_degrees'][node] = self.network.degree[node]
        
        # Save summary
        summary_file = self.output_dir / "fresh_network_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Network summary saved to {summary_file}")

        # Print summary
        print("\n" + "="*60)
        print("FRESH NETWORK ANALYSIS SUMMARY")
        print("="*60)
        print(f"Total genes analyzed: {summary['total_genes']}")
        print(f"Total eQTLs: {summary['total_eqtls']}")
        print(f"Total connections found: {summary['total_connections']}")
        print(f"Genes in network: {summary['genes_in_network']}")
        print("\nConnection types:")
        for conn_type, count in summary['connection_types'].items():
            print(f"  - {conn_type}: {count}")
        print("\nGene connectivity:")
        for gene, degree in sorted(summary['gene_degrees'].items(), key=lambda x: x[1], reverse=True):
            print(f"  - {gene}: {degree} connections")
        print("="*60)
        
        return summary
    
    def run_analysis(self):
        """Run complete fresh network analysis"""
        logger.info("Starting fresh real network analysis...")
        
        self.load_data()
        
        # Run different analysis types
        coexp_connections = self.analyze_coexpression()
        stat_connections = self.analyze_statistical_similarity()
        coloc_connections = self.analyze_coloc_evidence()
        
        # Save results
        connections_df = self.save_network_results()
        self.create_network_visualization()
        summary = self.generate_network_summary()
        
        print(f"\nFRESH CONNECTIONS GENERATED:")
        print(f"  - Co-expression: {coexp_connections} connections")
        print(f"  - Statistical Similarity: {stat_connections} connections")
        print(f"  - Coloc Evidence: {coloc_connections} connections")
        print("Fresh analysis completed!")
        
        return connections_df, summary

def main():
    print("FRESH REAL NETWORK ANALYSIS")
    print("="*60)
    print("Based on authentic GTEx eQTL data")
    print("="*60)
    
    analyzer = FreshRealNetworkAnalysis()
    connections_df, summary = analyzer.run_analysis()
    
    return connections_df, summary

if __name__ == "__main__":
    main()