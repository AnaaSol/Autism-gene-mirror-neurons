#!/usr/bin/env python3
"""
MEJORAS ESTADÍSTICAS PARA ANÁLISIS eQTL
Implementa correcciones FDR, intervalos de confianza, poder estadístico,
validación cruzada y análisis de sensibilidad
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import r2_score
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class StatisticalImprovements:
    
    def __init__(self, significance_level=0.05):
        self.alpha = significance_level
        self.results = {}
        
    def apply_multiple_testing_corrections(self, pvalues, methods=['fdr_bh', 'bonferroni', 'holm']):
        """
        Aplicar múltiples correcciones por comparaciones múltiples
        """
        corrections = {}
        
        for method in methods:
            if method == 'bonferroni':
                # Manual Bonferroni
                corrected_p = np.minimum(pvalues * len(pvalues), 1.0)
                significant = corrected_p < self.alpha
            else:
                # Use statsmodels
                significant, corrected_p, alpha_sidak, alpha_bonf = multipletests(
                    pvalues, method=method, alpha=self.alpha
                )
            
            corrections[method] = {
                'corrected_pvalues': corrected_p,
                'significant': significant,
                'n_significant': significant.sum(),
                'prop_significant': significant.mean()
            }
        
        return corrections
    
    def calculate_confidence_intervals(self, estimates, standard_errors, df=None, alpha=0.05):
        """
        Calcular intervalos de confianza para estimaciones
        """
        if df is None:
            # Use normal distribution for large samples
            t_crit = stats.norm.ppf(1 - alpha/2)
        else:
            # Use t-distribution
            t_crit = stats.t.ppf(1 - alpha/2, df)
        
        margin_error = t_crit * standard_errors
        
        ci_lower = estimates - margin_error
        ci_upper = estimates + margin_error
        
        return ci_lower, ci_upper
    
    def bootstrap_confidence_intervals(self, data, statistic_func, n_bootstrap=1000, alpha=0.05):
        """
        Intervalos de confianza por bootstrapping
        """
        bootstrap_stats = []
        
        for _ in range(n_bootstrap):
            # Resample with replacement
            bootstrap_sample = np.random.choice(data, size=len(data), replace=True)
            stat = statistic_func(bootstrap_sample)
            bootstrap_stats.append(stat)
        
        bootstrap_stats = np.array(bootstrap_stats)
        
        # Percentile method
        lower_percentile = (alpha/2) * 100
        upper_percentile = (1 - alpha/2) * 100
        
        ci_lower = np.percentile(bootstrap_stats, lower_percentile)
        ci_upper = np.percentile(bootstrap_stats, upper_percentile)
        
        return ci_lower, ci_upper, bootstrap_stats
    
    def power_analysis_correlation(self, effect_size, sample_size, alpha=0.05):
        """
        Análisis de poder para correlaciones
        """
        # Fisher z-transformation
        z_effect = 0.5 * np.log((1 + effect_size) / (1 - effect_size))
        se_z = 1 / np.sqrt(sample_size - 3)
        
        # Critical value
        z_alpha = stats.norm.ppf(1 - alpha/2)
        
        # Power calculation
        z_beta = abs(z_effect) / se_z - z_alpha
        power = stats.norm.cdf(z_beta)
        
        return power
    
    def power_analysis_regression(self, r_squared, sample_size, n_predictors=1, alpha=0.05):
        """
        Análisis de poder para regresión
        """
        # Effect size f²
        f_squared = r_squared / (1 - r_squared)
        
        # Non-centrality parameter
        lambda_nc = f_squared * sample_size
        
        # Critical F value
        df1 = n_predictors
        df2 = sample_size - n_predictors - 1
        f_crit = stats.f.ppf(1 - alpha, df1, df2)
        
        # Power using non-central F distribution
        power = 1 - stats.ncf.cdf(f_crit, df1, df2, lambda_nc)
        
        return power
    
    def sample_size_for_power(self, effect_size, desired_power=0.8, alpha=0.05, test_type='correlation'):
        """
        Calcular tamaño muestral necesario para poder deseado
        """
        if test_type == 'correlation':
            # Fisher z-transformation
            z_effect = 0.5 * np.log((1 + effect_size) / (1 - effect_size))
            z_alpha = stats.norm.ppf(1 - alpha/2)
            z_beta = stats.norm.ppf(desired_power)
            
            # Required sample size
            n_required = ((z_alpha + z_beta) / z_effect) ** 2 + 3
            
        elif test_type == 'regression':
            # Approximation for regression
            f_squared = effect_size / (1 - effect_size)
            z_alpha = stats.norm.ppf(1 - alpha/2)
            z_beta = stats.norm.ppf(desired_power)
            
            n_required = ((z_alpha + z_beta) ** 2) / f_squared + 3
        
        return int(np.ceil(n_required))
    
    def cross_validation_analysis(self, X, y, model_func, cv_folds=5, random_state=42):
        """
        Validación cruzada para evaluar robustez
        """
        kf = KFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
        
        cv_scores = []
        cv_predictions = []
        cv_residuals = []
        
        for train_idx, test_idx in kf.split(X):
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
            
            # Fit model
            model = model_func(X_train, y_train)
            
            # Predict
            y_pred = model.predict(X_test)
            
            # Calculate metrics
            r2 = r2_score(y_test, y_pred)
            cv_scores.append(r2)
            
            # Store predictions and residuals
            cv_predictions.extend(y_pred)
            cv_residuals.extend(y_test - y_pred)
        
        cv_results = {
            'scores': cv_scores,
            'mean_score': np.mean(cv_scores),
            'std_score': np.std(cv_scores),
            'predictions': cv_predictions,
            'residuals': cv_residuals
        }
        
        return cv_results
    
    def sensitivity_analysis_thresholds(self, pvalues, effect_sizes, thresholds=[0.05, 0.01, 0.001, 1e-4, 1e-5]):
        """
        Análisis de sensibilidad para diferentes umbrales de significancia
        """
        sensitivity_results = {}
        
        for threshold in thresholds:
            significant_mask = pvalues < threshold
            n_significant = significant_mask.sum()
            
            if n_significant > 0:
                mean_effect = effect_sizes[significant_mask].mean()
                median_effect = effect_sizes[significant_mask].median()
                min_effect = effect_sizes[significant_mask].min()
                max_effect = effect_sizes[significant_mask].max()
            else:
                mean_effect = median_effect = min_effect = max_effect = np.nan
            
            sensitivity_results[threshold] = {
                'n_significant': n_significant,
                'proportion': n_significant / len(pvalues),
                'mean_effect_size': mean_effect,
                'median_effect_size': median_effect,
                'min_effect_size': min_effect,
                'max_effect_size': max_effect
            }
        
        return sensitivity_results
    
    def outlier_analysis(self, data, method='iqr', threshold=1.5):
        """
        Análisis de outliers y su impacto
        """
        if method == 'iqr':
            Q1 = data.quantile(0.25)
            Q3 = data.quantile(0.75)
            IQR = Q3 - Q1
            
            lower_bound = Q1 - threshold * IQR
            upper_bound = Q3 + threshold * IQR
            
            outliers_mask = (data < lower_bound) | (data > upper_bound)
            
        elif method == 'zscore':
            z_scores = np.abs(stats.zscore(data))
            outliers_mask = z_scores > threshold
        
        outlier_results = {
            'outliers_mask': outliers_mask,
            'n_outliers': outliers_mask.sum(),
            'proportion_outliers': outliers_mask.mean(),
            'outlier_values': data[outliers_mask],
            'clean_data': data[~outliers_mask]
        }
        
        return outlier_results
    
    def effect_size_interpretation(self, effect_sizes, context='correlation'):
        """
        Interpretación de tamaños de efecto según estándares Cohen
        """
        if context == 'correlation':
            # Cohen's conventions for correlations
            small_threshold = 0.1
            medium_threshold = 0.3
            large_threshold = 0.5
        elif context == 'regression':
            # Cohen's f² for regression
            small_threshold = 0.02
            medium_threshold = 0.15
            large_threshold = 0.35
        
        abs_effects = np.abs(effect_sizes)
        
        interpretation = pd.cut(
            abs_effects,
            bins=[0, small_threshold, medium_threshold, large_threshold, np.inf],
            labels=['Negligible', 'Small', 'Medium', 'Large'],
            include_lowest=True
        )
        
        interpretation_summary = interpretation.value_counts()
        
        return interpretation, interpretation_summary
    
    def generate_statistical_report(self, analysis_results, output_file=None):
        """
        Generar reporte estadístico comprehensivo
        """
        report = []
        report.append("REPORTE DE MEJORAS ESTADÍSTICAS")
        report.append("=" * 50)
        report.append("")
        
        # Multiple testing corrections
        if 'multiple_testing' in analysis_results:
            report.append("1. CORRECCIONES POR COMPARACIONES MÚLTIPLES:")
            corrections = analysis_results['multiple_testing']
            
            for method, results in corrections.items():
                n_sig = results['n_significant']
                prop_sig = results['prop_significant'] * 100
                report.append(f"   {method.upper()}:")
                report.append(f"      Significativas: {n_sig} ({prop_sig:.1f}%)")
            report.append("")
        
        # Power analysis
        if 'power_analysis' in analysis_results:
            report.append("2. ANÁLISIS DE PODER ESTADÍSTICO:")
            power_results = analysis_results['power_analysis']
            
            report.append(f"   Poder promedio: {power_results['mean_power']:.3f}")
            report.append(f"   Poder mediano: {power_results['median_power']:.3f}")
            report.append(f"   Análisis con poder ≥0.8: {power_results['adequate_power']}")
            report.append("")
        
        # Effect sizes
        if 'effect_sizes' in analysis_results:
            report.append("3. INTERPRETACIÓN DE TAMAÑOS DE EFECTO:")
            effect_summary = analysis_results['effect_sizes']['summary']
            
            for category, count in effect_summary.items():
                report.append(f"   {category}: {count}")
            report.append("")
        
        # Cross-validation
        if 'cross_validation' in analysis_results:
            report.append("4. VALIDACIÓN CRUZADA:")
            cv_results = analysis_results['cross_validation']
            
            report.append(f"   R² promedio: {cv_results['mean_score']:.3f}")
            report.append(f"   R² std: {cv_results['std_score']:.3f}")
            report.append("")
        
        # Sensitivity analysis
        if 'sensitivity' in analysis_results:
            report.append("5. ANÁLISIS DE SENSIBILIDAD:")
            sens_results = analysis_results['sensitivity']
            
            for threshold, results in sens_results.items():
                n_sig = results['n_significant']
                report.append(f"   Umbral p<{threshold}: {n_sig} significativas")
            report.append("")
        
        report_text = "\n".join(report)
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(report_text)
        
        return report_text

def create_improved_visualizations(data, improvements_results, output_dir):
    """
    Crear visualizaciones mejoradas con estadísticas robustas
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Multiple testing corrections comparison
    if 'multiple_testing' in improvements_results:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        methods = list(improvements_results['multiple_testing'].keys())
        n_significant = [improvements_results['multiple_testing'][m]['n_significant'] 
                        for m in methods]
        
        bars = ax.bar(methods, n_significant, color=['#E74C3C', '#3498DB', '#2ECC71'])
        ax.set_ylabel('Número de Asociaciones Significativas')
        ax.set_title('Comparación de Correcciones por Comparaciones Múltiples')
        
        # Add value labels
        for bar, value in zip(bars, n_significant):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                   str(value), ha='center', va='bottom', fontweight='bold')
        
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(output_dir / "multiple_testing_comparison.png", dpi=300)
        plt.close()
    
    # 2. Power analysis visualization
    if 'power_analysis' in improvements_results:
        power_data = improvements_results['power_analysis']['individual_power']
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Power distribution
        ax1.hist(power_data, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        ax1.axvline(0.8, color='red', linestyle='--', label='Poder adecuado (0.8)')
        ax1.set_xlabel('Poder Estadístico')
        ax1.set_ylabel('Frecuencia')
        ax1.set_title('Distribución del Poder Estadístico')
        ax1.legend()
        
        # Power by gene (if available)
        if 'power_by_gene' in improvements_results['power_analysis']:
            power_by_gene = improvements_results['power_analysis']['power_by_gene']
            genes = list(power_by_gene.keys())
            powers = list(power_by_gene.values())
            
            bars = ax2.bar(genes, powers, color=['#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6'])
            ax2.axhline(0.8, color='red', linestyle='--', alpha=0.7, label='Poder adecuado')
            ax2.set_ylabel('Poder Estadístico Promedio')
            ax2.set_title('Poder Estadístico por Gen')
            ax2.set_xticklabels(genes, rotation=45)
            ax2.legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / "power_analysis.png", dpi=300)
        plt.close()
    
    # 3. Sensitivity analysis
    if 'sensitivity' in improvements_results:
        sens_data = improvements_results['sensitivity']
        
        thresholds = list(sens_data.keys())
        n_significant = [sens_data[t]['n_significant'] for t in thresholds]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot([-np.log10(t) for t in thresholds], n_significant, 
                marker='o', linewidth=2, markersize=8, color='#E74C3C')
        ax.set_xlabel('-log10(Umbral de Significancia)')
        ax.set_ylabel('Número de Asociaciones Significativas')
        ax.set_title('Análisis de Sensibilidad - Umbrales de Significancia')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / "sensitivity_analysis.png", dpi=300)
        plt.close()
    
    print(f"Visualizaciones mejoradas guardadas en: {output_dir}")

if __name__ == "__main__":
    print("MÓDULO DE MEJORAS ESTADÍSTICAS")
    print("Listo para implementar en análisis eQTL")
    print("=" * 50)