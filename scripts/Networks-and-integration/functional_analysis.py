#!/usr/bin/env python3
"""
Análisis funcional detallado de las variantes top
Predicción de impacto funcional y priorización para validación
"""

import pandas as pd
import numpy as np

class FunctionalVariantAnalysis:
    def __init__(self):
        # Datos de las variantes strong del análisis
        self.strong_variants = [
            {
                'gene': 'FOXP2',
                'snp': 'rsFOXP2_7_114294028',
                'chr': 7,
                'position': 114294028,
                'score': 12.786236,
                'distance_kb': 3.537,
                'tissue': 'Brain_Frontal_Cortex_BA9',
                'priority_rank': 1
            },
            {
                'gene': 'FOXP1', 
                'snp': 'rsFOXP1_3_71110207',
                'chr': 3,
                'position': 71110207,
                'score': 11.383794,
                'distance_kb': 5.544,
                'tissue': 'Brain_Cortex',
                'priority_rank': 2
            },
            {
                'gene': 'CNTNAP2',
                'snp': 'rsCNTNAP2_7_147031685',
                'chr': 7,
                'position': 147031685,
                'score': 10.711493,
                'distance_kb': 4.769,
                'tissue': 'Brain_Anterior_Cingulate_BA24',
                'priority_rank': 3
            },
            {
                'gene': 'FOXP2',
                'snp': 'rsFOXP2_7_114290452',
                'chr': 7,
                'position': 114290452,
                'score': 9.463692,
                'distance_kb': 0.039,  # ¡Directamente en el gen!
                'tissue': 'Brain_Anterior_Cingulate_BA24',
                'priority_rank': 4
            },
            {
                'gene': 'CACNA1C',
                'snp': 'rsCACNA1C_12_2638178',
                'chr': 12,
                'position': 2638178,
                'score': 9.816281,
                'distance_kb': 7.224,
                'tissue': 'Brain_Cerebellum',
                'priority_rank': 5
            }
        ]
        
        # Coordenadas de genes (TSS)
        self.gene_coordinates = {
            'FOXP2': {'chr': 7, 'tss': 114086327, 'end': 114693772},
            'FOXP1': {'chr': 3, 'tss': 71011832, 'end': 71585540},
            'CNTNAP2': {'chr': 7, 'tss': 145818786, 'end': 147891751},
            'CACNA1C': {'chr': 12, 'tss': 2162089, 'end': 2798436}
        }
    
    def predict_regulatory_impact(self, variant):
        """Predecir impacto regulatorio de una variante"""
        
        gene = variant['gene']
        distance = variant['distance_kb'] * 1000  # Convertir a bp
        
        # Clasificación por distancia al TSS
        if distance <= 1000:
            region_type = "Promoter_Proximal"
            regulatory_score = 0.9
        elif distance <= 10000:
            region_type = "Promoter_Distal"
            regulatory_score = 0.7
        elif distance <= 100000:
            region_type = "Enhancer_Proximal"
            regulatory_score = 0.5
        else:
            region_type = "Enhancer_Distal"
            regulatory_score = 0.3
        
        # Ajustar por gene-specific factors
        gene_weights = {
            'FOXP2': 1.2,  # Alto impacto conocido
            'FOXP1': 1.1,  # Relacionado con FOXP2
            'CNTNAP2': 1.0,  # Función sináptica
            'CACNA1C': 0.9   # Canal iónico
        }
        
        regulatory_score *= gene_weights.get(gene, 1.0)
        
        return {
            'region_type': region_type,
            'regulatory_score': min(regulatory_score, 1.0),
            'predicted_mechanism': self.predict_mechanism(region_type, gene)
        }
    
    def predict_mechanism(self, region_type, gene):
        """Predecir mecanismo molecular"""
        
        mechanisms = {
            "Promoter_Proximal": [
                "Alteración de sitios de unión de factores de transcripción",
                "Modificación de estructura de cromatina",
                "Cambio en iniciación de transcripción"
            ],
            "Promoter_Distal": [
                "Regulación por enhancers proximales",
                "Modulación de bucles de cromatina",
                "Efectos en modificaciones de histonas"
            ],
            "Enhancer_Proximal": [
                "Actividad de enhancers específicos de tejido",
                "Interacciones enhancer-promoter",
                "Regulación dependiente de contexto"
            ],
            "Enhancer_Distal": [
                "Regulación a larga distancia",
                "Efectos de super-enhancers",
                "Modulación de dominios topológicos"
            ]
        }
        
        # Mecanismos específicos por gen
        gene_specific = {
            'FOXP2': "Regulación de targets downstream: CNTNAP2, SRPX2, MET",
            'FOXP1': "Co-regulación con FOXP2 en desarrollo neural",
            'CNTNAP2': "Modulación de función sináptica y conectividad",
            'CACNA1C': "Alteración de señalización de calcio neuronal"
        }
        
        base_mechanisms = mechanisms.get(region_type, ["Mecanismo desconocido"])
        specific_mechanism = gene_specific.get(gene, "Función específica no definida")
        
        return {
            'general_mechanisms': base_mechanisms,
            'gene_specific_mechanism': specific_mechanism
        }
    
    def calculate_validation_priority(self, variant):
        """Calcular prioridad de validación experimental"""
        
        # Factores de priorización
        score_weight = variant['score'] / 12.79  # Normalizado al máximo
        
        # Bonus por proximidad al gen
        distance_bonus = max(0, 1 - (variant['distance_kb'] / 50))
        
        # Bonus por relevancia del tejido
        tissue_weights = {
            'Brain_Frontal_Cortex_BA9': 1.0,  # Máxima relevancia
            'Brain_Cortex': 0.9,
            'Brain_Anterior_Cingulate_BA24': 0.8,
            'Brain_Cerebellum': 0.7
        }
        tissue_weight = tissue_weights.get(variant['tissue'], 0.5)
        
        # Bonus por evidencia previa del gen
        gene_evidence = {
            'FOXP2': 1.0,  # Evidencia máxima
            'FOXP1': 0.8,
            'CNTNAP2': 0.7,
            'CACNA1C': 0.6
        }
        evidence_weight = gene_evidence.get(variant['gene'], 0.5)
        
        # Score final
        priority_score = (score_weight * 0.4 + 
                         distance_bonus * 0.3 + 
                         tissue_weight * 0.2 + 
                         evidence_weight * 0.1)
        
        return min(priority_score, 1.0)
    
    def generate_functional_report(self):
        """Generar reporte funcional completo"""
        
        print("=" * 60)
        print("ANÁLISIS FUNCIONAL DE VARIANTES TOP")
        print("=" * 60)
        print()
        
        analysis_results = []
        
        for variant in self.strong_variants:
            print(f"VARIANTE {variant['priority_rank']}: {variant['snp']}")
            print(f"   Gen: {variant['gene']}")
            print(f"   Posición: chr{variant['chr']}:{variant['position']:,}")
            print(f"   Co-localización Score: {variant['score']:.2f}")
            print(f"   Distancia al TSS: {variant['distance_kb']:.1f} kb")
            print(f"   Tejido: {variant['tissue']}")
            
            # Análisis regulatorio
            reg_analysis = self.predict_regulatory_impact(variant)
            print(f"   Región Regulatoria: {reg_analysis['region_type']}")
            print(f"   Score Regulatorio: {reg_analysis['regulatory_score']:.2f}")
            
            # Prioridad de validación
            priority = self.calculate_validation_priority(variant)
            print(f"   Prioridad Validación: {priority:.2f}")
            
            # Mecanismo predicho
            mech = reg_analysis['predicted_mechanism']
            print(f"   Mecanismo Específico: {mech['gene_specific_mechanism']}")
            
            print()
            
            # Guardar para análisis
            analysis_results.append({
                **variant,
                'regulatory_score': reg_analysis['regulatory_score'],
                'validation_priority': priority,
                'region_type': reg_analysis['region_type']
            })
        
        # Top 3 para validación experimental
        sorted_variants = sorted(analysis_results, 
                               key=lambda x: x['validation_priority'], 
                               reverse=True)
        
        print("TOP 3 VARIANTES PARA VALIDACIÓN EXPERIMENTAL:")
        print("=" * 50)
        
        for i, variant in enumerate(sorted_variants[:3], 1):
            print(f"{i}. {variant['snp']} ({variant['gene']})")
            print(f"   Priority Score: {variant['validation_priority']:.3f}")
            print(f"   Rationale: Score alto ({variant['score']:.1f}) + "
                  f"distancia corta ({variant['distance_kb']:.1f}kb) + "
                  f"tejido relevante")
            print()
        
        # Recomendaciones experimentales
        print("RECOMENDACIONES EXPERIMENTALES:")
        print("=" * 40)
        print("1. Ensayos de actividad promotora/enhancer")
        print("2. ChIP-seq para factores de transcripción relevantes")
        print("3. Análisis de expresión alelo-específica")
        print("4. Estudios de editing genético (CRISPR)")
        print("5. Análisis de cromatina (ATAC-seq, Hi-C)")
        print()
        
        return analysis_results
    
    def export_for_validation(self, output_file="/home/ana/Desktop/Autism-gene-mirror-neurons/resultados/functional_variants_for_validation.csv"):
        """Exportar variantes priorizadas para validación"""
        
        analysis_results = []
        
        for variant in self.strong_variants:
            reg_analysis = self.predict_regulatory_impact(variant)
            priority = self.calculate_validation_priority(variant)
            
            analysis_results.append({
                'variant_id': variant['snp'],
                'gene': variant['gene'],
                'chromosome': variant['chr'],
                'position': variant['position'],
                'colocalization_score': variant['score'],
                'distance_to_tss_kb': variant['distance_kb'],
                'brain_tissue': variant['tissue'],
                'regulatory_region': reg_analysis['region_type'],
                'regulatory_score': reg_analysis['regulatory_score'],
                'validation_priority': priority,
                'recommended_assays': 'Promoter activity, ChIP-seq, CRISPR editing',
                'predicted_mechanism': reg_analysis['predicted_mechanism']['gene_specific_mechanism']
            })
        
        df = pd.DataFrame(analysis_results)
        df = df.sort_values('validation_priority', ascending=False)
        
        # Guardar archivo
        df.to_csv(output_file, index=False)
        print(f"Variantes exportadas a: {output_file}")
        
        return df

# EJECUTAR ANÁLISIS
if __name__ == "__main__":
    analyzer = FunctionalVariantAnalysis()
    results = analyzer.generate_functional_report()
    df = analyzer.export_for_validation()
    
    print("RESUMEN ESTADÍSTICO:")
    print(f"Variantes analizadas: {len(results)}")
    print(f"Score promedio: {np.mean([r['score'] for r in results]):.2f}")
    print(f"Distancia promedio: {np.mean([r['distance_kb'] for r in results]):.1f} kb")
