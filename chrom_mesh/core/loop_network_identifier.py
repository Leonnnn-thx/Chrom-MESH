#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Chromatin Multi-Enhancer Loop Network Module
============================================

Constructs the multi-enhancer network targeting genes based on chromatin interactions.
Analyzes regulatory depth using BFS to classify genes into Deep and Shallow regulation.
"""

import pandas as pd
import numpy as np
import networkx as nx
from collections import defaultdict
from typing import Dict, List, Any
import logging
from pathlib import Path

logger = logging.getLogger("Chrom-MESH.Identifier")

class RegulatoryNetworkAnalyzer:
    """
    Analyze regulatory depth of genes based on 3D chromatin topology.
    """
    
    def __init__(
        self,
        loop_file: str,
        enhancer_file: str,
        output_dir: str = "."
    ):
        self.loop_file = Path(loop_file)
        self.enhancer_file = Path(enhancer_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.loops = None
        self.enhancers = None
        self.graph = None
        self.gene_layers = {}
        self.deep_genes = []
        self.shallow_genes = []
        self.col_map = {}  # Column name mapping
        
        logger.info("RegulatoryNetworkAnalyzer initialized")

    def _detect_separator(self, file_path: Path) -> str:
        with open(file_path, 'r') as f:
            first_line = f.readline()
            if ',' in first_line and first_line.count(',') > first_line.count('\t'):
                return ','
            return '\t'

    def _split_items(self, item_string: Any) -> List[str]:
        if pd.isna(item_string):
            return []
        s = str(item_string).strip()
        if not s or s.lower() == 'nan':
            return []
        return [x.strip() for x in s.split(',') if x.strip() and x.strip().lower() != 'nan']

    def _normalize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df.columns = [str(c).replace("\ufeff", "").strip() for c in df.columns]
        return df

    def _create_column_mapping(self):
        """Map standard internal names to actual CSV column names"""
        possible_mappings = {
            'intersection_type': [
            'intersection_type', 'Intersection_Type', 'IntersectionType',
            'interaction_type', 'Interaction_Type', 'InteractionType'
            ],
            'p_l_genes': ['P_L_genes', 'p_l_genes', 'PL_genes'],
            'p_r_genes': ['P_R_genes', 'p_r_genes', 'PR_genes'],
            'e_l_ids': ['E_L_IDs', 'E_L_ids', 'e_l_ids', 'EL_IDs'],
            'e_r_ids': ['E_R_IDs', 'E_R_ids', 'e_r_ids', 'ER_IDs'],
            'distance': ['Distance', 'distance', 'dist']
        }
        
        self.col_map = {}
        for key, candidates in possible_mappings.items():
            found = False
            for cand in candidates:
                if cand in self.loops.columns:
                    self.col_map[key] = cand
                    found = True
                    break
            if not found:
                self.col_map[key] = None
        
        # Validate critical columns
        required = ['intersection_type', 'p_l_genes', 'p_r_genes', 'e_l_ids', 'e_r_ids']
        missing = [k for k in required if self.col_map[k] is None]
        if missing:
            raise ValueError(f"Missing required columns in loop file: {missing}")

    def load_data(self):
        logger.info("Loading data...")
        
        sep = self._detect_separator(self.loop_file)
        self.loops = pd.read_csv(self.loop_file, sep=sep)
        self.loops = self._normalize_columns(self.loops)
        logger.info(f"Loaded {len(self.loops)} loops from {self.loop_file.name}")

        self._create_column_mapping()
        
        # 2. Load Enhancers
        e_sep = self._detect_separator(self.enhancer_file)
        try:
            self.enhancers = pd.read_csv(self.enhancer_file, sep=e_sep)
            self.enhancers = self._normalize_columns(self.enhancers)

            if 'enhancer_id' not in self.enhancers.columns and len(self.enhancers.columns) >= 4:
                self.enhancers = pd.read_csv(
                    self.enhancer_file, sep=e_sep, header=None, 
                    names=['chr', 'start', 'end', 'enhancer_id']
                )
        except Exception as e:
            logger.error(f"Error loading enhancers: {e}")
            raise

        logger.info(f"Loaded {len(self.enhancers)} enhancers")

    def build_network(self) -> nx.Graph:
        if self.loops is None:
            self.load_data()

        logger.info("Building regulatory network...")
        self.graph = nx.Graph()
        
        stats_counter = defaultdict(int)
        
        # --- Process Tier 1: Enhancer-Promoter ---
        tier1_mask = self.loops[self.col_map['intersection_type']] == 'Enhancer-Promoter'
        tier1_loops = self.loops[tier1_mask]
        
        for _, row in tier1_loops.iterrows():
            genes_L = self._split_items(row[self.col_map['p_l_genes']])
            enhs_R = self._split_items(row[self.col_map['e_r_ids']])
            genes_R = self._split_items(row[self.col_map['p_r_genes']])
            enhs_L = self._split_items(row[self.col_map['e_l_ids']])
            dist = row[self.col_map['distance']] if self.col_map['distance'] else 0
            
            for g in genes_L:
                pid = f"P:{g}"
                self.graph.add_node(pid, node_type='promoter', genes=g)
                for e in enhs_R:
                    eid = f"E:{e}"
                    self.graph.add_node(eid, node_type='enhancer', enhancer_id=e)
                    self.graph.add_edge(pid, eid, interaction_type='tier1', distance=dist)
                    stats_counter['tier1_edges'] += 1

            for g in genes_R:
                pid = f"P:{g}"
                self.graph.add_node(pid, node_type='promoter', genes=g)
                for e in enhs_L:
                    eid = f"E:{e}"
                    self.graph.add_node(eid, node_type='enhancer', enhancer_id=e)
                    self.graph.add_edge(pid, eid, interaction_type='tier1', distance=dist)
                    stats_counter['tier1_edges'] += 1

        # --- Process Tier 2: Enhancer-Enhancer ---
        tier2_mask = self.loops[self.col_map['intersection_type']] == 'Enhancer-Enhancer'
        tier2_loops = self.loops[tier2_mask]
        
        for _, row in tier2_loops.iterrows():
            enhs_L = self._split_items(row[self.col_map['e_l_ids']])
            enhs_R = self._split_items(row[self.col_map['e_r_ids']])
            dist = row[self.col_map['distance']] if self.col_map['distance'] else 0
            
            for e1 in enhs_L:
                eid1 = f"E:{e1}"
                if not self.graph.has_node(eid1):
                    self.graph.add_node(eid1, node_type='enhancer', enhancer_id=e1)
                
                for e2 in enhs_R:
                    eid2 = f"E:{e2}"
                    if not self.graph.has_node(eid2):
                        self.graph.add_node(eid2, node_type='enhancer', enhancer_id=e2)
                    
                    self.graph.add_edge(eid1, eid2, interaction_type='tier2', distance=dist)
                    stats_counter['tier2_edges'] += 1

        logger.info(f"Network built: {self.graph.number_of_nodes()} nodes, "
                    f"{self.graph.number_of_edges()} unique edges, "
                    f"Tier1 Edges: {stats_counter['tier1_edges']}, "
                    f"Tier2 Edges: {stats_counter['tier2_edges']}")
        return self.graph

    def analyze(self) -> pd.DataFrame:
        if self.graph is None:
            self.build_network()
            
        logger.info("Assigning regulatory layers via BFS...")
        
        promoter_nodes = [n for n, d in self.graph.nodes(data=True) if d.get('node_type') == 'promoter']
        
        results = []
        self.deep_genes = []
        self.shallow_genes = []
        self.gene_layers = {}

        for promoter in promoter_nodes:
            gene_name = self.graph.nodes[promoter]['genes']
            
            tier1_neighbors = {n for n in self.graph.neighbors(promoter) 
                               if self.graph.nodes[n].get('node_type') == 'enhancer'}
            
            tier2_neighbors = set()
            for t1_node in tier1_neighbors:
                next_neighbors = self.graph.neighbors(t1_node)
                for t2_node in next_neighbors:
                    if (t2_node != promoter and 
                        t2_node not in tier1_neighbors and 
                        self.graph.nodes[t2_node].get('node_type') == 'enhancer'):
                        tier2_neighbors.add(t2_node)
            
            t1_ids = sorted([self.graph.nodes[n].get('enhancer_id', n) for n in tier1_neighbors])
            t2_ids = sorted([self.graph.nodes[n].get('enhancer_id', n) for n in tier2_neighbors])
            
            t1_count = len(t1_ids)
            t2_count = len(t2_ids)
            
            reg_type = "Unclassified"
            if t1_count > 0:
                if t2_count > 0:
                    reg_type = "Deep"
                    self.deep_genes.append(gene_name)
                else:
                    reg_type = "Shallow"
                    self.shallow_genes.append(gene_name)
            
            self.gene_layers[gene_name] = {
                'promoter_node': promoter,
                'tier1': list(tier1_neighbors),
                'tier2': list(tier2_neighbors),
                'tier1_count': t1_count,
                'tier2_count': t2_count
            }
            
            results.append({
                "Gene": gene_name,
                "Regulation_Type": reg_type,
                "Tier1_Count": t1_count,
                "Tier2_Count": t2_count,
                "Total_Enhancers": t1_count + t2_count,
                "Tier1_Enhancers": ";".join(t1_ids),
                "Tier2_Enhancers": ";".join(t2_ids)
            })
            
        logger.info(f"Classified {len(self.deep_genes)} Deep and {len(self.shallow_genes)} Shallow genes.")
        
        return pd.DataFrame(results)

    # ==========================================
    # NEW METHODS ADDED HERE
    # ==========================================

    def create_summary_statistics(self, filename='regulatory_depth_summary.txt'):
        """
        Create comprehensive summary statistics
        """
        output_file = self.output_dir / filename
        logger.info("Creating summary statistics...")
        
        with open(output_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("REGULATORY DEPTH ANALYSIS - SUMMARY STATISTICS\n")
            f.write("(interaction_type based analysis)\n")
            f.write("=" * 70 + "\n\n")
            
            # Network statistics
            f.write("NETWORK STATISTICS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Total nodes: {self.graph.number_of_nodes()}\n")
            f.write(f"Total edges: {self.graph.number_of_edges()}\n")
            
            promoter_count = sum(
                1 for n, attr in self.graph.nodes(data=True)
                if attr.get('node_type') == 'promoter'
            )
            enhancer_count = sum(
                1 for n, attr in self.graph.nodes(data=True)
                if attr.get('node_type') == 'enhancer'
            )
            
            f.write(f"Promoter nodes: {promoter_count}\n")
            f.write(f"Enhancer nodes: {enhancer_count}\n\n")
            
            # Interaction type statistics
            f.write("INTERACTION TYPE STATISTICS\n")
            f.write("-" * 70 + "\n")
            # Fixed: Use col_map instead of col
            interaction_col = self.col_map['intersection_type']
            interaction_counts = self.loops[interaction_col].value_counts()
            for itype, count in interaction_counts.items():
                f.write(f"{itype}: {count}\n")
            f.write("\n")
            
            # Gene classification
            f.write("GENE CLASSIFICATION\n")
            f.write("-" * 70 + "\n")
            total_genes = len(self.gene_layers)
            f.write(f"Total genes analyzed: {total_genes}\n")
            if total_genes > 0:
                f.write(f"Deep Regulation genes: {len(self.deep_genes)} "
                       f"({len(self.deep_genes)/total_genes*100:.1f}%)\n")
                f.write(f"Shallow Regulation genes: {len(self.shallow_genes)} "
                       f"({len(self.shallow_genes)/total_genes*100:.1f}%)\n\n")
            
            # Regulatory layer statistics
            f.write("REGULATORY LAYER STATISTICS\n")
            f.write("-" * 70 + "\n")
            f.write("(All enhancers are sorted and deduplicated per gene)\n\n")
            
            if total_genes > 0:
                tier1_counts = [l['tier1_count'] for l in self.gene_layers.values()]
                tier2_counts = [l['tier2_count'] for l in self.gene_layers.values()]
                
                f.write("Tier 1 enhancers per gene:\n")
                f.write(f"  Mean: {np.mean(tier1_counts):.2f}\n")
                f.write(f"  Median: {np.median(tier1_counts):.0f}\n")
                f.write(f"  Std Dev: {np.std(tier1_counts):.2f}\n")
                f.write(f"  Range: {np.min(tier1_counts):.0f} - {np.max(tier1_counts):.0f}\n\n")
                
                f.write("Tier 2 enhancers per gene:\n")
                f.write(f"  Mean: {np.mean(tier2_counts):.2f}\n")
                f.write(f"  Median: {np.median(tier2_counts):.0f}\n")
                f.write(f"  Std Dev: {np.std(tier2_counts):.2f}\n")
                f.write(f"  Range: {np.min(tier2_counts):.0f} - {np.max(tier2_counts):.0f}\n\n")
            
            # Top genes by regulatory complexity
            f.write("TOP 10 GENES BY REGULATORY COMPLEXITY\n")
            f.write("-" * 70 + "\n")
            
            sorted_genes = sorted(
                self.gene_layers.items(),
                key=lambda x: x[1]['tier1_count'] + x[1]['tier2_count'],
                reverse=True
            )
            
            for i, (gene, layers) in enumerate(sorted_genes[:10], 1):
                reg_type = 'Deep' if gene in self.deep_genes else 'Shallow'
                f.write(f"{i:2d}. {gene:20s} - {reg_type:8s} - "
                       f"Tier1: {layers['tier1_count']:3d}, "
                       f"Tier2: {layers['tier2_count']:3d}, "
                       f"Total: {layers['tier1_count'] + layers['tier2_count']:3d}\n")
            
            f.write("\n" + "=" * 70 + "\n")
        
        logger.info(f"Saved summary statistics to: {output_file}")
        return self

    def create_gene_lists(self, 
                         deep_file='deep_regulation_genes.txt',
                         shallow_file='shallow_regulation_genes.txt'):
        """
        Create gene lists for Deep and Shallow regulation
        """
        logger.info("Creating gene lists...")
        
        deep_path = self.output_dir / deep_file
        shallow_path = self.output_dir / shallow_file
        
        # Save Deep genes
        with open(deep_path, 'w') as f:
            for gene in sorted(self.deep_genes):
                f.write(f"{gene}\n")
        logger.info(f"Saved {len(self.deep_genes)} Deep genes to: {deep_path}")
        
        # Save Shallow genes
        with open(shallow_path, 'w') as f:
            for gene in sorted(self.shallow_genes):
                f.write(f"{gene}\n")
        logger.info(f"Saved {len(self.shallow_genes)} Shallow genes to: {shallow_path}")
        
        return self

    def visualize(self):
        """Generate statistical visualizations"""
        try:
            import matplotlib.pyplot as plt
        except ImportError as e:
            raise ImportError("Visualization requires matplotlib.") from e

        if not self.gene_layers:
            logger.warning("No analysis results to visualize. Run analyze() first.")
            return

        logger.info("Generating visualizations...")
        
        tier1_counts = [d['tier1_count'] for d in self.gene_layers.values()]
        tier2_counts = [d['tier2_count'] for d in self.gene_layers.values()]
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Regulatory Depth Analysis', fontsize=16, fontweight='bold')
        
        # 1. Pie Chart
        ax = axes[0, 0]
        sizes = [len(self.deep_genes), len(self.shallow_genes)]
        labels = [f'Deep\n({sizes[0]})', f'Shallow\n({sizes[1]})']
        ax.pie(sizes, labels=labels, colors=['#FF6B6B', '#4ECDC4'], autopct='%1.1f%%', startangle=90)
        ax.set_title('Gene Classification')
        
        # 2. Tier 1 Dist
        ax = axes[0, 1]
        ax.hist(tier1_counts, bins=30, color='#FFD93D', edgecolor='black', alpha=0.7)
        ax.set_title('Tier 1 Enhancer Distribution')
        ax.set_xlabel('Count')
        
        # 3. Tier 2 Dist
        ax = axes[1, 0]
        ax.hist(tier2_counts, bins=30, color='#6BCB77', edgecolor='black', alpha=0.7)
        ax.set_title('Tier 2 Enhancer Distribution')
        ax.set_xlabel('Count')
        
        # 4. Scatter
        ax = axes[1, 1]
        colors = ['#FF6B6B' if g in self.deep_genes else '#4ECDC4' for g in self.gene_layers.keys()]
        ax.scatter(tier1_counts, tier2_counts, c=colors, alpha=0.6, edgecolors='k')
        ax.set_xlabel('Tier 1')
        ax.set_ylabel('Tier 2')
        ax.set_title('Tier 1 vs Tier 2 Complexity')
        
        plt.tight_layout()
        out_path = self.output_dir / "regulatory_statistics.png"
        plt.savefig(out_path, dpi=300)
        plt.close()
        logger.info(f"Saved visualization to {out_path}")

    def save_results(self, filename: str = "gene_regulatory_layers.csv"):
        """Save main analysis results to CSV and trigger sub-reports"""
        if not self.gene_layers:
            raise ValueError("No results to save. Run analyze() first.")
        
        # 1. Save Main CSV
        df = pd.DataFrame([
            {
                'Gene': k,
                'Regulation_Type': 'Deep' if k in self.deep_genes else ('Shallow' if k in self.shallow_genes else 'Other'),
                'Tier1_Count': v['tier1_count'],
                'Tier2_Count': v['tier2_count'],
                'Total_Enhancers': v['tier1_count'] + v['tier2_count'],
                'Tier1_Enhancers': ";".join(sorted([self.graph.nodes[n].get('enhancer_id', n) for n in v['tier1']])),
                'Tier2_Enhancers': ";".join(sorted([self.graph.nodes[n].get('enhancer_id', n) for n in v['tier2']]))
            }
            for k, v in self.gene_layers.items()
        ])
        
        df = df.sort_values(['Regulation_Type', 'Total_Enhancers'], ascending=[True, False])
        
        out_csv = self.output_dir / filename
        df.to_csv(out_csv, index=False)
        logger.info(f"Saved detailed table to {out_csv}")
        
        # 2. Trigger the new methods automatically
        self.create_gene_lists()
        self.create_summary_statistics()

    def get_statistics(self) -> Dict[str, Any]:
        if not self.gene_layers:
            return {}
        
        t1 = [v['tier1_count'] for v in self.gene_layers.values()]
        t2 = [v['tier2_count'] for v in self.gene_layers.values()]
        
        return {
            "total_genes": len(self.gene_layers),
            "deep_genes": len(self.deep_genes),
            "shallow_genes": len(self.shallow_genes),
            "avg_tier1": np.mean(t1),
            "avg_tier2": np.mean(t2)
        }


def run_identify(args):
    logger.info("🚀 Starting Regulatory Depth Analysis...")
    
    try:
        analyzer = RegulatoryNetworkAnalyzer(
            loop_file=args.loop_file,
            enhancer_file=args.enhancer_file,
            output_dir=args.output_dir
        )
        
        analyzer.load_data()
        analyzer.build_network()
        results_df = analyzer.analyze()
        
        # This will now save the CSV, the txt lists, and the summary stats
        analyzer.save_results() 
        
        if hasattr(args, 'viz') and args.viz:
            analyzer.visualize()
            
        stats = analyzer.get_statistics()
        logger.info("=== Analysis Summary ===")
        for k, v in stats.items():
            val = f"{v:.2f}" if isinstance(v, float) else v
            logger.info(f" - {k.replace('_', ' ').title()}: {val}")
            
        logger.info(f"✅ Analysis complete. Results saved to: {args.output_dir}")
        return results_df

    except Exception as e:
        logger.error(f"❌ Analysis failed: {str(e)}")
        raise  

def run_analyze(
    loop_file=None,
    enhancer_file=None,
    output_dir=".",
    save=True
) -> pd.DataFrame:
    if hasattr(loop_file, "__dict__") and enhancer_file is None:
        args = loop_file
        loop_file = getattr(args, "loops", None) or getattr(args, "loop_file", None) or getattr(args, "l", None)
        enhancer_file = getattr(args, "enhancers", None) or getattr(args, "enhancer_file", None) or getattr(args, "e", None)
        output_dir = getattr(args, "output", None) or getattr(args, "output_dir", None) or getattr(args, "o", None) or output_dir
        min_tier1 = getattr(args, "min_tier1", 1)
        min_tier2 = getattr(args, "min_tier2", 1)
        save = getattr(args, "save", save)
    else:
        min_tier1, min_tier2 = 1, 1

    if loop_file is None or enhancer_file is None:
        raise ValueError("run_analyze requires loop_file and enhancer_file.")

    output_dir = Path(output_dir)
    if output_dir.suffix in {".csv", ".tsv", ".txt"}:
        output_dir = output_dir.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    analyzer = RegulatoryNetworkAnalyzer(
        loop_file=str(loop_file),
        enhancer_file=str(enhancer_file),
        output_dir=str(output_dir)
    )
    analyzer.load_data()
    analyzer.build_network()
    df = analyzer.analyze()
    
    if "Tier1_Count" in df.columns and "Tier2_Count" in df.columns:
        df = df[(df["Tier1_Count"] >= min_tier1) | (df["Tier2_Count"] >= min_tier2)]

    if save:
        # Save main CSV, lists, and summary
        analyzer.save_results()
        print(f"[Chrom-MESH] Saved comprehensive results to: {output_dir}")

    return df

def main():
    import argparse
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(description="Chrom-MESH Identifier")
    parser.add_argument("-l", "--loop_file", required=True)
    parser.add_argument("-e", "--enhancer_file", required=True)
    parser.add_argument("-o", "--output_dir", default=".")
    parser.add_argument("--viz", action="store_true")

    args = parser.parse_args()
    run_identify(args)

if __name__ == "__main__":
    main()
