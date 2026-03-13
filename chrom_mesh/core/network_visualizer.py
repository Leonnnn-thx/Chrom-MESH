#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from typing import List, Optional

logger = logging.getLogger("Chrom-MESH.Visualize")

class NetworkVisualizer:
    """
    Extracts genes and their associated Tier 1/2 enhancers from analysis results.
    """
    
    def __init__(self, gene_layers_file: str, loop_file: str):
        self.gene_layers_file = gene_layers_file
        self.loop_file = loop_file
        self.gene_layers_df = None
        self.loops_df = None
        self.graph = nx.Graph()
        self.col_map = {}

    def load_data(self):
        logger.info("Loading data for visualization...")
        self.gene_layers_df = pd.read_csv(self.gene_layers_file)
        sep = ',' if str(self.loop_file).lower().endswith('.csv') else '\t'
        self.loops_df = pd.read_csv(self.loop_file, sep=sep)
        required_cols = {
            'itype': ['intersection_type', 'Intersection_Type', 'interaction_type', 'Interaction_Type'],
            'p_l': ['P_L_genes', 'p_l_genes'],
            'p_r': ['P_R_genes', 'p_r_genes'],
            'e_l': ['E_L_IDs', 'e_l_ids'],
            'e_r': ['E_R_IDs', 'e_r_ids']
        }
        for key, candidates in required_cols.items():
            for c in candidates:
                if c in self.loops_df.columns:
                    self.col_map[key] = c
                    break
        return self

    def _split_ids(self, val):
        if pd.isna(val) or str(val).lower() == 'nan': return []
        return [i.strip() for i in str(val).split(',')]

    def build_network(self, target_genes: List[str]):
        logger.info(f"Building network for {len(target_genes)} genes...")
        genes_set = set(target_genes)
        tier1_enhancers = set()

        # 1. Tier 1: Promoter-Enhancer
        t1_loops = self.loops_df[self.loops_df[self.col_map['itype']] == 'Enhancer-Promoter']
        for _, row in t1_loops.iterrows():
            genes = self._split_ids(row[self.col_map['p_l']]) + self._split_ids(row[self.col_map['p_r']])
            enhs = self._split_ids(row[self.col_map['e_l']]) + self._split_ids(row[self.col_map['e_r']])
            
            for g in genes:
                if g in genes_set:
                    p_id = f"P:{g}"
                    self.graph.add_node(p_id, node_type='promoter', label=g)
                    for e in enhs:
                        e_id = f"E:{e}"
                        self.graph.add_node(e_id, node_type='enhancer', tier='tier1')
                        self.graph.add_edge(p_id, e_id, itype='tier1')
                        tier1_enhancers.add(e_id)

        # 2. Tier 2: Enhancer-Enhancer 
        t2_loops = self.loops_df[self.loops_df[self.col_map['itype']] == 'Enhancer-Enhancer']
        for _, row in t2_loops.iterrows():
            e_ls = [f"E:{e}" for e in self._split_ids(row[self.col_map['e_l']])]
            e_rs = [f"E:{e}" for e in self._split_ids(row[self.col_map['e_r']])]
            
            if any(e in tier1_enhancers for e in e_ls + e_rs):
                for el in e_ls:
                    for er in e_rs:
                        if not self.graph.has_node(el): self.graph.add_node(el, node_type='enhancer', tier='tier2')
                        if not self.graph.has_node(er): self.graph.add_node(er, node_type='enhancer', tier='tier2')
                        self.graph.add_edge(el, er, itype='tier2')

    def draw(self, output_path: str, figsize=(12, 10)):
        if self.graph.number_of_nodes() == 0:
            logger.warning("No nodes found for the specified genes.")
            return

        plt.figure(figsize=figsize)
        pos = nx.spring_layout(self.graph, k=0.3, seed=42)
        
        nodes = self.graph.nodes(data=True)
        p_nodes = [n for n, d in nodes if d.get('node_type') == 'promoter']
        e_nodes = [n for n, d in nodes if d.get('node_type') == 'enhancer']
        
        nx.draw_networkx_nodes(self.graph, pos, nodelist=p_nodes, node_color='#FF7F50', node_size=700, label='Promoter')
        nx.draw_networkx_nodes(self.graph, pos, nodelist=e_nodes, node_color='#87CEEB', node_size=400, label='Enhancer')
        
        edges = self.graph.edges(data=True)
        t1_edges = [(u, v) for u, v, d in edges if d.get('itype') == 'tier1']
        t2_edges = [(u, v) for u, v, d in edges if d.get('itype') == 'tier2']
        
        nx.draw_networkx_edges(self.graph, pos, edgelist=t1_edges, width=1.5, alpha=0.6, edge_color='gray')
        nx.draw_networkx_edges(self.graph, pos, edgelist=t2_edges, width=1.0, alpha=0.4, edge_color='blue', style='dashed')

        labels = {n: d.get('label', '') for n, d in nodes if d.get('node_type') == 'promoter'}
        nx.draw_networkx_labels(self.graph, pos, labels, font_size=10, font_weight='bold')
        
        plt.title("Chrom-MESH: Regulatory Network Visualization", fontsize=15)
        plt.legend(scatterpoints=1)
        plt.axis('off')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Visualization saved to {output_path}")

    def draw(self, output_file: str, dpi: int = 300, **kwargs):
        """
        Draw and save network figure.
        """
        if self.graph.number_of_nodes() == 0:
            raise ValueError("Graph is empty. Please run build_network(...) first.")

        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(self.graph, seed=42)

        promoter_nodes = [n for n, d in self.graph.nodes(data=True) if d.get("node_type") == "promoter"]
        enhancer_nodes = [n for n, d in self.graph.nodes(data=True) if d.get("node_type") == "enhancer"]

        nx.draw_networkx_nodes(self.graph, pos, nodelist=promoter_nodes, node_color="#4C78A8", node_size=700, label="Promoter")
        nx.draw_networkx_nodes(self.graph, pos, nodelist=enhancer_nodes, node_color="#F58518", node_size=500, label="Enhancer")
        nx.draw_networkx_edges(self.graph, pos, alpha=0.6)

        labels = {n: d.get("label", n) for n, d in self.graph.nodes(data=True)}
        nx.draw_networkx_labels(self.graph, pos, labels=labels, font_size=8)

        plt.legend()
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
        plt.close()

        logger.info(f"Saved network figure to: {output_file}")
        print(f"[Chrom-MESH] Saved figure: {output_file}")



def run_visualize(args):
    target_genes = []

    if args.genes:
        target_genes = [g.strip() for g in args.genes.split(',') if g.strip()]
    elif args.gene_file:
        if not os.path.exists(args.gene_file):
            raise FileNotFoundError(f"Gene file not found: {args.gene_file}")
        with open(args.gene_file, 'r') as f:
            target_genes = [line.strip() for line in f if line.strip()]

    viz = NetworkVisualizer(args.gene_layers, args.loops)
    viz.load_data()

    if not target_genes:
        logger.info(f"No specific genes provided. Using top {args.top_n} genes from analysis.")
        if "Gene" in viz.gene_layers_df.columns:
            gene_col = "Gene"
        elif "gene" in viz.gene_layers_df.columns:
            gene_col = "gene"
        else:
            raise ValueError("No Gene/gene column found in gene_layers file.")
        target_genes = viz.gene_layers_df[gene_col].head(args.top_n).tolist()

    viz.build_network(target_genes)

    out_file = args.output if args.output else "regulatory_network.pdf"
    viz.draw(out_file, dpi=args.dpi)

