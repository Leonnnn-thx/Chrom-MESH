#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Loop Annotator Module
=====================
Annotates chromatin loops by classifying anchors as promoters/enhancers.
Based on specific logic for 7-column promoter files and 4-column enhancer files.
"""

import pandas as pd
import numpy as np
from typing import Optional, Tuple, Dict, List
import logging
from pathlib import Path
import os


logger = logging.getLogger(__name__)

class LoopAnnotator:
    """
    Annotate chromatin loops with enhancer-promoter interactions.
    
    This class identifies whether loop anchors overlap with promoter or
    enhancer regions and classifies loops into specific interaction types
    (pure_PP, pure_EE, pure_EP, mixed_EP).

    Parameters
    ----------
    loop_file : str or Path
        Path to BEDPE file containing chromatin loops.
    enhancer_file : str or Path
        Path to BED file containing enhancer regions (4 columns: chr, start, end, id).
    promoter_file : str or Path
        Path to BED file containing promoter regions (7 columns).
    output_file : str or Path
        Path to save the annotated results.
    
    Attributes
    ----------
    loops : pd.DataFrame
        Loaded loop data.
    enhancers_dict : Dict[str, pd.DataFrame]
        Enhancer data grouped by chromosome.
    promoters_dict : Dict[str, pd.DataFrame]
        Promoter data grouped by chromosome.
    """
    
    def __init__(
        self,
        loop_file: str,
        enhancer_file: str,
        promoter_file: str,
        output_file: str
    ):
        self.loop_file = Path(loop_file)
        self.enhancer_file = Path(enhancer_file)
        self.promoter_file = Path(promoter_file)
        self.output_file = Path(output_file)
        
        # Validate files
        self._validate_files()
        
        # Data containers
        self.loops = None
        self.enhancers_dict = {}
        self.promoters_dict = {}
        self.annotated_loops = None
        
        logger.info("LoopAnnotator initialized")
    
    def _validate_files(self):
        """Validate input files exist"""
        for file_path in [self.loop_file, self.enhancer_file, self.promoter_file]:
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {file_path}")
    
    def load_data(self):
        """Load all input files and preprocess them into dictionaries for fast lookup"""
        logger.info("Loading input files...")
        
        # 1. Load Loops (BEDPE)
        # Using specific columns as per the source logic
        try:
            self.loops = pd.read_csv(
                self.loop_file,
                sep='\t',
                header=None,
                comment='#',
                usecols=[0, 1, 2, 4, 5, 6],
                names=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'],
                dtype={'chr1': str, 'start1': int, 'end1': int, 'chr2': str, 'start2': int, 'end2': int}
            )
            logger.info(f"Loaded {len(self.loops)} loops")
        except Exception as e:
            logger.error(f"Failed to load loop file: {e}")
            raise

        # 2. Load Enhancers (4-column BED)
        # Grouping by chromosome immediately for performance
        try:
            enh_df = pd.read_csv(
                self.enhancer_file, 
                sep='\t', 
                header=None, 
                comment='#',
                usecols=[0, 1, 2, 3],
                names=['chr', 'start', 'end', 'id'],
                dtype={'chr': str, 'start': int, 'end': int, 'id': str}
            )
            self.enhancers_dict = {chrom: group.copy() for chrom, group in enh_df.groupby('chr')}
            logger.info(f"Loaded {len(enh_df)} enhancers across {len(self.enhancers_dict)} chromosomes")
        except Exception as e:
            logger.error(f"Failed to load enhancer file: {e}")
            raise

        # 3. Load Promoters (7-column BED)
        # Columns: chr, start, end, gene, gene_chr, gene_start, gene_end
        try:
            prom_df = pd.read_csv(
                self.promoter_file,
                sep='\t',
                header=None,
                comment='#',
                usecols=[0, 1, 2, 3, 4, 5, 6],
                names=['chr', 'start', 'end', 'gene', 'gene_chr', 'gene_start', 'gene_end'],
                dtype={
                    'chr': str, 'start': int, 'end': int, 'gene': str,
                    'gene_chr': str, 'gene_start': int, 'gene_end': int
                }
            )
            self.promoters_dict = {chrom: group.copy() for chrom, group in prom_df.groupby('chr')}
            logger.info(f"Loaded {len(prom_df)} promoters across {len(self.promoters_dict)} chromosomes")
        except Exception as e:
            logger.error(f"Failed to load promoter file: {e}")
            raise

    def _get_overlaps(self, chrom: str, start: int, end: int, ref_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """
        Helper to find overlaps for a specific anchor against a reference dictionary (promoters or enhancers).
        """
        if chrom not in ref_dict:
            return pd.DataFrame()
        
        ref_df = ref_dict[chrom]
        # Standard overlap logic: (StartA < EndB) and (EndA > StartB)
        hits = ref_df[
            (ref_df['start'] < end) & (ref_df['end'] > start)
        ]
        return hits

    def annotate(self) -> pd.DataFrame:
        """
        Annotate all loops and apply classification logic.
        
        Returns
        -------
        pd.DataFrame
            Annotated loops with classification columns.
        """
        if self.loops is None:
            self.load_data()
        
        logger.info("Annotating loops (this may take a while)...")
        
        # Prepare result containers
        # Using a list of dicts is often faster than appending to separate lists
        results_data = []
        total_loops = len(self.loops)
        
        for idx, loop in self.loops.iterrows():
            if (idx + 1) % 5000 == 0:
                logger.info(f"Processing loop {idx + 1} / {total_loops}...")

            # Define anchors
            anchors = {
                'L': {'chr': loop['chr1'], 'start': loop['start1'], 'end': loop['end1']},
                'R': {'chr': loop['chr2'], 'start': loop['start2'], 'end': loop['end2']}
            }
            
            row_result = {}
            
            # Process Left (L) and Right (R) anchors
            for side, coords in anchors.items():
                # 1. Promoter Overlaps
                prom_hits = self._get_overlaps(coords['chr'], coords['start'], coords['end'], self.promoters_dict)
                
                is_promoter = not prom_hits.empty
                row_result[f'P_{side}'] = is_promoter
                
                if is_promoter:
                    genes = prom_hits['gene'].dropna().unique()
                    row_result[f'P_{side}_genes'] = ','.join(genes)
                    # Create coordinate strings: gene_chr:gene_start-gene_end
                    coords_str = prom_hits.apply(
                        lambda r: f"{r['gene_chr']}:{r['gene_start']}-{r['gene_end']}", axis=1
                    )
                    row_result[f'P_{side}_coords'] = ','.join(coords_str)
                else:
                    row_result[f'P_{side}_genes'] = ''
                    row_result[f'P_{side}_coords'] = ''

                # 2. Enhancer Overlaps
                enh_hits = self._get_overlaps(coords['chr'], coords['start'], coords['end'], self.enhancers_dict)
                
                is_enhancer = not enh_hits.empty
                row_result[f'E_{side}'] = is_enhancer
                
                if is_enhancer:
                    ids = enh_hits['id'].dropna().unique()
                    row_result[f'E_{side}_IDs'] = ','.join(ids)
                    coords_str = enh_hits.apply(
                        lambda r: f"{r['chr']}:{r['start']}-{r['end']}", axis=1
                    )
                    row_result[f'E_{side}_coords'] = ','.join(coords_str)
                else:
                    row_result[f'E_{side}_IDs'] = ''
                    row_result[f'E_{side}_coords'] = ''
            
            results_data.append(row_result)
        
        # Merge results back to main DataFrame
        logger.info("Integrating annotation results...")
        results_df = pd.DataFrame(results_data)
        self.annotated_loops = pd.concat([self.loops.reset_index(drop=True), results_df], axis=1)
        
        # Apply Classification Logic
        self._classify_loops()
        
        return self.annotated_loops

    def _classify_loops(self):
        """
        Apply the specific logical rules to classify loops into pure_PP, pure_EE, etc.
        """
        logger.info("Calculating loop classifications...")
        df = self.annotated_loops
        
        # Boolean Logic Shortcuts
        PL, PR = df['P_L'], df['P_R']
        EL, ER = df['E_L'], df['E_R']
        
        # 1. Pure P-P: Both are P, neither are E
        df['pure_PP'] = PL & PR & (~EL) & (~ER)
        
        # 2. Pure E-E: Both are E, neither are P
        df['pure_EE'] = EL & ER & (~PL) & (~PR)
        
        # 3. Pure E-P: One is pure E, one is pure P
        df['pure_EP'] = ( (EL & ~PL) & (PR & ~ER) ) | \
                        ( (PL & ~EL) & (ER & ~PR) )
        
        # 4. Mixed E-P: Complex cases involving dual identity
        # Logic from source: (EL & PR & (PL or ER)) OR (ER & PL & (PR or EL))
        df['mixed_EP'] = (
            (EL & PR & (PL | ER)) |
            (ER & PL & (PR | EL))
        )
        
        self.annotated_loops = df

    def save_results(self):
        """Save annotated loops to the output file with specific column ordering"""
        if self.annotated_loops is None:
            raise ValueError("No annotated loops to save. Run annotate() first.")
        
        # Ensure output directory exists
        if self.output_file.parent:
            self.output_file.parent.mkdir(parents=True, exist_ok=True)
            
        # Define final column order
        final_column_order = [
            'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2',
            'P_L', 'P_R', 'E_L', 'E_R',
            'P_L_coords', 'P_R_coords',
            'E_L_coords', 'E_R_coords', 'P_L_genes', 'P_R_genes', 'E_L_IDs', 'E_R_IDs',
            'pure_PP', 'pure_EE', 'pure_EP', 'mixed_EP'
        ]
        
        # Ensure all columns exist (fill missing with empty string if necessary)
        for col in final_column_order:
            if col not in self.annotated_loops.columns:
                self.annotated_loops[col] = ''
                
        # Save
        self.annotated_loops[final_column_order].to_csv(
            self.output_file, 
            sep=',', 
            index=False, 
            na_rep=''
        )
        logger.info(f"Results saved to: {self.output_file}")

    def get_statistics(self) -> Dict[str, int]:
        """Return basic counts of classifications"""
        if self.annotated_loops is None:
            return {}
        
        stats = {
            "total_loops": len(self.annotated_loops),
            "pure_PP": self.annotated_loops['pure_PP'].sum(),
            "pure_EE": self.annotated_loops['pure_EE'].sum(),
            "pure_EP": self.annotated_loops['pure_EP'].sum(),
            "mixed_EP": self.annotated_loops['mixed_EP'].sum()
        }
        return stats


# --- Command Line Interface ---

def main():
    import argparse
    
    # Setup basic logging configuration for CLI usage
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    parser = argparse.ArgumentParser(
        description="Annotate chromatin loops (chrom-MESH core module).",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("-l", "--loop_file", required=True, help="Input BEDPE loop file")
    parser.add_argument("-e", "--enhancer_file", required=True, help="Input Enhancer BED (4-col)")
    parser.add_argument("-p", "--promoter_file", required=True, help="Input Promoter BED (7-col)")
    parser.add_argument("-o", "--output_file", required=True, help="Output CSV file path")

    args = parser.parse_args()

    # Instantiate and run
    try:
        annotator = LoopAnnotator(
            loop_file=args.loop_file,
            enhancer_file=args.enhancer_file,
            promoter_file=args.promoter_file,
            output_file=args.output_file
        )
        
        annotator.annotate()
        annotator.save_results()
        
        # Print stats
        stats = annotator.get_statistics()
        logger.info("=== Annotation Statistics ===")
        for k, v in stats.items():
            logger.info(f"{k}: {v}")
            
    except Exception as e:
        logger.error(f"Process failed: {e}")
        exit(1)

if __name__ == "__main__":
    main()
