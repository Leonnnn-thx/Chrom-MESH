#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Chrom-MESH Command Line Interface
=================================

Main entry point for all Chrom-MESH commands.
"""

import argparse
import sys
import logging
from pathlib import Path

from chrom_mesh.__version__ import __version__
from chrom_mesh.utils.logger import setup_logger


def main():
    """Main CLI entry point"""
    
    parser = argparse.ArgumentParser(
        prog="chrom-mesh",
        description="Chrom-MESH: Chromatin Multi-Enhancer Spatial Hierarchy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Annotate loops
  chrom-mesh annotate -l loops.bedpe -e enhancers.bed -p promoters.bed -o output.csv
  
  # Analyze regulatory depth
  chrom-mesh analyze -l annotated_loops.csv -e enhancers.bed -o gene_layers.csv
  
  # Visualize network
  chrom-mesh visualize -g gene_layers.csv -l annotated_loops.csv -G MYC,BCL6

For more information, visit: https://github.com/159357thx/Chrom-MESH
        """
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version=f"Chrom-MESH {__version__}"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    subparsers = parser.add_subparsers(
        dest="command",
        help="Available commands"
    )
    
    # Annotate command
    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Annotate chromatin loops"
    )
    annotate_parser.add_argument(
        "-l", "--loops",
        required=True,
        help="Input loop file (BEDPE format)"
    )
    annotate_parser.add_argument(
        "-e", "--enhancers",
        required=True,
        help="Enhancer BED file"
    )
    annotate_parser.add_argument(
        "-p", "--promoters",
        required=True,
        help="Promoter BED file"
    )
    annotate_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file"
    )

    
    # Analyze command
    analyze_parser = subparsers.add_parser(
        "analyze",
        help="Analyze regulatory depth"
    )
    analyze_parser.add_argument(
        "-l", "--loops",
        required=True,
        help="Annotated loop file (CSV format)"
    )
    analyze_parser.add_argument(
        "-e", "--enhancers",
        required=True,
        help="Enhancer BED file"
    )
    analyze_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file"
    )
    analyze_parser.add_argument(
        "--min-tier1",
        type=int,
        default=1,
        help="Minimum Tier 1 enhancers (default: 1)"
    )
    analyze_parser.add_argument(
        "--min-tier2",
        type=int,
        default=1,
        help="Minimum Tier 2 enhancers for deep regulation (default: 1)"
    )
    
    # Visualize command
    visualize_parser = subparsers.add_parser(
        "visualize",
        help="Visualize regulatory network"
    )
    visualize_parser.add_argument(
        "-g", "--gene-layers",
        required=True,
        help="Gene layers file (CSV format)"
    )
    visualize_parser.add_argument(
        "-l", "--loops",
        required=True,
        help="Annotated loop file (CSV format)"
    )
    visualize_parser.add_argument(
        "-G", "--genes",
        help="Comma-separated list of genes to visualize"
    )
    visualize_parser.add_argument(
        "-f", "--gene-file",
        help="File containing genes (one per line)"
    )
    visualize_parser.add_argument(
        "-o", "--output",
        default="regulatory_network.pdf",
        help="Output file (default: regulatory_network.pdf)"
    )
    visualize_parser.add_argument(
        "-n", "--top-n",
        type=int,
        default=10,
        help="Number of top genes to visualize (default: 10)"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    setup_logger(level=log_level)
    
    # Execute command
    if args.command == "annotate":
        from chrom_mesh.cli.annotate_loops import run_annotate
        run_annotate(args)
    elif args.command == "analyze":
        from chrom_mesh.cli.analysis_network import run_analyze
        run_analyze(args)
    elif args.command == "visualize":
        from chrom_mesh.cli.visualize_network import run_visualize
        run_visualize(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
