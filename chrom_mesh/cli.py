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

from chrom_mesh.__version__ import __version__
from chrom_mesh.utils.logger import setup_logger


def main():
    parser = argparse.ArgumentParser(
        prog="chrom-mesh",
        description="Chrom-MESH: Chromatin Multi-Enhancer Spatial Hierarchy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Identify enhancers from H3K27Ac + ATAC inputs
  chrom-mesh identify \\
    --h3k27ac-bed-dir ./h3k27ac_peaks \\
    --h3k27ac-bw-dir ./h3k27ac_bw \\
    --atac-bed-dir ./atac_bed \\
    --atac-bw-dir ./atac_bw \\
    -o ./results

  # Annotate loops
  chrom-mesh annotate -l loops.bedpe -e enhancers.bed -p promoters.bed -o output.csv

  # Analyze regulatory depth
  chrom-mesh analyze -l annotated_loops.csv -e enhancers.bed -o gene_layers.csv

  # Visualize network
  chrom-mesh visualize -g gene_layers.csv -l annotated_loops.csv -G MYC,BCL6
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
        help="Available commands",
        required=True
    )

    # -----------------------------
    # Identify enhancers (shell backend)
    # -----------------------------
    identify_parser = subparsers.add_parser(
        "identify",
        help="Identify active enhancers from H3K27Ac and ATAC-seq data"
    )

    # keep legacy bam arg for backward compatibility 
    identify_parser.add_argument(
        "-b", "--bam",
        help="(Legacy) Input BAM file; not used by Identify_AEs.sh pipeline"
    )

    identify_parser.add_argument(
        "--h3k27ac-bed-dir",
        required=True,
        help="H3K27Ac narrowPeak/BED directory (maps to -h)"
    )
    identify_parser.add_argument(
        "--h3k27ac-bw-dir",
        required=True,
        help="H3K27Ac BigWig directory (maps to -H)"
    )
    identify_parser.add_argument(
        "--atac-bed-dir",
        required=True,
        help="ATAC BED directory (maps to -a)"
    )
    identify_parser.add_argument(
        "--atac-bw-dir",
        required=True,
        help="ATAC BigWig directory (maps to -A)"
    )
    identify_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory (maps to -o)"
    )

    # tool/path overrides for script
    identify_parser.add_argument(
        "--bigwig-tool",
        default=None,
        help="Path or command name for bigWigAverageOverBed (maps to -b)"
    )
    identify_parser.add_argument(
        "--bigwig-info-tool",
        default=None,
        help="Path or command name for bigWigInfo (maps to --bigwig-info-tool)"
    )
    identify_parser.add_argument(
        "--bedtools",
        default=None,
        help="Path or command name for bedtools"
    )
    identify_parser.add_argument(
        "--python",
        dest="python_exec",
        default=None,
        help="Python executable used by shell pipeline"
    )

    # patterns
    identify_parser.add_argument(
        "--h3k27ac-pattern",
        default="*_peaks.narrowPeak",
        help="H3K27Ac file matching pattern (maps to -p)"
    )
    identify_parser.add_argument(
        "--atac-pattern",
        default="*.bed",
        help="ATAC BED file matching pattern (maps to -n)"
    )
    identify_parser.add_argument(
        "--h3k27ac-bw-pattern",
        default="*.bw",
        help="H3K27Ac BigWig matching pattern (maps to -w)"
    )
    identify_parser.add_argument(
        "--atac-bw-pattern",
        default="*.bw",
        help="ATAC BigWig matching pattern (maps to -W)"
    )

    # legacy options kept for compatibility (currently unused by shell script)
    identify_parser.add_argument("--workdir", help="(Legacy) Working directory")
    identify_parser.add_argument("--min-width", type=int, default=100, help="(Legacy) Minimum peak width")
    identify_parser.add_argument("--max-width", type=int, default=1500, help="(Legacy) Maximum peak width")

    # -----------------------------
    # Annotate chromatin loops
    # -----------------------------
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

    # -----------------------------
    # Analyze chromatin loop networks
    # -----------------------------
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
        help="Output file path"
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

    # -----------------------------
    # Visualize networks
    # -----------------------------
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
        help="Comma-separated gene list, e.g. MYC,BCL6"
    )
    visualize_parser.add_argument(
        "-f", "--gene-file",
        help="File containing genes (one per line)"
    )
    visualize_parser.add_argument(
        "-o", "--output",
        default="regulatory_network.pdf",
        help="Output directory"
    )
    visualize_parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure DPI (default: 300)"
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
    if args.command == "identify":
        from chrom_mesh.core import run_identify
        run_identify(args)
    elif args.command == "annotate":
        from chrom_mesh.core import run_annotate
        run_annotate(args)
    elif args.command == "analyze":
        from chrom_mesh.core import run_analyze
        run_analyze(args)
    elif args.command == "visualize":
        from chrom_mesh.core import run_visualize
        run_visualize(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
