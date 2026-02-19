#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Chrom-MESH: Chromatin Multi-Enhancer Spatial Hierarchy
=====================================================================

A comprehensive toolkit for analyzing chromatin multi-enhancer loop networks and gene regulation.

Main Components:
----------------
- EnhancerIdentifier: Identify active enhancers from H3K27Ac and ATAC-seq data
- LoopAnnotator: Annotate chromatin loops with enhancer, promoter data
- NetworkAnalyzer: Analyze regulatory depth and gene regulation patterns
- NetworkVisualizer: Visualize gene regulatory networks

Quick Start:
-----------
>>> from chrom_mesh import LoopAnnotator, NetworkAnalyzer, NetworkVisualizer
>>> 
>>> # Annotate loops
>>> annotator = LoopAnnotator(
...     loop_file="loops.bedpe",
...     enhancer_file="enhancers.bed",
...     promoter_file="promoters.bed"
... )
>>> annotated_loops = annotator.annotate()
>>> 
>>> # Analyze regulatory depth
>>> analyzer = NetworkAnalyzer(
...     loop_file="annotated_loops.csv",
...     enhancer_file="enhancers.bed"
... )
>>> results = analyzer.analyze()
>>> 
>>> # Visualize network
>>> visualizer = NetworkVisualizer(
...     gene_layers="gene_layers.csv",
...     loops="annotated_loops.csv"
... )
>>> visualizer.visualize(genes=["MYC", "BCL6"])

Modules:
--------
- core: Core analysis modules
- utils: Utility functions
- cli: Command-line interface
"""

from .__version__ import __version__, __author__, __email__

# Import main classes
from .core.loop_annotator import LoopAnnotator
from .core.loop_network_identifer import RegulatoryNetworkAnalyzer as NetworkAnalyzer
from .core.network_visualizer import RegulatoryNetworkVisualizer as NetworkVisualizer

# Import utility functions
from .utils.file_io import load_bed, load_bedpe, save_results
from .utils.logger import setup_logger

__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "LoopAnnotator",
    "NetworkAnalyzer",
    "NetworkVisualizer",
    "load_bed",
    "load_bedpe",
    "save_results",
    "setup_logger",
]

# Package metadata
__title__ = "Chrom-MESH"
__description__ = "A toolkit for analyzing 3D chromatin multi-enhancer loop networks and gene regulation"
__url__ = "https://github.com/159thx/Chrom-MESH"
__license__ = "MIT"
