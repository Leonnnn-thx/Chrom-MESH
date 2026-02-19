#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Chrom-MESH: Chromatin Multi-Enhancer Spatial Hierarchy
A comprehensive toolkit for analyzing chromatin multi-enhancer loop networks and gene regulation
"""

from setuptools import setup, find_packages
import os

# Read version
version = {}
with open("chrom_mesh/__version__.py") as fp:
    exec(fp.read(), version)

# Read README
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="chrom-mesh",
    version=version["__version__"],
    author="Tao Hongxu",
    author_email="kq_thx@163.com",
    description="A toolkit for analyzing 3D chromatin architecture and gene regulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/159thx/Chrom-MESH",
    packages=find_packages(exclude=["tests", "docs", "examples"]),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
            "mypy>=0.900",
            "sphinx>=4.0",
            "sphinx-rtd-theme>=1.0",
        ],
        "goatools": [
            "goatools>=1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "chrom-mesh=chrom_mesh.cli.main:main",
            "chrom-mesh-identify=chrom_mesh.cli.identify_enhancers:main",
            "chrom-mesh-annotate=chrom_mesh.cli.annotate_loops:main",
            "chrom-mesh-analyze=chrom_mesh.cli.analyze_depth:main",
            "chrom-mesh-visualize=chrom_mesh.cli.visualize_network:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords=[
        "bioinformatics",
        "genomics",
        "chromatin",
        "3D genome",
        "enhancer",
        "promoter",
        "gene regulation",
        "Hi-C",
        "ChIP-seq",
        "ATAC-seq",
    ],
)
