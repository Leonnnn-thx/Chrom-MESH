# Chrom-MESH

**Chromatin Multi-scale Enhancer-promoter Spatial Hierarchy**

A comprehensive toolkit for analyzing chromatin multi-enhancer loop networks and gene regulation.

## Features

-   🧬 **Enhancer Identification**: Identify active enhancers from H3K27Ac and ATAC-seq data
-   🔗 **Loop Annotation**: Annotate chromatin loops with enhancer/promoter data
-   📊 **Multi-Enhancer loop Networks**: Construct chromatin multi-enhancer loop networks and Classify genes into deep/shallow regulation categories
-   🎨 **Network Visualization**: Visualize gene regulatory networks

## Workflow

![](images/Chrom-MESH.jpg)

## Python Environment

``` bash
# Create a conda/mamba environment
mamba create -n chrom-mesh python=3.10 -y
mamba activate chrom-mesh

# Or using venv
python -m venv .venv
source .venv/bin/activate   # macOS/Linux
# .venv\Scripts\activate    # Windows PowerShell
```

## Requirements

-   Python: \>=3.10,\<3.13
-   OS: Linux/macOS (Windows via WSL recommended)
-   bedtools \>= 2.30
-   bigWigAverageOverBed (UCSC)
-   bigWigInfo (UCSC)

## Install Chrom-MESH

``` bash
# Clone the repository
git clone https://github.com/Leonnnn-thx/Chrom-MESH.git
cd Chrom-MESH

# Install the package in development mode
pip install -e .

# Verify installation
python -c "import chrom_mesh; print(chrom_mesh.__version__)"
chrom-mesh --help
python -m chrom_mesh.cli --help
```

## Dependency

``` bash
# Quick check of dependencies
which bedtools
which bigWigAverageOverBed
which bigWigInfo

# If not installed, install via conda/mamba 
mamba install -c bioconda bedtools ucsc-bigwigaverageoverbed ucsc-bigwiginfo -y

# Optional dependencies (for visualization)
pip install matplotlib seaborn
```

## Usage

### Quick Start with Test Data

Before using your own data, we strongly recommend running the test pipeline to understand the complete workflow.

``` bash
# Test for identifying enhancers 

chrom-mesh identify \
  --h3k27ac-bed-dir test_data/identify_enhancer/h3k27ac_bed \
  --h3k27ac-bw-dir  test_data/identify_enhancer/h3k27ac_bw \
  --atac-bed-dir    test_data/identify_enhancer/atac_bed \
  --atac-bw-dir     test_data/identify_enhancer/atac_bw \
  --output test_output/identify_enhancer

# output files in ./test_output/identify_enhancer
# Results of each steps are in ./test_output/identify_enhancer/analysis_summary.txt
```

``` bash
# Test for annotating chromatin loops

chrom-mesh annotate \
  -l test_data/annotate_loops/loops.chr1_1_40000000.bedpe \
  -e test_data/annotate_loops/enhancers.chr1_1_40000000.bed \
  -p test_data/annotate_loops/promoters.chr1_1_40000000.bed \
  -o test_output/annotate_loops/loops.chr1_1_40000000.annotated.csv
  
# output files in ./test_output/annotate_loops
```

``` bash
# Test for constructing chromatin multi-enhancer loop networks

chrom-mesh analyze \
    -l ./test_data/analyze_networks/loops.chr1_1_40000000.annotated.csv \
    -e ./test_data/analyze_networks/enhancers.chr1_1_40000000.bed \
    -o ./test_output/analyze_networks

# output files in ./test_output/analyze_networks/gene_regulatory_depth.csv
```

``` bash
# Test for visualization
chrom-mesh visualize \
  -g ./test_data/visualize/gene_regulatory_depth.csv \
  -l ./test_data/visualize/loops.chr1_1_40000000.annotated.csv \
  -G Tcea1,Rb1cc1 \
  -o ./test_output/visualize/network_Tcea1_Rb1cc1.pdf
```

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/Badgerliu/EP_annotator/blob/main/LICENSE) file for details.

The MIT License allows for:

-   ✅ Commercial use

-   ✅ Modification

-   ✅ Distribution

-   ✅ Private use

Only requires:

-   📄 License and copyright notice

## Contact

For questions and support, please contact: kq_thx\@163.com

## Development Status

This project is currently under active development.
