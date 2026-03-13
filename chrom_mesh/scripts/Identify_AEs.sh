#!/usr/bin/env bash
# Author: Tao Hongxu
# Date: 2025-12-02 
# Description: Identify Active Enhancers (AEs) by integrating H3K27Ac CUT&Tag and ATAC-seq data.
# Usage: ./Identify_AEs.sh [options]

set -euo pipefail
shopt -s nullglob

# ------------------------------
# Help
# ------------------------------
show_help() {
    cat << EOF
Usage: ${0##*/} [OPTIONS]

Pipeline for identifying active enhancers

Required arguments:
    -h, --h3k27ac-bed DIR      H3K27Ac narrowPeak file directory
    -H, --h3k27ac-bw DIR       H3K27Ac BigWig file directory
    -a, --atac-bed DIR         ATAC-seq NFR BED file directory
    -A, --atac-bw DIR          ATAC-seq BigWig file directory
    -o, --output DIR           Output directory

Optional arguments:
    -b, --bigwig-tool PATH     bigWigAverageOverBed executable/path (default: auto from PATH: bigWigAverageOverBed)
        --bigwig-info-tool PATH bigWigInfo executable/path (default: auto from PATH: bigWigInfo)
        --bedtools PATH        bedtools executable/path (default: auto from PATH: bedtools)
        --python PATH          python executable/path (default: auto from PATH: python3)

    -p, --h3k27ac-pattern STR  H3K27Ac file matching pattern (default: *_peaks.narrowPeak)
    -n, --atac-pattern STR     ATAC BED file matching pattern (default: *.bed)
    -w, --h3k27ac-bw-pattern   H3K27Ac BigWig matching pattern (default: *.bw)
    -W, --atac-bw-pattern      ATAC BigWig matching pattern (default: *.bw)

    --help                     Show this help message

Example:
    ${0##*/} \\
        -h ./narrowpeak \\
        -H ./bw \\
        -a /path/to/ATAC \\
        -A /path/to/ATAC/bw \\
        -o ./results

EOF
    exit 0
}

# ------------------------------
# Default parameters
# ------------------------------
BIGWIG_AVG_EXEC="bigWigAverageOverBed"
BIGWIG_INFO_EXEC="bigWigInfo"
BEDTOOLS_EXEC="bedtools"
PYTHON_EXEC="python3"

H3K27AC_PATTERN="*_peaks.narrowPeak"
ATAC_BED_PATTERN="*.bed"
H3K27AC_BW_PATTERN="*.bw"
ATAC_BW_PATTERN="*.bw"

H3K27AC_BED_DIR=""
H3K27AC_BW_DIR=""
ATAC_BED_DIR=""
ATAC_BW_DIR=""
OUTPUT_DIR=""

# ------------------------------
# Utility functions
# ------------------------------
resolve_executable() {
    local exe_input="$1"
    local friendly="$2"

    if [[ "$exe_input" == */* ]]; then
        if [[ -x "$exe_input" ]]; then
            echo "$exe_input"
            return 0
        fi
        echo "Error: ${friendly} is not executable: $exe_input" >&2
        exit 1
    fi

    local resolved
    resolved="$(command -v "$exe_input" 2>/dev/null || true)"
    if [[ -n "$resolved" && -x "$resolved" ]]; then
        echo "$resolved"
        return 0
    fi

    echo "Error: ${friendly} not found in PATH: $exe_input" >&2
    exit 1
}

# ------------------------------
# Parse arguments
# ------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--h3k27ac-bed)
            H3K27AC_BED_DIR="$2"; shift 2 ;;
        -H|--h3k27ac-bw)
            H3K27AC_BW_DIR="$2"; shift 2 ;;
        -a|--atac-bed)
            ATAC_BED_DIR="$2"; shift 2 ;;
        -A|--atac-bw)
            ATAC_BW_DIR="$2"; shift 2 ;;
        -o|--output)
            OUTPUT_DIR="$2"; shift 2 ;;
        -b|--bigwig-tool)
            BIGWIG_AVG_EXEC="$2"; shift 2 ;;
        --bigwig-info-tool)
            BIGWIG_INFO_EXEC="$2"; shift 2 ;;
        --bedtools)
            BEDTOOLS_EXEC="$2"; shift 2 ;;
        --python)
            PYTHON_EXEC="$2"; shift 2 ;;
        -p|--h3k27ac-pattern)
            H3K27AC_PATTERN="$2"; shift 2 ;;
        -n|--atac-pattern)
            ATAC_BED_PATTERN="$2"; shift 2 ;;
        -w|--h3k27ac-bw-pattern)
            H3K27AC_BW_PATTERN="$2"; shift 2 ;;
        -W|--atac-bw-pattern)
            ATAC_BW_PATTERN="$2"; shift 2 ;;
        --help)
            show_help ;;
        *)
            echo "Error: Unknown argument '$1'" >&2
            echo "Use --help to see help information" >&2
            exit 1 ;;
    esac
done

# ------------------------------
# Validate required args
# ------------------------------
if [[ -z "$H3K27AC_BED_DIR" || -z "$H3K27AC_BW_DIR" || -z "$ATAC_BED_DIR" || -z "$ATAC_BW_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing required arguments" >&2
    echo "Use --help to see help information" >&2
    exit 1
fi

for dir in "$H3K27AC_BED_DIR" "$H3K27AC_BW_DIR" "$ATAC_BED_DIR" "$ATAC_BW_DIR"; do
    if [[ ! -d "$dir" ]]; then
        echo "Error: Directory does not exist: $dir" >&2
        exit 1
    fi
done

# ------------------------------
# Resolve executables
# ------------------------------
BIGWIG_AVG_EXEC="$(resolve_executable "$BIGWIG_AVG_EXEC" "bigWigAverageOverBed")"
BIGWIG_INFO_EXEC="$(resolve_executable "$BIGWIG_INFO_EXEC" "bigWigInfo")"
BEDTOOLS_EXEC="$(resolve_executable "$BEDTOOLS_EXEC" "bedtools")"
PYTHON_EXEC="$(resolve_executable "$PYTHON_EXEC" "python3")"

AWK_EXEC="$(resolve_executable "awk" "awk")"
SORT_EXEC="$(resolve_executable "sort" "sort")"
CAT_EXEC="$(resolve_executable "cat" "cat")"
UNIQ_EXEC="$(resolve_executable "uniq" "uniq")"
TEE_EXEC="$(resolve_executable "tee" "tee")"
WC_EXEC="$(resolve_executable "wc" "wc")"
DATE_EXEC="$(resolve_executable "date" "date")"
BASENAME_EXEC="$(resolve_executable "basename" "basename")"

# ------------------------------
# Create output dirs
# ------------------------------
PROCESS_H3K27AC="${OUTPUT_DIR}/process/H3K27Ac"
PROCESS_ATAC="${OUTPUT_DIR}/process/ATAC"
SCORE_DIR="${OUTPUT_DIR}/scoring_results"
LOG_DIR="${OUTPUT_DIR}/logs"
QC_DIR="${OUTPUT_DIR}/quality_control"
PLOT_DIR="${OUTPUT_DIR}/plots"
FINAL_DIR="${OUTPUT_DIR}/final_results"

mkdir -p "${PROCESS_H3K27AC}" "${PROCESS_ATAC}" "${SCORE_DIR}" "${LOG_DIR}" \
         "${QC_DIR}" "${PLOT_DIR}" "${FINAL_DIR}"

# ------------------------------
# Logging
# ------------------------------
LOG_FILE="${LOG_DIR}/pipeline_$(${DATE_EXEC} +%Y%m%d_%H%M%S).log"
exec > >("${TEE_EXEC}" -a "${LOG_FILE}")
exec 2>&1

start_time="$(${DATE_EXEC} +%s)"

echo "=========================================================="
echo "Active Enhancer Identification Pipeline"
echo "=========================================================="
echo "Started at: $(${DATE_EXEC})"
echo ""
echo "Configuration:"
echo "  H3K27Ac BED directory: ${H3K27AC_BED_DIR}"
echo "  H3K27Ac BigWig directory: ${H3K27AC_BW_DIR}"
echo "  ATAC BED directory: ${ATAC_BED_DIR}"
echo "  ATAC BigWig directory: ${ATAC_BW_DIR}"
echo "  Output directory: ${OUTPUT_DIR}"
echo ""
echo "Resolved tools:"
echo "  bigWigAverageOverBed: ${BIGWIG_AVG_EXEC}"
echo "  bigWigInfo:           ${BIGWIG_INFO_EXEC}"
echo "  bedtools:             ${BEDTOOLS_EXEC}"
echo "  python:               ${PYTHON_EXEC}"
echo "=========================================================="
echo ""

# ------------------------------
# Step 1: Process H3K27Ac
# ------------------------------
echo "[Step 1] Processing H3K27Ac narrowPeak files"
echo "----------------------------------------------------------"

H3K27AC_PEAK_FILES=( "${H3K27AC_BED_DIR}"/${H3K27AC_PATTERN} )
if [[ ${#H3K27AC_PEAK_FILES[@]} -eq 0 ]]; then
    echo "Error: No files matching ${H3K27AC_PATTERN} found in ${H3K27AC_BED_DIR}" >&2
    exit 1
fi

echo "Found ${#H3K27AC_PEAK_FILES[@]} H3K27Ac peak files"

H3K27AC_BED_FILES=()
for peak_file in "${H3K27AC_PEAK_FILES[@]}"; do
    [[ -e "$peak_file" ]] || continue

    basename_file="$(${BASENAME_EXEC} "${peak_file}")"
    sample_name="${basename_file%_peaks.narrowPeak}"
    sample_name="${sample_name%.narrowPeak}"
    output_bed="${PROCESS_H3K27AC}/${sample_name}.bed"

    echo "  Processing: ${basename_file} -> ${sample_name}.bed"
    "${AWK_EXEC}" 'BEGIN{FS=OFS="\t"} {print $1, $2, $3}' "${peak_file}" > "${output_bed}"

    H3K27AC_BED_FILES+=( "${output_bed}" )
done

echo "  Generated ${#H3K27AC_BED_FILES[@]} BED files"

H3K27AC_MERGED_BED="${PROCESS_H3K27AC}/H3K27Ac_merged.bed"
H3K27AC_FINAL_BED="${PROCESS_H3K27AC}/H3K27Ac_final_peaks.bed"

echo "  Merging all H3K27Ac peaks..."
"${CAT_EXEC}" "${H3K27AC_BED_FILES[@]}" | \
    "${SORT_EXEC}" -k1,1V -k2,2n | \
    "${BEDTOOLS_EXEC}" merge -i stdin > "${H3K27AC_MERGED_BED}"

"${AWK_EXEC}" 'BEGIN{OFS="\t"} {print $1, $2, $3, "H3K27Ac_peak_" NR}' \
    "${H3K27AC_MERGED_BED}" > "${H3K27AC_FINAL_BED}"

peak_count="$("${WC_EXEC}" -l < "${H3K27AC_FINAL_BED}")"
echo "  Final H3K27Ac peaks count: ${peak_count}"
echo "  Output file: ${H3K27AC_FINAL_BED}"
echo ""

# ------------------------------
# Step 2: Process ATAC
# ------------------------------
echo "[Step 2] Processing ATAC-seq NFR BED files"
echo "----------------------------------------------------------"

ATAC_BED_FILES=( "${ATAC_BED_DIR}"/${ATAC_BED_PATTERN} )
if [[ ${#ATAC_BED_FILES[@]} -eq 0 ]]; then
    echo "Error: No files matching ${ATAC_BED_PATTERN} found in ${ATAC_BED_DIR}" >&2
    exit 1
fi

echo "Found ${#ATAC_BED_FILES[@]} ATAC BED files"

ATAC_PROCESSED_FILES=()
for bed_file in "${ATAC_BED_FILES[@]}"; do
    [[ -e "$bed_file" ]] || continue

    basename_file="$(${BASENAME_EXEC} "${bed_file}")"
    sample_name="${basename_file%.bed}"
    output_bed="${PROCESS_ATAC}/${sample_name}.bed"

    echo "  Processing: ${basename_file}"
    "${AWK_EXEC}" 'BEGIN{FS=OFS="\t"} {print $1, $2, $3}' "${bed_file}" > "${output_bed}"

    ATAC_PROCESSED_FILES+=( "${output_bed}" )
done

echo "  Processed ${#ATAC_PROCESSED_FILES[@]} BED files"

ATAC_MERGED_BED="${PROCESS_ATAC}/ATAC_NFR_merged.bed"
ATAC_FINAL_BED="${PROCESS_ATAC}/ATAC_NFR_final_peaks.bed"

echo "  Merging all ATAC NFR peaks..."
"${CAT_EXEC}" "${ATAC_PROCESSED_FILES[@]}" | \
    "${SORT_EXEC}" -k1,1V -k2,2n | \
    "${BEDTOOLS_EXEC}" merge -i stdin > "${ATAC_MERGED_BED}"

"${AWK_EXEC}" 'BEGIN{OFS="\t"} {print $1, $2, $3, "NFR_" NR}' \
    "${ATAC_MERGED_BED}" > "${ATAC_FINAL_BED}"

nfr_count="$("${WC_EXEC}" -l < "${ATAC_FINAL_BED}")"
echo "  Final NFR peaks count: ${nfr_count}"
echo "  Output file: ${ATAC_FINAL_BED}"
echo ""

# ------------------------------
# Step 3: Collect valid BigWigs
# ------------------------------
echo "[Step 3] Collecting BigWig files"
echo "----------------------------------------------------------"

H3K27AC_BW_FILES=()
for bw_file in "${H3K27AC_BW_DIR}"/${H3K27AC_BW_PATTERN}; do
    [[ -e "$bw_file" ]] || continue
    if "${BIGWIG_INFO_EXEC}" "$bw_file" &> /dev/null; then
        H3K27AC_BW_FILES+=( "$bw_file" )
        echo "  ✓ H3K27Ac BigWig: $(${BASENAME_EXEC} "$bw_file")"
    else
        echo "  ✗ Skipping corrupted file: $(${BASENAME_EXEC} "$bw_file")"
    fi
done

if [[ ${#H3K27AC_BW_FILES[@]} -eq 0 ]]; then
    echo "Error: No valid H3K27Ac BigWig files found" >&2
    exit 1
fi
echo "  Found ${#H3K27AC_BW_FILES[@]} valid H3K27Ac BigWig files"

ATAC_BW_FILES=()
for bw_file in "${ATAC_BW_DIR}"/${ATAC_BW_PATTERN}; do
    [[ -e "$bw_file" ]] || continue
    if "${BIGWIG_INFO_EXEC}" "$bw_file" &> /dev/null; then
        ATAC_BW_FILES+=( "$bw_file" )
        echo "  ✓ ATAC BigWig: $(${BASENAME_EXEC} "$bw_file")"
    else
        echo "  ✗ Skipping corrupted file: $(${BASENAME_EXEC} "$bw_file")"
    fi
done

if [[ ${#ATAC_BW_FILES[@]} -eq 0 ]]; then
    echo "Error: No valid ATAC BigWig files found" >&2
    exit 1
fi
echo "  Found ${#ATAC_BW_FILES[@]} valid ATAC BigWig files"
echo ""

# ------------------------------
# Step 4: Score Peaks with BigWig signals
# ------------------------------
echo "[Step 4] Scoring peaks with BigWig signals"
echo "----------------------------------------------------------"

PYTHON_SCORING_SCRIPT=$(cat <<'PYTHON_EOF'
import os
import sys
import subprocess
import pandas as pd
import shutil

def run_bigWigAverageOverBed(executable, bigwig_file, bed_file, output_file):
    try:
        command = [executable, bigwig_file, bed_file, output_file]
        subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"    ✓ Successfully quantified: {os.path.basename(bigwig_file)}")
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to run bigWigAverageOverBed - {bigwig_file}", file=sys.stderr)
        print(f"Error message: {e.stderr}", file=sys.stderr)
        sys.exit(1)

def process_and_combine_scores(executable, bed_file, output_dir, output_prefix, bigwig_files):
    print(f"  Processing: {output_prefix}")
    print(f"    BED file: {os.path.basename(bed_file)}")
    print(f"    Number of BigWig files: {len(bigwig_files)}")

    temp_dir = os.path.join(output_dir, f"temp_{output_prefix}")
    os.makedirs(temp_dir, exist_ok=True)

    all_scores_df = None
    col_names = ["name", "size", "covered", "sum", "mean0", "mean"]

    for bigwig_file in bigwig_files:
        sample_name = os.path.basename(bigwig_file).replace(".bw", "").replace(".bigWig", "")
        temp_output_file = os.path.join(temp_dir, f"{sample_name}_scores.tab")

        run_bigWigAverageOverBed(executable, bigwig_file, bed_file, temp_output_file)

        try:
            scores_df = pd.read_csv(temp_output_file, sep='\t', header=None, names=col_names)
            current_scores = scores_df[['name', 'mean']].rename(columns={'mean': sample_name})

            if all_scores_df is None:
                all_scores_df = current_scores
            else:
                all_scores_df = pd.merge(all_scores_df, current_scores, on='name', how='outer')

        except pd.errors.EmptyDataError:
            print(f"    ⚠ Warning: Output from {bigwig_file} is empty, skipping", file=sys.stderr)
        except Exception as e:
            print(f"    ✗ Error: Failed to process {temp_output_file}: {e}", file=sys.stderr)
            sys.exit(1)

    if all_scores_df is not None:
        final_output_path = os.path.join(output_dir, f"combined_scores_{output_prefix}.tsv")
        all_scores_df.to_csv(final_output_path, sep='\t', index=False, float_format='%.5f')
        print(f"    ✓ Saved results: {os.path.basename(final_output_path)}")
        print(f"    ✓ Number of peaks: {len(all_scores_df)}")
        print(f"    ✓ Number of samples: {len(all_scores_df.columns) - 1}")
    else:
        print(f"    ✗ Warning: No results generated for {output_prefix}", file=sys.stderr)

    shutil.rmtree(temp_dir)

def main():
    if len(sys.argv) < 6:
        print("Usage: python3 script.py <executable> <bed_file> <output_dir> <output_prefix> <bigwig_file1> [bigwig_file2...]", file=sys.stderr)
        sys.exit(1)

    executable = sys.argv[1]
    bed_file = sys.argv[2]
    output_dir = sys.argv[3]
    output_prefix = sys.argv[4]
    bigwig_files = sys.argv[5:]

    if not os.path.exists(bed_file):
        print(f"Error: BED file does not exist: {bed_file}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(output_dir):
        print(f"Error: Output directory does not exist: {output_dir}", file=sys.stderr)
        sys.exit(1)

    process_and_combine_scores(executable, bed_file, output_dir, output_prefix, bigwig_files)

if __name__ == "__main__":
    main()
PYTHON_EOF
)

echo "Scoring H3K27Ac peaks..."
"${PYTHON_EXEC}" -c "$PYTHON_SCORING_SCRIPT" \
    "${BIGWIG_AVG_EXEC}" \
    "${H3K27AC_FINAL_BED}" \
    "${SCORE_DIR}" \
    "H3K27Ac" \
    "${H3K27AC_BW_FILES[@]}"
echo ""

echo "Scoring ATAC NFR peaks..."
"${PYTHON_EXEC}" -c "$PYTHON_SCORING_SCRIPT" \
    "${BIGWIG_AVG_EXEC}" \
    "${ATAC_FINAL_BED}" \
    "${SCORE_DIR}" \
    "ATAC_NFR" \
    "${ATAC_BW_FILES[@]}"
echo ""

# ------------------------------
# Step 5: H3K27Ac QC
# ------------------------------
echo "[Step 5] H3K27Ac Peak Quality Control and Filtering"
echo "----------------------------------------------------------"

PYTHON_H3K27AC_QC=$(cat <<'PYTHON_EOF'
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def filter_h3k27ac_peaks(score_file, bed_file, output_dir, threshold=1, fold_change_max=3):
    print(f"  Reading signal matrix: {score_file}")
    df_scores = pd.read_csv(score_file, sep='\t')

    print(f"  Reading BED file: {bed_file}")
    df_bed = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'name'])

    df = pd.merge(df_bed, df_scores, on='name', how='inner')
    total_peaks = len(df)
    print(f"  Total peaks: {total_peaks}")

    signal_cols = [col for col in df.columns if col not in ['chr', 'start', 'end', 'name']]

    if len(signal_cols) < 2:
        print("  Warning: Less than 2 replicates, skipping fold change filtering")
        consistent_peaks = df[df[signal_cols].min(axis=1) > threshold]
    else:
        print(f"  Calculating Fold Change (using {len(signal_cols)} replicates)...")
        df["FoldChange"] = df[signal_cols].max(axis=1) / df[signal_cols].min(axis=1)

        print(f"  Applying filtering criteria:")
        print(f"    - All replicate signals > {threshold}")
        print(f"    - Fold Change < {fold_change_max}")

        consistent_peaks = df[
            (df[signal_cols].min(axis=1) > threshold) &
            (df["FoldChange"] < fold_change_max)
        ]

    num_consistent = len(consistent_peaks)
    percentage = (num_consistent / total_peaks * 100) if total_peaks > 0 else 0.0

    print(f"\n  Filtering results:")
    print(f"    Total peaks: {total_peaks}")
    print(f"    Consistent peaks: {num_consistent}")
    print(f"    Retention rate: {percentage:.2f}%")

    output_bed = f"{output_dir}/filtered_H3K27Ac_peaks.bed"
    consistent_peaks[['chr', 'start', 'end', 'name']].to_csv(output_bed, sep='\t', header=False, index=False)
    print(f"  Saved filtered peaks: {output_bed}")

    output_full = f"{output_dir}/filtered_H3K27Ac_peaks_with_signals.tsv"
    consistent_peaks.to_csv(output_full, sep='\t', index=False)
    print(f"  Saved complete data: {output_full}")

    plot_file = f"{output_dir}/H3K27Ac_QC_plots.pdf"
    print(f"  Generating visualization plots: {plot_file}")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    if len(signal_cols) >= 2:
        ax = axes[0]
        ax.scatter(df[signal_cols[0]], df[signal_cols[1]], alpha=0.3, s=10, label='All Peaks', color='gray')
        ax.scatter(consistent_peaks[signal_cols[0]], consistent_peaks[signal_cols[1]], alpha=0.5, s=10, label='Consistent Peaks', color='red')
        ax.set_xlabel(signal_cols[0])
        ax.set_ylabel(signal_cols[1])
        ax.set_title('Replicate Consistency')
        ax.legend()
        ax.set_xscale('log')
        ax.set_yscale('log')

    ax = axes[1]
    counts = [total_peaks, num_consistent]
    labels = ['Total Peaks', 'Consistent Peaks']
    colors = ['skyblue', 'green']
    bars = ax.bar(labels, counts, color=colors)
    ax.set_ylabel('Number of Peaks')
    ax.set_title('Peak Filtering Summary')
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        pct = (count / total_peaks * 100) if total_peaks > 0 else 0
        ax.text(bar.get_x() + bar.get_width()/2., height, f'{count}\n({pct:.1f}%)', ha='center', va='bottom')

    ax = axes[2]
    for col in signal_cols:
        ax.hist(np.log10(df[col] + 1), bins=50, alpha=0.5, label=col)
    ax.axvline(np.log10(threshold + 1), color='red', linestyle='--', label=f'Threshold={threshold}')
    ax.set_xlabel('log10(Signal + 1)')
    ax.set_ylabel('Frequency')
    ax.set_title('Signal Distribution')
    ax.legend()

    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <score_file> <bed_file> <output_dir>")
        sys.exit(1)

    filter_h3k27ac_peaks(sys.argv[1], sys.argv[2], sys.argv[3])
PYTHON_EOF
)

"${PYTHON_EXEC}" -c "$PYTHON_H3K27AC_QC" \
    "${SCORE_DIR}/combined_scores_H3K27Ac.tsv" \
    "${H3K27AC_FINAL_BED}" \
    "${QC_DIR}"

FILTERED_H3K27AC_BED="${QC_DIR}/filtered_H3K27Ac_peaks.bed"
echo ""

# ------------------------------
# Step 6: ATAC NFR QC
# ------------------------------
echo "[Step 6] ATAC NFR Quality Control and Filtering"
echo "----------------------------------------------------------"

PYTHON_NFR_QC=$(cat <<'PYTHON_EOF'
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def filter_nfr_peaks(score_file, bed_file, output_dir, signal_threshold=0.5, fold_change_max=5, max_width=1500):
    print(f"  Reading signal matrix: {score_file}")
    df_scores = pd.read_csv(score_file, sep='\t')

    print(f"  Reading BED file: {bed_file}")
    df_bed = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'name'])
    df_bed['Width'] = df_bed['end'] - df_bed['start']

    df = pd.merge(df_bed, df_scores, on='name', how='inner')
    total_peaks = len(df)
    print(f"  Total NFRs: {total_peaks}")

    signal_cols = [col for col in df.columns if col not in ['chr', 'start', 'end', 'name', 'Width']]
    print(f"  Using {len(signal_cols)} replicates for filtering")

    df["FoldChange"] = df[signal_cols].max(axis=1) / df[signal_cols].min(axis=1)

    print("  Applying filtering criteria:")
    print(f"    - All replicate signals > {signal_threshold}")
    print(f"    - Fold Change < {fold_change_max}")
    print(f"    - Width <= {max_width} bp")

    consistent_peaks = df[
        (df[signal_cols].min(axis=1) > signal_threshold) &
        (df["FoldChange"] < fold_change_max)
    ]
    num_after_signal = len(consistent_peaks)

    filtered_peaks = consistent_peaks[consistent_peaks['Width'] <= max_width]
    num_final = len(filtered_peaks)

    print(f"\n  Filtering results:")
    print(f"    Total NFRs: {total_peaks}")
    p1 = (num_after_signal / total_peaks * 100) if total_peaks > 0 else 0
    p2 = (num_final / total_peaks * 100) if total_peaks > 0 else 0
    print(f"    After signal/consistency filtering: {num_after_signal} ({p1:.2f}%)")
    print(f"    After width filtering: {num_final} ({p2:.2f}%)")

    print(f"\n  Width statistics:")
    print(f"    Mean width: {df['Width'].mean():.1f} bp")
    print(f"    Median width: {df['Width'].median():.1f} bp")
    print(f"    Maximum width: {df['Width'].max()} bp")
    if len(filtered_peaks) > 0:
        print(f"    Mean width after filtering: {filtered_peaks['Width'].mean():.1f} bp")
    else:
        print(f"    Mean width after filtering: NA")

    output_bed = f"{output_dir}/HEPM_NFR_merged_fixed_no_abnormal.bed"
    filtered_peaks[['chr', 'start', 'end', 'name']].to_csv(output_bed, sep='\t', header=False, index=False)
    print(f"\n  Saved filtered NFRs: {output_bed}")

    output_full = f"{output_dir}/filtered_NFR_peaks_with_signals.tsv"
    filtered_peaks.to_csv(output_full, sep='\t', index=False)
    print(f"  Saved complete data: {output_full}")

    plot_file = f"{output_dir}/NFR_QC_plots.pdf"
    print(f"  Generating visualization plots: {plot_file}")

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    if len(signal_cols) >= 2:
        ax = axes[0, 0]
        ax.scatter(df[signal_cols[0]], df[signal_cols[1]], alpha=0.3, s=10, label='All NFRs', color='gray')
        ax.scatter(filtered_peaks[signal_cols[0]], filtered_peaks[signal_cols[1]], alpha=0.5, s=10, label='Filtered NFRs', color='blue')
        ax.set_xlabel(signal_cols[0]); ax.set_ylabel(signal_cols[1])
        ax.set_title('NFR Replicate Consistency')
        ax.legend()
        ax.set_xscale('log'); ax.set_yscale('log')

    ax = axes[0, 1]
    counts = [total_peaks, num_after_signal, num_final]
    labels = ['Total', 'After Signal\nFiltering', 'After Width\nFiltering']
    colors = ['skyblue', 'orange', 'green']
    bars = ax.bar(labels, counts, color=colors)
    ax.set_ylabel('Number of NFRs')
    ax.set_title('NFR Filtering Summary')
    for bar, count in zip(bars, counts):
        h = bar.get_height()
        pct = (count / total_peaks * 100) if total_peaks > 0 else 0
        ax.text(bar.get_x() + bar.get_width()/2., h, f'{count}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=9)

    ax = axes[1, 0]
    ax.hist(df['Width'], bins=50, alpha=0.5, label='Before', color='gray')
    ax.hist(filtered_peaks['Width'], bins=50, alpha=0.7, label='After', color='green')
    ax.axvline(max_width, color='red', linestyle='--', label=f'Max Width={max_width}')
    ax.set_xlabel('NFR Width (bp)')
    ax.set_ylabel('Frequency')
    ax.set_title('NFR Width Distribution')
    ax.legend()

    ax = axes[1, 1]
    for col in signal_cols:
        ax.hist(np.log10(df[col] + 1), bins=50, alpha=0.5, label=col)
    ax.axvline(np.log10(signal_threshold + 1), color='red', linestyle='--', label=f'Threshold={signal_threshold}')
    ax.set_xlabel('log10(Signal + 1)')
    ax.set_ylabel('Frequency')
    ax.set_title('NFR Signal Distribution')
    ax.legend()

    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <score_file> <bed_file> <output_dir>")
        sys.exit(1)

    filter_nfr_peaks(sys.argv[1], sys.argv[2], sys.argv[3])
PYTHON_EOF
)

"${PYTHON_EXEC}" -c "$PYTHON_NFR_QC" \
    "${SCORE_DIR}/combined_scores_ATAC_NFR.tsv" \
    "${ATAC_FINAL_BED}" \
    "${QC_DIR}"

FILTERED_NFR_BED="${QC_DIR}/HEPM_NFR_merged_fixed_no_abnormal.bed"
echo ""

# ------------------------------
# Step 7: Generate flanked regions
# ------------------------------
echo "[Step 7] Generating H3K27Ac Flanked Regions"
echo "----------------------------------------------------------"

PYTHON_FLANKED=$(cat <<'PYTHON_EOF'
import sys

def generate_flanked_regions(input_file, output_file, min_size=100, max_size=1500):
    print(f"  Reading input file: {input_file}")
    print(f"  Size range: {min_size}-{max_size} bp")

    filtered_regions = []
    try:
        with open(input_file, 'r') as infile:
            prev_chr, prev_end = None, None
            for line in infile:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                chr_, start, end = parts[0], int(parts[1]), int(parts[2])

                if prev_chr == chr_:
                    gap_size = start - prev_end
                    if min_size <= gap_size <= max_size:
                        filtered_regions.append([chr_, prev_end, start])

                prev_chr, prev_end = chr_, end

        with open(output_file, 'w') as out:
            for region in filtered_regions:
                out.write('\t'.join(map(str, region)) + '\n')

        print(f"  Found {len(filtered_regions)} flanked regions")
        print(f"  Saved to: {output_file}")
        return len(filtered_regions)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_bed> <output_bed> <min_size> <max_size>")
        sys.exit(1)

    generate_flanked_regions(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
PYTHON_EOF
)

H3K27AC_FLANKED_BED="${QC_DIR}/H3K27Ac_Flanked_region.bed"
"${PYTHON_EXEC}" -c "$PYTHON_FLANKED" \
    "${FILTERED_H3K27AC_BED}" \
    "${H3K27AC_FLANKED_BED}" \
    100 \
    1500
echo ""

# ------------------------------
# Step 8: Overlap analysis
# ------------------------------
echo "[Step 8] Finding NFR overlaps with H3K27Ac"
echo "----------------------------------------------------------"

NFR_OVERLAP_PEAKS="${FINAL_DIR}/HEPM_NFR_overlapping_H3K27ac_peaks.bed"
echo "  Calculating NFR overlap with H3K27Ac peaks..."
"${BEDTOOLS_EXEC}" intersect -wa -a "${FILTERED_NFR_BED}" -b "${FILTERED_H3K27AC_BED}" > "${NFR_OVERLAP_PEAKS}"
overlap_peaks_count="$("${WC_EXEC}" -l < "${NFR_OVERLAP_PEAKS}")"
echo "  Found ${overlap_peaks_count} NFRs overlapping with H3K27Ac peaks"

NFR_OVERLAP_FLANKED="${FINAL_DIR}/HEPM_NFR_overlapping_H3K27ac_flanked.bed"
echo "  Calculating NFR overlap with H3K27Ac flanked regions..."
"${BEDTOOLS_EXEC}" intersect -wa -a "${FILTERED_NFR_BED}" -b "${H3K27AC_FLANKED_BED}" > "${NFR_OVERLAP_FLANKED}"
overlap_flanked_count="$("${WC_EXEC}" -l < "${NFR_OVERLAP_FLANKED}")"
echo "  Found ${overlap_flanked_count} NFRs overlapping with H3K27Ac flanked regions"

NFR_COMBINED="${FINAL_DIR}/HEPM_NFR_combined_overlaps.bed"
echo "  Merging both overlap types..."
"${CAT_EXEC}" "${NFR_OVERLAP_PEAKS}" "${NFR_OVERLAP_FLANKED}" > "${NFR_COMBINED}"

NFR_FINAL="${FINAL_DIR}/HEPM_NFR_combined_overlaps_sorted.bed"
echo "  Sorting and removing duplicates..."
"${SORT_EXEC}" -k1,1 -k2,2n "${NFR_COMBINED}" | "${UNIQ_EXEC}" > "${NFR_FINAL}"

final_count="$("${WC_EXEC}" -l < "${NFR_FINAL}")"
echo "  Final NFR count: ${final_count}"
echo "  Output file: ${NFR_FINAL}"
echo ""

# ------------------------------
# Step 9: Final summary
# ------------------------------
echo "[Step 9] Generating Final Summary"
echo "----------------------------------------------------------"

SUMMARY_FILE="${OUTPUT_DIR}/analysis_summary.txt"

cat > "${SUMMARY_FILE}" << EOF
========================================================
Active Enhancer Identification Pipeline - Complete Analysis Summary
========================================================
Analysis time: $(${DATE_EXEC})

Input Data:
---------------------------------------------------------
H3K27Ac Peaks directory: ${H3K27AC_BED_DIR}
H3K27Ac BigWig directory: ${H3K27AC_BW_DIR}
ATAC BED directory: ${ATAC_BED_DIR}
ATAC BigWig directory: ${ATAC_BW_DIR}

Output directory: ${OUTPUT_DIR}

Resolved tools:
---------------------------------------------------------
bigWigAverageOverBed: ${BIGWIG_AVG_EXEC}
bigWigInfo: ${BIGWIG_INFO_EXEC}
bedtools: ${BEDTOOLS_EXEC}
python: ${PYTHON_EXEC}

Processing Results:
---------------------------------------------------------
H3K27Ac:
  - Number of input peak files: ${#H3K27AC_PEAK_FILES[@]}
  - Number of BigWig files: ${#H3K27AC_BW_FILES[@]}
  - Merged peaks count: ${peak_count}
  - QC-filtered peaks count: $("${WC_EXEC}" -l < "${FILTERED_H3K27AC_BED}")
  - Flanked regions count: $("${WC_EXEC}" -l < "${H3K27AC_FLANKED_BED}")

ATAC-seq NFR:
  - Number of input BED files: ${#ATAC_BED_FILES[@]}
  - Number of BigWig files: ${#ATAC_BW_FILES[@]}
  - Merged NFR count: ${nfr_count}
  - QC-filtered NFR count: $("${WC_EXEC}" -l < "${FILTERED_NFR_BED}")

Overlap Analysis:
  - NFR overlapping with H3K27Ac peaks: ${overlap_peaks_count}
  - NFR overlapping with H3K27Ac flanked regions: ${overlap_flanked_count}
  - Final merged NFR count: ${final_count}

Key Output Files:
---------------------------------------------------------
QC-filtered files:
  - ${FILTERED_H3K27AC_BED}
  - ${FILTERED_NFR_BED}
  - ${H3K27AC_FLANKED_BED}

Final results:
  - ${NFR_FINAL}

Visualization plots:
  - ${QC_DIR}/H3K27Ac_QC_plots.pdf
  - ${QC_DIR}/NFR_QC_plots.pdf

Signal matrices:
  - ${SCORE_DIR}/combined_scores_H3K27Ac.tsv
  - ${SCORE_DIR}/combined_scores_ATAC_NFR.tsv
========================================================
EOF

"${CAT_EXEC}" "${SUMMARY_FILE}"

end_time="$(${DATE_EXEC} +%s)"
runtime=$((end_time - start_time))
minutes=$((runtime / 60))
seconds=$((runtime % 60))

echo ""
echo "=========================================================="
echo "Pipeline Completed!"
echo "=========================================================="
echo "Completion time: $(${DATE_EXEC})"
echo "Total runtime: ${minutes} minutes ${seconds} seconds"
echo "Log file: ${LOG_FILE}"
echo "Summary file: ${SUMMARY_FILE}"
echo "=========================================================="
