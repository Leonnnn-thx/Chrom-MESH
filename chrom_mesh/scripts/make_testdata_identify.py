#!/usr/bin/env python3
from pathlib import Path
import shutil
import random
import subprocess
import tempfile

# ====== 你需要改这4个源目录 ======
SRC_H3K27AC_BED = Path("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/1.E145_Mes_AEs_Pros/H3K27ac_CUTTag")
SRC_H3K27AC_BW  = Path("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/1.E145_Mes_AEs_Pros/H3K27ac_CUTTag")
SRC_ATAC_BED    = Path("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/1.E145_Mes_AEs_Pros/ATAC")
SRC_ATAC_BW     = Path("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/1.E145_Mes_AEs_Pros/ATAC")

OUT_ROOT = Path("test_data/identify_enhancer")

MAX_SAMPLES_PER_TYPE = 3
BED_MAX_LINES = 1000000
RANDOM_SEED = 42

TARGET_CHR = "chr1"
TARGET_START = 1
TARGET_END = 10_000_000

def pick_files(src: Path, suffixes):
    files = [p for p in src.iterdir() if p.is_file() and p.suffix.lower() in suffixes]
    files.sort()
    random.shuffle(files)
    return files[:MAX_SAMPLES_PER_TYPE]

def has_cmd(cmd: str) -> bool:
    return shutil.which(cmd) is not None

def run_cmd(cmd):
    subprocess.run(cmd, check=True)

def build_chrom_sizes_from_bw(in_bw: Path, out_sizes: Path):
    
    p = subprocess.run(
        ["bigWigInfo", "-chroms", str(in_bw)],
        check=True, capture_output=True, text=True
    )
    lines = []
    for line in p.stdout.splitlines():
        parts = line.strip().split()
        if len(parts) >= 3 and parts[2].isdigit():
            chrom = parts[0]
            size = parts[2]
            if not chrom.endswith(":"):  
                lines.append(f"{chrom}\t{size}\n")
    out_sizes.write_text("".join(lines), encoding="utf-8")

def crop_bw_files(src: Path, dst: Path):
    if not (has_cmd("bigWigToBedGraph") and has_cmd("bedGraphToBigWig") and has_cmd("bigWigInfo")):
        raise RuntimeError("缺少 UCSC 工具：bigWigToBedGraph / bedGraphToBigWig / bigWigInfo")

    dst.mkdir(parents=True, exist_ok=True)
    picked = pick_files(src, {".bw", ".bigwig"})
    out_files = []

    for f in picked:
        stem = f.stem
        out_bw = dst / f"{stem}.{TARGET_CHR}_{TARGET_START}_{TARGET_END}.bw"

        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            bg = td / "region.bedGraph"
            bg_sorted = td / "region.sorted.bedGraph"
            sizes = td / "chrom.sizes"

            run_cmd([
                "bigWigToBedGraph",
                f"-chrom={TARGET_CHR}",
                f"-start={TARGET_START}",
                f"-end={TARGET_END}",
                str(f),
                str(bg)
            ])

            build_chrom_sizes_from_bw(f, sizes)

        
            with bg.open("r", encoding="utf-8", errors="ignore") as fin:
                rows = [ln for ln in fin if ln.strip()]
            rows.sort(key=lambda x: (x.split()[0], int(x.split()[1])))
            with bg_sorted.open("w", encoding="utf-8") as fout:
                fout.writelines(rows)

            run_cmd(["bedGraphToBigWig", str(bg_sorted), str(sizes), str(out_bw)])
            out_files.append(out_bw)

    return picked, out_files

def truncate_and_filter_bed_files(src: Path, dst: Path, mark_h3k27ac: bool = False):
    dst.mkdir(parents=True, exist_ok=True)
    picked = pick_files(src, [".bed", ".narrowpeak"])
    out_files = []

    if not picked:
        print(f"[WARN] No BED/narrowPeak files found in: {src}")

    for f in picked:
        base = f"{f.stem}.{TARGET_CHR}_{TARGET_START}_{TARGET_END}"
        if mark_h3k27ac:
            out = dst / f"{base}_peaks.narrowPeak"
        else:
            out = dst / f"{base}.bed"

        n_written = 0
        n_total = 0
        n_badcol = 0
        n_badcoord = 0
        n_wrong_chr = 0
        n_no_overlap = 0

        with f.open("r", encoding="utf-8", errors="ignore") as fin, \
             out.open("w", encoding="utf-8") as fout:
            for line in fin:
                s = line.strip()
                if not s or s.startswith(("#", "track", "browser")):
                    continue

                n_total += 1
                parts = s.split()  
                if len(parts) < 3:
                    n_badcol += 1
                    continue

                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    n_badcoord += 1
                    continue

                if chrom != TARGET_CHR:
                    n_wrong_chr += 1
                    continue

                if end <= TARGET_START or start >= TARGET_END:
                    n_no_overlap += 1
                    continue

            
                new_start = max(start, TARGET_START)
                new_end = min(end, TARGET_END)
                if new_start >= new_end:
                    n_badcoord += 1
                    continue

                parts[1] = str(new_start)
                parts[2] = str(new_end)
                fout.write("\t".join(parts) + "\n")
                n_written += 1

                if n_written >= BED_MAX_LINES:
                    break

        out_files.append(out)
        print(
            f"[BED] {f.name}: total={n_total}, kept={n_written}, "
            f"badcol={n_badcol}, badcoord={n_badcoord}, wrong_chr={n_wrong_chr}, no_overlap={n_no_overlap}"
        )

    return picked, out_files

def main():
    random.seed(RANDOM_SEED)

    out_h3_bed = OUT_ROOT / "h3k27ac_bed"
    out_h3_bw  = OUT_ROOT / "h3k27ac_bw"
    out_atac_bed = OUT_ROOT / "atac_bed"
    out_atac_bw  = OUT_ROOT / "atac_bw"

    print("Generating test data under:", OUT_ROOT.resolve())
    print(f"Target region: {TARGET_CHR}:{TARGET_START}-{TARGET_END}")

    h3_bed_src, h3_bed_out = truncate_and_filter_bed_files(SRC_H3K27AC_BED, out_h3_bed, mark_h3k27ac=True)
    h3_bw_src,  h3_bw_out  = crop_bw_files(SRC_H3K27AC_BW, out_h3_bw)
    a_bed_src,  a_bed_out  = truncate_and_filter_bed_files(SRC_ATAC_BED, out_atac_bed, mark_h3k27ac=False)
    a_bw_src,   a_bw_out   = crop_bw_files(SRC_ATAC_BW, out_atac_bw)

    print(f"H3K27ac BED: {len(h3_bed_src)} selected, {len(h3_bed_out)} generated")
    print(f"H3K27ac BW : {len(h3_bw_src)} selected, {len(h3_bw_out)} generated")
    print(f"ATAC BED   : {len(a_bed_src)} selected, {len(a_bed_out)} generated")
    print(f"ATAC BW    : {len(a_bw_src)} selected, {len(a_bw_out)} generated")
    print("Done.")

if __name__ == "__main__":
    main()
