#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse

TARGET_CHR = "chr1"
TARGET_START = 1
TARGET_END = 40_000_000

OUT_DIR = Path("./test_data/annotate_loops")


def overlap(start: int, end: int, qstart: int, qend: int) -> bool:
    # 区间 [start,end) 与 [qstart,qend) 是否有重叠
    return not (end <= qstart or start >= qend)


def normalize_chr(s: str) -> str:
    return str(s).strip()


def parse_int(x: str):
    try:
        return int(x)
    except Exception:
        return None


def filter_enhancer(in_file: Path, out_file: Path):
    """
    enhancer: 前4列 = chr, start, end, id
    """
    total, kept = 0, 0
    with in_file.open("r", encoding="utf-8", errors="ignore") as fin, \
         out_file.open("w", encoding="utf-8") as fout:
        for line in fin:
            s = line.strip()
            if not s or s.startswith(("#", "track", "browser")):
                continue
            total += 1
            parts = s.split()  # 兼容空格/tab
            if len(parts) < 4:
                continue

            chrom = normalize_chr(parts[0])
            st = parse_int(parts[1])
            ed = parse_int(parts[2])
            if st is None or ed is None or st >= ed:
                continue
            if chrom != TARGET_CHR:
                continue
            if not overlap(st, ed, TARGET_START, TARGET_END):
                continue

            # 裁剪到目标区间
            parts[1] = str(max(st, TARGET_START))
            parts[2] = str(min(ed, TARGET_END))
            fout.write("\t".join(parts[:4]) + "\n")
            kept += 1

    print(f"[enhancer] total={total}, kept={kept}, out={out_file}")


def filter_promoter(in_file: Path, out_file: Path):
    """
    promoter: 前7列 = chr, start, end, gene, gene_chr, gene_start, gene_end
    """
    total, kept = 0, 0
    with in_file.open("r", encoding="utf-8", errors="ignore") as fin, \
         out_file.open("w", encoding="utf-8") as fout:
        for line in fin:
            s = line.strip()
            if not s or s.startswith(("#", "track", "browser")):
                continue
            total += 1
            parts = s.split()
            if len(parts) < 7:
                continue

            chrom = normalize_chr(parts[0])
            st = parse_int(parts[1])
            ed = parse_int(parts[2])

            gst = parse_int(parts[5])
            ged = parse_int(parts[6])

            if st is None or ed is None or st >= ed:
                continue
            if gst is None or ged is None or gst >= ged:
                continue
            if chrom != TARGET_CHR:
                continue
            if not overlap(st, ed, TARGET_START, TARGET_END):
                continue

            # 裁剪 promoter 区间（前三列）
            parts[1] = str(max(st, TARGET_START))
            parts[2] = str(min(ed, TARGET_END))
            fout.write("\t".join(parts[:7]) + "\n")
            kept += 1

    print(f"[promoter] total={total}, kept={kept}, out={out_file}")


def detect_loop_layout(parts):
    """
    返回 (chr1_i, start1_i, end1_i, chr2_i, start2_i, end2_i)
    优先标准布局 [0,1,2,3,4,5]，否则尝试备用 [0,1,2,4,5,6]
    """
    # 标准
    if len(parts) >= 6:
        c1, s1, e1, c2, s2, e2 = 0, 1, 2, 3, 4, 5
        if parse_int(parts[s1]) is not None and parse_int(parts[e1]) is not None \
           and parse_int(parts[s2]) is not None and parse_int(parts[e2]) is not None:
            return (c1, s1, e1, c2, s2, e2)

    # 备用
    if len(parts) >= 7:
        c1, s1, e1, c2, s2, e2 = 0, 1, 2, 4, 5, 6
        if parse_int(parts[s1]) is not None and parse_int(parts[e1]) is not None \
           and parse_int(parts[s2]) is not None and parse_int(parts[e2]) is not None:
            return (c1, s1, e1, c2, s2, e2)

    return None


def filter_loops(in_file: Path, out_file: Path):
    """
    loop BEDPE: 至少6列
    保留：两端都在 chr1，且各自与目标区间重叠
    """
    total, kept = 0, 0
    with in_file.open("r", encoding="utf-8", errors="ignore") as fin, \
         out_file.open("w", encoding="utf-8") as fout:
        for line in fin:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            total += 1
            parts = s.split()
            layout = detect_loop_layout(parts)
            if layout is None:
                continue

            c1i, s1i, e1i, c2i, s2i, e2i = layout
            chr1 = normalize_chr(parts[c1i])
            st1 = parse_int(parts[s1i]); ed1 = parse_int(parts[e1i])
            chr2 = normalize_chr(parts[c2i])
            st2 = parse_int(parts[s2i]); ed2 = parse_int(parts[e2i])

            if None in (st1, ed1, st2, ed2):
                continue
            if st1 >= ed1 or st2 >= ed2:
                continue

            if chr1 != TARGET_CHR or chr2 != TARGET_CHR:
                continue
            if not overlap(st1, ed1, TARGET_START, TARGET_END):
                continue
            if not overlap(st2, ed2, TARGET_START, TARGET_END):
                continue

            # 裁剪两端坐标并输出标准6列
            st1 = max(st1, TARGET_START); ed1 = min(ed1, TARGET_END)
            st2 = max(st2, TARGET_START); ed2 = min(ed2, TARGET_END)
            if st1 >= ed1 or st2 >= ed2:
                continue

            fout.write(f"{chr1}\t{st1}\t{ed1}\t{chr2}\t{st2}\t{ed2}\n")
            kept += 1

    print(f"[loop] total={total}, kept={kept}, out={out_file}")


def main():
    ap = argparse.ArgumentParser(
        description="Extract demo annotate files in chr1:1-10000000 to ./test_data/annotate_loops"
    )
    ap.add_argument("--promoter", required=True, help="Path to promoter file (7 columns)")
    ap.add_argument("--enhancer", required=True, help="Path to enhancer file (4 columns)")
    ap.add_argument("--loop", required=True, help="Path to loop BEDPE file (>=6 columns)")
    args = ap.parse_args()

    promoter = Path(args.promoter)
    enhancer = Path(args.enhancer)
    loop = Path(args.loop)

    for fp in [promoter, enhancer, loop]:
        if not fp.exists():
            raise FileNotFoundError(f"Not found: {fp}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    out_promoter = OUT_DIR / f"promoters.{TARGET_CHR}_{TARGET_START}_{TARGET_END}.bed"
    out_enhancer = OUT_DIR / f"enhancers.{TARGET_CHR}_{TARGET_START}_{TARGET_END}.bed"
    out_loop = OUT_DIR / f"loops.{TARGET_CHR}_{TARGET_START}_{TARGET_END}.bedpe"

    print(f"Target region: {TARGET_CHR}:{TARGET_START}-{TARGET_END}")
    print(f"Output dir: {OUT_DIR.resolve()}")

    filter_promoter(promoter, out_promoter)
    filter_enhancer(enhancer, out_enhancer)
    filter_loops(loop, out_loop)

    print("Done.")


if __name__ == "__main__":
    main()
