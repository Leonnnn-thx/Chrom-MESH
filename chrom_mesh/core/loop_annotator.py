#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Loop Annotator Module
=====================
Annotates chromatin loops by classifying anchors as promoters/enhancers.
Based on specific logic for 7-column promoter files and 4-column enhancer files.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict

import pandas as pd

logger = logging.getLogger(__name__)

OUTPUT_COLUMNS = [
    "chr1", "start1", "end1", "chr2", "start2", "end2",
    "P_L", "P_R", "E_L", "E_R",
    "P_L_coords", "P_R_coords",
    "E_L_coords", "E_R_coords",
    "P_L_genes", "P_R_genes",
    "E_L_IDs", "E_R_IDs",
    "pure_PP", "pure_EE", "pure_EP", "mixed_EP",
    "intersection_type",
]


class LoopAnnotator:
    def __init__(
        self,
        loop_file: str | Path,
        enhancer_file: str | Path,
        promoter_file: str | Path,
        output_file: str | Path,
    ) -> None:
        self.loop_file = Path(loop_file)
        self.enhancer_file = Path(enhancer_file)
        self.promoter_file = Path(promoter_file)
        self.output_file = Path(output_file)

        self._validate_files_exist()

        self.loops: pd.DataFrame | None = None
        self.enhancers_dict: Dict[str, pd.DataFrame] = {}
        self.promoters_dict: Dict[str, pd.DataFrame] = {}
        self.annotated_loops: pd.DataFrame | None = None

        logger.info("LoopAnnotator initialized")

    def _validate_files_exist(self) -> None:
        for fp in [self.loop_file, self.enhancer_file, self.promoter_file]:
            if not fp.exists():
                raise FileNotFoundError(f"File not found: {fp}")

    @staticmethod
    def _normalize_chr(series: pd.Series) -> pd.Series:
        return series.astype(str).str.strip()

    def _load_loops_bedpe(self) -> pd.DataFrame:
        raw = pd.read_csv(
            self.loop_file,
            sep="\t",
            header=None,
            comment="#",
            dtype=str,
            engine="c",
        )

        if raw.shape[1] < 6:
            raise ValueError(f"Loop BEDPE has too few columns ({raw.shape[1]}). Need at least 6.")

        use_standard = False
        col3 = raw.iloc[:, 3].astype(str).str.lower() if raw.shape[1] >= 4 else pd.Series([], dtype=str)
        if len(col3) > 0:
            like_chr = col3.str.match(r"^(chr)?([0-9]+|x|y|m|mt)$", na=False).mean()
            use_standard = like_chr > 0.5

        if use_standard:
            cols = [0, 1, 2, 3, 4, 5]
            logger.info("Detected standard BEDPE layout: [0,1,2,3,4,5]")
        else:
            if raw.shape[1] < 7:
                raise ValueError("BEDPE format ambiguous and <7 columns; cannot map alternative layout [0,1,2,4,5,6].")
            cols = [0, 1, 2, 4, 5, 6]
            logger.warning("Using alternative BEDPE layout: [0,1,2,4,5,6]. Please verify input format.")

        loops = raw.iloc[:, cols].copy()
        loops.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2"]

        loops["chr1"] = self._normalize_chr(loops["chr1"])
        loops["chr2"] = self._normalize_chr(loops["chr2"])

        for c in ["start1", "end1", "start2", "end2"]:
            loops[c] = pd.to_numeric(loops[c], errors="coerce")

        bad = loops[["start1", "end1", "start2", "end2"]].isna().any(axis=1).sum()
        if bad > 0:
            logger.warning("Dropping %d loop rows due to non-numeric coordinates", bad)
            loops = loops.dropna(subset=["start1", "end1", "start2", "end2"]).copy()

        loops[["start1", "end1", "start2", "end2"]] = loops[["start1", "end1", "start2", "end2"]].astype(int)

        invalid = ((loops["start1"] >= loops["end1"]) | (loops["start2"] >= loops["end2"])).sum()
        if invalid > 0:
            logger.warning("Dropping %d loop rows with invalid intervals (start >= end)", invalid)
            loops = loops[(loops["start1"] < loops["end1"]) & (loops["start2"] < loops["end2"])].copy()

        loops.reset_index(drop=True, inplace=True)
        logger.info("Loaded %d loops", len(loops))
        return loops

    def _load_enhancers(self) -> Dict[str, pd.DataFrame]:
        enh = pd.read_csv(
            self.enhancer_file,
            sep="\t",
            header=None,
            comment="#",
            usecols=[0, 1, 2, 3],
            names=["chr", "start", "end", "id"],
            dtype={"chr": str, "start": "Int64", "end": "Int64", "id": str},
            engine="c",
        )

        enh["chr"] = self._normalize_chr(enh["chr"])
        enh["start"] = pd.to_numeric(enh["start"], errors="coerce")
        enh["end"] = pd.to_numeric(enh["end"], errors="coerce")

        bad = enh[["start", "end"]].isna().any(axis=1).sum()
        if bad > 0:
            logger.warning("Dropping %d enhancer rows due to bad coordinates", bad)
            enh = enh.dropna(subset=["start", "end"]).copy()

        enh[["start", "end"]] = enh[["start", "end"]].astype(int)
        invalid = (enh["start"] >= enh["end"]).sum()
        if invalid > 0:
            logger.warning("Dropping %d enhancer rows with invalid intervals", invalid)
            enh = enh[enh["start"] < enh["end"]].copy()

        grouped = {c: g.sort_values("start").reset_index(drop=True) for c, g in enh.groupby("chr", sort=False)}
        logger.info("Loaded %d enhancers across %d chromosomes", len(enh), len(grouped))
        return grouped

    def _load_promoters(self) -> Dict[str, pd.DataFrame]:
        prom = pd.read_csv(
            self.promoter_file,
            sep="\t",
            header=None,
            comment="#",
            usecols=[0, 1, 2, 3, 4, 5, 6],
            names=["chr", "start", "end", "gene", "gene_chr", "gene_start", "gene_end"],
            dtype={
                "chr": str,
                "start": "Int64",
                "end": "Int64",
                "gene": str,
                "gene_chr": str,
                "gene_start": "Int64",
                "gene_end": "Int64",
            },
            engine="c",
        )

        prom["chr"] = self._normalize_chr(prom["chr"])
        prom["gene_chr"] = self._normalize_chr(prom["gene_chr"])
        for c in ["start", "end", "gene_start", "gene_end"]:
            prom[c] = pd.to_numeric(prom[c], errors="coerce")

        bad = prom[["start", "end", "gene_start", "gene_end"]].isna().any(axis=1).sum()
        if bad > 0:
            logger.warning("Dropping %d promoter rows due to bad coordinates", bad)
            prom = prom.dropna(subset=["start", "end", "gene_start", "gene_end"]).copy()

        prom[["start", "end", "gene_start", "gene_end"]] = prom[["start", "end", "gene_start", "gene_end"]].astype(int)
        invalid = (prom["start"] >= prom["end"]).sum()
        if invalid > 0:
            logger.warning("Dropping %d promoter rows with invalid intervals", invalid)
            prom = prom[prom["start"] < prom["end"]].copy()

        grouped = {c: g.sort_values("start").reset_index(drop=True) for c, g in prom.groupby("chr", sort=False)}
        logger.info("Loaded %d promoters across %d chromosomes", len(prom), len(grouped))
        return grouped

    def load_data(self) -> None:
        logger.info("Loading input files...")
        self.loops = self._load_loops_bedpe()
        self.enhancers_dict = self._load_enhancers()
        self.promoters_dict = self._load_promoters()

    @staticmethod
    def _overlap_filter(df: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
        return df[(df["start"] < end) & (df["end"] > start)]

    def _get_overlaps(self, chrom: str, start: int, end: int, ref_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        ref = ref_dict.get(chrom)
        if ref is None or ref.empty:
            return pd.DataFrame()
        return self._overlap_filter(ref, start, end)

    @staticmethod
    def _join_unique(values: pd.Series) -> str:
        vals = [str(v).strip() for v in values.dropna().tolist() if str(v).strip() != ""]
        return ",".join(sorted(set(vals))) if vals else ""

    def annotate(self) -> pd.DataFrame:
        if self.loops is None:
            self.load_data()

        assert self.loops is not None
        logger.info("Annotating loops (may take time for large datasets)...")

        total = len(self.loops)
        step = max(5000, total // 20) if total > 0 else 5000
        records = []

        for i, loop in self.loops.iterrows():
            if (i + 1) % step == 0:
                logger.info("Processing %d/%d loops...", i + 1, total)

            anchors = {
                "L": {"chr": loop["chr1"], "start": int(loop["start1"]), "end": int(loop["end1"])},
                "R": {"chr": loop["chr2"], "start": int(loop["start2"]), "end": int(loop["end2"])},
            }

            row = {}
            for side, a in anchors.items():
                p_hits = self._get_overlaps(a["chr"], a["start"], a["end"], self.promoters_dict)
                row[f"P_{side}"] = not p_hits.empty
                row[f"P_{side}_genes"] = self._join_unique(p_hits["gene"]) if not p_hits.empty else ""
                row[f"P_{side}_coords"] = self._join_unique(
                    p_hits.apply(lambda r: f"{r['gene_chr']}:{int(r['gene_start'])}-{int(r['gene_end'])}", axis=1)
                ) if not p_hits.empty else ""

                e_hits = self._get_overlaps(a["chr"], a["start"], a["end"], self.enhancers_dict)
                row[f"E_{side}"] = not e_hits.empty
                row[f"E_{side}_IDs"] = self._join_unique(e_hits["id"]) if not e_hits.empty else ""
                row[f"E_{side}_coords"] = self._join_unique(
                    e_hits.apply(lambda r: f"{r['chr']}:{int(r['start'])}-{int(r['end'])}", axis=1)
                ) if not e_hits.empty else ""

            records.append(row)

        self.annotated_loops = pd.concat([self.loops.reset_index(drop=True), pd.DataFrame.from_records(records)], axis=1)
        self._classify_loops()
        self._finalize_schema()
        return self.annotated_loops

    def _classify_loops(self) -> None:
        if self.annotated_loops is None:
            raise ValueError("annotated_loops is None. Run annotate() first.")

        df = self.annotated_loops
        for c in ["P_L", "P_R", "E_L", "E_R"]:
            df[c] = df[c].fillna(False).astype(bool)

        PL, PR = df["P_L"], df["P_R"]
        EL, ER = df["E_L"], df["E_R"]

        df["pure_PP"] = PL & PR & (~EL) & (~ER)
        df["pure_EE"] = EL & ER & (~PL) & (~PR)
        df["pure_EP"] = ((EL & ~PL) & (PR & ~ER)) | ((PL & ~EL) & (ER & ~PR))
        df["mixed_EP"] = (EL & PR & (PL | ER)) | (ER & PL & (PR | EL))

        for c in ["pure_PP", "pure_EE", "pure_EP", "mixed_EP"]:
            df[c] = df[c].astype(bool)

        df["intersection_type"] = "Other"

        m_pp = df["pure_PP"]
        m_ee = df["pure_EE"]
        m_ep = df["pure_EP"] | df["mixed_EP"]

        df.loc[m_pp, "intersection_type"] = "Promoter-Promoter"
        df.loc[m_ee, "intersection_type"] = "Enhancer-Enhancer"
        df.loc[m_ep, "intersection_type"] = "Enhancer-Promoter"

        unassigned = df["intersection_type"].eq("Other")


        U_L = ~PL & ~EL
        U_R = ~PR & ~ER
        m_pu = ((PL & U_R) | (PR & U_L)) & unassigned
        m_eu = ((EL & ~PL & U_R) | (ER & ~PR & U_L)) & unassigned
        m_uu = (U_L & U_R) & unassigned


        df.loc[m_pu, "intersection_type"] = "Promoter-Undefined"
        df.loc[m_eu, "intersection_type"] = "Enhancer-Undefined"
        df.loc[m_uu, "intersection_type"] = "Undefined-Undefined"

        self.annotated_loops = df

    def _finalize_schema(self) -> None:
        if self.annotated_loops is None:
            raise ValueError("annotated_loops is None. Run annotate() first.")

        df = self.annotated_loops

        for c in OUTPUT_COLUMNS:
            if c not in df.columns:
                df[c] = False if c in {"P_L", "P_R", "E_L", "E_R", "pure_PP", "pure_EE", "pure_EP", "mixed_EP"} else ""

        for c in ["P_L", "P_R", "E_L", "E_R", "pure_PP", "pure_EE", "pure_EP", "mixed_EP"]:
            df[c] = df[c].fillna(False).astype(bool)

        for c in ["P_L_coords", "P_R_coords", "E_L_coords", "E_R_coords", "P_L_genes", "P_R_genes", "E_L_IDs", "E_R_IDs"]:
            df[c] = df[c].fillna("").astype(str)

        df["intersection_type"] = df["intersection_type"].fillna("Other").astype(str)

        self.annotated_loops = df

    def get_statistics(self) -> Dict[str, int]:
        if self.annotated_loops is None:
            return {}
        df = self.annotated_loops

        stats = {
            "total_loops": int(len(df)),
            "pure_PP": int(df["pure_PP"].sum()),
            "pure_EE": int(df["pure_EE"].sum()),
            "pure_EP": int(df["pure_EP"].sum()),
            "mixed_EP": int(df["mixed_EP"].sum()),
        }

        type_counts = df["intersection_type"].value_counts(dropna=False).to_dict()
        for k, v in type_counts.items():
            stats[f"intersection_type::{k}"] = int(v)

        return stats

    def save_results(self) -> None:
        if self.annotated_loops is None:
            raise ValueError("No annotated loops to save. Run annotate() first.")

        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        out = self.annotated_loops.copy()
        for c in OUTPUT_COLUMNS:
            if c not in out.columns:
                out[c] = ""
        out[OUTPUT_COLUMNS].to_csv(self.output_file, sep=",", index=False, na_rep="")
        logger.info("Results saved to: %s", self.output_file)

    def get_statistics(self) -> Dict[str, int]:
        if self.annotated_loops is None:
            return {}
        df = self.annotated_loops
        return {
            "total_loops": int(len(df)),
            "pure_PP": int(df["pure_PP"].sum()),
            "pure_EE": int(df["pure_EE"].sum()),
            "pure_EP": int(df["pure_EP"].sum()),
            "mixed_EP": int(df["mixed_EP"].sum()),
        }


def run_annotate(args) -> Dict[str, int]:
    """
    Expected argparse fields from cli.py:
    - args.loops
    - args.enhancers
    - args.promoters
    - args.output
    """
    annotator = LoopAnnotator(
        loop_file=args.loops,
        enhancer_file=args.enhancers,
        promoter_file=args.promoters,
        output_file=args.output,
    )
    annotator.annotate()
    annotator.save_results()
    stats = annotator.get_statistics()

    logger.info("=== Annotation Statistics ===")
    for k, v in stats.items():
        logger.info("%s: %s", k, v)

    return stats
