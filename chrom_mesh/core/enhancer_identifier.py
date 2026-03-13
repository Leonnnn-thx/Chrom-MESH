#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Chrom-MESH Enhancer Identifier Wrapper
======================================

Python wrapper for calling the shell pipeline:
    chrom_mesh/scripts/Identify_AEs.sh
"""

from __future__ import annotations

import os
import subprocess
import logging
from pathlib import Path
from typing import List, Optional
import importlib.resources as pkg_resources


logger = logging.getLogger(__name__)


def _append_if_value(cmd: List[str], flag: str, value: Optional[str]) -> None:
    if value is not None and str(value).strip() != "":
        cmd.extend([flag, str(value)])


def _resolve_script_path() -> Path:
    try:
        script_ref = pkg_resources.files("chrom_mesh.scripts").joinpath("Identify_AEs.sh")
        with pkg_resources.as_file(script_ref) as script_path:
            return Path(script_path)
    except Exception:
        here = Path(__file__).resolve()
        fallback = here.parent.parent / "scripts" / "Identify_AEs.sh"
        return fallback


def run_identify(args) -> None:
    """
    Required args (from cli.py):
      - h3k27ac_bed_dir
      - h3k27ac_bw_dir
      - atac_bed_dir
      - atac_bw_dir
      - output

    Optional args:
      - bigwig_tool
      - bigwig_info_tool
      - bedtools
      - python_exec
      - h3k27ac_pattern
      - atac_pattern
      - h3k27ac_bw_pattern
      - atac_bw_pattern
    """
    script_path = _resolve_script_path()

    if not script_path.exists():
        raise FileNotFoundError(f"Identify_AEs.sh not found: {script_path}")

    # Ensure executable bit is present (best effort).
    try:
        current_mode = script_path.stat().st_mode
        script_path.chmod(current_mode | 0o111)
    except Exception:
        pass

    cmd: List[str] = ["bash", str(script_path)]

    # Required mappings
    cmd.extend(["-h", str(args.h3k27ac_bed_dir)])
    cmd.extend(["-H", str(args.h3k27ac_bw_dir)])
    cmd.extend(["-a", str(args.atac_bed_dir)])
    cmd.extend(["-A", str(args.atac_bw_dir)])
    cmd.extend(["-o", str(args.output)])

    # Optional mappings to shell script
    _append_if_value(cmd, "-b", getattr(args, "bigwig_tool", None))
    _append_if_value(cmd, "--bigwig-info-tool", getattr(args, "bigwig_info_tool", None))
    _append_if_value(cmd, "--bedtools", getattr(args, "bedtools", None))
    _append_if_value(cmd, "--python", getattr(args, "python_exec", None))

    _append_if_value(cmd, "-p", getattr(args, "h3k27ac_pattern", None))
    _append_if_value(cmd, "-n", getattr(args, "atac_pattern", None))
    _append_if_value(cmd, "-w", getattr(args, "h3k27ac_bw_pattern", None))
    _append_if_value(cmd, "-W", getattr(args, "atac_bw_pattern", None))

    # Optional legacy args are ignored intentionally:
    # bam, workdir, min_width, max_width

    logger.info("Running enhancer identification pipeline...")
    logger.debug("Command: %s", " ".join(cmd))

    env = os.environ.copy()
    env.setdefault("LC_ALL", "C")
    env.setdefault("LANG", "C")

    try:
        subprocess.run(cmd, check=True, env=env)
    except subprocess.CalledProcessError as e:
        logger.error("Enhancer identification failed with exit code %s", e.returncode)
        raise SystemExit(e.returncode) from e

    logger.info("Enhancer identification completed successfully.")
