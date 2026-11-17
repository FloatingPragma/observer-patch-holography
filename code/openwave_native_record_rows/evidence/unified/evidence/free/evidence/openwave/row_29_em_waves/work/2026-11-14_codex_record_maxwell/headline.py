#!/usr/bin/env python3
"""Regenerate the finite native-record Maxwell component diagnostics."""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).resolve().parent / "experiment/probe.py"),
               run_name="__main__")
