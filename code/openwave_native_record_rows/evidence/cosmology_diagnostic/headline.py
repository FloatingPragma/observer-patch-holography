#!/usr/bin/env python3
import runpy
from pathlib import Path
import sys
sys.dont_write_bytecode = True
shared = Path(__file__).resolve().parents[3] / 'shared/operational_gravity_20261114'
sys.path.insert(0, str(shared))
runpy.run_path(str(shared/'headline.py'), run_name='__main__')
