#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys
raise SystemExit(subprocess.call([sys.executable,str(Path(__file__).with_name("check_boost.py"))]))
