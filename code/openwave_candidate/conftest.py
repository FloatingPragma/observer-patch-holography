"""Expose the portable package to repository-wide test collection."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "candidate"))
