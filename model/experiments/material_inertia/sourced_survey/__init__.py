"""Expose the survey's shared station instrument from the repository source tree."""

from pathlib import Path
import sys

# The clocks are a Python package under code/, alongside the top-level model/.
code_root = str(Path(__file__).resolve().parents[4] / "code")
if code_root not in sys.path:
    sys.path.append(code_root)
