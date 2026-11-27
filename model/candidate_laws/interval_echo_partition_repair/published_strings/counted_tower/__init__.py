"""Counted small-ball entropy contrasts on one published-string survey."""
from pathlib import Path
import sys

# The modules here share the ten-frame convention of the matched Einstein theorem under code/.
# The model CLI and the model test job do not set PYTHONPATH.
CODE = Path(__file__).resolve().parents[5] / "code"
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
