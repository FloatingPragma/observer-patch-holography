"""Import the portable candidate during RER-wide pytest collection.

The repository's pytest configuration takes precedence over the candidate's
nested pyproject.toml. Keep the authenticated portable files unchanged and
make their checked-in package available to the enclosing test harness.
"""
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parent / "candidate"))
