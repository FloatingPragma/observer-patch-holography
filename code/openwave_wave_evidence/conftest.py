"""Keep sealed historical programs out of repository-wide pytest collection.

Their original test modules have repeated filenames and private-Git fixtures.
The documented reproduce.py tests command runs all 21 amended checks against
the bundled source, with the historical bytes preserved.
"""

collect_ignore = ["evidence"]
