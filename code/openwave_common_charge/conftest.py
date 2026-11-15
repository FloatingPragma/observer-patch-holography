"""Collection policy for the sealed OpenWave snapshot.

``workspace/`` holds immutable input bytes that ``replay.py`` checks against
SNAPSHOT_MANIFEST.json as an exact file set. They are not modules of the
current Python package: collecting them would import frozen copies of current
tests under clashing module names and write unmanifested ``__pycache__`` files
into the snapshot. ``test_custody.py`` and ``replay.py`` verify and execute it.
"""

collect_ignore = ["workspace"]
