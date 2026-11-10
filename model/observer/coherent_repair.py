"""Named coherent observer, separate from the committed native record law.

Replay with ``python -m model.observer.coherent_repair --check``.
"""
from model.coherent_seam_repair import PROGRAM_ID, Observer, State, Session
from model.coherent_seam_repair.__main__ import main

__all__ = ['PROGRAM_ID', 'Observer', 'State', 'Session']

if __name__ == '__main__':
    main()
