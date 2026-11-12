"""Bounded source words for complete prepared lepton cycles."""
import os

for _name in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ.setdefault(_name,'1')
