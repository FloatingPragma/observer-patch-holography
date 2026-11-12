"""Exact outcome-conditioned compression of the complete lepton archive."""
import os
for _key in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ.setdefault(_key,'1')
