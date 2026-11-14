"""Public row-1 readout: native observer-only admission precedes decoding.

The raw count helper alone is not an admission interface.
"""
from pathlib import Path
import sys
import topology
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(ROOT/'reverse-engineering-reality/code/openwave_row1_integrity'),
              str(ROOT/'openwave-tuning/next/CHAIN3/row1_integrity_audit')]
import integrity
from observer_only import admit_record_only


def admit_and_decode(claim,G,region,retries):
    verdict=admit_record_only(claim,G,region,retries,producer=integrity)
    if not verdict['admitted']:return verdict
    degree=topology.raw_pointer_degree(claim['meter'],integrity.meter.CYCLE)
    charge=None if degree is None else -degree
    if charge!=verdict['readout']:
        return {'admitted':False,'reasons':['Validated native pointer and raw spatial decode disagree'],
                'compatible_pointer_initials':[]}
    return {**verdict,'raw_spatial_pointer_degree':degree}
