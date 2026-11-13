from fractions import Fraction as F
from math import comb
import json
import numpy as np
from scipy.linalg import eigh
from .certificate import calculate,singlet_dimension
from .complement import local_rows,three_cell_singlet_dimension
from .source import projected_source
from .exact import certified_levels
from .chain_bound import fusion


def test_complete_exterior_and_fusion_dimensions():
    rows=local_rows()
    assert len(rows)==97
    for n in range(13):
        assert sum((p+1)*(q+1)*(p+q+2)//2*(s+1)*(i+1)*m
                   for N,p,q,s,i,m in rows if N==n)==comb(12,n)
    assert three_cell_singlet_dimension()==singlet_dimension(3)==572572
    assert singlet_dimension(4)==176729280
    assert dict(fusion((1,0),(0,1)))=={(0,0):1,(1,1):1}
    assert dict(fusion((1,0),(1,0)))=={(0,1):1,(2,0):1}


def test_exact_all_multiplets_match_full_compressed_operator():
    B=F(calculate()['three_body_coefficient'])
    for cells in(3,4):
        source=projected_source(cells,B,True)
        reference_Iz=F(-1,2)if cells%2 else F(0)
        predicted=[]
        for row in certified_levels(cells,B):
            spin,isospin=F(row['spin']),F(row['isospin'])
            for lo,hi in row['source_energy_intervals']:
                for twice_iz in range(-int(2*isospin),int(2*isospin)+1,2):
                    from .source import PARAMETERS
                    shift=-PARAMETERS['delta']*(F(twice_iz,2)-reference_Iz)
                    predicted.extend([float((lo+hi)/2+shift)]*int(2*spin+1))
        actual=eigh(source['H'],eigvals_only=True)
        assert len(predicted)==source['dimension']
        assert np.max(abs(np.sort(predicted)-actual))<1e-8


def test_full_source_intervals_failure_records_and_heldouts():
    result=calculate()
    masses=result['mass_comparisons']
    lo,hi=map(F,masses['triton']['source_GeV_interval'])
    assert lo<F(masses['triton']['target_GeV'])<hi
    for name in('helion','alpha'):
        lo,hi=map(F,masses[name]['source_GeV_interval'])
        assert not lo<=F(masses[name]['target_GeV'])<=hi
    for graph in result['complete_graphs']:
        records=graph['complete_original_preparation_and_spectral_records']
        assert F(records['preparation_success'])+F(records['preparation_failure'])==1
        assert F(graph['complement_floor_at_reference_Iz'])>F(graph['P_upper'])
        assert F(graph['uniform_low_energy_error_GeV'])<F(1,50000)
        assert sum(r['rank']*(2*F(r['spin'])+1)*(2*F(r['isospin'])+1)
                   for r in graph['low_multiplets'])==graph['complete_low_magnetic_dimension']
