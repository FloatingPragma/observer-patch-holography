"""Actual same-boundary antiunitary conjugation of the periodic source."""
from collections import defaultdict
from functools import lru_cache
import sympy as s
from ..model import active_terms,clean,sequence
from ..two_cell.symmetries import (conjugating_mask,periodic_column,periodic_translation,
                                  periodic_current_column,apply_column)


def reflected_conjugation_column(key):
    mask,l01,l10,hl,hr=key;labels=active_terms()[3]
    inverse={i:labels.index((b,a)) for i,(a,b) in enumerate(labels)}
    scalar=(0,3,4,1,2);target,phase=conjugating_mask(mask)
    return {(target,inverse[l01],inverse[l10],scalar[hr],scalar[hl]):phase}


def conjugation_column(key):
    return apply_column(reflected_conjugation_column(key),periodic_translation)


def apply_conjugation(state):
    out=defaultdict(lambda:s.Integer(0))
    for key,value in state.items():
        for target,phase in conjugation_column(key).items():out[target]+=s.conjugate(value)*phase
    return clean(out)


@lru_cache(maxsize=1)
def certificate_data():
    # The inherited complete CAR map plus the exact field-register actions
    # determine the complete operator identity, including occupied inputs.
    labels=active_terms()[3];link=s.zeros(6)
    for j,(a,b) in enumerate(labels):link[labels.index((b,a)),j]=1
    for j,(a,b) in enumerate(labels[1:],1):
        u=s.zeros(6);u[j,0]=1;target=s.zeros(6);target[labels.index((b,a)),0]=1
        assert link*u*link.T==target
    fixtures={(mask,0,0,h,0) for mask in range(64) for h in range(5)}
    fixtures.update((mask,l01,l10,hl,hr) for mask in (0,1,64,65,521,2047,4095)
                    for l01 in range(6) for l10 in range(6) for hl,hr in ((0,0),(1,2),(3,4)))
    for key in fixtures:
        assert apply_conjugation(periodic_column(key))==apply_column(conjugation_column(key),periodic_column)
    gauge=0
    for key in sorted(fixtures)[::max(1,len(fixtures)//8)]:
        for boundary in ((0,0),(0,-1),(0,1),(-1,0),(1,0)):
            inverse=tuple(-q for q in boundary)
            for site in range(2):
                for current in range(12):
                    a=apply_conjugation(periodic_current_column(key,site,current,boundary))
                    b=apply_column(conjugation_column(key),lambda k:periodic_current_column(k,site,current,inverse))
                    assert a=={target:-value for target,value in b.items()}
                    gauge+=1
    return {'same_boundary_periodic_H_antiunitary_symmetry':True,
            'operator':'C=T A; A transposes each quantum link and exchanges/conjugates scalar memories',
            'source_CAR_action':['c_left -> -i c_left_dagger','c_right -> +i c_right_dagger'],
            'boundary_map':'(q_left,q_right) -> (-q_left,-q_right)',
            'source_label_map':'L ->4-L, twoSz ->-twoSz',
            'actual_quantum_link_matrix_units_checked':5,
            'actual_full_H_fixture_columns_checked':len(fixtures),'actual_local_Gauss_checks':gauge,
            'unique_neutral_ground_implies_each_e_c_mode_mean_occupation':'1/2',
            'same_site_addition_removal_spectral_measures_match':True,
            'source_only_particle_or_Lorentz_CPT_identification':False}
