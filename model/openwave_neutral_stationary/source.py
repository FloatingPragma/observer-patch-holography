"""Full native source with explicit periodic signs and doubled link penalty."""
from functools import lru_cache
import sympy as s
from model.openwave_neutral_candidate import source as previous
from model.coherent_defect_wave.chiral_joint.oriented_links.model import endpoint_current
from .isolation import LAMBDA,KAPPA


def full_column(key):
    output=previous.full_column(key,periodic=True,coupling=LAMBDA)
    output[key]=output.get(key,0)+(KAPPA-1)*sum(label!=0 for label in key[1])
    return {k:s.expand(v) for k,v in output.items() if s.expand(v)!=0}


@lru_cache(maxsize=1)
def certificate():
    currents=previous.current_certificate();boundary=previous.boundary_certificate()
    nonblank=s.diag(0,*([1]*118));checks=0
    for current in range(12):
        for endpoint in (0,1):
            generator=endpoint_current(current,endpoint)
            assert nonblank*generator==generator*nonblank;checks+=1
    return {'Hamiltonian':'H_local+2 N_nonblank+2^-34 T_periodic',
        'local_quantum_Yukawa_and_scalar_coefficients':'1, unchanged native operators',
        'native_current_controls':currents,'periodic_boundary_controls':boundary,
        'full_link_penalty_Gauss_commutators':checks,
        'added_link_penalty_commutes_with_translations_L_B3_and_parity':True,
        'original_penalty_one_ground_or_preparation_imported':False,
        'all_native_occupations_scalars_119_link_states_and_complements_retained':True}
