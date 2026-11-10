from pathlib import Path
import sys,json,hashlib
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[2]
import numpy as np
from scipy import sparse
from model.cube_eta_protection import cube,repair,source
from model.field_stress.energy import GaussianInteger as GI,commutator
from model.field_symmetry import symmetry
from model.protected_particle_readout import bounds,circuits
from functools import lru_cache

@lru_cache(None)
def certificate():
    ops=cube.operators(); eta=ops.eta_plus
    terms=tuple(term for group in cube.GEOMETRY.groups for term in group)
    assert len(terms)==24 and len(ops.bonds)==12
    pairs=[]
    for n in range(0,len(terms),2):
     left,right=terms[n:n+2]
     assert left[0]//2==right[0]//2 and left[1]//2==right[1]//2
     assert len({left[0],left[1],right[0],right[1]})==4
     # Each ordered adjacent pair consists of commuting, disjoint even CAR terms.
     matrix=np.zeros((16,16),complex)
     for p,q,z in (left,right):matrix[p,q]+=16*z;matrix[q,p]+=16*z.conjugate()
     bond=symmetry._lift(matrix,16)
     assert commutator(bond,eta).is_zero()
     pairs.append({'terms':[[p,q,[int((16*z).real),int((16*z).imag)]] for p,q,z in (left,right)],'exact_commutes_with_eta_plus':True})
    D=np.array([sum((mask>>(2*j))&3==3 for j in range(8)) for mask in range(65536)],dtype=np.int64)
    double=GI.from_exact(sparse.diags(D,dtype=np.int64))
    assert (commutator(double,eta)-eta).is_zero()
    vac=np.zeros(65536,dtype=np.int64);vac[0]=1
    # Local pair creation carries canonical CAR sign; its origin support is0,1.
    a0=symmetry._annihilator(16,0);a1=symmetry._annihilator(16,1)
    q=(a0.adjoint()@a1.adjoint()).real
    A=eta.real-q
    parent_R=eta.real@(eta.real@vac)
    parent_FR=A@(A@vac)-4*A@(q@vac)
    R=np.zeros(65536,dtype=np.int64)
    for mask,value in enumerate(source.reference_values()):R[sum(3<<(2*j) for j in range(8) if mask>>j&1)]=value
    FR=R*np.array([(-1)**((mask&3)==3) for mask in range(65536)],dtype=np.int64)
    assert np.array_equal(eta.real@parent_R,6*R)
    assert np.array_equal(eta.real@parent_FR,6*FR)
    a,b=repair.exact_embedding()
    assert (eta.real.T@a).nnz==0 and (eta.real.T@b).nnz==0
    rc=repair.certificate();bc=bounds.source_certificate()
    e=F(bc['preparation_complete_state_error'])+F(bc['pointer_preparation_error'])+2*F(bc['one_field_segment_source_error'])
    retained=json.loads((ROOT/'model/protected_readout_numerics/retained.json').read_text())
    p=retained['executed_native_readout']['terminal_native_repair_probability']
    value={'schema':'oph.protected_readout_kernel_audit.v1','exact_template_bond_pairs':pairs,
     'all_repair_columns_killed_by_eta_minus':True,'all1820_eta_image_columns_annihilated':rc['all_native_outcomes_annihilate_protected_image'],
     'preparation_integer_identities':['eta_plus parent_R=6 R_integer','eta_plus parent_FR=6 FR_integer'],
     'contact_identity':'[D,eta_plus]=eta_plus',
     'exact_target_Strang_template_preserves_protected_eta_image':True,
     'exact_target_template_all_terminal_native_probabilities':0,
     'literal_finite_control_word_terminal_probability_is_forced_zero':False,
     'literal_complete_source_error_without_Strang_term':str(e),
     'literal_all_terminal_native_probability_upper':str(e*e),
     'literal_all_terminal_native_probability_upper_decimal':float(e*e),
     'retained_binary64_terminal_probability':p,
     'backend_arithmetic_error_included_in_literal_control_certificate':False,
     'retained_scientific_numbers_or_original_comparison_modified':False,
     'sources':{str(x):hashlib.sha256((ROOT/x).read_bytes()).hexdigest() for x in [Path('model/protected_readout_numerics/retained.json'),Path('model/protected_particle_readout/packet.py'),Path('model/protected_particle_readout/execute.py'),Path('model/protected_particle_readout/experiment.py'),Path('model/protected_particle_readout/circuits.py'),Path('model/protected_particle_readout/bounds.py'),Path('model/cube_eta_protection/repair.py'),Path('model/cube_eta_protection/source.py'),Path('model/observer/source_scalar/compiler.py')]}}
    return value
