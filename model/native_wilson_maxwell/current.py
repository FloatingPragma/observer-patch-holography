"""Exact Dirac current, scalar-band obstruction and Wilson carrier controls."""
from pathlib import Path
import sys
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT))
from model.native_block_dirac.algebra import spinor_frames,pauli


def matrices():
    alpha=tuple(s.BlockMatrix([[s.zeros(2),v],[v,s.zeros(2)]]).as_explicit() for v in pauli())
    beta=s.diag(1,1,-1,-1)
    for a in alpha:
        assert a*a==s.eye(4) and a*beta+beta*a==s.zeros(4)
    return alpha,beta


def transverse_countercontrol():
    alpha,_=matrices()
    k=s.Rational(3,4)
    omega,plus,_,_=spinor_frames((0,0,k),s.Integer(1))
    _,minus,_,_=spinor_frames((0,0,-k),s.Integer(1))
    density=(plus.H*minus).applyfunc(s.simplify)
    current=(plus.H*alpha[0]*minus).applyfunc(s.simplify)
    longitudinal=(plus.H*alpha[2]*minus).applyfunc(s.simplify)
    assert omega==s.Rational(5,4)
    assert density==s.Rational(4,5)*s.eye(2)
    assert current==s.Matrix([[0,s.Rational(3,5)],[-s.Rational(3,5),0]])
    assert longitudinal==s.zeros(2)
    # Frechet derivative of sqrt(m^2+P^2) in the transverse plane wave.
    p=s.Matrix([0,0,-k]);pout=-p;epsilon=s.Matrix([1,0,0])
    scalar_vertex=(p+pout).dot(epsilon)/(2*omega)
    assert scalar_vertex==0
    return {'mass_unit':'actual selected positive native excitation m',
        'incoming_momentum_over_m':['0','0','-3/4'],
        'outgoing_momentum_over_m':['0','0','3/4'],
        'transverse_polarization':['1','0','0'],
        'omega_over_m':'5/4',
        'Dirac_density_matrix':[[str(x) for x in row] for row in density.tolist()],
        'Dirac_transverse_current_matrix':[[str(x) for x in row] for row in current.tolist()],
        'Dirac_longitudinal_current_matrix':[[str(x) for x in row] for row in longitudinal.tolist()],
        'scalar_square_root_minimal_coupling_vertex':'0',
        'canonical_particle_number_density_matrix':'I2',
        'nonzero_transverse_spin_flip_amplitude':'3/5',
        'charge_or_free_dispersion_equality_fixes_current':False,
        'countercontrol_persists_at_fixed_physical_momentum_in_refinement':True}


def wilson_controls():
    alpha,beta=matrices();z=s.symbols('z0:3',real=True)
    mass,spacing=s.symbols('m a',positive=True)
    momenta=[s.sin(x)/spacing for x in z]
    M=mass+sum(1-s.cos(x) for x in z)/spacing
    h=M*beta+sum((p*A for p,A in zip(momenta,alpha)),s.zeros(4))
    omega2=M*M+sum(p*p for p in momenta)
    assert (h*h-omega2*s.eye(4)).applyfunc(s.simplify)==s.zeros(4)
    bond=-(beta+s.I*alpha[0])/2
    singular_square=bond.H*bond
    assert singular_square*singular_square==singular_square and s.trace(singular_square)==2
    complete=s.BlockMatrix([[s.zeros(4),bond],[bond.H,s.zeros(4)]]).as_explicit()
    assert complete**3==complete and s.trace(complete**2)==4 and s.trace(complete)==0
    # Off-diagonal occupation effects in the canonical local charge frame.
    # They are absent in an onsite particle-plus/antiparticle-minus density.
    _,U,V,_=spinor_frames((0,0,s.Rational(3,4)),s.Integer(1))
    _,U0,V0,_=spinor_frames((0,0,0),s.Integer(1))
    pair_density=(U.H*V0).applyfunc(s.simplify)
    pair_current=(U0.H*alpha[0]*V0).applyfunc(s.simplify)
    assert pair_density!=s.zeros(2)
    assert pair_current==pauli()[0]
    return {'Wilson_r':'1','mass_parameter':'m>0, selected actual native G1 source excitation',
        'one_body_symbol':'sum_j alpha_j sin(a p_j)/a + beta*(m+sum_j(1-cos(a p_j))/a)',
        'full_four_by_four_square_is_Omega_squared_identity':True,
        'free_positive_energy_lower':'Omega_W(p)>=m',
        'complete_one_bond_one_body_spectrum':{'-1/a':2,'0':4,'+1/a':2},
        'full_Fock_bond_norm':'2/a, every occupation retained',
        'full_finite_range_weighted_matter_interaction_bound':'24 lambda_spatial/a on a cubic lattice',
        'canonical_local_charge':'rho_x=sum_(c=1..4) Psi_xc^dagger Psi_xc - 2',
        'canonical_local_charge_spectrum':[-2,-1,0,1,2],
        'total_charge':'sum_x rho_x=N_plus-N_minus',
        'local_density_equals_original_native_onsite_N_plus_minus_N_minus':False,
        'pair_density_nonzero_control':[[str(x) for x in row] for row in pair_density.tolist()],
        'zero_momentum_pair_current_matrix':[[str(x) for x in row] for row in pair_current.tolist()],
        'interacting_current_retains_pair_creation_and_annihilation':True,
        'locality_weight_requires_actual_local_tensor_link_operators':True,
        'transverse_total_occupation_cutoff_automatically_satisfies_that_locality_domain':False}


def calculate():
    return {'scalar_band_current_countercontrol':transverse_countercontrol(),
            'candidate_Wilson_current_interface':wilson_controls()}


if __name__=='__main__':
    import json
    print(json.dumps(calculate(),indent=2))
