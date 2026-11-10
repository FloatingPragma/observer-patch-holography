"""Reproducible finite free-field checks and source custody."""
from itertools import product
import hashlib
from pathlib import Path
import sympy as s

from model.particle_common_contract.provenance import import_closure
from ..spatial.dirac_route import clifford,symbol as floquet_symbol
from .spectrum import (PROGRAM_ID,source_mass,symbol,energy,projectors,
                       continuum_bound,corners)
from .fock import act,exterior_matrix,rest_instrument
from .source import certificate as source_certificate
from .position import certificate as position_certificate

ROOT=Path(__file__).resolve().parents[3]


def car_checks():
    size=20;masks=(0,(1<<size)-1,0x55555,0xAAAAA,1,0x391AB)
    checks=0
    def composed(mask,first,second):
        one=act(mask,*first,size)
        if one is None:return {}
        two=act(one[0],*second,size)
        return {} if two is None else {two[0]:one[1]*two[1]}
    def add(left,right):
        return {key:left.get(key,0)+right.get(key,0) for key in left.keys()|right.keys()
                if left.get(key,0)+right.get(key,0)}
    for mask,i,j in product(masks,range(size),range(size)):
        mixed=add(composed(mask,(j,True),(i,False)),composed(mask,(i,False),(j,True)))
        assert mixed==({mask:1} if i==j else {})
        same=add(composed(mask,(j,True),(i,True)),composed(mask,(i,True),(j,True)))
        assert same=={}
        checks+=2
    return checks


def algebra_certificate():
    alpha,beta=clifford();p=s.symbols('p0:3',real=True);M=s.symbols('M',real=True)
    h=sum((p[i]*alpha[i] for i in range(3)),s.zeros(4))+M*beta
    assert (h*h-(sum(x*x for x in p)+M*M)*s.eye(4)).applyfunc(s.expand)==s.zeros(4)
    conjugator=alpha[1]*beta
    assert conjugator*conjugator.conjugate()==s.eye(4)
    assert all(conjugator*x.conjugate()*conjugator.H==x for x in alpha)
    assert conjugator*beta.conjugate()*conjugator.H==-beta
    momenta=((0,0,0),(s.pi,0,0),(s.pi,s.pi/2,-s.pi),(2*s.pi,0,0))
    for p in momenta:
        a=s.Rational(1,2);matrix=symbol(a,p);e=energy(a,p);plus,minus=projectors(a,p)
        clean=lambda z:z.applyfunc(s.simplify)
        assert matrix.H==matrix and clean(matrix*matrix-e*e*s.eye(4))==s.zeros(4)
        assert clean(plus*plus-plus)==clean(minus*minus-minus)==s.zeros(4)
        assert clean(plus*minus)==s.zeros(4) and clean(plus+minus)==s.eye(4)
        assert s.simplify(s.trace(plus))==s.simplify(s.trace(minus))==2
        assert clean(conjugator*matrix.conjugate()*conjugator.H+symbol(a,tuple(-x for x in p)))==s.zeros(4)
    # Every degree of the four-mode exterior lift is complete and unitary.
    a=s.pi/12;k=(s.pi/2,0,s.pi);u=floquet_symbol(k,2*a)
    dimensions=[]
    for number in range(5):
        lifted=exterior_matrix(u,number)
        assert (lifted.H*lifted-s.eye(lifted.rows)).applyfunc(s.simplify)==s.zeros(lifted.rows)
        dimensions.append(lifted.rows)
    assert sum(dimensions)==16
    low=[];edge=[]
    mass_coin=s.cos(2*a)*s.eye(4)-s.I*s.sin(2*a)*beta
    for bits in product((0,1),repeat=3):
        assert (floquet_symbol(tuple(s.pi*b for b in bits),2*a)-(-1)**sum(bits)*mass_coin).applyfunc(s.simplify)==s.zeros(4)
        (low if sum(bits)%2==0 else edge).append(list(bits))
    return {'exact_symbol_projector_charge_conjugation_points':len(momenta),
            'Clifford_mass_shell_identity_exact':True,
            'charge_conjugation_matrix':[[str(v) for v in row] for row in conjugator.tolist()],
            'particle_hole_charges':[-1,1],
            'charge_is_source_number_relative_to_filled_sea':True,
            'all_four_mode_exterior_lift_dimensions':dimensions,
            'Floquet_even_corner_small_phase_branches':low,
            'Floquet_odd_corner_zone_boundary_branches':edge,
            'Wilson_clock_is_Floquet_logarithm':False,
            'twenty_mode_CAR_column_checks':car_checks()}


def continuum_certificate():
    import numpy as np
    from scipy.linalg import expm
    alpha,beta=clifford();p=(s.Rational(1,5),s.Rational(-1,4),s.Rational(1,3));window=s.Rational(1,3)
    hd=2*beta+sum((p[i]*alpha[i] for i in range(3)),s.zeros(4))
    hd=np.array(hd,dtype=complex);ed=float(s.sqrt(sum(x*x for x in p)+4))
    alphas=[np.array(matrix,dtype=complex) for matrix in alpha];bn=np.array(beta,dtype=complex)
    pm=(np.eye(4)-hd/ed)/2;rows=[]
    for a in (s.Rational(1,4),s.Rational(1,8),s.Rational(1,16)):
        # Numerical controls evaluate the proved symbol directly. Symbolic
        # simplification of transcendental test momenta supplies no extra proof.
        av=float(a);pv=np.array([float(x) for x in p]);k=av*pv
        effective_mass=2+float(np.sum(1-np.cos(k)))/av
        ha=sum((np.sin(k[i])*alphas[i]/av for i in range(3)),np.zeros((4,4),complex))+effective_mass*bn
        ea=float(np.sqrt(np.sum(np.sin(k)**2)/av**2+effective_mass**2))
        pa=(np.eye(4)-ha/ea)/2
        bound=continuum_bound(a,window)
        errors=(float(np.linalg.norm(ha-hd,2)),float(np.linalg.norm(expm(-1j*ha)-expm(-1j*hd),2)),
                float(np.linalg.norm(pa-pm,2)))
        assert errors[0]<float(bound) and errors[1]<float(bound) and errors[2]<float(bound/2)
        rows.append({'spacing':str(a),'symbol_and_unit_time_bound':str(bound),
                     'ground_projector_bound':str(bound/2),'floating_diagnostic_errors':[round(x,12) for x in errors]})
    return {'exact_window_symbol_bound':'b(a,P)=a^2*P^3/2+3*a*P^2/2',
            'exact_window_energy_bound':'b(a,P)',
            'exact_window_dynamics_bound':'abs(t)*b(a,P)',
            'exact_window_ground_covariance_bound':'b(a,P)/2',
            'normal_ordered_2n_field_Wick_error':'n*n!*max_two_point_error',
            'limiting_mass_shell':'E^2=|p|^2+4',
            'one_finite_energy_momentum_valley_per_face_channel':True,
            'numerical_controls_do_not_replace_analytic_bounds':rows,
            'Wilson_corners_at_spacing_one_half':corners(s.Rational(1,2))}


def hashes():
    here=Path(__file__).parent
    paths={str(p.relative_to(ROOT)) for p in here.iterdir()
           if p.suffix in ('.py','.md') and not p.name.startswith('test_')}
    paths.update(('model/requirements.txt','model/coherent_defect_wave/mobile_flux/PROOF.md',
                  'model/coherent_defect_wave/mobile_flux/ROUTING_COMPLETION.md',
                  'model/coherent_defect_wave/spatial/ROUTE.md',
                  'model/continuum_free_limit/README.md',
                  'model/coherent_defect_wave/physical_attachment/PROOF.md',
                  'evidence/native_geometric_source_20260925/spec.json',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(paths)) if (ROOT/name).is_file()}


def certificate():
    return {'schema':'oph.coherent_defect_wave.many_body_route.v1','program':PROGRAM_ID,
            'source':source_certificate(),'algebra':algebra_certificate(),
            'literal_spatial_action':position_certificate(),
            'complete_rest_instrument':rest_instrument(),'continuum':continuum_certificate(),
            'reference':{'state':'pure filled negative Wilson bands',
                         'finite_cube_sea_occupation':'10*V','finite_cube_sea_energy':'-10*sum_p E_a(p)',
                         'unique_positive_excitation_ground':True,
                         'full_Fock_gap':'2','relative_charge_zero_gap':'4'},
            'scope':{'complete_declared_graded_free_source':True,
                     'source_encoded_particle_hole_readouts':True,
                     'controlled_massive_free_Dirac_field_and_vacuum_limit':True,
                     'CAR_statistics_are_declared_architecture':True,
                     'native_gauge_currents_nontrivial_on_empty_native_link_code':False,
                     'five_faces_identified_as_five_measured_species':False,
                     'bosonic_collision_law_or_native_interaction_automatically_transferred':False,
                     'interacting_relativistic_limit_or_SM_or_GeV':False},
            'input_sha256':hashes()}
