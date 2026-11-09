"""An explicit law perturbation preserving complete conditional pixel closure."""
from fractions import Fraction as F
import hashlib
import json

from model.charged_lepton_response import p_dependence as shape
from model.native_flavor_orbit.intervals import I, sqrt_rational, pi_interval
from model.mcpr_native_memory_bridge import certificate as elementary
from model.pixel_mass_closure import packet as owner, memory

ROOT=shape.ROOT
SHIFT=F(1,10**14)
SOURCE_INPUTS=(
 'model/pixel_mass_closure/results.json','model/pixel_mass_closure/PROOF.md',
 'model/pixel_mass_closure/closure.py','model/pixel_mass_closure/memory.py',
 'model/charged_lepton_response/p_dependence.py',
 'code/particles/leptons/derive_charged_face_incidence_conditional_theorem.py')


def box(row):
    if type(row)is not dict or set(row)!= {'lo','hi'}:raise ValueError('owning exact interval required')
    lo,hi=F(row['lo']),F(row['hi'])
    if lo>hi:raise ValueError('ordered interval required')
    return I(lo,hi)


def audit_row(row,shift=SHIFT):
    if type(shift) not in (F,int,str):raise ValueError('exact rational shift required')
    shift=F(shift)
    if not 0<shift<=F(1,10000):raise ValueError('positive attenuation perturbation in proved domain required')
    if not row['strict_self_map'] or not row['conditional_unique_root_in_domain']:
        raise ValueError('original certified root required')
    d=box(row['uniqueness_domain_alpha']);image=box(row['fixed_point_alpha'])
    lipschitz=F(row['lipschitz_upper']);derivative=box(row['derivative'])
    if not max(abs(derivative.lo),abs(derivative.hi))<=lipschitz<F(734,100000):
        raise ValueError('actual contraction certificate differs')
    margin=min(image.lo-d.lo,d.hi-image.hi)
    if not margin>0:raise ValueError('strict centered image required')
    checks=row['implicit_and_mass_enclosures'];u=box(checks['alpha_u'])
    for entry,sign in ((checks['alpha_u_verification'],'R_u_sign_definite'),
                       (checks['mz_verification'],'h_m_sign_definite')):
        if not entry['endpoint_signs_verified'] or not entry[sign]:raise ValueError('implicit source branch unverified')
    if box(checks['alpha_inv']).lo<=101:raise ValueError('original full response lower bound required')
    phi=(1+sqrt_rational(5))/2
    # The pi enclosure is rational; this square-root endpoint extension is outward.
    pi=pi_interval();root_pi=I(sqrt_rational(pi.lo).lo,sqrt_rational(pi.hi).hi)
    pixel=phi+root_pi*d
    evaluated=shape.evaluate(pixel,u)
    lower=min(x.lo for x in evaluated['positive_roots'])
    if lower<=F(1,100):raise ValueError('uniform positive root margin failed')
    chi=evaluated['chi'];shifted=chi+shift
    if not shifted.lo>chi.hi:raise ValueError('perturbed attenuation is not separated on full root boxes')
    difference=shift/25
    if not difference<margin:raise ValueError('perturbed map not certified as a strict self-map')
    q0=(1+elementary.exponential(-2*chi))/3
    qs=(1+elementary.exponential(-2*shifted))/3
    if not qs.hi<q0.lo:raise ValueError('changed dimensionless shape is unresolved')
    error=600*shift
    if not error<1:raise ValueError('mass enlargement outside rational bound')
    original=[box(x) for x in checks['masses_over_v']]
    expanded=[I(x.lo*(1-error),x.hi/(1-error)) for x in original]
    inverse=memory.certify({'implicit_and_mass_enclosures':{'masses_over_v':[
        {'lo':str(x.lo),'hi':str(x.hi)} for x in expanded]}})
    if len(inverse['inverse_branches'])!=2:raise ArithmeticError('complete controller branches missing')
    return {'design':row['design'],'shift':str(shift),
        'original_alpha_domain':d.strings(),'original_image':image.strings(),
        'original_lipschitz_upper':str(lipschitz),'strict_image_margin':str(margin),
        'whole_domain_u':u.strings(),'whole_domain_positive_root_lower':str(lower),
        'original_response_lower':str(box(checks['alpha_inv']).lo),
        'log_mass_derivative_upper':'600','response_derivative_upper':'400',
        'closure_uniform_perturbation_upper':str(difference),
        'perturbed_alpha_image':(image+I(-difference,difference)).strings(),
        'distance_to_original_root_upper':str(difference/(1-lipschitz)),
        'original_chi_interval':chi.strings(),'perturbed_chi_interval':shifted.strings(),
        'original_Q_interval':q0.strings(),'perturbed_Q_interval':qs.strings(),
        'Q_intervals_strictly_disjoint':True,
        'perturbed_fixed_point_existence_certified':True,
        'perturbed_numeric_root_computed':False,'explicit_shift_root_uniqueness_certified':False,
        'perturbed_masses_over_v':[x.strings() for x in expanded],
        'perturbed_complete_memory_inverse':inverse,
        'source_program_is_identical_to_original':False,
        'physical_law_selected_or_empirical_fit':False}


def certificate():
    pins=owner.source_bindings();retained=shape.authenticate(owner.OUTPUT,pins)
    rows=[audit_row(r) for r in retained['closures']]
    if [r['design'] for r in rows]!=['mssm_beta4_proxy','sm_color_beta7']:
        raise ValueError('both actual conditional scale designs required')
    if pins!=owner.source_bindings():raise RuntimeError('source changed during audit')
    return {'schema':'oph.capacity_nonselection.v1','rows':rows,
        'parent_receipt_sha256':hashlib.sha256(owner.OUTPUT.read_bytes()).hexdigest(),
        'parent_source_pin_count':len(pins),'original_source_law_or_receipts_modified':False,
        'parent_full_root_producer_cold_replayed_here':False,
        'old_fixed_law_unique_root_disputed':False,
        'pixel_self_consistency_selects_attenuation_law':False,
        'additional_derivative_radius_for_unique_parametric_family_quantified':False,
        'perturbed_native_tape_executed':False,'physical_particle_assignments':[],
        'experimental_mass_inputs':[]}
