"""Exact family-frame witness to the independently configured charged source."""
from functools import lru_cache
import numpy as np
import sympy as s
from model.standard_model_source_map import koide_source as configured
from model.lepton_cycle import instrument
from model.lepton_register_round import carrier,rounds
from . import analyzer


def _clean(matrix):
    return matrix.applyfunc(s.simplify)


def _strings(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


@lru_cache(maxsize=1)
def certificate():
    # Extract the existing source polynomial from its own executable values.
    # No native effect enters these four configured-source evaluations.
    samples={z:configured.exact_source(z)['C_exclusive'] for z in (s.S.One,-s.S.One,s.I,-s.I)}
    constant=_clean((samples[1]+samples[-1])/2)
    cosine=_clean((samples[1]-samples[-1])/2)
    sine=_clean((samples[s.I]-samples[-s.I])/2)
    u,v=s.symbols('u v',real=True)
    c=constant+u*cosine+v*sine
    z=u+s.I*v
    prescribed=s.sqrt(configured.WEIGHTS[0])*s.eye(3)+s.sqrt(configured.WEIGHTS[1])*(z*configured.P+s.conjugate(z)*configured.P**2)
    if _clean(c-prescribed)!=s.zeros(3):raise ArithmeticError('configured matrix polynomial changed')
    omega=-s.Rational(1,2)+s.sqrt(3)*s.I/2
    fourier=_clean(s.Matrix(3,3,lambda j,k:omega**(-j*k)/s.sqrt(3)))
    if _clean(fourier.H*fourier)!=s.eye(3) or _clean(fourier.H*configured.P*fourier-s.diag(1,omega,omega**2))!=s.zeros(3):
        raise ArithmeticError('literal Fourier frame failed')
    cs=(u,-u/2-s.sqrt(3)*v/2,-u/2+s.sqrt(3)*v/2)
    roots=[s.sqrt(s.Rational(5,108))+2*s.sqrt(s.Rational(5,216))*x for x in cs]
    native=s.diag(*(x*x for x in roots))
    if _clean(fourier.H*c*fourier-s.sqrt(3)*s.diag(*roots))!=s.zeros(3):
        raise ArithmeticError('native amplitudes do not intertwine the existing configured C')
    if _clean(fourier*native*fourier.H-c.H*c/3)!=s.zeros(3):
        raise ArithmeticError('native effect does not intertwine the configured charged coefficient')
    rotation,geometric_frame=configured._native_rotation()
    if _clean(rotation*geometric_frame-geometric_frame*configured.P)!=s.zeros(3):
        raise ArithmeticError('existing native geometric family frame failed')
    native_to_geometric=_clean(geometric_frame*fourier)
    if _clean(native_to_geometric.H*native_to_geometric)!=s.eye(3):
        raise ArithmeticError('complete geometric frame is not unitary')
    theta=s.Rational(2,9)
    selected=c.subs({u:s.cos(theta),v:s.sin(theta)})
    selected_y=_clean(selected.H*selected/3)
    _,values,_=instrument.analyzer_data(float(analyzer.PHASE))
    f=np.asarray(fourier,complex);yn=np.asarray(selected_y.evalf(18),complex)
    error=float(np.linalg.norm(f@np.diag(values)@f.conj().T-yn))
    if error>2e-13:raise ArithmeticError('actual native analyzer effect misses the literal family map')
    walk=carrier.carrier_walk();double=rounds.register_branches(walk)
    neighbours=carrier.neighbours()[rounds.PORT]
    other=min(p for p in range(12) if p!=rounds.PORT and p not in neighbours)
    pair=rounds.native_pair_branches(walk,(rounds.PORT,other))
    if any(double[k]/pair[k]!=s.Rational(5,4) for k in ('a2','b2','b2_minus','trace')):
        raise ArithmeticError('the prepared/native-pair factor is not 5/4')
    negative={
        'extra_square_Y_squared':float(np.linalg.norm(yn@yn-yn)),
        'native_pair_substitution_or_lost_five_fourths':float(np.linalg.norm((4/5)*yn-yn)),
        'conditional_active_probability_division':float(np.linalg.norm(yn/float(s.Rational(5,18))-yn)),
        'wrong_conjugate_Fourier_mode_assignment':float(np.linalg.norm(f.conj()@np.diag(values)@f.T-yn))}
    if min(negative.values())<1e-3:raise ArithmeticError('a forbidden identification escaped the negative controls')
    return {'existing_source':'standard_model_source_map.koide_source.exact_source',
        'existing_source_polynomial_extracted_at_exact_phases':['1','-1','I','-I'],
        'selected_phase_exact':'exp(2 I/9)',
        'existing_rational_phase_API_used_outside_domain':False,
        'extension':'the exact existing polynomial in Re(z), Im(z), evaluated at cos(2/9), sin(2/9); no matrix is defined from the native effect',
        'configured_C_at_selected_phase_exact':_strings(selected),
        'character_to_configured_cyclic_frame_exact':_strings(fourier),
        'native_geometric_to_cyclic_frame_exact':_strings(geometric_frame),
        'character_to_native_geometric_frame_exact':_strings(native_to_geometric),
        'exact_amplitude_identity':'F† C_ex F=sqrt(3) diag(lambda_k)',
        'exact_effect_identity':'F E_native F†=C_ex† C_ex/3',
        'exact_geometric_identity':'(O F) E_native (O F)†=O (C_ex† C_ex/3) O†',
        'actual_native_effect_numeric_residual':error,
        'prepared_double_trace_exact':str(double['trace']),
        'distinct_native_pair_trace_exact':str(pair['trace']),
        'prepared_to_native_pair_factor_exact':'5/4',
        'negative_control_operator_distances':negative,
        'admission':{'same_gauge_trivial_three_family_representation':True,
            'literal_effect_to_existing_configured_action_frame_map':True,
            'full_CAR_dynamics_supplied_by_this_compiler':True,
            'configured_eight_flag_and_native_3600_address_complete_source_laws_identical':False,
            'native_pair_preparation_substituted':False,'extra_square_of_effect_used':False,
            'phase_selected_by_passive_source':False,
            'canonical_physical_Higgs_residue_or_weak_scale_identified':False,
            'observed_interacting_lepton_dictionary_proved':False}}
