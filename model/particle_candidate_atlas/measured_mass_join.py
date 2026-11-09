"""Necessary matching budgets, preserving the four configured source boxes."""
from fractions import Fraction as F
from math import prod
from pathlib import Path
import hashlib
import json

from model.charged_lepton_response import p_dependence as source
from model.native_flavor_orbit.intervals import I, rounded_out
from model.mcpr_native_memory_bridge import certificate as elementary

ROOT=source.ROOT
BRANCHES=('forward','measured_endpoint','mssm_beta4_proxy','sm_color_beta7')
SOURCE_INPUTS=source.SOURCE_INPUTS+(
    'model/charged_lepton_response/MASS_RESPONSE_BRIDGE.md',
    'model/pixel_mass_closure/closure.py',
    'model/mcpr_native_memory_bridge/certificate.py',)


def rational(x):
    if type(x) not in (str,int,F):raise ValueError('exact rational required')
    try:return F(x.numerator,x.denominator) if type(x)is F else F(x)
    except (ValueError,ZeroDivisionError) as exc:raise ValueError('exact rational required') from exc


def box(x):
    if type(x) not in (list,tuple) or len(x)!=2:raise ValueError('two exact endpoints required')
    a,b=map(rational,x)
    if not 0<a<=b:raise ValueError('positive ordered endpoints required')
    return a,b


def cube_bounds(value,bits=100):
    q=rational(value)
    if q<0 or type(bits)is not int or not 1<=bits<=256:raise ValueError('nonnegative value and bounded integer precision required')
    scale=1<<bits;n=(q.numerator*scale**3)//q.denominator
    lo,hi=0,1<<((n.bit_length()+2)//3)
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**3<=n:lo=mid
        else:hi=mid
    lower=F(lo,scale)
    upper=lower if lower**3==q else F(lo+1,scale)
    if not lower**3<=q<=upper**3:raise ArithmeticError('integer cube enclosure failed')
    return lower,upper


def higgs_unit_masses(row):
    u=I(*box(row['source']['u']))
    squares=[I(*box(v)) for v in row['source']['shape_squares']]
    if len(squares)!=3:raise ValueError('three positive source shape squares required')
    kappa=(u/50)/(1+u/31+u*u/310)
    cube=elementary.exponential(3*kappa)/(2*6**14*prod(squares))
    scale=I(cube_bounds(cube.lo)[0],cube_bounds(cube.hi)[1])
    return [rounded_out(scale*x,28) for x in squares],rounded_out(scale,28)


def _inputs(ratios,references):
    r,t=[box(x) for x in ratios],[box(x) for x in references]
    if len(r)!=3 or len(t)!=3 or r[0]!=(F(1),F(1)) or t[0]!=(F(1),F(1)):
        raise ValueError('three positive ratios with the same exact unit denominator required')
    if any(a[1]>=b[0] for a,b in zip(r,r[1:])) or any(a[1]>=b[0] for a,b in zip(t,t[1:])):
        raise ValueError('separated ascending ratios required')
    return r,t


def necessary_budgets(ratios,references):
    r,t=_inputs(ratios,references)
    dlo=max(F(1),max(a/u for (a,b),(l,u) in zip(t,r))/min(b/l for (a,b),(l,u) in zip(t,r)))
    dhi=max(F(1),max(a/l for (a,b),(l,u) in zip(t,r))/min(b/u for (a,b),(l,u) in zip(t,r)))
    elo=max(F(0),max((t[j][0]*r[i][0]-t[i][1]*r[j][1])/(t[j][0]+t[i][1]) for i in range(3) for j in range(3)))
    ehi=max(F(0),max((t[j][0]*r[i][1]-t[i][1]*r[j][0])/(t[j][0]+t[i][1]) for i in range(3) for j in range(3)))
    return {'minimal_component_multiplier_spread_interval':[str(dlo),str(dhi)],
        'necessary_relative_matching_error_lower':str((dlo-1)/(dlo+1)),
        'necessary_equal_L_R_kinetic_epsilon_lower':str((dlo-1)/(dlo+1)),
        'necessary_one_sided_kinetic_epsilon_lower':str((dlo*dlo-1)/(dlo*dlo+1)),
        'minimal_additive_E_over_source_light_mass_interval':[str(elo),str(ehi)],
        'a_nonzero_correction_is_necessary_in_this_constant_matching_class':dlo>1,
        'actual_source_correction_upper_bound_available':False}


def additive_witness(ratios,references):
    """Sharpness control in the allowed matrix class; no source correction claim."""
    r=[rational(v) for v in ratios]
    b=necessary_budgets([[x,x] for x in r],references)
    eta=F(b['minimal_additive_E_over_source_light_mass_interval'][0])
    windows=[box(x) for x in references]
    lo=max((x-eta)/u for x,(l,u) in zip(r,windows))
    hi=min((x+eta)/l for x,(l,u) in zip(r,windows))
    if lo>hi or hi<=0:raise ArithmeticError('additive feasibility interval failed')
    scale=(max(F(0),lo)+hi)/2
    shifted=[min(max(x,scale*l),scale*u) for x,(l,u) in zip(r,windows)]
    if shifted[0]!=scale or max(abs(x-y) for x,y in zip(r,shifted))>eta:
        raise ArithmeticError('attaining diagonal control failed')
    return {'eta':str(eta),'electron_multiplier':str(scale),
        'shifted_singular_values':list(map(str,shifted)),
        'allowed_matrix_control_is_computed_source_self_energy':False}


def summarize(capacity,reference):
    if [r['branch'] for r in capacity['branches']]!=list(BRANCHES):raise ValueError('all four source branches required')
    if (capacity['source_parameters_fitted_by_this_certificate'] is not False
        or capacity['physical_interacting_mass_predictions']!=[]):raise ValueError('source scope changed')
    refs=[[F(1),F(1)]]
    for name in ('muon_electron_ratio','tau_electron_ratio'):
        c,e=map(rational,(reference[name]['value'],reference[name]['standard_uncertainty']))
        if not 0<e<c:raise ValueError('positive reference estimate and uncertainty required')
        refs.append([c-e,c+e])
    electron=rational(reference['electron_MeV']['value'])/1000
    rows=[]
    for row in capacity['branches']:
        masses,scale=higgs_unit_masses(row)
        ratios=row['source']['charged_action_ratios']
        budgets=necessary_budgets(ratios,refs)
        eta_lo,eta_hi=map(F,budgets['minimal_additive_E_over_source_light_mass_interval'])
        qlo,qhi=box(row['source']['positive_root_Koide_ratio'])
        balance_d=max(F(1),F(2,3)/qhi)
        rows.append({'branch':row['branch'],'ancestry':row['ancestry'],
            'certified_source_coordinate_boxes':{k:row['source'][k] for k in ('P','u','chi','delta','rho')},
            'masses_over_declared_Higgs_amplitude':[x.strings() for x in masses],
            'declared_common_amplitude_A':scale.strings(),
            'original_action_ratios':ratios,
            'source_positive_root_Koide_interval':[str(qlo),str(qhi)],
            'optional_electron_central_unit_GeV':[(I(*box(r))*electron).strings() for r in ratios],
            'optional_GeV_column_propagates_electron_uncertainty':False,
            'matching_budgets':budgets,
            'necessary_additive_norm_over_Higgs_interval':[str(eta_lo*masses[0].lo),str(eta_hi*masses[0].hi)],
            'additional_exact_balance_hypothesis':{
                'assumed_for_reference_comparison':False,
                'necessary_Q_distortion_lower':str(balance_d),
                'necessary_relative_matching_error_lower':str((balance_d-1)/(balance_d+1))},
            'complete_memory_inverse_branches':row['admitted_complete_memory_inverse_branches'],
            'physical_mass_predictions':[]})
    return {'schema':'oph.measured_mass_join.v1','branches':rows,
        'reference_ratio_windows':[list(map(str,r)) for r in refs],
        'reference_window_kind':'quoted CODATA2022 standard uncertainty; no joint significance claim',
        'optional_electron_central_energy_unit_GeV':str(electron),
        'source_parameters_altered_to_match_reference':False,
        'computed_interacting_correction':False,'physical_particle_map':[],
        'proved_physical_correction_upper_bounds':[],
        'generic_positive_family_CAR_and_free_field_map_is_conditional':True}


def certificate():
    capacity=source.certificate()
    path=ROOT/source.REFERENCE;raw=path.read_bytes();reference=json.loads(raw)
    result=summarize(capacity,reference)
    result['authenticated_parent_receipts_and_reference']=capacity['authenticated_receipts_and_reference']
    result['calibration_reference_sha256']=hashlib.sha256(raw).hexdigest()
    result['capacity_certificate_sha256']=hashlib.sha256((json.dumps(capacity,sort_keys=True,indent=2)+'\n').encode()).hexdigest()
    return result
