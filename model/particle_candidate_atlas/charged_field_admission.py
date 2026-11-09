"""Generic positive-family full-CAR/free-field attachment and scoped admissions."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
import hashlib
from pathlib import Path

import numpy as np
import sympy as s
from model.common_source import native_matter_attachment as matter
from model.common_source import native_family_record_memory as memory
from model.common_source import native_record_yukawa_field as field
from model.standard_model_source_map import feedback_limit as feedback
from model.standard_model_source_map import algebra as source_algebra
from model.standard_model_source_map import fermion_attachment
from model.charged_lepton_response import relativistic_mass

ROOT=Path(matter.__file__).resolve().parents[2]
SOURCE_INPUTS=(
    'model/standard_model_source_map/particle_map.json',
    'model/common_source/native_matter_attachment.py',
    'model/common_source/native_family_record_memory.py',
    'model/common_source/native_record_yukawa_field.py',
    'model/standard_model_source_map/FEEDBACK_LIMIT.md',
    'model/standard_model_source_map/FERMION_ATTACHMENT.md',
    'model/standard_model_source_map/KINETIC_ATTACHMENT.md',
    'model/standard_model_source_map/GAUGE_ATTACHMENT.md',
    'model/mcpr_native_memory_bridge/PROOF.md',
    'model/charged_lepton_response/RELATIVISTIC_MASS.md',
    'model/charged_lepton_response/MASS_RESPONSE_BRIDGE.md',
)

def rational(value):
    if type(value) not in (int,str,F):raise ValueError('an exact rational is required')
    try:return F(value.numerator,value.denominator) if type(value)is F else F(value)
    except (ValueError,ZeroDivisionError) as exc:raise ValueError('an exact rational is required') from exc

def authenticated_parent():
    path=ROOT/'model/standard_model_source_map/particle_map.json'
    raw=path.read_bytes();parent=json.loads(raw)
    for name,digest in parent['source_bindings'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('parent source pin differs: '+name)
    if parent['source_weights']!=list(source_algebra.source_components()):
        raise ValueError('parent occupation order or current weights differ')
    return parent,{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(raw).hexdigest(),
                   'source_pin_count':len(parent['source_bindings'])}

def clean(a):
    return s.ImmutableSparseMatrix(s.Matrix(a).applyfunc(lambda x:s.factor(s.radsimp(s.expand(x)))))

def positive_family(value):
    y=s.ImmutableMatrix(value)
    if y.shape!=(3,3) or any(v.has(s.Float) or v.free_symbols or v.is_finite is not True for v in y):
        raise ValueError('an exact finite three-by-three family matrix is required')
    if clean(y-y.H)!=s.zeros(3):raise ValueError('the family coefficient must be Hermitian')
    if any(s.simplify(y[:j,:j].det()).is_positive is not True for j in (1,2,3)):
        raise ValueError('positive definiteness must be proved by exact principal minors')
    return y

@lru_cache(maxsize=1)
def _source():
    rep=matter.representation();ids=list(rep['selected'])
    weak=rep['operators']['weak_C2'].extract(ids,ids)
    left=s.diag(*(int(weak[i,i]!=0) for i in range(15)))
    lines=[r for r in matter.invariant_higgs_lines() if r['neutral_rank']==2]
    if len(lines)!=1:raise ArithmeticError('unique source charged-Higgs line required')
    d=clean(left*lines[0]['components'][1]*(s.eye(15)-left))
    if d.rank()!=1 or clean(d*d)!=s.zeros(15) or s.trace(d.H*d)!=1:
        raise ArithmeticError('source charged partial-isometry normalization differs')
    pl,pr=clean(d*d.H),clean(d.H*d)
    if clean(pl*pr)!=s.zeros(15) or clean(pl*pl-pl)!=s.zeros(15):
        raise ArithmeticError('source charged supports are not orthogonal')
    return d,pl,pr

def family_kernel(value):
    y=positive_family(value);d,pl,pr=_source()
    directed=s.kronecker_product(d,y.conjugate())
    k=clean(directed+directed.T)
    expected=s.kronecker_product(pl,(y*y.H).conjugate())+s.kronecker_product(pr.T,y.H*y)
    if clean(k*k.H-expected)!=s.zeros(45):raise ArithmeticError('generic doubled singular spectrum failed')
    if k!=k.T:raise ArithmeticError('all-left mass kernel is not symmetric')
    return k

def physical_controller_budget(steps,time,p_upper,*,cells=1,gain_upper=1,higgs_norm=1):
    if type(cells)is not int or cells<1:raise ValueError('positive finite cell count required')
    p,g,v,t=map(rational,(p_upper,gain_upper,higgs_norm,time))
    if not 0<=p<1 or min(g,v,t)<0:raise ValueError('admitted probability and nonnegative bounds required')
    joint=1-(1-p)**cells;norm=2*cells*g*v
    b=feedback.averaging_bound(steps,t,joint,feedback_norm=norm)
    return {k:str(v) if isinstance(v,F) else v for k,v in b.items()}|{
        'joint_nonreset_probability':str(joint),'full_Fock_pairing_norm_upper':str(norm),
        'physical_CAR_modes':90*cells,'volume_uniform':False,
        'native_attempt_clock_identified':False,'original_path_or_failures_postselected':False}

def native_C_control():
    source=memory.architecture();T=s.zeros(12);weights=[s.S.Zero]*12;count=0
    for old in range(12):
        for cc,cr,ec,er,_,w in source['records']:
            c,e=cc==0 and cr==1,ec==0 and er in(4,5)
            dest=memory._destination(old,c,e,2,source['gC'],source['gE'])
            T[dest,old]+=w
            if c:weights[dest]+=source['stationary'][old]*w
            count+=1
    mean=sum((w*p for w,p in zip(weights,source['projectors'],strict=True)),s.zeros(3))
    y=positive_family(clean(mean));k=family_kernel(y)
    if T!=source['T'] or count!=1536 or s.trace(y)!=s.Rational(1,12):
        raise ArithmeticError('complete native C-duty source mean differs')
    other=feedback.record_model(F(1,12),F(1,6),'C')['mean']
    if clean(y-other)!=s.zeros(3):raise ArithmeticError('fine/coarse source means differ')
    report=field.dispersion_certificate(np.array(k.evalf(17),complex),speed=1)
    if report['dispersion_spectral_error']>1e-12 or report['Hermiticity_error']>1e-12:
        raise ArithmeticError('actual180-dimensional native free symbol failed')
    return y,k,{'fine_records_all_old_inputs':count,'native_fine_records_per_old_input':128,
        'source_p':'1/12','source_q':'1/6','gain':'1','duty':'C',
        'unconditional_trace':'1/12','physical_pairing_norm':'1/6',
        'Y':[[str(x) for x in row] for row in y.tolist()],
        'exact_positive_principal_minors':[str(s.simplify(y[:j,:j].det())) for j in (1,2,3)],
        'full_native_180_symbol_diagnostic':report,
        'diagnostic_is_outward_roundoff_proof':False,'coherent_Koide_balance_assumed':False}

def spectator_census(k,parent):
    weights=parent['source_weights'];d,pl,pr=_source()
    if weights!=list(source_algebra.source_components()):raise ValueError('source weight order differs')
    zero_internal=[i for i in range(15) if pl[i,i]==0 and pr[i,i]==0]
    if len(zero_internal)!=13:raise ArithmeticError('charged support count differs')
    for i in zero_internal:
        for f in range(3):
            if any(k[3*i+f,j]!=0 for j in range(45)):raise ArithmeticError('a spectator acquired a mass')
    rows=[{'source_mask':weights[i]['occupation_mask'],'Q':weights[i]['Q'],
           'family_copies':3} for i in zero_internal]
    charged=sum(3 for r in rows if F(r['Q'])!=0);neutral=sum(3 for r in rows if F(r['Q'])==0)
    if (charged,neutral)!=(36,3):raise ArithmeticError('source spectator charge census differs')
    # A color-singlet triple of down-left zero modes has total charge -1.
    down=[i for i in zero_internal if weights[i]['Q']=='-1/3' and weights[i]['weak_C2']=='3/4']
    if len(down)!=3:raise ArithmeticError('three source color components required')
    modes=[fermion_attachment.mode_index(0,i,0) for i in down]
    mask=sum(1<<m for m in modes)
    rep=matter.representation();q=rep['operators']['T3']+rep['operators']['q']/6
    charges=[q[rep['selected'][i],rep['selected'][i]] for i in down]
    if sum(charges)!=-1 or mask.bit_count()!=3:raise ArithmeticError('spectator charge control failed')
    # Literal native ports mix internal directions. Recover the pure color
    # basis through the owning exact twelve-port control map before lifting.
    color=tuple(matter.exterior_derivative(matter.spin_lift(g)).extract(
        list(rep['selected']),list(rep['selected'])) for g in matter.canonical_currents()[:8])
    for g in color:
        full=s.kronecker_product(s.eye(2),g,s.eye(3))
        if fermion_attachment.one_body_action(full,mask):
            raise ArithmeticError('zero-mode triple is not a global color singlet')
    return {'massless_quark_internal_coefficients':charged,'massless_neutral_internal_coefficients':neutral,
            'total_zero_internal_coefficients':39,'rows':rows,
            'free_zero_momentum_control':{'physical_CAR_modes':90,'occupied_mode_indices':modes,
                'particle_number':3,'electric_charge':'-1','global_color_generators_annihilate':8,
                'color_basis_admitted_as_combinations_of_all12_native_ports':True,
                'energy_above_a_free_ground_with_empty_spectator_zero_modes':'0',
                'local_Gauss_constraint_or_interacting_hadron_claim':False},
            'whole_electric_sector_threshold_equals_each_charged_family_mass':False}

def default_action_exclusions(parent):
    data=parent['finite_source_attachments']['source_action_masses']
    rows=[]
    for policy in data['policies']:
        c=next(r for r in policy['channels'] if r['channel']=='charged_lepton')
        box=list(map(F,c['Koide_diagnostic']['positive_root_Koide_ratio_interval']))
        if box[0]>box[1] or box[0]<=F(2,3)<=box[1]:
            raise ValueError('owning exact Koide exclusion absent')
        rows.append({'duty':policy['policy'],'action_mass_ratios':c['ratios_to_smallest_interval'],
            'positive_root_Koide_interval':list(map(str,box)),'two_thirds_excluded':True})
    if [r['duty'] for r in rows]!=['always','C','E','CE']:raise ValueError('all four duties required')
    c=rows[1]['action_mass_ratios']
    reference=(F('206.7682827')-F('.0000046'),F('206.7682827')+F('.0000046'))
    ratio=tuple(map(F,c[1]))
    if not ratio[1]<reference[0]:raise ValueError('literal C-duty reference exclusion failed')
    return {'current_gated_actions':rows,'literal_C_duty_muon_ratio_disjoint_from_reference':True,
        'reference_muon_electron_CODATA2022_interval':list(map(str,reference)),
        'reference_url':'https://physics.nist.gov/cuu/pdf/all.pdf',
        'scope':'direct canonical-free action identification under one common unit only; no exclusion of uncomputed interactions',
        'distinct_coherent_default_ratios':['1','1','34+24*sqrt(2)'],
        'coherent_default_is_current_gated_C_action':False}

def admission_rows(parent):
    expected=('electron','muon','tau','up','charm','top','down','strange','bottom',
              'electron_neutrino','muon_neutrino','tau_neutrino','photon','gluon','W_plus','W_minus','Z','Higgs')
    if [r['target'] for r in parent['rows']]!=list(expected):raise ValueError('exact18-target interface required')
    rows=[]
    for i,name in enumerate(expected):
        green=['exact internal-current interface']
        red=[];amber=['natural-particle identity','physical clock and selected interacting state','source-defect creation map']
        if i<12:green+=['configured odd CAR statistics','specified free kinetic spinor structure']
        if i<3:
            green+=['general positive-Y full-Fock pairing','three free Dirac families per valley','same-clock free dispersion']
            red=[{'claim':'literal default coherent source has three unequal charged masses',
                  'reason':'z=1 gives an exact repeated pair; a common unit cannot remove it'},
                 {'claim':'literal current-gated C action directly gives observed charged ratios under one common unit',
                  'reason':'exact owning ratio boxes give approximately1:119.900581:896.516443'},
                 {'claim':'charged-only family rest mass is the unrestricted electric-sector gap',
                  'reason':'same free action has massless charged spectators'}]
            amber+=['interacting charged threshold' if i==0 else 'unstable resonance and weak-decay width',
                    'kinetic/self-energy/observable-matching error bounds','source phase/controller selection or calibration ancestry']
        elif i<9:
            red=[{'claim':'this charged-only free spectrum supplies positive quark masses',
                  'reason':'the corresponding36 internal coefficients vanish exactly'}]
            amber+=['quark family mass/mixing action','gauge-invariant short-distance or hadronic matching']
        elif i<12:
            red=[{'claim':'this charged-only free action has nonzero neutrino mass splittings',
                  'reason':'all three neutral internal coefficients vanish exactly'}]
            amber+=['neutral mass mechanism','flavor-to-mass mixing; no individual flavor mass assignment']
        else:
            amber+=['propagating gauge/Higgs dynamics and selected vacuum',
                    'physical polarization or scalar response','pole/threshold/resonance observable']
        rows.append({'target':name,'green_proved_in_specified_model':green,
            'red_exactly_excluded_assertions':red,'amber_unproved_attachments':amber,
            'physical_particle_identified':False,'physical_mass':None,'defect_link':None})
    return rows

def certificate():
    y,k,control=native_C_control()
    complex_control=s.Matrix([[2,s.I,0],[-s.I,3,1],[0,1,4]])
    family_kernel(complex_control)
    clifford=relativistic_mass.clifford_certificate()
    parent,binding=authenticated_parent()
    return {'generic_positive_family_theorem':True,
        'whole_Fock_pairing_map':'V(Y), real-linear; norm2tr(Y) for positiveY',
        'source_C_duty_control':control,'complex_positive_family_control':[[str(x) for x in r] for r in complex_control.tolist()],
        'exact_Clifford_certificate':clifford,'spectators':spectator_census(k,parent),
        'authenticated_parent':binding,'distinct_default_action_exclusions':default_action_exclusions(parent),
        'retained_valleys':[{'bits':list(b),'winding':-(-1)**sum(b)} for b in product((0,1),repeat=3)],
        'global_scoped_obstructions':[
            {'claim':'the unmodified sine regulator has only one Weyl valley',
             'reason':'all8 zeros are retained, with4 positive and4 negative windings'},
            {'claim':'a nilpotent quantum raising link equals an identity background on a nonzero invariant space',
             'reason':'nilpotency forbids eigenvalue1; a background approximation needs separate state/error data'},
            {'claim':'preparing an excitation of the fixed free Hamiltonian creates a new kinetic winding',
             'reason':'the specified Hamiltonian symbol and its winding are unchanged by the state preparation'}],
        'configured_family_count':3,'families_selected_by_primitive_repair':False,
        'MCPR_attachment':'conditional complete C-duty inverse plus configured full-CAR feedback; no coherent-analyzer equivalence assumed',
        'controller_budget_example':physical_controller_budget(10**6,'1','1/12',cells=1),
        'admission_rows':admission_rows(parent),'physical_particle_map':[],
        'physical_mass_predictions':[],'interacting_pole_or_threshold_proved':False,
        'new_native_full_Fock_pulse_tape_executed':False}
