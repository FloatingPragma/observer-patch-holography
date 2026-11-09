"""Authenticate finite candidate inputs and derive an exact comparison atlas."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/'docs/research/OPH_PARTICLE_CANDIDATE_ATLAS.svg'
AUDIT = HERE/'audit.json'
CHILD = 'model/charged_lepton_response/results.json'
MAP = 'model/charged_lepton_response/candidate_map.json'
PARENT = 'model/standard_model_source_map/particle_map.json'
REFS = 'model/particle_candidate_atlas/reference_values.json'
CERTIFICATES = {'full_CAR_realization', 'phase_and_action_readout',
    'source_phase_selection_audit', 'free_relativistic_field',
    'conditional_mass_response', 'calibrated_sanity_test',
    'conditional_P_dependence', 'candidate_particle_map'}
TARGETS = ('electron','muon','tau','up','charm','top','down','strange','bottom',
           'electron_neutrino','muon_neutrino','tau_neutrino',
           'photon','gluon','W_plus','W_minus','Z','Higgs')
BRANCHES = ('forward','measured_endpoint','mssm_beta4_proxy','sm_color_beta7')


def canonical(obj):
    return json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)+'\n'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(path):
    def pairs(items):
        d={}
        for key,value in items:
            if key in d: raise ValueError('duplicate JSON key: '+key)
            d[key]=value
        return d
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))


def authenticate(root, pins):
    if not isinstance(pins,dict) or not pins: raise ValueError('missing input bindings')
    for name, expected in pins.items():
        path=(root/name).resolve()
        if not path.is_relative_to(root.resolve()) or digest(path)!=expected:
            raise ValueError('stale or escaped atlas input: '+name)


def load_inputs(root=ROOT):
    child=load_json(root/CHILD)
    authenticate(root,child.get('source_bindings'))
    candidate=load_json(root/MAP)
    parent=load_json(root/PARENT)
    authenticate(root,parent.get('source_bindings'))
    refs=load_json(root/REFS)
    from tools.build_completion_program import integration_source_paths
    from . import (memory_diffusion, charged_field_admission, measured_mass_join,
                   interacting_chiral, capacity_nonselection, conditional_standard_model,
                   calibrated_triplet)
    from model.common_source.chiral_candidate import action
    direct={MAP,PARENT,REFS,
        'code/particles/leptons/derive_charged_face_incidence_conditional_theorem.py',
        'model/native_energy_measurement/PROOF.md',
        'model/particle_memory_response/PROOF.md',
        'model/compact_defect_pair/results.json',
        'model/compact_defect_pair/PROOF.md',
        'model/compact_defect_pair/certification.py',
        'model/mcpr_native_memory_bridge/PROOF.md',
        'model/standard_model_source_map/FEEDBACK_LIMIT.md',
        'model/charged_lepton_response/RELATIVISTIC_MASS.md'}
    for module in (memory_diffusion, charged_field_admission, measured_mass_join,
                   interacting_chiral, capacity_nonselection, conditional_standard_model,
                   calibrated_triplet):
        direct.update(module.SOURCE_INPUTS)
    direct.update(action.SOURCE_FILES)
    direct.update(str(p.relative_to(root)) for p in (root/'model/particle_candidate_atlas').iterdir()
                  if p.is_file() and p.suffix in ('.py','.md'))
    paths=integration_source_paths(root,{
        'domain_and_proof':'model/particle_candidate_atlas/README.md',
        'source_receipts':[CHILD], 'independent_checks':sorted(direct)})
    pins={name:digest(root/name) for name in sorted(paths)}
    return child,candidate,parent,refs,pins


def box(xs):
    if not isinstance(xs,(list,tuple)) or len(xs)!=2: raise ValueError('interval pair required')
    lo,hi=map(F,xs)
    if lo>hi: raise ValueError('reversed interval')
    return lo,hi


def decimal(q, digits=30):
    q=F(q)
    with localcontext() as ctx:
        ctx.prec=digits
        return str(Decimal(q.numerator)/Decimal(q.denominator))


def number(xs, places=9):
    lo,hi=box(xs)
    with localcontext() as ctx:
        ctx.prec=max(60,places+30)
        mid=(Decimal(lo.numerator)/Decimal(lo.denominator)+Decimal(hi.numerator)/Decimal(hi.denominator))/2
        return f'{mid:.{places}f}'


def compare(xs, center, uncertainty):
    lo,hi=box(xs);r,u=F(center),F(uncertainty)
    if u<=0: raise ValueError('positive reference uncertainty required')
    return {'source_interval':[str(lo),str(hi)],'numerical_width':str(hi-lo),
        'center_decimal':decimal((lo+hi)/2),
        'reference_center':str(r),'reference_standard_uncertainty':str(u),
        'residual':[str(lo-r),str(hi-r)],
        'residual_over_quoted_uncertainty':[str((lo-r)/u),str((hi-r)/u)],
        'relative_residual_ppm':[str((lo/r-1)*10**6),str((hi/r-1)*10**6)],
        'inside_quoted_window':r-u<=lo<=hi<=r+u,
        'joint_statistical_significance':False}


def spectrum_crosscheck(cal,p):
    """Independent numerical diagnostic; rational certificates retain proof authority."""
    import mpmath as mp
    with mp.workdps(100):
        def q(x):
            f=F(x);return mp.mpf(f.numerator)/f.denominator
        def middle(xs):return q(sum(box(xs))/2)
        def inside(xs,value):return q(xs[0])<=value<=q(xs[1])
        phase=cal['phase'];theta=mp.atan2(q(phase['imaginary']),q(phase['real']))
        roots=sorted(mp.sqrt(mp.mpf(10)/72)+2*mp.sqrt(mp.mpf(5)/72)*mp.cos(theta+2*mp.pi*j/3) for j in range(3))
        cases=[('calibrated_coherent',roots,cal['comparison']['source_action_ratios'])]
        for row in p['branches']:
            source=row['source'];P,u=middle(source['P']),middle(source['u'])
            chi=u*u/(512*(1-u/77))
            delta=mp.mpf(2)/9+(P-(1+mp.sqrt(5))/2)**2/(21*mp.pi*(1-u/27-u*u/135))
            rho=mp.sqrt(2)*mp.exp(-chi)
            roots=sorted(1+rho*mp.cos(delta+2*mp.pi*j/3) for j in range(3))
            cases.append((row['branch'],roots,source['charged_action_ratios']))
        for name,roots,ratios in cases:
            if not all(inside(interval,(root/roots[0])**2) for root,interval in zip(roots,ratios)):
                raise ValueError('independent Fourier spectrum outside source enclosure: '+name)
    return {'branches_checked':5,'decimal_precision':100,'all_inside_rational_enclosures':True,
            'diagnostic_replaces_rigorous_source_enclosures':False}


def project(child,candidate,parent,refs,pins):
    """Pure projection, with every physical/nonphysical boundary admitted."""
    if set(child.get('certificates',{}))!=CERTIFICATES: raise ValueError('eight child certificates required')
    if child['certificates']['candidate_particle_map']!=candidate:
        raise ValueError('candidate export differs from authenticated child')
    for obj in (child,candidate,parent):
        if obj.get('physical_particle_map')!=[]: raise ValueError('physical assignment not admitted')
    if child.get('physical_mass_predictions')!=[] or candidate.get('physical_mass_predictions')!=[]:
        raise ValueError('physical mass claim not admitted')
    if child.get('source_selected_phase') is not None: raise ValueError('phase selection not admitted')
    if child.get('calibrated_action_comparison_is_independent_prediction',False):
        raise ValueError('calibrated comparison cannot be called independent prediction')
    if [r['target'] for r in candidate['rows']]!=list(TARGETS): raise ValueError('exact eighteen targets required')
    nulls=('physical_mass','physical_spin','physical_lifetime','physical_width',
           'identified_topological_defect','topological_invariant','defect_creation_map',
           'physical_particle_or_resonance_map')
    for row in candidate['rows']:
        if row['physical_particle_identified'] or any(row[k] is not None for k in nulls):
            raise ValueError('candidate promoted to physical object')
    defects=candidate['separate_defect_inventory']
    if not defects['supplied'] or len(defects['entries'])!=11 or defects['species_links']:
        raise ValueError('separate eleven-object inventory required')
    if len({r['id'] for r in defects['entries']})!=11: raise ValueError('duplicate protected object')
    weights=parent['source_weights']
    if len(weights)!=15 or len({r['occupation_mask'] for r in weights})!=15:
        raise ValueError('fifteen distinct source weights required')
    groups={}
    for row in weights:
        if F(row['Q'])!=F(row['T3'])+F(row['Y']) or 6*F(row['Y'])!=F(row['q6Y']):
            raise ValueError('weight charge identity fails')
        groups.setdefault(row['q6Y'],[]).append(row['occupation_mask'])
    if {k:len(v) for k,v in groups.items()}!={'1':6,'-4':3,'2':3,'-3':2,'6':1}:
        raise ValueError('five-irrep multiplicities differ')
    by_mask={r['occupation_mask']:r for r in weights}
    for row in candidate['rows']:
        for component in row['internal_components'].values():
            for mask in component['source_masks']:
                if mask not in by_mask or any(F(component[k])!=F(by_mask[mask][k]) for k in ('Q','T3','q6Y')):
                    raise ValueError('candidate internal component differs from source weights')
    if {r['target'] for r in refs['masses']}!=set(TARGETS) or len(refs['masses'])!=18:
        raise ValueError('eighteen reference channels required')
    if any(r['unit'] not in ('GeV',None) for r in refs['masses']):
        raise ValueError('all reference energies must use GeV')
    cal=child['certificates']['calibrated_sanity_test']
    if cal['phase']['tau_used_to_choose_phase'] or cal['phase_calibration_is_source_derivation']:
        raise ValueError('calibration ancestry differs')
    if cal['exact_positive_root_Koide_ratio']!='2/3': raise ValueError('balanced branch altered')
    p=child['certificates']['conditional_P_dependence']
    if [b['branch'] for b in p['branches']]!=list(BRANCHES): raise ValueError('four distinct P branches required')
    if p['physical_interacting_mass_predictions'] or p['source_parameters_fitted_by_this_certificate']:
        raise ValueError('P branch scope differs')
    reference=cal['reference']
    for key in ('muon_electron_ratio','tau_electron_ratio'):
        ref=refs['ratios'][key.replace('_ratio','')+'_CODATA2022']
        if any(reference[key][k]!=ref[k] for k in ('value','standard_uncertainty')):
            raise ValueError('reference editions require separate comparison')
    tau_pdg=next(r for r in refs['masses'] if r['target']=='tau')
    electron=F(reference['electron_MeV']['value'])/1000
    descriptions=[('calibrated_coherent',cal['comparison']['source_action_ratios'],
                   'mu/e calibrated; tau held out',['2/3','2/3'])]
    for b in p['branches']:
        if b['admitted_complete_memory_inverse_branches']!=2: raise ValueError('complete inverse branch count differs')
        if box(b['source']['positive_root_Koide_ratio'])[1]>=F(2,3): raise ValueError('P action is not balanced')
        descriptions.append((b['branch'],b['source']['charged_action_ratios'],b['ancestry'],
                             b['source']['positive_root_Koide_ratio']))
    rows=[]
    for name,ratios,ancestry,q in descriptions:
        if len(ratios)!=3 or box(ratios[0])!=(F(1),F(1)): raise ValueError('complete normalized triple required')
        if not all(box(ratios[i])[1]<box(ratios[i+1])[0] for i in (0,1)):
            raise ValueError('ordered positive triple required')
        rr={}
        for i,key in [(1,'muon_electron_ratio'),(2,'tau_electron_ratio')]:
            r=reference[key];rr[key]=compare(ratios[i],r['value'],r['standard_uncertainty'])
        gev=[str(v*electron) for v in box(ratios[2])]
        rr['tau_GeV_CODATA2022']=compare(gev,F(reference['tau_MeV']['value'])/1000,F(reference['tau_MeV']['standard_uncertainty'])/1000)
        rr['tau_GeV_PDG2026']=compare(gev,tau_pdg['value'],tau_pdg['uncertainty'])
        rows.append({'branch':name,'ratios':ratios,'ancestry':ancestry,'positive_root_Koide_ratio':q,
                     'comparisons':rr,'optional_tau_GeV':gev,
                     'optional_charged_GeV':[[str(v*electron) for v in box(r)] for r in ratios],
                     'electron_unit_is_supplied':True,
                     'optional_GeV_values_propagate_electron_uncertainty':False})
    return {'schema':'oph.particle_candidate_atlas.v1','input_sha256':pins,
        'source_binding_count':len(pins),'candidate_rows':candidate['rows'],
        'source_weight_groups':groups,'protected_objects':defects['entries'],
        'reference':refs,'lepton_branches':rows,
        'phase_radians_approx':cal['phase']['phase_radians_approx'],
        'physical_particle_map':[],'physical_mass_predictions':[],
        'source_experiments_replayed_by_atlas':False,
        'independent_fourier_spectrum_diagnostic':spectrum_crosscheck(cal,p),
        'canonical_child_export_consumed':False,
        'comparison_has_physical_theory_error_budget':False,
        'supplied_decimal_coordinate_error_bound':None,
        'proof_links':[
            'model/standard_model_source_map/FERMION_ATTACHMENT.md',
            'model/charged_lepton_response/CAR_KOIDE.md',
            'model/charged_lepton_response/RELATIVISTIC_MASS.md',
            'model/charged_lepton_response/SOURCE_PHASE_SELECTION.md',
            'model/charged_lepton_response/P_DEPENDENCE.md',
            'model/charged_lepton_response/MASS_RESPONSE_BRIDGE.md',
            'model/charged_lepton_response/CANDIDATE_MAP.md',
            'model/particle_candidate_atlas/COMPATIBILITY.md']}


def build(root=ROOT):
    args=load_inputs(root)
    result=project(*args)
    result['finite_joins']=science_projection(args[0],args[2],tuple(sorted(args[-1].items())))
    if load_inputs(root)[-1]!=args[-1]:
        raise RuntimeError('atlas source inventory or bytes changed during construction')
    result['canonical_child_export_consumed']=True
    from .graphic import render
    svg=render(result)
    return result,svg


def science_projection(child,parent,source_key=None):
    if source_key is None:
        return _calculate_science(child,parent)
    return json.loads(_science_text(canonical(child),canonical(parent),source_key))


@lru_cache(maxsize=1)
def _science_text(child_text,parent_text,source_key):
    # The complete input hash map is part of the cache key; callers get fresh
    # JSON trees, never mutable cached certificate objects.
    return canonical(_calculate_science(json.loads(child_text),json.loads(parent_text)))


def _calculate_science(child,parent):
    """Bounded arithmetic plus explicitly identified retained source theorems.

    Complete native experiments and large operator proofs are owned by their
    scientific modules. This projection does not relabel them as executed here.
    """
    from . import (charged_field_admission as field, measured_mass_join as masses,
                   capacity_nonselection as capacity, defect_field_join as defect,
                   native_response as response, attenuation, interacting_chiral,
                   conditional_standard_model as benchmark, calibrated_triplet)
    cal=child['certificates']['calibrated_sanity_test']
    p=child['certificates']['conditional_P_dependence']
    sm=benchmark.certificate()
    fitted=calibrated_triplet.certificate()
    if (not fitted['all_three_masses_used_as_inputs']
            or fitted['held_out_mass_predictions']!=0
            or fitted['source_selected_controller']
            or fitted['physical_mass_matching_proved']
            or fitted['scale_is_selected_physical_Higgs_vev']
            or fitted['native_tape_for_this_fitted_point_executed']
            or fitted['complete_retained_flags']!=16
            or not fitted['all_output_errors_less_than_1e_minus_36_GeV']
            or fitted['reference_central_values_are_exact_true_masses']):
        raise ValueError('three-input controller calibration promoted to prediction')
    if ([r['target'] for r in sm['target_rows']]!=list(TARGETS)
            or sm['physical_OPH_assignments']!=0
            or sm['source_selected_renormalized_parameters']
            or sm['source_selected_vacuum_and_clock']
            or sm['naive_eight_valleys_removed']
            or sm['interacting_continuum_existence_proved']
            or len(sm['neutrino_mass_fields'])!=3
            or any(r['numerical_mass'] is not None for r in sm['neutrino_mass_fields'])
            or any(r['source_predicted_physical_mass'] is not None for r in sm['target_rows'])):
        raise ValueError('conditional benchmark promoted to source selection')
    mass=masses.summarize(p,cal['reference'])
    cap=capacity.certificate()
    if (cap['pixel_self_consistency_selects_attenuation_law']
            or cap['old_fixed_law_unique_root_disputed']
            or len(cap['rows'])!=2):
        raise ValueError('capacity nonselection scope differs')
    for row in cap['rows']:
        if (F(row['shift'])!=F(1,10**14)
                or not row['perturbed_fixed_point_existence_certified']
                or not row['Q_intervals_strictly_disjoint']
                or row['explicit_shift_root_uniqueness_certified']
                or row['source_program_is_identical_to_original']):
            raise ValueError('capacity perturbation domain differs')
    compact={'source':defect.source_receipt(), 'spatial_clusters':defect.parity_schur(),
             'charged_groups':defect.phase_zero_branches(),
             'original_core_weight':defect.localization_readout()}
    if (compact['charged_groups']['complete_positive_charged_Nambu_rows']!=24
            or compact['charged_groups']['static_gap_floor']!='1/40000'
            or compact['charged_groups']['three_observed_particle_species_assigned']
            or [r['positive_Nambu_rows'] for r in compact['charged_groups']['rows']]!=[8,8,4,4]
            or F(compact['spatial_clusters']['leading_positive_energy'])!=F(1,256)
            or F(compact['spatial_clusters']['parity_resolvent_error'])>=F(1,25000)):
        raise ValueError('compact theorem domain changed')
    weights=response.c_duty_weights('1/12','1/6')
    theta_cos=F(99,101)
    native={'p':'1/12','q':'1/6','tangent_half_angle':'1/10',
        'all_four_coarse_record_weights':list(map(str,weights)),
        'logical_purity_loss_lower':str(response.purity_floor('1/12','1/6','1/10')),
        'full_CAR_vacuum_purity_loss_lower':str(2*F(str(weights[0]))*F(str(weights[2]))*(1-theta_cos**4)),
        'scope':'finite unread pure-kick channel; no full-archive or weak-limit exclusion'}
    retained=[
      {'id':'complete_attenuator','claim':'16 retained analyzer flags; configured r and phase; full-CAR lift has 1440 modes',
       'proof':'model/particle_candidate_atlas/BARE_EFFECTIVE_CONTROLLER.md',
       'code':'model/particle_candidate_atlas/attenuation.py'},
      {'id':'generic_positive_Y','claim':'full-CAR pairing norm 2tr(Y); complete feedback and specified Gaussian free-field limit',
       'proof':'model/particle_candidate_atlas/CHARGED_FIELD_ADMISSION.md',
       'code':'model/particle_candidate_atlas/charged_field_admission.py'},
      {'id':'quantum_Higgs','claim':'rank-one family vertex norm2; all12 finite Gauss commutators vanish; arbitrary common H0',
       'proof':'model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md',
       'code':'model/particle_candidate_atlas/interacting_chiral.py'},
      {'id':'native_diffusion','claim':'all4 duties at default p,q leak outside the C3 circulant algebra; 1536 original rows per duty',
       'proof':'model/particle_candidate_atlas/COMPATIBILITY.md',
       'code':'model/particle_candidate_atlas/memory_diffusion.py'},
      {'id':'native_purity','claim':'finite C-duty unread channel is mixed on an exact test state; full archive retained',
       'proof':'model/particle_candidate_atlas/NATIVE_RESPONSE.md',
       'code':'model/particle_candidate_atlas/native_response.py'},
    ]
    return {'retained_theorem_evidence':retained,
        'executed_bounded_controls':['exact matching floors and Higgs-unit intervals',
            'conditional SM anomaly, quadratic action, electroweak mass and field census identities',
            'all-three-input inverse controller and rational output error intervals',
            'capacity selfmap perturbation and complete inverse inequalities',
            'compact source pin checks and rational Schur/gap/original-weight bounds',
            'four original C-duty weights and purity lower bounds',
            'finite quantum-link identity obstruction and complete weak-kick budgets',
            'direct default-action exclusions from authenticated source coefficients'],
        'complete_native_experiments_replayed_here':False,
        'large_operator_certificates_replayed_here':False,
        'conditional_standard_model':sm,
        'all_three_input_controller':fitted,
        'measured_mass_join':mass,'capacity_nonselection':cap,'compact_defect_join':compact,
        'native_unread_response':native,
        'quantum_link_boundary':interacting_chiral.identity_link_control(),
        'complete_attenuator_budget':attenuation.budget(1024,'1/16'),
        'positive_Y_controller_budget':field.physical_controller_budget(10**6,'1','1/12',cells=1),
        'distinct_default_action_exclusions':field.default_action_exclusions(parent),
        'target_admissions':field.admission_rows(parent),
        'physical_particle_map':[],'physical_mass_predictions':[]}


def write():
    audit,svg=build()
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    OUTPUT.write_text(svg)
    AUDIT.write_text(canonical(audit))


def check():
    audit,svg=build()
    if OUTPUT.read_text()!=svg or AUDIT.read_text()!=canonical(audit):
        raise ValueError('atlas SVG or exact comparison audit differs from authenticated inputs')
