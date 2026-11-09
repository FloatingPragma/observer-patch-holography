"""Candidate source interfaces for eighteen targets, without physical assignments."""
from fractions import Fraction as F
import hashlib
import json

from tools import strict_json
from model.standard_model_source_map import algebra
from model.charged_lepton_response import calibrated_sanity

ROOT = calibrated_sanity.REFERENCE.parents[2]
PARENT = 'model/standard_model_source_map/particle_map.json'
DEFECTS = 'model/protected_defect_catalog/catalog.json'
SOURCE_INPUTS = (PARENT, DEFECTS, 'model/mcpr_native_memory_bridge/results.json',
                 'model/pixel_mass_closure/results.json',
                 'model/charged_lepton_response/calibration_reference.json')
TARGETS = ('electron','muon','tau','up','charm','top','down','strange','bottom',
           'electron_neutrino','muon_neutrino','tau_neutrino',
           'photon','gluon','W_plus','W_minus','Z','Higgs')
P_BRANCHES = ('forward','measured_endpoint','mssm_beta4_proxy','sm_color_beta7')
GROUPS = ((TARGETS[:3],('charged_lepton_left','charged_lepton_right_conjugate')),
          (TARGETS[3:6],('up_left','up_right_conjugate')),
          (TARGETS[6:9],('down_left','down_right_conjugate')),
          (TARGETS[9:12],('neutrino_left',)))


def _copy(value):
    return json.loads(json.dumps(value,sort_keys=True,allow_nan=False))


def _box(pair, *, positive=True):
    if type(pair) is not list or len(pair)!=2 or any(type(x) is not str for x in pair):
        raise ValueError('two exact rational endpoint strings required')
    a,b=map(F,pair)
    if a>b or positive and a<=0:
        raise ValueError('ordered positive source interval required')
    return a,b


def _triple(rows, *, ratios=False):
    if type(rows) is not list or len(rows)!=3:
        raise ValueError('exactly three source intervals required')
    intervals=[_box(x) for x in rows]
    if any(a[1]>=b[0] for a,b in zip(intervals,intervals[1:])):
        raise ValueError('strict ascending source modes required')
    if ratios and intervals[0]!=(F(1),F(1)):
        raise ValueError('source ratios must retain the smallest mode as one')
    return _copy(rows)


def calibrated_summary():
    """Exact chosen-phase intervals, without repeating the native or free-field check."""
    reference=strict_json.load(calibrated_sanity.REFERENCE)
    phase=calibrated_sanity.choose_phase(reference['muon_electron_ratio']['value'])
    return {'phase':phase,'comparison':calibrated_sanity.compare(reference,phase),
            'phase_calibration_is_source_derivation':False,
            'physical_particle_assignments':[],
            'full_native_or_free_field_replay_performed_here':False}


def _branches(calibrated,p_dependence):
    if calibrated.get('physical_particle_assignments')!=[] or calibrated.get(
            'phase_calibration_is_source_derivation') is not False:
        raise ValueError('calibration must not become a physical assignment or derivation')
    phase=calibrated['phase'];comparison=calibrated['comparison']
    x,y=F(phase['real']),F(phase['imaginary'])
    if x*x+y*y!=1 or y<=0 or phase.get('tau_used_to_choose_phase') is not False or phase.get(
            'phase_derived_from_repair') is not False:
        raise ValueError('one-input configured phase ancestry changed')
    if comparison.get('comparison_proves_physical_lepton_identity') is not False or comparison.get(
            'physical_loop_correction_bounds_supplied') is not False or comparison.get(
            'tau_reference_consumed_by_source_phase_or_kernel') is not False:
        raise ValueError('unproved physical or held-out calibration promotion')
    coherent={'kind':'calibrated_coherent','source_operator':'Y=C(z)^2/3',
              'source_state':'specified filled-negative-energy Gaussian reference for the free field',
              'source_clock':'declared weak-pulse and quadratic field clock; native overhead comparison separate',
              'input_ancestry':'muon/electron ratio selects z; electron scale is optional; tau is held out',
              'phase':{k:phase[k] for k in ('real','imaginary','half_angle_tangent',
                       'calibration_input_muon_electron_ratio')},
              'ascending_action_coefficients':_triple(comparison['source_action_coefficients']),
              'ascending_action_ratios':_triple(comparison['source_action_ratios'],ratios=True),
              'coefficient_unit':'dimensionless Y; free energies use |v_H| times its singular values',
              'exact_positive_root_Koide_ratio':'2/3',
              'source_phase_derived':False,'physical_particle_assignment_proved':False,
              'free_branch_theorem_owner':'model/charged_lepton_response/RELATIVISTIC_MASS.md',
              'full_CAR_owner':'model/charged_lepton_response/CAR_KOIDE.md'}
    if (p_dependence.get('source_parameters_fitted_by_this_certificate') is not False or
        p_dependence.get('experimental_values_consumed_by_source_evaluation') is not False or
        p_dependence.get('physical_interacting_mass_predictions')!=[] or
        p_dependence.get('primitive_finite_assumptions_force_this_capacity_law') is not False):
        raise ValueError('conditional pixel source ancestry or physical scope changed')
    rows=p_dependence['branches']
    if len(rows)!=4 or {r['branch'] for r in rows}!=set(P_BRANCHES):
        raise ValueError('all four distinct retained capacity branches required')
    variants=[]
    for name in P_BRANCHES:
        row=next(r for r in rows if r['branch']==name);source=row['source']
        if row.get('admitted_complete_memory_inverse_branches')!=2:
            raise ValueError('both complete inverse-memory branches required')
        for key in ('P','u','rho','delta','chi','positive_root_Koide_ratio'):
            _box(source[key])
        if F(source['positive_root_Koide_ratio'][1])>=F(2,3):
            raise ValueError('MCPR deformation must not be relabeled exact Koide balance')
        variants.append({'branch':name,'ancestry':row['ancestry'],
                         'source_coordinates':{k:_copy(source[k]) for k in ('P','u','delta','rho','chi')},
                         'ascending_shape_squares':_triple(source['shape_squares']),
                         'ascending_action_ratios':_triple(source['charged_action_ratios'],ratios=True),
                         'positive_root_Koide_ratio':_copy(source['positive_root_Koide_ratio']),
                         'complete_inverse_memory_branch_count':2})
    mcpr={'kind':'P_MCPR','source_operator':'g Z_C, a positive sum of P2,P3,P10 in the complete reset2 memory',
          'source_projector_components':['P2','P3','P10'],
          'individual_projector_is_a_mass_eigenmode':False,
          'source_state':'stationary complete twelve-state memory for the coefficient realization',
          'source_clock':'declared correlated weak-pulse clock; physical clock unproved',
          'input_ancestry':'declared historically target-informed MCPR law and separately specified capacity branch',
          'variants':variants,
          'coefficient_scale':'t_j=s_j*(exp(3*kappa)/(2*6^14*product(s)))^(1/3); kappa=(u/50)/(1+u/31+u^2/310)',
          'shape_square_is_absolute_mass_coefficient':False,
          'same_as_calibrated_coherent_action':False,'physical_particle_assignment_proved':False,
          'native_memory_owner':'model/mcpr_native_memory_bridge/PROOF.md',
          'capacity_owner':'model/pixel_mass_closure/PROOF.md'}
    return [coherent,mcpr]


def _defect_inventory(catalogue):
    if catalogue is None:
        return {'supplied':False,'entries':[],'species_links':[]}
    if catalogue.get('physical_particle_map')!=[] or catalogue.get('physical_mass_prediction_count')!=0:
        raise ValueError('existing defect inventory contains a changed physical admission')
    entries=catalogue['entries']
    if len(entries)!=11 or len({r['id'] for r in entries})!=11:
        raise ValueError('the separate eleven-object inventory is required')
    if any(r.get('identified_Standard_Model_species')!=[] or r.get('physical_mass') is not None or
           r.get('same_source_relativistic_particle_attachment_proved') is not False for r in entries):
        raise ValueError('a defect-to-species bridge was introduced without proof')
    return {'supplied':True,'entries':[
        {k:_copy(r[k]) for k in ('id','kind','source','protected_object','topology')}
        for r in entries], 'species_links':[],
        'exhaustive_over_all_possible_sources_or_defects':False}


def summarize(parent_map, calibrated, p_dependence, optional_defectinventory=None):
    """Pure candidate projection; no species is admitted by a matching name."""
    if parent_map.get('physical_particle_map')!=[] or parent_map.get('physical_mass_prediction_count')!=0:
        raise ValueError('the parent physical-admission boundary changed')
    parents=parent_map['rows']
    if len(parents)!=18 or {r['target'] for r in parents}!=set(TARGETS):
        raise ValueError('all eighteen distinct target channels required')
    components=algebra.match_components(parent_map['source_weights'])
    if components!=parent_map['complete_internal_component_match']:
        raise ValueError('the complete fifteen-weight dictionary changed')
    gauge=parent_map['gauge_and_Higgs_readout']
    if (gauge['independent_color_current_directions']!=8 or
        gauge['weak_charged_direction_charges']!=['1','-1'] or
        gauge['Higgs_vacuum_selected_by_source'] is not False or
        gauge['adjoint_directions_identified_with_propagating_bosons'] is not False):
        raise ValueError('internal gauge directions must not be promoted to physical bosons')
    by_target={r['target']:r for r in parents}
    expected={name:keys for names,keys in GROUPS for name in names}
    rows=[]
    for target in TARGETS:
        parent=by_target[target]
        if (any(parent.get(k) is not None for k in ('identified_OPH_particle_source',
                'identified_topological_defect','predicted_mass')) or
            any(parent.get(k) is not False for k in ('physical_particle_identified',
                'physical_pole_or_resonance_proved'))):
            raise ValueError('an unproved parent source/particle promotion was supplied')
        kind=('fermion_species' if target in TARGETS[:9] else
              'neutrino_flavor_channel' if target in TARGETS[9:12] else 'boson_species')
        if parent['target_kind']!=kind:
            raise ValueError('flavor channels, fermions and boson targets must remain distinct')
        keys=expected.get(target,())
        internal={k:components[k] for k in keys}
        if parent['internal_components']!=internal:
            raise ValueError('a target has the wrong exact internal-current components')
        row={'target':target,'target_kind':parent['target_kind'],
             'status':'candidate interface; physical identity unproved',
             'internal_components':_copy(internal),
             'candidate_name_is_a_hypothesis':True,
             'physical_particle_identified':False,'physical_mass':None,
             'physical_spin':None,'physical_lifetime':None,'physical_width':None,
             'topological_invariant':None,'identified_topological_defect':None,
             'defect_creation_map':None,'physical_particle_or_resonance_map':None,
             'missing_defect_bridge':'No source topology invariant, creation map and spectral-sector map connect this channel to a defect.',
             'candidate_sorted_charged_mode_index':None,'charged_source_branch_kinds':[],
             'source_mass_coefficient':None}
        if target in TARGETS[:3]:
            j=TARGETS.index(target)
            row.update({'candidate_sorted_charged_mode_index':j,
                'charged_source_branch_kinds':['calibrated_coherent','P_MCPR'],
                'source_object':'jth ascending spectral subspace of the selected charged-family action',
                'source_operator_interface':'odd CAR creator in the charged left/right-conjugate components, with the Higgs-line pairing',
                'source_mass_coefficient':'branch-specific intervals and common scales in charged_source_branches; never a physical mass',
                'missing_physical_links':['source/phase/controller selection or explicit calibration ancestry',
                    'selected interacting vacuum and chiral/mirror treatment',
                    'common physical clock and controlled kinetic/mass corrections',
                    'infrared threshold and observable match' if j==0 else 'unstable-particle resonance, width and observable match',
                    'independent test of this sorted-family species hypothesis']})
        elif target in TARGETS[3:9]:
            row.update({'source_object':'configured color-triplet and conjugate internal CAR component; family coordinate unspecified',
                'source_operator_interface':'same-current odd CAR creators and supplied Weyl/Yukawa menu',
                'missing_physical_links':['choose and justify one family action and mass/mixing basis',
                    'interacting gauge vacuum, chiral regulator and physical clock',
                    'interacting color-sector treatment and the matching convention for a quark mass',
                    'gauge-invariant hadronic or short-distance observable comparison']})
        elif target in TARGETS[9:12]:
            row.update({'source_object':'configured neutral component of the left weak doublet',
                'source_operator_interface':'odd CAR flavor channel; no right-neutrino mode in M15',
                'missing_physical_links':['neutrino mass mechanism and Majorana/Dirac assignment',
                    'flavor-to-mass mixing and oscillation response',
                    'common source state, clock and physical observable map']})
        else:
            objects={'photon':'Q=T3+Y direction conditional on a selected neutral Higgs vacuum',
                'gluon':'eight exact color-current directions; one adjoint target channel',
                'W_plus':'weak raising current T+ with [Q,T+]=T+',
                'W_minus':'weak lowering current T- with [Q,T-]=-T-',
                'Z':'neutral span(T3,Y), with a mixing coefficient to be selected',
                'Higgs':'radial response of a supplied neutral Higgs doublet direction'}
            missing={
                'photon':['selected neutral Higgs/gauge vacuum and unbroken physical Q',
                    'massless transverse radiation sector and physical clock',
                    'gauge-invariant electromagnetic observable comparison'],
                'gluon':['dynamical color gauge action and interacting vacuum',
                    'gauge-invariant color-sector observables; no colored free pole inferred',
                    'physical clock and independent strong-interaction comparison'],
                'W_plus':['selected broken-phase gauge action and Higgs vacuum',
                    'physical polarizations, charged resonance and decay width',
                    'common clock, current normalization and independent comparison'],
                'W_minus':['selected broken-phase gauge action and Higgs vacuum',
                    'physical polarizations, conjugate charged resonance and decay width',
                    'common clock, current normalization and independent comparison'],
                'Z':['selected neutral mixing, broken-phase action and vacuum',
                    'physical polarizations, neutral resonance and decay width',
                    'common clock and independent neutral-current comparison'],
                'Higgs':['selected stable vacuum and radial scalar kinetic response',
                    'scalar resonance, interactions and decay width',
                    'common clock and independent scalar-observable comparison']}
            row.update({'source_object':objects[target],
                'source_operator_interface':'finite current/doublet interface; propagation and physical spin are not inferred',
                'color_components':8 if target=='gluon' else 1,
                'missing_physical_links':missing[target]})
        rows.append(row)
    return {'schema':'oph.charged_lepton_response.candidate_map.v1',
            'target_count':18,'source_internal_weight_count':15,
            'rows':rows,'charged_source_branches':_branches(calibrated,p_dependence),
            'separate_defect_inventory':_defect_inventory(optional_defectinventory),
            'target_dictionary_reference':parent_map['target_dictionary_reference'],
            'target_names_enter_source_weight_calculation':False,
            'physical_particle_map':[],'physical_mass_predictions':[],
            'all_eighteen_targets_physically_identified':False,
            'topological_defect_species_assignment_derived':False,
            'chosen_action_candidates_exhaust_all_OPH_sources':False,
            'parent_scientific_experiments_replayed_by_projection':False}


def certificate():
    from . import p_dependence
    hashes=lambda: {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_INPUTS}
    before=hashes()
    parent=strict_json.load(ROOT/PARENT);defects=strict_json.load(ROOT/DEFECTS)
    result=summarize(parent,calibrated_summary(),p_dependence.certificate(),defects)
    if hashes()!=before:
        raise RuntimeError('a retained input changed during candidate projection')
    result['direct_input_sha256']=before
    return result
