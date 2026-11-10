"""A source-bound candidate dictionary, with no fitted or identified masses."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from model.particle_common_contract import contract, provenance
from model.standard_model_source_map import algebra
from model.contact_scattering.collision import collision_certificate, execute as execute_contact
from . import finite_sector_ratios

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
INPUTS=(
    'model/protected_defect_catalog/catalog.json',
    'model/standard_model_source_map/particle_map.json',
    'model/particle_atlas_horizontal/families.json',
    'model/particle_primary_regional_fermion/results.json',
    'model/particle_primary_composite_spin/results.json',
    'model/particle_primary_reduced_sector/results.json',
    'model/particle_primary_loop_fusion/results.json',
    'model/moving_contact_scattering/results.json',
    'model/scattering_packet_readout/results.json',
)
SOURCE_INPUTS=tuple('model/standard_model_source_map/candidates/'+n for n in
    ('__init__.py','__main__.py','model.py','test_model.py','PROOF.md'))+INPUTS+finite_sector_ratios.SOURCE_INPUTS+(
        'model/contact_scattering/README.md','model/contact_scattering/collision.py',
        'model/contact_scattering/test_collision.py','model/moving_contact_scattering/PROOF.md',
        'model/scattering_packet_readout/PROOF.md')
REGIONAL_ID='particle_primary_regional_fermion:regional_endpoint_irrep'
PAIR_ID='particle_primary_regional_fermion:internal_pair_probe_irreps'
LOOP_ID='particle_primary_reduced_sector:isolated_A1u_cluster'
STATUSES=('possible','unsupported','excluded')
REFERENCES={
    'gauge_and_Higgs':'https://pdg.lbl.gov/2025/tables/rpp2025-sum-gauge-higgs-bosons.pdf',
    'leptons':'https://pdg.lbl.gov/2025/tables/rpp2025-sum-leptons.pdf',
    'quarks':'https://pdg.lbl.gov/2025/tables/rpp2025-sum-quarks.pdf',
}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticated_inputs():
    values={};pins={};cache={}
    for path in INPUTS:
        file=ROOT/path;value=json.loads(file.read_text());values[path]=value;pins[path]=sha(file)
        for name,expected in value.get('source_sha256',value.get('source_bindings',{})).items():
            p=(ROOT/name).resolve()
            if not p.is_relative_to(ROOT) or not p.is_file():raise ValueError('missing candidate ancestry: '+name)
            if name not in cache:cache[name]=sha(p)
            actual=cache[name]
            if actual!=expected:raise ValueError('stale candidate ancestry: '+name)
            pins[name]=expected
        def validate_binding(item):
            if isinstance(item,dict):
                if isinstance(item.get('result'),dict) and 'result_sha256' in item:
                    if contract.digest(item['result'])!=item['result_sha256']:
                        raise ValueError('candidate source-result binding differs')
                for child in item.values():validate_binding(child)
            elif isinstance(item,list):
                for child in item:validate_binding(child)
        validate_binding(value)
    return values,pins


def candidate_rows(values):
    old=values[INPUTS[1]];components=algebra.match_components()
    if components!=old['complete_internal_component_match']:
        raise ValueError('live exact internal weights differ from the frozen map')
    targets={t['id']:t for lane in values[INPUTS[2]]['lanes'] for t in lane['targets']}
    regional=values[INPUTS[3]]['members'][0]['result']
    if regional['endpoint_charge_increment']!=6 or regional['regional_endpoint_irrep']!='G1':
        raise ValueError('the primary charged regional endpoint changed')
    source=[]
    for name,c in components.items():
        weak=2 if c['q6Y'] in ('-3','1') else 1
        source.append({'id':name,'id_origin':'existing complete_internal_component_match key',
            'kind':'native_internal_component','protected_defect_id':None,
            'source_domain':'declared M15 exterior register; same native gauge representation in the primary Wilson source',
            'exact':{'source_masks':c['source_masks'],'SU3':c['color_representation'],
                     'SU2_dimension':weak,'Q6':c['q6Y'],'Y':str(F(c['q6Y'])/6),
                     'T3':c['T3'],'comparison_Q':c['Q']},
            'native_spin_statement':'primary Dirac spin fiber; regional Gauss-complete G1 records explicitly attached only for Q6=6',
            'physical_Lorentz_spin':None,'physical_flavor':None,
            'same_source_particle_sector_identified':False})
    source.extend([
        {'id':REGIONAL_ID,'id_origin':'existing regional receipt field',
         'kind':'actual_charged_regional_CAR_factor','protected_defect_id':None,
         'source_domain':'original primary member zero; complement and reference retained',
         'exact':{'SU3':'1','SU2_dimension':1,'Q6':6,'Y':'1','comparison_Q':'1',
                  'source_cube_spin':'G1','regional_CAR_parity':1,'complete_pair_CAR_parity':0},
         'internal_component_id':'charged_lepton_right_conjugate',
         'best_guess':'charged-lepton-singlet-type component; conjugate endpoint has opposite charge',
         'possible_named_targets':['electron','muon','tau'],'flavor_assignment':None,
         'source_readout':'Q6=±6 and chosen-axis source-spin increments±1/2; each of four raw joint rows has original pair-energy weight>1/16',
         'mass_readout_domain':'complete pair transfer, not an endpoint mass',
         'comparison_hypotheses':['Q=T3+Q6/6 is the physical current dictionary','finite G1 factor becomes a relativistic spin-one-half field','a chiral, flavor and dressed charged-threshold attachment'],
         'status':'possible','physical_particle_identified':False},
        {'id':PAIR_ID,'id_origin':'existing regional receipt field','kind':'complete_compensated_pair_probe',
         'protected_defect_id':None,'exact':{'total_Q6_increment':0,'CAR_parity':0,'internal_probe_irreps':{'A1':1,'T1':3}},
         'best_guess':'pair apparatus with retained charged endpoints; no single-photon assignment',
         'possible_named_targets':[],'status':'unsupported','physical_particle_identified':False},
    ])
    spin=values[INPUTS[4]]['selected_member_bound']['result']
    for family in spin['multiplet_families']:
        for signed in family['signed_multiplets']:
            m=signed['multiplet'];identifier='particle_primary_composite_spin:'+family['label']+':'+','.join(map(str,m['delta']))
            composite_hypothesis={
                'weak_pair_charge_six':'neutral lepton-composite or sterile-like channel from LL e^c; not an identified active neutrino',
                'color_triple':'neutral-antibaryon (antineutron-type) composite channel from u^c(d^c)^2, with the opposite increment the neutron-type conjugate',
                'color_determinant_charge_six_hole':'neutral even composite channel; no selected named-particle hypothesis',
            }[family['label']]
            if family['label']=='color_triple' and m['delta'][0]<0:
                composite_hypothesis='neutral-baryon (neutron-type) conjugate channel; physical baryon identity unproved'
            if family['label']=='weak_pair_charge_six' and m['delta'][1]<0:
                composite_hypothesis='conjugate neutral-lepton-composite or sterile-like channel; physical antiparticle identity unproved'
            source.append({'id':identifier,'id_origin':'existing multiplet family label and signed native increment',
                'kind':'actual_local_gauge_singlet_composite','protected_defect_id':None,
                'exact':{'native_occupation_increment':m['delta'],'SU3':'1','SU2_dimension':1,'Q6':0,
                         'CAR_parity':m['CAR_parity'],'source_cube_irreps':m['proper_binary_cube_irreps']},
                'best_guess':composite_hypothesis,
                'hypothesis_scope':'composite-family analogy conditional on physical current, baryon/lepton number and Lorentz-spin dictionaries; sign reverses the conjugate assignment',
                'possible_named_targets':[],'status':'unsupported','physical_particle_identified':False,
                'reason':'Gauge-singlet and finite cube labels do not establish active weak interactions, Lorentz spin or a neutrino mass/flavor field.'})
    line=values[INPUTS[5]]['selected_member_bound']['result'];spectrum=line['interacting_source_spectrum']
    source.append({'id':LOOP_ID,'id_origin':'existing reduced-sector receipt field',
        'kind':'actual_isolated_gauge_transition','protected_defect_id':None,
        'exact':{'CAR_parity':0,'native_gauge_representation':'singlet','source_cube_irrep':'A1u',
                 'transition_interval_source_units':spectrum['gauge_ground_to_A1u_energy_gap_interval'],
                 'transition_matrix_element_squared_lower':spectrum['actual_ground_to_odd_matrix_element_squared_strict_lower']},
        'best_guess':'finite neutral collective gauge line; no named photon, gluon, Z or Higgs assignment',
        'possible_named_targets':[],'status':'unsupported','physical_particle_identified':False,
        'reason':'This rare empty-matter-sector line has no physical momentum or Lorentz-spin dictionary; its original weight vanishes under the selected refinement.'})
    target_rows=[]
    boson_sources={'photon':['Q=T3+Y'],'gluon':['color_current_1..8'],
        'W_plus':['T_plus'],'W_minus':['T_minus'],'Z':['span(T3,Y)'],'Higgs':['declared_Higgs_doublet_radial_direction']}
    for row in old['rows']:
        name=row['target'];props=targets[name];ids=list(row['internal_components'])
        if ids:
            if 'neutrino' in name:guess='active-neutrino-type component of the native lepton doublet; flavor and mass mixing unresolved'
            elif name in ('electron','muon','tau'):guess='charged-lepton family components; no e/mu/tau assignment'
            elif name in ('up','charm','top'):guess='up-type quark family components; no generation assignment'
            else:guess='down-type quark family components; no generation assignment'
            missing=['physical chiral current and Lorentz-spin dictionary','one common interacting dressed particle/threshold or resonance sector','source-derived flavor and mass mixing','dispersion and quantitative complete detector matching']
            actual=[REGIONAL_ID] if name in ('electron','muon','tau') else []
        else:
            guess={'photon':'unbroken neutral gauge direction, conditional on its vacuum and mixing',
                   'gluon':'color-adjoint gauge-field candidate from eight native current directions',
                   'W_plus':'charged weak raising-current direction','W_minus':'charged weak lowering-current direction',
                   'Z':'neutral weak/hypercharge combination, with mixing unspecified',
                   'Higgs':'radial scalar of the declared Higgs doublet, with vacuum unspecified'}[name]
            ids=boson_sources[name];actual=[]
            missing=['source-selected field and vacuum/mixing attachment','physical spin-one or scalar propagation','complete particle/threshold or resonance detector','source-derived masses or proved masslessness']
        target_rows.append({'target':name,'comparison_status':'possible',
            'status_domain':'internal representation or current direction only',
            'benchmark':{k:props[k] for k in ('symbol','spin','charge','color','antiparticle')},
            'source_component_ids':ids,'actual_primary_excitation_ids':actual,'best_guess':guess,
            'source_flavor_index':None,'source_derived_mass_GeV':None,
            'physical_Lorentz_spin_proved':False,'physical_particle_identified':False,
            'protected_defect_id':None,'missing':missing})
    protected=[]
    for row in values[INPUTS[0]]['entries']:
        protected.append({'id':row['id'],'id_origin':'existing protected_defect_catalog entry',
            'kind':row['kind'],'source_domain':row['source'],'protected_object':row['protected_object'],
            'topology':row['topology'],'candidate_named_targets':[],
            'comparison_status':'excluded' if row['id']=='bare_native_source' else 'unsupported',
            'reason':('The complete source-only branch commutant is scalar in this specified bare domain.'
                      if row['id']=='bare_native_source' else 'This source charge or texture has no same-source Standard Model current/spin/detector attachment.'),
            'evidence':row['evidence']})
    exclusions=[
        {'source_id':REGIONAL_ID,'target':'photon','status':'excluded',
         'domain':'direct charge- and CAR-preserving component map under Q=T3+Q6/6',
         'reason':'The charged odd endpoint has comparison Q=+1 and G1 central sign−1; the photon benchmark is neutral and bosonic.'},
        {'source_id':PAIR_ID,'target':'photon','status':'unsupported',
         'domain':'identification of the full compensated pair as one propagating photon',
         'reason':'Even parity and zero total charge do not supply a massless transverse field, its dispersion or a single-particle detector.'},
        {'source_id':'particle_primary_composite_spin:weak_pair_charge_six','target':'active_neutrino','status':'excluded',
         'domain':'direct map of the local singlet to an unbroken lepton-doublet component with the same weak currents',
         'reason':'Weak singlet and weak doublet have different current action. A separate Higgs-dressed/current dictionary could evade this direct-map exclusion and is unproved.'},
        {'source_id':'particle_primary_composite_spin:color_triple','target':'quark','status':'excluded',
         'domain':'direct color-current-preserving elementary quark assignment',
         'reason':'The composite is a color singlet; the quark benchmark is a color triplet.'},
        {'source_id':LOOP_ID,'target':'Higgs','status':'unsupported',
         'domain':'named physical identity from one finite A1u gauge line',
         'reason':'A cubic total-state irrep does not determine intrinsic Lorentz spin. No radial-Higgs current, vacuum or physical pole dictionary is supplied.'},
        {'source_id':'charged_lepton_right_conjugate','target':'three_lepton_generations','status':'excluded',
         'domain':'deriving three flavors by repeating one internal component',
         'reason':'One Q6=6 weight and one regional spinor do not distinguish electron, muon and tau.'},
    ]
    return source,target_rows,protected,exclusions


def validate(value):
    if value['schema']!='oph.particle_candidate_correspondence.v1':raise ValueError('candidate schema differs')
    if value['physical_particle_map'] or value['mass_predictions_GeV'] or value['fit_parameters_used']:
        raise ValueError('a candidate dictionary cannot introduce physical assignments or fitted masses')
    targets=value['target_rows'];ids=[r['target'] for r in targets]
    expected={'electron','muon','tau','up','charm','top','down','strange','bottom',
              'electron_neutrino','muon_neutrino','tau_neutrino','photon','gluon','W_plus','W_minus','Z','Higgs'}
    if len(ids)!=18 or set(ids)!=expected:raise ValueError('all eighteen distinct target channels are required')
    for row in targets:
        if (row['comparison_status'] not in STATUSES or row['source_flavor_index'] is not None
                or row['source_derived_mass_GeV'] is not None or row['physical_particle_identified']
                or row['physical_Lorentz_spin_proved'] or row['protected_defect_id'] is not None):
            raise ValueError('a speculative row was promoted beyond its source evidence')
    if len({r['id'] for r in value['protected_catalogue_rows']})!=11:raise ValueError('existing protected catalogue IDs were lost')
    by={r['id']:r for r in value['source_excitation_rows']}
    if len(by)!=len(value['source_excitation_rows']):raise ValueError('duplicate source excitation ID')
    if by[REGIONAL_ID]['exact']['complete_pair_CAR_parity']!=0 or by[REGIONAL_ID]['exact']['regional_CAR_parity']!=1:
        raise ValueError('regional and global parity were conflated')
    if by[REGIONAL_ID]['flavor_assignment'] is not None or by[LOOP_ID]['possible_named_targets']:
        raise ValueError('unsupported flavor or gauge-line species assignment')
    if any(r['protected_defect_id'] is not None for r in by.values()):raise ValueError('a source excitation was silently relabeled as a protected topological defect')
    if len(value['internal_gauge_multiplets'])!=5:raise ValueError('internal multiplet table changed')
    for row in value['internal_gauge_multiplets']:
        e=by[row['component_ids'][0]]['exact']
        if any(row[k]!=e[k] for k in ('SU3','SU2_dimension','Y','Q6')):
            raise ValueError('multiplet display differs from the exact component weight')
    ratio=value['finite_sector_ratio_control']
    if contract.digest(ratio['result'])!=ratio['result_sha256'] or ratio['result']['physical_mass_table_GeV']:
        raise ValueError('finite threshold comparison was promoted or altered')
    collision=value['collision_status']
    if collision['moving_two_packet_collision_executed'] or collision['physical_scattering_cross_sections_calculated']:
        raise ValueError('static source fusion was promoted to physical scattering')
    if (not collision['separate_contact_extension']['finite_interaction_echo_executed']
            or not collision['separate_contact_extension']['small_32_row_native_control_executed']
            or collision['separate_contact_extension']['fine_spacing_packet_experiment_executed']):
        raise ValueError('executed contact controls and unexecuted fine-grid experiment were conflated')
    separate=collision['separate_contact_extension']
    if (separate['executed_interacting_echo']['source_resources']['native_pulses']!=22531
            or not separate['executed_interacting_echo']['all_outputs_retained']
            or F(separate['executed_free_echo']['pointer_probability_upper'])>F(1,10**24)):
        raise ValueError('executed contact echo changed its complete source control')


def contact_execution_controls():
    """Replay literal words; export stable admitted enclosures, not float snapshots."""
    result=[]
    for interaction in (1,0):
        report=execute_contact(2,interaction=interaction)['report']
        probability=report['pointer_probabilities'][0]
        if interaction and not .0469<probability<.047:
            raise ArithmeticError('executed contact interaction left its admitted interval')
        if not interaction and not probability<1e-24:
            raise ArithmeticError('executed free echo left its admitted interval')
        result.append({'interaction':interaction,'model_duration':'4','axis':2,
            'source_resources':{k:v for k,v in report['source_resources'].items() if not isinstance(v,float)},
            'pointer_probability_interval':['469/10000','47/1000'] if interaction else ['0','1/1000000000000000000000000'],
            'pointer_probability_upper':'47/1000' if interaction else '1/1000000000000000000000000',
            'floating_admission_tolerance':'1/50000000',
            'reference_Gram_error_within_declared_tolerance':True,
            'complete_vector_error_within_ideal_bound_plus_tolerance':True,
            'all_outputs_retained':report['all_outputs_retained'],
            'physical_hardware_error_guarantee':False})
    return result


def build():
    values,pins=authenticated_inputs();source,targets,protected,exclusions=candidate_rows(values)
    paths=provenance.definition_inputs(SOURCE_INPUTS)
    pins.update(provenance.definition_hashes(paths))
    ratio=finite_sector_ratios.certificate();pins.update(ratio['source_sha256'])
    fusion=values[INPUTS[6]]['selected_member_bound']['result']
    echo,free_echo=contact_execution_controls()
    by={r['id']:r for r in source};multiplets=[]
    for label,ids in (('Q_L',('up_left','down_left')),('u^c',('up_right_conjugate',)),
                      ('d^c',('down_right_conjugate',)),('L_L',('neutrino_left','charged_lepton_left')),
                      ('e^c',('charged_lepton_right_conjugate',))):
        e=by[ids[0]]['exact']
        multiplets.append({'comparison_label':label,'component_ids':list(ids),
                          'SU3':e['SU3'],'SU2_dimension':e['SU2_dimension'],
                          'Q6':e['Q6'],'Y':e['Y'],
                          'component_T3':[by[i]['exact']['T3'] for i in ids],
                          'component_Q':[by[i]['exact']['comparison_Q'] for i in ids],
                          'convention':'all-left Weyl comparison; physical chirality not proved'})
    value={'schema':'oph.particle_candidate_correspondence.v1',
        'map_kind':'best-guess internal correspondence with explicit source domains',
        'status_legend':{'possible':'exact internal/current match motivates this hypothesis; physical species unproved',
                         'unsupported':'available source facts do not select this named particle',
                         'excluded':'the precisely stated direct assignment violates an exact source property'},
        'source_excitation_rows':source,'target_rows':targets,'protected_catalogue_rows':protected,
        'internal_gauge_multiplets':multiplets,
        'rejection_controls':exclusions,'benchmark_references':REFERENCES,
        'finite_sector_ratio_control':ratio,
        'collision_status':{
            'source_receipt':INPUTS[6],
            'actual_static_loop_fusion_energy_separation_proved':fusion['actual_static_loop_fusion_energy_separation_proved'],
            'adjacent_fused_cluster_rank':fusion['full_complement_fusion_bounds']['interacting_fusion_cluster_rank'],
            'adjacent_fused_energy_interval':fusion['full_complement_fusion_bounds']['interacting_fusion_cluster_interval'],
            'opposite_disjoint_energy_window':fusion['full_complement_fusion_bounds']['disjoint_final_energy_window'],
            'each_complete_pair_row_original_weight_lower':fusion['each_complete_pair_row_original_weight_strict_lower_factored'],
            'source_domain':'complete sequential static Wilson-loop creation in the rare empty-matter reducing sector',
            'moving_two_packet_collision_executed':False,
            'moving_collision_false_flag_domain':'primary source and separate fine-grid packet program; finite contact echo and small native controls were executed',
            'physical_scattering_cross_sections_calculated':False,
            'separate_contact_extension':{
                'source_domain':'supplied Weyl/contact Hamiltonian and preparation; not the selected primary gauge source',
                'finite_interaction_echo_executed':True,
                'four_field_modes':4,'native_pulses':22531,
                'executed_interacting_echo':echo,'executed_free_echo':free_echo,
                'numerical_scope':'literal source execution with declared floating admission; not a sampled CERN histogram',
                'exact_ideal_word_signal_certificate':collision_certificate(),
                'small_32_row_native_control_executed':True,
                'small_native_controls':values[INPUTS[8]]['executed_small_native_controls'],
                'fine_spacing_packet_experiment_executed':values[INPUTS[8]]['finite_observation']['fine_spacing_native_experiment_executed'],
                'finite_complete_scattering_observation_theorem':values[INPUTS[8]]['scope']['finite_complete_scattering_observation_theorem'],
                'physical_particle_identification':values[INPUTS[8]]['scope']['physical_particle_identification'],
                'relativistic_interacting_QFT_constructed':values[INPUTS[8]]['scope']['relativistic_interacting_QFT_constructed'],
                'evidence':['model/contact_scattering/README.md',INPUTS[7],INPUTS[8]]},
            'remaining':['prepare separated moving incoming states in the physical sector','control asymptotic or finite detector scattering error','compute inclusive complete outgoing records and compare physical channels']},
        'physical_particle_map':[],'mass_predictions_GeV':[],'fit_parameters_used':False,
        'repeated_family_targets_are_independent_source_flavors':False,
        'gauge_weights_define_Lorentz_spin':False,'source_sha256':pins}
    validate(value)
    if provenance.definition_hashes(pins)!=pins:raise ValueError('candidate inputs changed during export')
    return value


def markdown(value):
    validate(value)
    text=['# Candidate source-to-particle correspondence','',
          'This is a best-guess comparison dictionary. Exact source charges and finite spin labels support family hypotheses; none of its rows is an identified physical particle or a fitted mass. Source excitation IDs and protected-defect catalogue IDs are separate.','',
          '## Exact internal gauge weights','',
          'The five patterns use the all-left-Weyl comparison convention, Y=Q6/6 and Q=T3+Y. These are internal representation weights, not a proof of physical chirality or Lorentz spin.','',
          '| Component | SU(3) | SU(2) dimension | Y | Q6 | T3 components | Q components |',
          '|---|---|---|---|---|---|---|']
    for row in value['internal_gauge_multiplets']:
        text+=['| '+row['comparison_label']+' | '+row['SU3']+' | '+str(row['SU2_dimension'])+' | '+row['Y']+' | '+row['Q6']+' | '+', '.join(row['component_T3'])+' | '+', '.join(row['component_Q'])+' |']
    text+=['','## Named elementary comparison targets','',
          '| Named target | Best source lead | Benchmark spin; charge; color | Status |',
          '|---|---|---|---|']
    for row in value['target_rows']:
        b=row['benchmark'];text.append('| '+row['target']+' | '+', '.join(row['source_component_ids'])+' | '+b['spin']+'; '+b['charge']+'; '+b['color']+' | Possible internal match; physical attachment unproved |')
    text+=['','The electron, muon and tau rows reuse the same source component IDs. The up-type, down-type and neutrino rows likewise share internal patterns. No flavor or generation is inferred from repetition. Neutrino flavor labels are not individual mass eigenstates.','',
           '## The strongest actual charged lead','',
           'The primary Q6=6 regional endpoint has the exact weak/color-singlet weight and G1 source spinor transformation of a charged-lepton-singlet-type candidate. Under the all-left-Weyl convention it is the charged-lepton conjugate component; source conjugation reverses its charge. Its opposite endpoint and Wilson link remain in the full experiment. The complete pair is even. Coherent charge and spin counters give original energy records, but their pair transfer is not a single endpoint mass. A charged lepton also needs its other chiral component, physical flavor, interactions and a controlled dressed threshold or resonance.','',
           '## Other source excitations','',
           'The neutral LL e^c trion motivates a neutral-lepton-composite or sterile-like hypothesis; it has no established active-neutrino current dictionary. The color-singlet u^c(d^c)^2 trion motivates a neutral antibaryon/neutron-type composite hypothesis, with its opposite increment conjugate. Baryon/lepton number, physical spin and named composite identity remain hypotheses. These analogies sit outside the eighteen elementary targets. The isolated A1u gauge transition near 63 source units is a rare neutral collective mode, with original probability retained. Its cubic label and gauge origin do not identify a photon, gluon, Z or Higgs. Native gauge-current directions motivate boson candidates, but do not supply their massless or massive propagation.','',
           '## Finite energies, ratios and collisions','',
           'Exact common-ground bounds put every ordered ratio of the three neutral sector thresholds below 21.16 and the two odd-family ratios below 15.31. Thus a direct common-scale identification of those member-zero thresholds with the charged-lepton hierarchy fails the measured muon/electron ratio 206.7682827(46). This excludes that particular assignment, not other excitations, capacities or source dictionaries. See the [exact rational bounds](FINITE_SECTOR_RATIOS.md). No GeV mass is inferred.','',
           'Two adjacent outward Wilson loops cancel a shared edge and create a six-edge loop in an isolated 44-dimensional energy cluster [1311/14,1335/14]. Opposite loops have a complete energy window [124,128]. Every recorded pair retains original weight greater than c_g/40; the rare sector factor is not divided out. This is a computed static fusion response, not an executed collision of moving particles or a CERN cross section. [Full source fusion proof](../../particle_primary_loop_fusion/PROOF.md).','',
           'A separate supplied Weyl/contact extension has an executed four-mode interaction echo at c=1/8, U=1, t=4: 22,531 native pulses, numerical pointer signal about 0.04691534, and an ideal-word signal greater than 0.0361998125 versus the exact free value zero. The rational bound excludes floating and hardware errors; literal numerical execution uses its declared tolerance. The moving-packet route has a complete finite observation theorem and executed 32-row five-qubit controls; its fine-grid incoming-packet program is not executed. These are different source/apparatus domains, not primary-source or physical CERN scattering. [Contact echo](../../contact_scattering/README.md); [moving construction](../../moving_contact_scattering/PROOF.md); [packet readout and controls](../../scattering_packet_readout/PROOF.md).','',
           '## Assignments ruled out or unsupported','']
    for row in value['rejection_controls']:
        text+=['- **'+row['status'].capitalize()+': '+row['source_id']+' → '+row['target']+'.** '+row['domain']+'. '+row['reason']]
    text+=['','## Protected source objects','',
           'The existing catalogue retains its exact IDs and source domains. None has a proved map to a named Standard Model particle. The bare complete source-only commutant excludes nontrivial protected charges in that specific domain; memory, extended sources and prescribed textures have different domains.','']
    for row in value['protected_catalogue_rows']:
        text+=['- `'+row['id']+'`: '+row['kind']+'. '+row['protected_object']+'.']
    text+=['','## Evidence and next physical calculations','',
           'A physical identification needs a common source state and clock, physical current/spin attachment, controlled spatial propagation and an appropriate charged threshold, inclusive sector or resonance. Pair costs, variational sector bounds and internal weights have different observable meanings. Source-selected flavor, scattering and decay records are required before comparing a complete Standard Model spectrum.','',
           '[Exact internal representation](../PROOF.md); [charged regional readout](../../particle_primary_regional_fermion/PROOF.md); [isolated gauge line](../../particle_primary_reduced_sector/PROOF.md); [protected catalogue](../../protected_defect_catalog/CATALOG.md).','',
           'Benchmark quantum numbers use the [PDG lepton tables]('+REFERENCES['leptons']+'), [quark tables]('+REFERENCES['quarks']+') and [gauge/Higgs tables]('+REFERENCES['gauge_and_Higgs']+'). They are comparison data, not source inputs or predictions.','']
    return '\n'.join(text)


def main():
    import argparse
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');args=parser.parse_args()
    value=build();outputs={HERE/'correspondence.json':json.dumps(value,indent=2,sort_keys=True,ensure_ascii=False)+'\n',HERE/'CANDIDATES.md':markdown(value)}
    for path,text in outputs.items():
        if args.write:path.write_text(text)
        elif not path.is_file() or path.read_text()!=text:raise ValueError('candidate export differs: '+str(path))
    print('PASS: scoped candidate source IDs and 18 target channels; no identified particle or fitted mass')


if __name__=='__main__':main()
