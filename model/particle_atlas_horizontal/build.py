"""Reproducible horizontal SVG, pan/zoom reader and bounded illustration audit."""
import argparse
import hashlib
from html import escape
import json
from pathlib import Path
from urllib.parse import unquote
import xml.etree.ElementTree as ET

from . import drawing, measurement, endpoints, source_connections

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SVG = ROOT/'docs/research/OPH_PARTICLE_ATLAS_V2.svg'
HTML = ROOT/'docs/research/OPH_PARTICLE_ATLAS_V2.html'
AUDIT = HERE/'audit.json'
V1 = ROOT/'model/particle_candidate_atlas/audit.json'
OUTPUTS = ROOT/'model/particle_observer_outputs/results.json'
DEPENDENCIES = ROOT/'model/particle_joint_reconstruction/dependencies.json'
CANDIDATES = ROOT/'model/standard_model_source_map/candidates/correspondence.json'
IDENTITY_TESTS = ROOT/'model/particle_candidate_confirmation/results.json'
PRIMARY_PROOFS = tuple('model/'+owner+'/PROOF.md' for owner in (
    'particle_primary_field_algebra', 'particle_filtered_probe',
    'particle_primary_stress', 'particle_positive_window',
    'particle_primary_survival', 'particle_primary_mass_readout',
    'particle_primary_spatial_algebra', 'particle_primary_capacity',
    'particle_compact_mass_probe', 'particle_primary_momentum',
    'particle_primary_channel_labels', 'particle_primary_spectral_scale',
    'particle_primary_endpoint_bounds', 'particle_primary_variational_readout',
    'particle_primary_joint_spectroscopy', 'particle_primary_occupation_bound',
    'particle_primary_residual_refinement', 'particle_primary_response_moments',
    'particle_primary_pair_channels', 'particle_primary_refinement_response',
    'particle_primary_gravity_premises', 'particle_primary_branch_resolution',
    'particle_primary_energy_tightness', 'particle_primary_origin_momentum',
    'particle_primary_open_coframe', 'particle_primary_charged_string',
    'particle_primary_odd_residual', 'particle_primary_momentum_transfer',
    'particle_primary_local_fluctuation', 'particle_primary_spectral_resolution',
    'particle_primary_charged_readout', 'particle_primary_staggered_response',
    'particle_primary_current_tail', 'particle_primary_sector_threshold',
    'particle_primary_populated_momentum', 'particle_primary_weighted_survival',
    'particle_primary_charge_conjugation', 'particle_primary_composite_spin',
    'particle_primary_composite_transport', 'particle_primary_selected_odd_response',
    'particle_primary_covariant_stress', 'particle_primary_selected_spin_response',
    'particle_primary_odd_resolvent', 'particle_primary_odd_low_energy',
    'particle_primary_geometric_cut', 'particle_primary_coherent_odd_response',
    'particle_primary_original_shell_moments', 'particle_primary_original_branch', 'particle_primary_critical_refinement', 'particle_primary_boundary_information',
    'particle_primary_odd_sum_rules', 'particle_primary_branch_certificate', 'particle_primary_work_spectrum', 'particle_primary_boundary_correlation',
    'particle_primary_reduced_sector', 'particle_primary_regional_fermion', 'particle_primary_weighted_spectral_tail', 'particle_primary_energy_response',
    'particle_primary_loop_fusion', 'particle_primary_equal_work_control'))
PRIMARY_RESPONSE_RECEIPTS = tuple('model/'+owner+'/results.json' for owner in (
    'particle_primary_branch_resolution', 'particle_primary_energy_tightness',
    'particle_primary_origin_momentum', 'particle_primary_open_coframe',
    'particle_primary_charged_string', 'particle_primary_odd_residual',
    'particle_primary_momentum_transfer', 'particle_primary_local_fluctuation',
    'particle_primary_spectral_resolution', 'particle_primary_charged_readout',
    'particle_primary_staggered_response', 'particle_primary_current_tail',
    'particle_primary_sector_threshold', 'particle_primary_populated_momentum',
    'particle_primary_weighted_survival', 'particle_primary_charge_conjugation', 'particle_primary_composite_spin',
    'particle_primary_composite_transport', 'particle_primary_selected_odd_response',
    'particle_primary_covariant_stress', 'particle_primary_selected_spin_response',
    'particle_primary_odd_resolvent', 'particle_primary_odd_low_energy',
    'particle_primary_geometric_cut', 'particle_primary_coherent_odd_response',
    'particle_primary_original_shell_moments', 'particle_primary_original_branch', 'particle_primary_critical_refinement', 'particle_primary_boundary_information',
    'particle_primary_odd_sum_rules', 'particle_primary_branch_certificate', 'particle_primary_work_spectrum', 'particle_primary_boundary_correlation',
    'particle_primary_reduced_sector', 'particle_primary_regional_fermion', 'particle_primary_weighted_spectral_tail', 'particle_primary_energy_response',
    'particle_primary_loop_fusion', 'particle_primary_equal_work_control'))
CONTENT = ('story.json','families.json','measurement_content.json','theorems.json','field_nature.json')


def canonical(value):
    return json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+'\n'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def links(value):
    if isinstance(value,dict):
        if 'path' in value and 'label' in value:yield value
        for v in value.values():yield from links(v)
    elif isinstance(value,list):
        for v in value:yield from links(v)


def scientific_receipt_records():
    """Authenticate frozen consequence receipts without replaying their experiments."""
    records=[];paths=set(PRIMARY_RESPONSE_RECEIPTS)
    def bindings(value):
        if isinstance(value,dict):
            if 'result_sha256' in value and isinstance(value.get('result'),dict):
                result=value['result']
                digest=hashlib.sha256(json.dumps(result,sort_keys=True,
                    separators=(',',':'),allow_nan=False).encode()).hexdigest()
                if digest!=value['result_sha256'] or not result.get('member_id'):
                    raise ValueError('primary response result binding differs')
                yield {'member_id':result.get('member_id'),
                       'result_sha256':value['result_sha256']}
            for child in value.values():yield from bindings(child)
        elif isinstance(value,list):
            for child in value:yield from bindings(child)
    for name in PRIMARY_RESPONSE_RECEIPTS:
        value=json.loads((ROOT/name).read_text())
        pins=value['source_sha256']
        for path,expected in pins.items():
            resolved=(ROOT/path).resolve()
            if not resolved.is_relative_to(ROOT) or not resolved.is_file() or sha(resolved)!=expected:
                raise ValueError('stale primary response input: '+path)
        paths.update(pins)
        # Validate every embedded copy before deduplicating exact bindings.
        # A repeated parity binding is one result, not an additional experiment.
        verified=list(bindings(value));distinct=[];seen=set()
        for row in verified:
            key=(row['member_id'],row['result_sha256'])
            if key not in seen:distinct.append(row);seen.add(key)
        records.append({'path':name,'schema':value['schema'],
                        'receipt_sha256':sha(ROOT/name),'source_pin_count':len(pins),
                        'member_result_bindings':distinct,
                        'scientific_experiments_replayed_by_atlas':False})
    return records,paths


def load():
    data={'v1':json.loads(V1.read_text()),'measurement':measurement.certificate(),
          'endpoints':endpoints.certificate(), 'double_path':measurement.double_path_certificate()}
    for name in CONTENT:data[Path(name).stem]=json.loads((HERE/name).read_text())
    data['joint_dependencies']=json.loads(DEPENDENCIES.read_text())
    data['observer_outputs']=json.loads(OUTPUTS.read_text())
    data['primary_response_receipt_records'],response_paths=scientific_receipt_records()
    data['candidate_correspondence']=json.loads(CANDIDATES.read_text())
    candidate=data['candidate_correspondence']
    data['candidate_identity_tests']=json.loads(IDENTITY_TESTS.read_text())
    from model.particle_candidate_confirmation.model import read_inputs as identity_inputs
    data['_candidate_identity_authenticated']=identity_inputs()
    for path,expected in data['candidate_identity_tests']['source_sha256'].items():
        resolved=(ROOT/path).resolve()
        if not resolved.is_relative_to(ROOT) or not resolved.is_file() or sha(resolved)!=expected:
            raise ValueError('stale candidate identity input: '+path)
    from model.standard_model_source_map.candidates.model import validate as validate_candidates
    validate_candidates(candidate)
    for path,expected in candidate['source_sha256'].items():
        resolved=(ROOT/path).resolve()
        if not resolved.is_relative_to(ROOT) or not resolved.is_file() or sha(resolved)!=expected:
            raise ValueError('stale candidate map input: '+path)
    for rel,expected in data['observer_outputs']['input_sha256'].items():
        if sha(ROOT/rel)!=expected:raise ValueError('stale observer-output input: '+rel)
    pins=data['v1']['input_sha256']
    for rel,expected in pins.items():
        path=(ROOT/rel).resolve()
        if not path.is_relative_to(ROOT) or sha(path)!=expected:
            raise ValueError('stale v1 scientific input: '+rel)
    if data['v1']['physical_mass_predictions'] or data['v1']['physical_particle_map']:
        raise ValueError('v1 physical scope changed; editorial review required')
    validate_content(data)
    paths=set(pins)|set(data['observer_outputs']['input_sha256'])|{str(V1.relative_to(ROOT)),str(OUTPUTS.relative_to(ROOT))}
    paths.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                 if p.is_file() and (p.suffix in ('.py','.md','.js') or p.name in CONTENT))
    paths.update({str(DEPENDENCIES.relative_to(ROOT)),
                  'model/particle_joint_reconstruction/DEPENDENCIES.md',
                  data['joint_dependencies']['conditional_candidate']['proof'],
                  data['joint_dependencies']['conditional_candidate']['producer']})
    paths.update(PRIMARY_PROOFS)
    paths.update(response_paths)
    paths.update(candidate['source_sha256'])
    paths.update(data['candidate_identity_tests']['source_sha256'])
    paths.update((str(IDENTITY_TESTS.relative_to(ROOT)),
                  'model/particle_candidate_confirmation/RESULTS.md'))
    paths.update((str(CANDIDATES.relative_to(ROOT)),
                  'model/standard_model_source_map/candidates/CANDIDATES.md'))
    for node in data['joint_dependencies']['nodes']:
        for evidence in node['evidence']:
            name=evidence['path']
            path=(ROOT/name).resolve()
            if not path.is_relative_to(ROOT) or not path.is_file():
                raise ValueError('missing or escaped dependency evidence: '+name)
            paths.add(name)
    for link in links([data[n.removesuffix('.json')] for n in CONTENT]):
        name=unquote(link['path'].split('#')[0])
        if name.startswith('https://'):continue
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError('missing or escaped reference: '+name)
        paths.add(name)
    data['input_sha256']={p:sha(ROOT/p) for p in sorted(paths)}
    return data


def validate_content(data):
    source_connections.validate(data['joint_dependencies'])
    ids=[t for lane in data['families']['lanes'] for t in lane['target_ids']]
    expected={r['target'] for r in data['v1']['candidate_rows']}
    if len(ids)!=18 or len(set(ids))!=18 or set(ids)!=expected:
        raise ValueError('complete elementary target coverage failed')
    from model.particle_candidate_confirmation.model import validate as validate_identity
    identity=data['candidate_identity_tests']
    validate_identity(identity,authenticated=data.get('_candidate_identity_authenticated'))
    comparison={r['target']:r for r in data['candidate_correspondence']['target_rows']}
    if {r['target'] for r in identity['target_rows']}!=expected:
        raise ValueError('candidate identity tests lost an elementary target')
    for row in identity['target_rows']:
        target=comparison[row['target']]
        if any(row[key]!=target[key] for key in ('source_component_ids','actual_primary_excitation_ids')):
            raise ValueError('candidate identity tests changed the comparison source assignment')
    window=identity['calculated_source_results']['charged']['original_compensated_dressed_energy_window']
    if (not window['same_original_source_state_reference_clock']
            or window['energy_window_is_a_single_endpoint_or_rest_mass']
            or window['physical_charged_lepton_identified']
            or window['complete_probe_outcomes']!=window['complete_success_outcomes']+1
            or window['additional_clock_loss_subtracted']!=0):
        raise ValueError('dressed pair response changed its complete source-record domain')
    named=data['field_nature']['families']
    if len(named)!=18 or {r['target'] for r in named}!=expected:
        raise ValueError('named-field coverage failed')
    if data['field_nature']['physical_particle_assignments'] or any(
            r['source_predicted_physical_mass'] is not None for r in named):
        raise ValueError('named-field explanation promoted a physical prediction')
    if data['observer_outputs']['physical_particle_assignments'] or data['observer_outputs']['source_predicted_GeV_masses']:
        raise ValueError('finite output explanation promoted a physical prediction')
    expected_stages=('loop','charge','local_bridge','texture','mode','field','apparatus','record')
    if tuple(s['id'] for s in data['story']['stages'])!=expected_stages:
        raise ValueError('hero identity chain changed')
    edges=data['story']['edges']
    if [(r['from'],r['to']) for r in edges]!=list(zip(expected_stages,expected_stages[1:])):
        raise ValueError('hero edges must cover the same sequence')
    if any(edges[i]['status']!='open' for i in (1,2,4)):
        raise ValueError('unproved physical bridge promoted')
    for lane in data['families']['lanes']:
        if len(lane['stages'])!=4 or any(t['source_physical_identification'] for t in lane['targets']):
            raise ValueError('family physical scope changed')
    for row in data['families']['calibration_paths']:
        if not all(row.get(k) for k in ('input','removal_strategy','success_condition','links')):
            raise ValueError('calibration path lost its required input or test')


def reader(svg,data):
    full_links=[];seen=set()
    extra=[{'label':label,'path':'model/particle_observer_outputs/'+name} for name,label in (
        ('HOLONOMY_TRANSPORT.md','Finite holonomy birth, work exchange and regional protection'),
        ('CONSENSUS_RECORDS.md','Record copying, observer meanings and finite event gluing'),
        ('BARE_EFFECTIVE_BOUNDARY.md','Same-clock spectral coefficient inverse'),
        ('CONTROLLER_READOUT.md','Controller attenuation and phase-class readout'),
        ('GAUGE_INPUT.md','Operational source couplings and conditional gauge ray'))]
    extra += [{'label':label, 'path':'model/'+owner+'/PROOF.md'} for owner,label in (
        ('particle_defect_attachment','Protected memory and holonomy coupled to complete field records'),
        ('particle_field_causality','Same-field causal influence and complete local-history bounds'),
        ('particle_gravity_response','Same-action energy transfer, work and geometric premises'))]
    extra += [{'label':label, 'path':'model/'+owner+'/PROOF.md'} for owner,label in (
        ('particle_scalar_ground','Full-carrier scalar ground, singlet preparation and complete occupation records'),
        ('particle_regional_geometry','Induced regional modular pair interaction and entropy response'),
        ('particle_field_preparation','Finite gauge-word preparation and same-clock Wilson-gap bounds'),
        ('particle_ground_charged','Complete charged spectra about the full scalar ground'),
        ('particle_prepared_field_readout','One prepared spatial field with all 61 local detector outcomes'),
        ('particle_local_gauge_state','Mixed gauge-state family with fixed-parameter local roughness bounds'),
        ('particle_local_wilson_density','Local Wilson low-mode weight and one-body sign approximation'),
        ('particle_local_history_family','Complete polynomial field records at fixed nonzero time and growing volume'),
        ('particle_sign_response_boundary','Native inverse-gap electric response and the spectral shortcut boundary'),
        ('particle_smooth_sign_field','Smooth interacting action with complete local quantum records'),
        ('particle_joint_translation_readout','Joint field/reference translations, position backaction and finite uncertainty'),
        ('particle_physical_causality','Refinement-uniform physical causal bound for complete local records'),
        ('particle_field_state_response','Complete interacting thermal preparation and flux records'),
        ('particle_geometric_response','Same-action spatial strain, work and pressure with declared geometry'),
        ('particle_spatial_ground_state','Gap-free full spatial low-energy preparation'),
        ('particle_ground_spectral_records','Charge-conjugate thresholds and complete source-energy records'),
        ('particle_gauge_characteristics','Native trace-dual gauge and matter characteristics'),
        ('particle_common_contract','One selected source, state, reference and clock'),
        ('particle_gauss_ground_state','Physical Gauss preparation on the full spatial carrier'),
        ('particle_charged_residue','Original charged-threshold weights'),
        ('particle_spatial_modular_response','Physical regional support and modular response'),
        ('particle_link_cutoff_limit','Fixed-volume electric cutoff and complete-record convergence'),
        ('particle_physical_boundary_state','Physical preparation with boundary charges and complete records'),
        ('particle_primary_locality','Complete local boundary records on the selected physical source'),
        ('particle_primary_field_algebra','Common filtered observable algebra, positive energy and finite pointer limits'),
        ('particle_primary_survival','Actual primary finite contrast and band; shrinking times and nonuniform weight'),
        ('particle_primary_mass_readout','Original transfer-to-excitation weights, rest-energy conditions and declared SI units'),
        ('particle_primary_spatial_algebra','Continuous regular spatial action; selected-state invariance and momentum conditions'),
        ('particle_primary_capacity','Source capacity in nats: ambient and physical bounds, references and execution cost'),
        ('particle_compact_mass_probe','Four-outcome readout of the same finite filter; two apparatus qubits and explicit service costs'),
        ('particle_primary_momentum','Complete spatial Fourier records, open-aperture error and the conserved-momentum conditions'),
        ('particle_primary_channel_labels','Five conserved occupation labels and a finite cubic scalar response with original weight'),
        ('particle_primary_spectral_scale','Exact decoupled source and local perturbative coefficients; no primary endpoint gap inference'),
        ('particle_primary_endpoint_bounds','Actual interacting ground lower bound and the unproved branch-isolation condition'),
        ('particle_primary_variational_readout','Eleven physical endpoint trial vectors, a ground upper bound and the full-space residual'),
        ('particle_primary_joint_spectroscopy','Complete initial/final energy and Fourier records with authenticated ground-energy intervals'),
        ('particle_primary_occupation_bound','Occupation-sensitive electric bound and substantial original neutral probe weight'),
        ('particle_primary_residual_refinement','Exact local Haar and CAR moments control the full omitted-space residual'),
        ('particle_primary_response_moments','Actual neutral pair energy band with original probability and complete finite clocks'),
        ('particle_primary_pair_channels','Populated cubic/parity channels and complete spatial pair records; no particle-spin assignment'),
        ('particle_primary_refinement_response','Scalable populated collective response and controlled finite-source defect density'),
        ('particle_primary_gravity_premises','Polynomial modular work, source-clock comparison and explicit geometric premises'),
        ('particle_primary_branch_resolution','Concentrated original pair bands, fixed-time correlation sign and the unresolved odd-sector trial'),
        ('particle_primary_energy_tightness','Exact collective moments and narrower finite energy bands; fixed-energy survival unproved'),
        ('particle_primary_origin_momentum','Retained common origins and commuting energy/momentum on the relative-observation domain'),
        ('particle_primary_open_coframe','Boundary-completed spatial strain with exact flat primary action and complete mean records'),
        ('particle_primary_charged_string','Opposite physical hypercharge endpoints with original selected-state energy weight'),
        ('particle_primary_odd_residual','An actual odd excitation enclosure and a separate coherent original-state response'),
        ('particle_primary_momentum_transfer','Complete coherent-offset records and joint energy/momentum sinc kernels on the relative-observation domain'),
        ('particle_primary_local_fluctuation','Bounded local unitary fluctuations with original energy tightness and an allowed elastic limit'),
        ('particle_primary_spectral_resolution','Prepared physical A1u energy window with original probability above one quarter and retained displaced state'),
        ('particle_primary_charged_readout','Coherent five-state meters join opposite charge transfers to original complete energy records'),
        ('particle_primary_staggered_response','Actual Wilson cancellation and populated coherent A2g energy bands on the selected state'),
        ('particle_primary_current_tail','Energy-form control joins unfiltered phase probes and filtered fields on common positive-energy time dynamics'),
        ('particle_primary_sector_threshold','Positive finite occupation-sector thresholds and complete prepared energy windows with compensating memory'),
        ('particle_primary_populated_momentum','Populated scanned spatial/energy records with actual Wilson phase dependence and original union weights'),
        ('particle_primary_weighted_survival','A necessary normalized observation strength for fixed positive-energy records and its critical-Weyl boundary'),
        ('particle_primary_charge_conjugation','Exact full-source particle/hole symmetry, opposite-sector spectra and conjugate complete detector laws'),
        ('particle_primary_composite_spin','Actual odd gauge-singlet G1 and H multiplets, native currents, graded exchange and complete prepared energy records'),
        ('particle_primary_composite_transport','Literal full-source third-order composite transport across complete six-state spin multiplets'),
        ('particle_primary_selected_odd_response','Original selected-state odd creation weight, complete clocks and a uniform cutoff-scale escaping energy band'),
        ('particle_primary_covariant_stress','Cube-covariant source energy and current, exact boundary improvement and conjugate regional response'),
        ('particle_primary_selected_spin_response','Original G1/H probe-energy records, exact twelve-mode CAR residuals and coherent native-number meters'),
        ('particle_primary_odd_resolvent','Actual coefficient-one full sector complement resolvent, prepared window and finite original final-energy probability'),
        ('particle_primary_odd_low_energy','Complete finite odd filters, exact selector energy laws and necessary dressing depth'),
        ('particle_primary_geometric_cut','Covariant induced regional state, actual boundary-representation entropy and contained energy'),
        ('particle_primary_coherent_odd_response','Positive G1/H energy records in every complete coherent Walsh spatial row on the original state'),
        ('particle_primary_original_shell_moments','Exact original even-density spectral moments, native hopping signs and complete finite-clock contrast'),
        ('particle_primary_original_branch','Polynomial detailed balance, original total-cube/current branches and non-single-line odd transfer measures'),
        ('particle_primary_critical_refinement','Controlled three-dimensional phase-work limit and source-independent profile-work ratios'),
        ('particle_primary_boundary_information','Original boundary-center coherence, complete recording work and regional information identity'),
        ('particle_primary_odd_sum_rules','Exact native chirality parity pairing, odd spectral sum rules and original graded propagation'),
        ('particle_primary_branch_certificate','Admitted43-coordinate sector upper, original low-energy records and exact symmetry reduction'),
        ('particle_primary_work_spectrum','Original modular fluctuation and work records, stronger moment bounds and fixed-energy survival boundary'),
        ('particle_primary_boundary_correlation','Actual boundary-interaction bound on original conditional correlations and center coherence'),
        ('particle_filtered_probe','Complete finite filtering apparatus and original Ramsey responses'),
        ('particle_primary_stress','Primary local energy currents and source-marginal entropy response'),
        ('particle_positive_window','Original positive-energy window weight from complete records'),
        ('particle_yang_mills_gap_contract','Same-clock local correlations and a continuum spectral-gap criterion'))]
    extra += [{'label':label, 'path':'model/particle_joint_reconstruction/'+name} for name,label in (
        ('DEPENDENCIES.md','Mass-first dependency graph with deferred pure-gauge and gravity branches'),
        ('CANDIDATE.md','Conditional chiral gauge-Higgs candidate and calibrated spectrum conventions'))]
    extra.append({'label':'All-time gap control from conserved energy and its refinement boundary',
                  'path':'model/particle_field_preparation/REFINEMENT.md'})
    extra += [{'label':'Best-guess named-particle map and scoped source IDs','path':'model/standard_model_source_map/candidates/CANDIDATES.md'},
              {'label':'Identity tests for all 18 named comparison targets','path':'model/particle_candidate_confirmation/RESULTS.md'},
              {'label':'Native mass coefficients, scalar dressing and complete charged records','path':'model/particle_candidate_charged/PROOF.md'},
              {'label':'Neutral composite current selection and coherent occupation-increment records','path':'model/particle_candidate_composites/PROOF.md'},
              {'label':'Gauge-loop scalar/vector response and exact physical current exclusions','path':'model/particle_candidate_bosons/PROOF.md'},
              {'label':'Exact finite-sector ratio exclusion','path':'model/standard_model_source_map/candidates/FINITE_SECTOR_RATIOS.md'},
              {'label':'Static Wilson-loop fusion on the primary source','path':'model/particle_primary_loop_fusion/PROOF.md'},
              {'label':'Same-source exact work with an elastic spectral limit','path':'model/particle_primary_equal_work_control/PROOF.md'}]
    for item in links([data[k] for k in ('story','families','measurement_content','theorems','field_nature')]+[extra]):
        href=item['path']+('#'+item['anchor'] if item.get('anchor') else '')
        if href in seen:continue
        seen.add(href)
        url=href if href.startswith('https://') else '../../'+href
        full_links.append(f'<li><a href="{escape(url,quote=True)}" target="_blank">{escape(item["label"])}</a></li>')
    return '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>OPH particle atlas v2: illustrated reader</title>
<style>*{box-sizing:border-box}body{margin:0;height:100vh;display:flex;flex-direction:column;overflow:hidden;background:#fbf8f0;color:#142c3c;font:15px Arial,sans-serif}header{min-height:76px;flex-shrink:0;background:#142c3c;color:white;display:flex;gap:10px;align-items:center;padding:12px 18px;flex-wrap:wrap;position:relative;z-index:3}header strong{font:22px Georgia,serif;margin-right:12px}button,a.download{border:1px solid #789ba7;background:#244454;color:white;border-radius:6px;padding:9px 13px;cursor:pointer;text-decoration:none;font-size:14px}button:hover,a.download:hover{background:#397078}#hint{font-size:12px;color:#b9d2db;margin-left:auto}#viewport{flex:1;min-height:0;overflow:hidden;touch-action:none;cursor:grab;position:relative}#art{transform-origin:0 0;position:absolute;left:0;top:0}#art svg{display:block;max-width:none}#zoom{min-width:48px;text-align:center}#refs{display:none;position:fixed;right:0;top:0;width:min(560px,90vw);height:100vh;overflow:auto;background:#fffaf0;border-left:1px solid #ccd5d5;z-index:4;padding:22px}#refs a{color:#285b9f}#refs li{margin:11px 0}#refs h2{font:26px Georgia,serif}#nav{position:fixed;bottom:20px;right:22px;z-index:3;display:flex;gap:8px}#nav button{box-shadow:0 3px 12px #0002}@media(max-width:850px){header{min-height:118px;max-height:40vh;overflow:auto;align-content:flex-start}#hint{display:none}}@media print{header,#nav,#refs{display:none}#viewport{overflow:visible;height:auto}#art{position:static;transform:none!important}#art svg{width:100%;height:auto}}</style>
<header><strong>Simulation-generated observables</strong><button data-row="origin">Defect origins</button><button data-row="story">Derivation chain</button><button data-row="observers">Observers & meaning</button><button data-row="experiment">The experiment</button><button data-row="theorems">Theorem bridges</button><button data-row="endpoints">Endpoints & slits</button><button data-row="families">Every family</button><button data-row="particles">Mass references</button><button data-row="bestguess">Best-guess map</button><button data-row="comparisons">Masses & collisions</button><button data-row="fieldnature">Photon, Higgs & matter</button><button data-row="requirements">Ten dependencies</button><button data-row="readouts">Actual field readouts</button><button data-row="contract">One program</button><button data-row="ground">Ground states &amp; cones</button><button data-row="refinement">Refinement & pressure</button><button data-row="calibrations">Remove calibrations</button><button id="overview">Overview</button><button id="references">Proofs & papers</button><a class="download" href="OPH_PARTICLE_ATLAS_V2.svg" download>SVG ↓</a><span style="flex-basis:100%;font-size:13px;color:#f2d79a">Named fields: photon · Higgs · W and Z · gluons · electron, muon, tau · quarks · neutrinos. Gold boxes show comparison particles; arrows show the construction and its proof boundaries.</span><span id="hint">Drag to pan · wheel to zoom · click blue proof links</span></header>
<main id="viewport" aria-label="Draggable and zoomable particle atlas"><div id="art">'''+svg+'''</div></main><aside id="refs"><button id="closeRefs">Close</button><h2>Exact proofs and paper sections</h2><p>Full references behind the illustrated steps. Supplied hypotheses and unproved bridges are labeled in the spread.</p><ul>'''+''.join(full_links)+'''</ul></aside>
<nav id="nav"><button id="left" aria-label="Pan left">←</button><button id="right" aria-label="Pan right">→</button><button id="out" aria-label="Zoom out">−</button><button id="in" aria-label="Zoom in">+</button><button id="zoom">100%</button></nav>
<script>
const viewport=document.querySelector('#viewport'),art=document.querySelector('#art');let scale=.25,x=0,y=0,drag=null;
function draw(){art.style.transform=`translate(${x}px,${y}px) scale(${scale})`;document.querySelector('#zoom').textContent=Math.round(scale*100)+'%'}
function zoom(f,cx=viewport.clientWidth/2,cy=viewport.clientHeight/2){const b=art.querySelector('svg').viewBox.baseVal;const floor=Math.min(.035,.97*viewport.clientWidth/b.width,.97*viewport.clientHeight/b.height);const next=Math.min(1.8,Math.max(floor,scale*f));x=cx-(cx-x)*next/scale;y=cy-(cy-y)*next/scale;scale=next;draw()}
function row(name){const places={origin:[0,0],observers:[38000,0],story:[10320,680],experiment:[10320,2320],theorems:[10320,4500],families:[10320,6150],particles:[26800,6150],calibrations:[32000,6150],endpoints:[10320,11500],fieldnature:[10320,13900],requirements:[10320,25400],readouts:[10320,29200],refinement:[10320,33800],ground:[10320,38950],contract:[10320,43850],bestguess:[0,30400],comparisons:[38000,30400]};scale=Math.min(.43,viewport.clientWidth/(['origin','observers','bestguess','comparisons'].includes(name)?9900:5600));x=35-places[name][0]*scale;y=35-places[name][1]*scale;draw()}
viewport.addEventListener('wheel',e=>{e.preventDefault();const r=viewport.getBoundingClientRect();zoom(Math.exp(-e.deltaY*.0015),e.clientX-r.left,e.clientY-r.top)},{passive:false});
viewport.addEventListener('pointerdown',e=>{if(e.target.closest('a'))return;drag={px:e.clientX,py:e.clientY,x,y};viewport.setPointerCapture(e.pointerId);viewport.style.cursor='grabbing'});
viewport.addEventListener('pointermove',e=>{if(!drag)return;x=drag.x+e.clientX-drag.px;y=drag.y+e.clientY-drag.py;draw()});
function stop(){drag=null;viewport.style.cursor='grab'}viewport.addEventListener('pointerup',stop);viewport.addEventListener('pointercancel',stop);
document.querySelectorAll('[data-row]').forEach(b=>b.onclick=()=>row(b.dataset.row));document.querySelector('#in').onclick=()=>zoom(1.3);document.querySelector('#out').onclick=()=>zoom(1/1.3);document.querySelector('#left').onclick=()=>{x+=viewport.clientWidth*.72;draw()};document.querySelector('#right').onclick=()=>{x-=viewport.clientWidth*.72;draw()};
document.querySelector('#overview').onclick=()=>{const b=art.querySelector('svg').viewBox.baseVal;scale=Math.min(viewport.clientWidth/b.width,viewport.clientHeight/b.height)*.97;x=(viewport.clientWidth-b.width*scale)/2;y=8;draw()};
document.querySelector('#references').onclick=()=>document.querySelector('#refs').style.display='block';document.querySelector('#closeRefs').onclick=()=>document.querySelector('#refs').style.display='none';
document.addEventListener('keydown',e=>{if(e.key==='ArrowRight'){x-=viewport.clientWidth*.5;draw()}if(e.key==='ArrowLeft'){x+=viewport.clientWidth*.5;draw()}if(e.key==='Escape')document.querySelector('#refs').style.display='none'});row('origin');
</script></html>\n'''


def build():
    data=load();svg=drawing.render(data);html=reader(svg,data)
    parsed=ET.fromstring(svg)
    for node in parsed.iter('{http://www.w3.org/2000/svg}a'):
        href=node.get('href','')
        if not href:raise ValueError('empty graphic link')
        if href.startswith(('https://','#')):continue
        path=(SVG.parent/unquote(href.split('#')[0])).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():raise ValueError('broken graphic link: '+href)
    audit={
        'schema':'oph.horizontal_particle_atlas.v2',
        'input_sha256':data['input_sha256'],
        'source_scientific_experiments_replayed':False,
        'source_physical_identifications':[],
        'reference_data_changed':False,
        'hero_stage_ids':[s['id'] for s in data['story']['stages']],
        'family_targets':{r['id']:r['target_ids'] for r in data['families']['lanes']},
        'measurement_illustration':data['measurement'],
        'primitive_endpoint_illustration':data['endpoints'],
        'double_path_illustration':data['double_path'],
        'finite_theorem_roles':data['theorems'],
        'named_field_meanings':data['field_nature'],
        'finite_observer_output_connections':data['observer_outputs']['finite_output_connections'],
        'primary_response_receipt_records':data['primary_response_receipt_records'],
        'candidate_correspondence':{
            'path':str(CANDIDATES.relative_to(ROOT)),
            'sha256':sha(CANDIDATES),
            'source_excitation_ids':[r['id'] for r in data['candidate_correspondence']['source_excitation_rows']],
            'named_targets':[r['target'] for r in data['candidate_correspondence']['target_rows']],
            'protected_catalogue_ids':[r['id'] for r in data['candidate_correspondence']['protected_catalogue_rows']],
            'finite_sector_ratio_result_sha256':data['candidate_correspondence']['finite_sector_ratio_control']['result_sha256'],
            'physical_particle_map':[], 'mass_predictions_GeV':[]},
        'candidate_identity_tests':{
            'path':str(IDENTITY_TESTS.relative_to(ROOT)),
            'sha256':sha(IDENTITY_TESTS),
            'schema':data['candidate_identity_tests']['schema'],
            'target_count':len(data['candidate_identity_tests']['target_rows']),
            'targets':[{key:row[key] for key in ('target','family','test_receipt','supporting_receipts','physical_identity_status')}
                       for row in data['candidate_identity_tests']['target_rows']],
            'original_dressed_pair_window':{key:data['candidate_identity_tests']['calculated_source_results']['charged']
                ['original_compensated_dressed_energy_window'][key] for key in (
                    'original_dressed_recorded_simple_strict_lower','reported_transfer_interval',
                    'complete_success_outcomes','complete_probe_outcomes',
                    'energy_window_is_a_single_endpoint_or_rest_mass')},
            'physical_particle_map':[], 'mass_predictions_GeV':[],
            'source_scientific_experiments_replayed_by_atlas':False},
        'common_construction_specification_resolved':True,
        'finite_physical_preparation_including_boundary_charges_resolved':True,
        'finite_primary_local_records_and_uniform_influence_constructed':True,
        'unified_dependency_status':{
            node['id']:node['status'] for node in data['joint_dependencies']['nodes']},
        'physical_continuum_and_spectrum_resolved':False,
        'conditional_candidate_discharges_primary_physical_dependencies':False,
        'full_theory_positive_global_gap_claimed':False,
        'canonical_joint_dependency_map':str(DEPENDENCIES.relative_to(ROOT)),
        'active_target':data['joint_dependencies']['active_target'],
        'deferred_branches':data['joint_dependencies']['deferred_branches'],
        'comparison_requirement_domains':data['observer_outputs']['comparison_requirement_domains'],
        'layout':{'width':parsed.get('width'),'height':parsed.get('height')},
        'svg_sha256':hashlib.sha256(svg.encode()).hexdigest(),
        'html_sha256':hashlib.sha256(html.encode()).hexdigest(),
    }
    if any(sha(ROOT/p)!=h for p,h in data['input_sha256'].items()):
        raise ValueError('inputs changed during illustration build')
    return svg,html,canonical(audit)


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    outputs=build()
    for path,text in zip((SVG,HTML,AUDIT),outputs):
        if a.write:path.write_text(text)
        elif not path.is_file() or path.read_text()!=text:raise ValueError('stale horizontal atlas: '+str(path))
    print('PASS: horizontal SVG, standalone reader and bounded audit')


if __name__=='__main__':main()
