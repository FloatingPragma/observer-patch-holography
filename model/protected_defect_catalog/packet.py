"""Catalogue protected objects with their exact completeness domains."""
from pathlib import Path
import hashlib
import json

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from model.completion import defect_algebra
from model.defect_protection import protection
from . import pauli, charge

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENTS = tuple('model/'+name+'/results.json' for name in (
    'defect_protection', 'local_exchange_charge', 'common_spin_protection',
    'cube_eta_protection', 'compact_defect_pair', 'defect_localization_readout'))


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def bindings():
    for name in PARENTS:
        parent = strict_json.load(ROOT/name)
        pins = parent.get('source_bindings', parent.get('source_sha256'))
        if not isinstance(pins, dict) or not pins:
            raise ValueError('catalogue parent has no source custody: '+name)
        for path, digest in pins.items():
            target = (ROOT/path).resolve()
            if not target.is_relative_to(ROOT) or hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError('stale catalogue parent: '+path)
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/protected_defect_catalog/PROOF.md',
        'source_receipts': list(PARENTS),
        'independent_checks': ['model/protected_defect_catalog/README.md']})
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def rows():
    return [
        {'id':'bare_native_source', 'kind':'exhaustive_exclusion',
         'source':'complete all-seam C/E on H66 or H78; source-only charges',
         'protected_object':'only scalar branch-commuting charges',
         'completeness':'entire branch commutant, every quantum basis',
         'preparation_readout':'native complete occupation/current records',
         'propagation_energy':'no nontrivial charge supplied by this classification',
         'topology':'no protected source-only reducing charge in this domain',
         'boundary':'history, coherent memory, recovery, partial seam support and refinement are outside this domain',
         'evidence':'model/completion/defect_algebra.py'},
        {'id':'ordered_exchange_charge', 'kind':'protected_memory_charge',
         'source':'H66 with a coherent orientation bit; complete ordered-cover repair',
         'protected_object':'Q=I66 tensor X, two sharp eigenvalue sectors',
         'completeness':'entire all-seam commutant span{I,Q}',
         'preparation_readout':'coherent memory preparation; complete charge pointer; local qutrit realization retains the same Q',
         'propagation_energy':'exact pair-configuration transport; no relativistic rest-energy map',
         'topology':'Z2 character of the ordered exchange cover',
         'boundary':'sharp charge is not strictly localized with reference exterior; local nonsharp excess is available; memory dephasing destroys protection',
         'evidence':'model/defect_protection/README.md',
         'same_object_realizations':['model/local_exchange_charge/README.md','model/protected_reference/PROOF.md']},
        {'id':'selected_spin_charge', 'kind':'conserved_symmetry_charge',
         'source':'actual eight-site selected spin cube, J=1, B=1/4; inherited two-flavor repair',
         'protected_object':'integer U1 spin charge Q; global X parity is exp(i pi Q); local Z creates two sharp Q=+/-1 channels',
         'completeness':'all 65536 Pauli-product words plus exact full charge spectral resolution; not all conserved operators',
         'preparation_readout':'complete native local-Z programs and ideal 256-outcome X-basis charge instrument; full native charge tape not executed',
         'propagation_energy':'positive finite same-H_B energy; contrast bound through ten repairs for T<=1/16 under the stated source scaling',
         'topology':'global internal U1 symmetry, not a winding number or identified electric charge',
         'boundary':'no nonidentity proper-support Pauli charge; coherent controls can change parity; no permanent localization or particle pole',
         'evidence':'model/defect_localization_readout/PROOF.md'},
        {'id':'single_occupancy_spin_memory', 'kind':'repair_protected_algebra',
         'source':'eight-site one-particle-per-site spin code with complete C/E repair',
         'protected_object':'full M256 operator algebra under repair alone',
         'completeness':'entire operator algebra of this code under the specified repair',
         'preparation_readout':'declared spin control and detector services',
         'propagation_energy':'exchange evolution acts nontrivially within the code',
         'topology':'no topological label supplied by identity repair',
         'boundary':'repair protection of every state does not make every state an energy eigenstate or a distinct species',
         'evidence':'model/defect_localization_readout/PROOF.md'},
        {'id':'common_hidden_spin_factor', 'kind':'joint_noiseless_subsystem',
         'source':'same six-site Weyl/contact field and complete charged-block repair',
         'protected_object':'at N=4, S=1, a qutrit with 105 multiplicity states, total dimension315',
         'completeness':'exhibited common spin factor; maximal joint noiseless algebra not classified',
         'preparation_readout':'complete ideal joined histories with arbitrary spin/reference correlations',
         'propagation_energy':'multiplicity evolves and localized carrier leaks; repair need not conserve field energy',
         'topology':'internal spin factor, no topological boundary invariant',
         'boundary':'global factor preservation does not identify localized matter, Lorentz spin or a physical spin apparatus',
         'evidence':'model/common_spin_protection/README.md'},
        {'id':'cube_eta_range', 'kind':'joint_protected_sector',
         'source':'original eight-site Weyl/contact field with the stated four-dimer repair',
         'protected_object':'rank1820 image of eta+ from N4 to N6; eta>=2 sector',
         'completeness':'entire image of this injective map; all protected sectors not classified',
         'preparation_readout':'source-prepared eta reference and complete local response detector',
         'propagation_energy':'same interacting action intertwines with an energy shift U; one probe has eight inelastic lines plus elastic weight',
         'topology':'eta symmetry representation, no spatial winding invariant',
         'boundary':'eight response lines are not eight species; protection uses the complete outside outcome for this repair domain',
         'evidence':'model/cube_eta_protection/README.md'},
        {'id':'compact_opposite_texture_cores', 'kind':'specified_topological_background',
         'source':'prescribed two-core texture in the declared cubic two-spin field',
         'protected_object':'opposite shell degrees +1 and -1, constant exterior, four localized band modes',
         'completeness':'specified texture family and certified isolated band, not all native defects',
         'preparation_readout':'declared texture and transported band controls; autonomous source formation unproved',
         'propagation_energy':'uniform finite texture/vacuum-energy bounds and controlled band transport',
         'topology':'genuine shell winding of the prescribed order-parameter field',
         'boundary':'winding does not force zero modes, four generations or physical rest masses; stability needs the stated dynamics',
         'evidence':'model/compact_defect_pair/PROOF.md'},
    ]


def build():
    pins = bindings()
    bare = defect_algebra.report()
    cover = {law:protection.algebra_certificate(law) for law in ('C','E')}
    bare_counts = {law:{domain:{key:bare['C_E'][law][domain][key] for key in (
        'dimension','generated_star_algebra_complex_dimension','commutant_complex_dimension')}
        for domain in ('full78','binary66')} for law in ('C','E')}
    cover_counts = {law:{key:record[key] for key in (
        'joint_dimension','generated_star_algebra_complex_dimension','commutant_complex_dimension',
        'complete_branch_X_commutation','ordered_cover_reachable_states','raw_maps_sha256')}
        for law,record in cover.items()}
    entries = rows()
    for entry in entries:
        entry['identified_Standard_Model_species'] = []
        entry['physical_mass'] = None
        entry['same_source_relativistic_particle_attachment_proved'] = False
    if bindings() != pins:
        raise RuntimeError('catalogue source changed during construction')
    return {'schema':'oph.protected_defect_catalog.v1', 'source_bindings':pins,
            'bare_complete_branch_classification':bare_counts,
            'ordered_cover_complete_branch_classification':cover_counts,
            'selected_spin_complete_Pauli_classification':pauli.census(),
            'same_source_charge_resolved_creation':charge.report(),
            'entries':entries,
            'catalogue_is_exhaustive_over_all_OPH_sources_or_defects':False,
            'full_parent_scientific_replays_performed':False,
            'parent_source_bindings_authenticated':True,
            'replayed_classifications':['bare all-seam branch commutant',
                                       'ordered-cover all-seam branch commutant',
                                       'actual selected-cube Pauli-product charges'],
            'physical_particle_map':[], 'physical_mass_prediction_count':0}


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('catalogue classification, scope or source custody changed')


def render(value):
    lines = ['# Protected sectors and defect candidates', '',
             'The catalogue records separate source architectures. It is not a list of physical species.',
             'Completeness applies only to the explicit domains in the table; there is no exhaustive census over all OPH realizations.', '',
             '| Object | Protecting source and object | Completeness | Particle attachment |',
             '| --- | --- | --- | --- |']
    for row in value['entries']:
        lines.append('| '+row['id']+' | '+row['protected_object']+' | '+row['completeness']+' | Unproved |')
    for row in value['entries']:
        lines += ['', '## '+row['id'].replace('_',' '), '',
                  '**Source:** '+row['source']+'.', '',
                  '**Preparation and readout:** '+row['preparation_readout']+'.', '',
                  '**Propagation and energy:** '+row['propagation_energy']+'.', '',
                  '**Topology:** '+row['topology']+'.', '',
                  '**Boundary:** '+row['boundary']+'.', '',
                  '[Owning proof or construction](../../'+row['evidence']+').']
    lines += ['', 'Source family, state, clock, physical spin/statistics, dispersion or resonance,',
              'detector coupling and experimental units must be joined before a catalogue',
              'entry can fill the [Standard Model correspondence](../standard_model_source_map/PARTICLE_MAP.md).', '']
    return '\n'.join(lines)
