"""Bind the executed scoped result to the checked public entry points."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def seal():
    result=json.loads((HERE/'RESULT.json').read_text())
    for filename,digest in result['source_sha256'].items():
        if sha(HERE/filename)!=digest:
            raise ValueError('integration source changed during run: '+filename)
    contract={
        'schema':'oph.chain.checked_integration_contract.v1',
        'entry_point':'python3 openwave-tuning/next/CHAIN/integration/run.py',
        'admission':'PASS_SCOPED_FINITE_COMPOSITIONS',
        'historical_source_snapshot':result['snapshot_execution'],
        'snapshot_manifest':'SNAPSHOT_MANIFEST.json','snapshot_manifest_sha256':sha(HERE/'SNAPSHOT_MANIFEST.json'),
        'native_process_to_common_physics_complete':False,
        'all31_physical_criteria_complete':False,
        'closed_quantum_links':result['transfer_gates']['closed_quantum_links'],
        'closed_classical_links':result['transfer_gates']['closed_classical_links'],
        'missing_transfers':result['transfer_gates']['missing_transfers'],
        'raw_helper_contract':{
            'raw_helpers_untrusted_input_safe':False,
            'clock_states':'Exact finite Hermitian positive matrices, trace in[0,1], shape384x384; trajectory trace1. Number-sector helpers require exact number-block diagonal input. Full unitary/native helpers retain cross-number coherences.',
            'clock_parameters':'Finite nonnegative time/rates/scales; generating marks in[0,1]. Trajectory requires0<rate<=uniform_rate and strictly ordered interior events with valid marks.',
            'force_parameters':'Finite nonnegative horizon and2x288 nonnegative real rates with positive total; two integer positions0..95,288 ordered distinct seams,valid target seam and literal endpoint function.',
            'theorem_domain':'Exact mathematical inputs satisfy these preconditions. Frozen independent scientific audits cover their declared valid members.',
            'frozen_source_files_preserved':True},
        'checked_facade_contract':{
            'file':'strict_facade.py','api':['StrictClock.evolve','StrictClock.native','StrictClock.semigroup','StrictClock.count_split','StrictClock.trajectory','StrictParticle.measure','StrictParticle.energy_coefficients','first_event'],
            'validates_before_raw_calls':True,
            'state_absolute_representation_tolerance':2e-12,'number_block_entry_tolerance':0.,
            'no_projection_clipping_symmetrization_or_renormalization':True,
            'compression_domain':'Exactly zero off-number-block entries required. The raw helper per-entry1e-13 tolerance is outside the checked compressed API domain; its accumulated coherence loss is demonstrated by a positive-matrix mutation.',
            'all_zero_classical_rates':'Explicitly refused because the raw exposure formula divides by total rate.',
            'source_identity':'first_event validates mathematical domains; classical_witness separately binds regenerated seam order and canonical native function. Quantum joins require an actual regenerated H/preparation/complete-clock-output witness.',
            'valid_input_map_equality':'Each facade returns the unchanged raw call with identical argument objects. Frozen comparator numerical max errors are zero in RESULT.json.',
            'positive_tolerance_boundary':'Tolerance permits numerical representations within the declared spectral error; no interval-positive certificate is claimed for near-zero eigenvalues.'},
        'audit_closure':{
            'independent_audits':{name:{k:v for k,v in entry.items() if k!='contract'} for name,entry in result['independent_audits'].items()},
            'raw_admission_findings_closed_for_checked_facade':True,
            'raw_admission_findings_changed_frozen_helpers':False,
            'excluded_physical_transfers_admitted':False,
            'guard_and_join_rejections':result['mutation_count']},
        'coherent_proper_code_transfer':'Declared reference-preserving Gauss-code load/export and routing; no fine field or boundary measurement. Spectral energy apparatus is grouped Luders with no native compiler claimed.',
        'quantum_coarsening':'Executed linear coarsening of all reservoir/count/regional labels preserves every stop family and348 energy poststate blocks; compared against frozen576x348 coefficients.',
        'field_bounds':{k:result['transfer_gates'][k] for k in ['preparation_radius_all384_labels','M6_complete_clock_TV_upper','M42_complete_clock_TV_upper_exact','M42_complete_clock_TV_upper','M42_extra_terminal_MIX_TV_upper','M42_particle_transfer']},
        'source_clock_map':'SOURCE_CLOCK_MAP.json','external_source_map':'EXTERNAL_SOURCE_MAP.json',
        'per_row_gate':'ALL31_GATE.json','row1_scope':result.get('row1_scope'),
        'input_pins':'INPUT_PINS.json','input_pins_sha256':sha(HERE/'INPUT_PINS.json'),
        'result':'RESULT.json','result_sha256':sha(HERE/'RESULT.json'),
        'proof':'PROOF.md','proof_sha256':sha(HERE/'PROOF.md'),
        'source_sha256':result['source_sha256'],
        'producer_files_modified':[],'official_rows_changed':[]}
    (HERE/'CONTRACT.json').write_text(json.dumps(contract,indent=2)+'\n')
    files={p.name:sha(p) for p in sorted(HERE.glob('*')) if p.is_file() and p.suffix in ('.py','.json','.md') and p.name!='MANIFEST.json'}
    (HERE/'MANIFEST.json').write_text(json.dumps({'schema':'oph.chain.integration_manifest.v1','files':files,'official_rows_changed':[]},indent=2)+'\n')


if __name__=='__main__':seal()
