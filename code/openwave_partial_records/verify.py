"""Verify the OpenWave partial-records evidence and recompute its results.

Checks every archived file in evidence/ against its original sealed manifest
(EVIDENCE.json), binds the sealed receipt fields of each row and recomputes the
results from the carrier graph, the native rule and the archived records
(recompute.py).

    python3 code/openwave_partial_records/verify.py
    python3 code/openwave_partial_records/verify.py --workspace DIR

With --workspace, DIR is a copy of the OpenWave workspace that contains the
original packets; every original sealed leaf is then checked as well.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / 'evidence'
BINDINGS = HERE / 'EVIDENCE.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(provider, name):
    return json.loads((EVIDENCE / provider / name).read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def manifest_leaves(path):
    """Leaves of a sealed manifest: sha256sum text or a JSON {'files': {...}} map."""
    if path.suffix == '.json':
        return dict(json.loads(path.read_text())['files'])
    leaves = {}
    for line in path.read_text().splitlines():
        digest, name = line.split('  ', 1)
        leaves[name] = digest
    return leaves


def providers():
    return {p['id']: p for p in json.loads(BINDINGS.read_text())['providers']}


def custody(workspace=None):
    bound = set()
    for pid, provider in providers().items():
        folder = EVIDENCE / pid
        manifest = folder / provider['manifest']
        require(sha(manifest) == provider['manifest_sha256'], f'{pid}: archived manifest differs')
        leaves = manifest_leaves(manifest)
        require(len(leaves) == provider['sealed_leaf_count'], f'{pid}: manifest leaf count')
        bound.add(manifest)
        for name, digest in provider['archived_files'].items():
            require(leaves.get(name) == digest, f'{pid}/{name}: not a leaf of the sealed manifest')
            require(sha(folder / name) == digest, f'{pid}/{name}: archived bytes differ')
            bound.add(folder / name)
        if workspace is not None:
            original = Path(workspace) / provider['original_path']
            require(sha(original / provider['manifest']) == provider['manifest_sha256'],
                    f'{pid}: original manifest differs')
            for name, digest in leaves.items():
                require(sha(original / name) == digest, f'{pid}/{name}: original leaf differs')
    on_disk = {p for p in EVIDENCE.rglob('*') if p.is_file()}
    require(on_disk == bound, 'evidence/ holds exactly the bound files')
    return len(bound)


def receipts():
    bindings = providers()

    # Rows 3 and 4: compact internal clock and its array replay.
    clock = read('compact_internal_clock', 'RECEIPT.json')
    sci = clock['scientific']
    src, cyc, clk, loc = sci['source'], sci['cycle'], sci['clock'], sci['localization']
    require(src['lambda'] == '209/210' and src['D_over_Var'] == '1/210', 'clock source eigenvalue')
    require(src['native_calls'] == 303660 and src['charged_native_columns'] == 302400
            and src['full_arc_dimension'] == 524288 and src['involution_columns'] == 524288, 'clock columns')
    require(cyc['triangle'] == [0, 1, 5] and cyc['walk_eigenvalue'] == 1
            and cyc['integer_divergence_zero'] and cyc['native_reversal_sign'] == -1, 'triangle reference')
    require(abs(clk['theta'] - math.acos(209 / 210)) < 1e-15
            and clk['internal_reference_frequency'] == clk['theta']
            and abs(clk['return_probability_frequency'] - 2 * clk['theta']) < 1e-15, 'theta')
    require(clk['all_position_event_labels_retained'] == 524288 and clk['witness_arc'] == 10243
            and abs(clk['witness_coefficients'][1]) > 1e-3, 'theta harmonic witness')
    require(sci['diagnostics']['full_record_harmonic_residual'] < 1e-12
            and sci['diagnostics']['walk_state_residual'] < 1e-12, 'harmonic residual')
    require(loc['all_integer_times_proved'] and loc['exact_integer_invariant_span_checked']
            and loc['within_original_twelve_port_cell'] and loc['position_envelope'] == list(range(12)),
            'one-cell confinement')
    review = read('compact_internal_clock_root_review', 'REPLAY.json')
    arrays = read('compact_internal_clock_root_review', 'FULL_ARRAY_REPLAY.json')
    require(review['exit'] == 0 and review['leaves'] == 36
            and review['manifest_sha256'] == bindings['compact_internal_clock']['manifest_sha256'], 'clock replay')
    require(arrays['all_arrays_byte_identical'] and arrays['artifacts'] == clock['artifacts'], 'clock arrays')

    # Row 6: C/CAR field embedding on the compact source.
    spin = read('compact_field_spin_join__existing_audit', 'RECEIPT.json')
    require(spin['C_complete_instrument_intertwiner'] and spin['L6_source_energy_unchanged'], 'C/CAR join')
    require(spin['native_calls'] == 302400 and spin['single_sector_columns_per_law'] == 604800
            and spin['CAR_C_and_E_columns'] == 1209600, 'CAR columns')
    review = read('compact_field_spin_root_review', 'REPLAY.json')
    require(review['exit'] == 0 and review['leaves'] == 15 and review['manifest_sha256']
            == bindings['compact_field_spin_join__existing_audit']['manifest_sha256'], 'spin replay')

    # Row 8: contact, contact energy and the propagating ring field.
    contact = read('neutral_source_completion', 'RECEIPT.json')
    energy = contact['contact_energy']
    require(contact['all_registered_component_checks_pass'] and contact['contact_branches'] == 2520
            and contact['literal_native_calls'] == 6400, 'contact branches')
    require(energy['preparation'] == [0, 1] and energy['mean_energy_decimal_diagnostic'] == 11.051928582384567
            and len(energy['bernstein_energy_pointer_probabilities']) == 17
            and energy['post_contact_energy_distribution'] == {'zero': 1}, 'contact energy pointer')
    review = read('neutral_boundary_root_review', 'REPLAY.json')
    require(review['exit'] == 0 and review['leaves'] == 24 and review['manifest_sha256']
            == bindings['neutral_source_completion']['manifest_sha256'], 'contact replay')
    ring = read('neutral_radiation_join', 'RECEIPT.json')
    field = ring['field']
    require(ring['all_checks_pass'] and all(ring['checks'].values()), 'ring gates')
    require(field['wave_equalities'] == 147840 and field['native_calls'] == 86751, 'ring wave equalities')
    require(not any(field['errors'].values()) and not any(field['source_errors'].values()), 'ring defects')
    flights = field['flights']
    require(len(flights) == 360 and all(
        f['mode'] == 'native' and f['received'] and f['final_field_norm'] == 0
        and [t['kind'] for t in f['transfers']] == ['emit', 'absorb']
        and f['transfers'][0]['change'] == [-2, 2] and f['transfers'][1]['change'] == [2, -2]
        for f in flights), 'contact flights launched and absorbed')

    # Row 20: residence records.
    res = read('row20_residence', 'RESULT.json')
    fit, samples = res['fit'], res['samples']
    require(samples['native'] == 5040 and samples['native_currents'] == 13213995
            and samples['native_censored'] == 0 and samples['control_censored'] == 0, 'residence samples')
    require(res['record_checks']['all_pass'] and res['record_checks']['native_branches'] == 3168
            and res['record_checks']['saved_tapes'] == 24, 'residence record checks')
    require(round(fit['exponent'], 6) == 1.035144
            and [round(x, 6) for x in fit['profile95_approx']] == [0.983376, 1.086373], 'residence fit')
    require(res['raw_sha256']['records.csv'] == sha(EVIDENCE / 'row20_residence' / 'records.csv'),
            'residence records are the analysed records')

    # Row 31: first-hit feedback well and the controlled Coulomb limit.
    probe = read('row31_probe_binding', 'RECEIPT.json')
    replay = read('row31_probe_binding', 'REPLAY.json')
    probe_leaves = manifest_leaves(EVIDENCE / 'row31_probe_binding' / 'MANIFEST.sha256')
    require(probe['actual_native_calls'] == 305157 and probe['all_registered_checks'] == 'pass'
            and not probe['supplied_gamma'] and not probe['supplied_numeric_potential'], 'probe controller')
    require(replay['byte_identical'] and replay['sha256']['RECEIPT.json'] == probe_leaves['RECEIPT.json']
            and replay['sha256']['RECORDS.json.gz'] == probe_leaves['RECORDS.json.gz'], 'probe replay')
    wells = probe['finite_wells']
    require(len(wells) == 3 and all(w['lowest_squared_energy'] < w['native_lambda']
                                    and abs(w['feedback_disabled_threshold'] - w['native_lambda']) < 1e-9
                                    for w in wells), 'attractive wells')
    proof = (EVIDENCE / 'row31_coulomb_limit' / 'PROOF.md').read_text()
    require('e_n = -(lambda c_*)²/(4n²)' in proof and 'Summing 2l+1 gives n²' in proof, 'shell theorem')
    unified = manifest_leaves(EVIDENCE / 'unified_wave_submission' / 'MANIFEST.sha256')
    base = 'evidence/orbital/evidence/openwave/row_31_orbital_levels/work/'
    require(unified[base + '2026-11-14_codex_coulomb_limit/PROOF.md']
            == bindings['row31_coulomb_limit']['archived_files']['PROOF.md']
            and unified[base + '2026-11-14_codex_coulomb_limit/MANIFEST.sha256']
            == bindings['row31_coulomb_limit']['manifest_sha256']
            and unified[base + '2026-11-14_codex_probe_binding/MANIFEST.sha256']
            == bindings['row31_probe_binding']['manifest_sha256'], 'row 31 packets in the unified submission')

    return {
        'row_3': {'native_calls': src['native_calls'], 'lambda': src['lambda'], 'theta': clk['theta'],
                  'outcomes_retained': clk['all_position_event_labels_retained'],
                  'witness_coefficients': clk['witness_coefficients'], 'arrays_byte_identical': True},
        'row_4': {'position_envelope': loc['position_envelope'], 'all_integer_times': True},
        'row_6': {'native_calls': spin['native_calls'], 'columns_per_law': spin['single_sector_columns_per_law']},
        'row_8': {'contact_branches': contact['contact_branches'], 'native_calls': contact['literal_native_calls'],
                  'mean_contact_energy': energy['mean_energy_decimal_diagnostic'],
                  'pointer_outcomes': len(energy['bernstein_energy_pointer_probabilities']),
                  'contact_flights_absorbed': len(flights), 'wave_equalities': field['wave_equalities']},
        'row_20': {'trajectories': samples['native'], 'currents': samples['native_currents'],
                   'exponent': fit['exponent'], 'profile95': fit['profile95_approx'],
                   'amplitude_b': fit['coefficients_a_b'][1]},
        'row_31': {'native_calls': probe['actual_native_calls'], 'first_hit': probe['source_hit_probabilities'],
                   'lowest_squared_energies': [w['lowest_squared_energy'] for w in wells],
                   'band_thresholds': [w['native_lambda'] for w in wells]},
    }


def public(result):
    """Drop floating diagnostics (keys starting with '_') from a recompute result."""
    return {k: v for k, v in result.items() if not k.startswith('_')}


def verify(workspace=None):
    spec = importlib.util.spec_from_file_location('partial_records_recompute', HERE / 'recompute.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    files = custody(workspace)
    sealed = receipts()
    graph, seams = module.load_graph()
    energy = read('neutral_source_completion', 'RECEIPT.json')['contact_energy']
    clock = module.clock(graph, seams)
    require(abs(clock['_theta_float'] - sealed['row_3']['theta']) < 1e-15, 'recomputed theta')
    require(all(abs(x - y) < 1e-14 for x, y in zip(clock['_witness'], sealed['row_3']['witness_coefficients'])),
            'recomputed witness coefficients')
    require(clock['position_envelope'] == sealed['row_4']['position_envelope'], 'recomputed envelope')
    orbital = module.orbital(seams, read('row31_probe_binding', 'RECEIPT.json'))
    require(all(abs(x - y) < 1e-8 for x, y in zip(orbital['_wells'], sealed['row_31']['lowest_squared_energies'])),
            'recomputed wells')
    recomputed = {'clock': public(clock), 'spin': public(module.spin(graph, seams)),
                  'contact': public(module.contact(seams, energy)),
                  'residence': public(module.residence(read('row20_residence', 'RESULT.json')['fit'])),
                  'orbital': public(orbital)}
    return {'archived_files_verified': files, 'providers': len(providers()), 'sealed': sealed,
            'recomputed': recomputed}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--workspace', help='copy of the OpenWave workspace with the original packets')
    args = parser.parse_args()
    print(json.dumps(verify(args.workspace), indent=2))
