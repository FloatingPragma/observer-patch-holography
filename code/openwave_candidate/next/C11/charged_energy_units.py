"""One explicit calibration of the actual C10 charged spectral response.

This reuses complete already executed reservoir/energy records. The displayed
GeV energies are calibrated source excitations, not predicted particle masses.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]


def generate(anchor_gev=.000511):
    if isinstance(anchor_gev, bool) or not isinstance(anchor_gev, (int, float)) or not math.isfinite(anchor_gev) or anchor_gev <= 0:
        raise ValueError('Positive finite declared GeV anchor required')
    joint_path = BASE / 'next/C10/JOINT_SOURCE.json'
    contract_path = BASE / 'next/C10/CANONICAL_SOURCE_CONTRACT.json'
    joint = json.loads(joint_path.read_text()); contract = json.loads(contract_path.read_text())
    for document in (joint, contract):
        for relative, expected in document['source_sha256'].items():
            path = (BASE / relative).resolve()
            if not path.is_relative_to(BASE) or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError('Changed interacting source: ' + relative)
    reference = joint['references']['global_allN']; ground = reference['reference_energy']
    response = reference['complete_reservoir_response']
    addition = next(r for r in response if r['source_charge_change'] == 1)
    # Selection is frozen by observed probability, not by target physical masses.
    peak = max(addition['complete_source_number_energy_records'], key=lambda r: (r['probability'], -r['energy']))
    source_anchor = peak['energy'] - ground
    if source_anchor <= 0:
        raise ValueError('Selected charged excitation is not positive')
    scale = anchor_gev / source_anchor
    records = []; peaks = []
    for branch in response:
        bins = [{**row, 'source_excitation': row['energy'] - ground,
                 'calibrated_excitation_GeV': scale * (row['energy'] - ground)}
                for row in branch['complete_source_number_energy_records']]
        records.append({'reservoir_initial': branch['reservoir_initial'], 'reservoir_final': branch['reservoir_final'],
                        'source_charge_change': branch['source_charge_change'],
                        'original_unconditional_probability': branch['probability'],
                        'all_number_energy_bins': bins,
                        'calibrated_unnormalized_apparatus_work_GeV': scale * branch['unnormalized_apparatus_work']})
        dominant = max(bins, key=lambda r: (r['probability'], -r['energy']))
        peaks.append({'source_charge_change': branch['source_charge_change'],
                      'reservoir_initial': branch['reservoir_initial'], 'reservoir_final': branch['reservoir_final'],
                      'dominant_source_number': dominant['source_number'],
                      'dominant_source_excitation': dominant['source_excitation'],
                      'unconditional_peak_probability': dominant['probability'],
                      'calibrated_peak_GeV': dominant['calibrated_excitation_GeV']})
    remove = next(p for p in peaks if p['source_charge_change'] == -1)
    add = next(p for p in peaks if p['source_charge_change'] == 1)
    checks = {'complete_four_reservoir_branches': len(records) == 4,
              'all_original_energy_groups_retained': all(len(a['all_number_energy_bins']) == len(b['complete_source_number_energy_records']) for a, b in zip(records, response)),
              'no_probability_conditioning_or_filter': all(a['probability'] == b['probability'] for ra, rb in zip(records, response) for a, b in zip(ra['all_number_energy_bins'], rb['complete_source_number_energy_records'])),
              'anchor_is_calibration_identity': math.isclose(add['calibrated_peak_GeV'], anchor_gev, abs_tol=1e-18, rel_tol=1e-14),
              'full_work_uses_one_common_scale': math.isclose(sum(r['calibrated_unnormalized_apparatus_work_GeV'] for r in records), scale * reference['mean_reservoir_apparatus_work'], abs_tol=1e-17, rel_tol=1e-13)}
    paths = [Path(__file__), joint_path, contract_path]
    return {'schema': 'oph.c11.charged_energy_units.v1', 'construction_id': contract['construction_id'],
            'checks': checks, 'source_reference_energy': ground,
            'calibration': {'selection_rule': 'largest unconditional probability in source addition response; lower energy breaks an exact tie',
                            'source_anchor_excitation': source_anchor, 'supplied_anchor_GeV': anchor_gev,
                            'GeV_per_source_energy': scale, 'empirical_anchors_used': 1,
                            'anchor_agreement_is_prediction': False,
                            'interpretation': 'The chosen0.000511GeV value is a supplied display calibration. The model does not identify this excitation with an electron or derive the value.'},
            'dominant_lines': peaks,
            'scale_free_removal_to_addition_peak_ratio': remove['dominant_source_excitation'] / source_anchor,
            'complete_response': records,
            'scope': 'All interacting source-number/energy bins, zeros and original unconditional weights retained. One positive multiplicative calibration is used; no energy offset or per-species fit. Physical species, rest momentum, radiative matching and charge conjugation are not established.',
            'physical_particle_masses_predicted': [], 'full_criterion_met': False, 'official_rows_changed': [],
            'source_sha256': {str(p.relative_to(BASE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true'); mode.add_argument('--check', action='store_true')
    args = ap.parse_args(); data = generate(); path = HERE / 'CHARGED_ENERGY_UNITS.json'
    if not all(v is True for v in data['checks'].values()):
        raise ValueError(data['checks'])
    if args.write:
        path.write_text(json.dumps(data, indent=2) + '\n')
        table = '\n'.join(f"| {p['source_charge_change']:+d} | {p['reservoir_initial']} → {p['reservoir_final']} | {p['dominant_source_excitation']:.10g} | {p['calibrated_peak_GeV']:.10g} | {p['unconditional_peak_probability']:.10g} |" for p in data['dominant_lines'])
        path.with_suffix('.md').write_text('# Calibrated charged source energies\n\nThe complete C10 interacting reservoir response is displayed with one supplied anchor: the most probable addition line is assigned0.000511GeV. This is a unit calibration, not a predicted electron mass. All groups, zeros and original probabilities remain in the JSON; no source is rerun.\n\n| Charge change | Reservoir | Source excitation | Calibrated GeV | Unconditional peak probability |\n| --- | --- | --- | --- | --- |\n' + table + '\n\nThe scale-free removal/addition peak ratio is ' + str(data['scale_free_removal_to_addition_peak_ratio']) + '. It reflects this fixed-background source and apparatus, with no automatic particle/antiparticle identification.\n')
    elif json.loads(path.read_text()) != data:
        raise ValueError('Calibrated complete response differs')
    print(json.dumps({'checks': data['checks'], 'dominant_lines': data['dominant_lines']}))
