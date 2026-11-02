#!/usr/bin/env python3
"""Build a proof-bound mass/energy catalogue without assigning physical species."""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction as F
import hashlib
import io
import json
from math import gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = 'tools/particle_mass_table_sources.json'
JSON_OUTPUT = 'model/completion/particle_mass_table.json'
MARKDOWN_OUTPUT = 'model/completion/PARTICLE_MASS_TABLE.md'
INPUTS = (
    'model/native_mass_channels/README.md', 'model/native_mass_channels/algebra.py',
    'model/native_mass_channels/results.json',
    'model/ballistic_free_boundary/README.md', 'model/ballistic_free_boundary/bounds.py',
    'model/ballistic_free_boundary/results.json',
    'model/spin_canonical_energy/README.md', 'model/spin_canonical_energy/results.json',
    'model/centered_particle_spectrum/README.md', 'model/centered_particle_spectrum/source.py',
    'model/centered_particle_spectrum/certificate.py', 'model/centered_particle_spectrum/witness.npz',
    'model/centered_particle_spectrum/results.json', 'model/centered_particle_spectrum/lines.csv',
    'tools/read_ordered_spin_energies.py', 'model/completion/ordered_spin_energy_readout.json',
    'model/native_ordered_ground/source.py', 'model/native_ordered_ground/certificate.py',
    'model/native_ordered_ground/witness.npz', 'model/native_ordered_ground/results.json',
)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def loads(data):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate JSON key')
            result[key] = value
        return result
    def bad(_):
        raise ValueError('nonfinite JSON constant')
    return json.loads(data, object_pairs_hook=pairs, parse_constant=bad)


def exact(value):
    if type(value) is int:
        return F(value)
    if (type(value) is not F or type(value.numerator) is not int
            or type(value.denominator) is not int or value.denominator <= 0
            or gcd(value.numerator, value.denominator) != 1):
        raise ValueError('canonical exact rational required')
    return F(value.numerator, value.denominator)


def designed_masses(coefficients):
    """Optional exact design, deliberately not a prediction of selected order.

    The default designs have rational vector norm. General coefficients are
    represented by their exact squared norm, without floating square roots.
    """
    if type(coefficients) not in (tuple, list) or len(coefficients) != 4:
        raise ValueError('four exact coefficients required')
    m = tuple(exact(x) for x in coefficients)
    r2 = sum((x*x for x in m[1:]), F())
    n, d = isqrt(r2.numerator), isqrt(r2.denominator)
    rational_norm = n*n == r2.numerator and d*d == r2.denominator
    branches = None
    ratio = None
    signed_levels = None
    if rational_norm:
        r = F(n, d)
        values = sorted((abs(m[0]-r), abs(m[0]+r)))
        branches = [str(x) for x in values]
        ratio = str(values[1]/values[0]) if values[0] else None
        multiplicities = {}
        for energy in values:
            for sign in (-1, 1):
                level = sign*energy
                multiplicities[level] = multiplicities.get(level, 0)+4
        signed_levels = [{'energy': str(e), 'multiplicity': count}
                         for e, count in sorted(multiplicities.items())]
    return {'coefficients': [str(x) for x in m], 'vector_norm_squared': str(r2),
            'rest_energy_branches': branches, 'heavy_over_light_ratio': ratio,
            'signed_zero_momentum_levels': signed_levels,
            'general_formula': '|m0-|m|| and |m0+|m||',
            'coefficient_status': 'externally supplied quadratic design',
            'is_physical_prediction': False, 'physical_species': None,
            'calibration_to_eV': None}


def arbitrary_ratio_design(ratio):
    q = exact(ratio)
    if q < 1:
        raise ValueError('ordered positive mass ratio must be at least one')
    return designed_masses(((q+1)/2, (q-1)/2, F(), F()))


def arbitrary_pair_design(larger, smaller):
    a, b = exact(larger), exact(smaller)
    if not a >= b >= 0:
        raise ValueError('require larger>=smaller>=0')
    return designed_masses(((a+b)/2, (a-b)/2, F(), F()))


def interval(values, positive=False):
    if type(values) is not list or len(values) != 2 or any(type(x) is not str for x in values):
        raise ValueError('two exact string interval endpoints required')
    lo, hi = map(F, values)
    if lo > hi or (positive and lo <= 0):
        raise ValueError('ordered positive interval required')
    return lo, hi


def decimal_interval(lo, hi, places=12):
    if lo > hi:
        raise ValueError('ordered interval required')
    scale = 10**places
    lower = (lo*scale).numerator//(lo*scale).denominator
    upper = -((-hi*scale).numerator//(-hi*scale).denominator)
    def display(n):
        sign = '-' if n < 0 else ''
        whole, part = divmod(abs(n), scale)
        return sign+str(whole)+'.'+str(part).zfill(places)
    return [display(lower), display(upper)]


def checked_snapshot(root):
    root = Path(root)
    manifest_bytes = (root/MANIFEST).read_bytes()
    expected = loads(manifest_bytes)
    if (type(expected) is not dict or set(expected) != {'schema', 'reviewed_sources'}
            or expected['schema'] != 'oph.particle_mass_table.sources.v1'
            or set(expected['reviewed_sources']) != set(INPUTS)):
        raise ValueError('complete reviewed mass-table source manifest required')
    snapshot = {name: (root/name).read_bytes() for name in INPUTS}
    for name, data in snapshot.items():
        digest = expected['reviewed_sources'][name]
        if type(digest) is not str or hashlib.sha256(data).hexdigest() != digest:
            raise ValueError('reviewed mass-table input changed: '+name)
    snapshot[MANIFEST] = manifest_bytes
    snapshot['tools/build_particle_mass_table.py'] = (root/'tools/build_particle_mass_table.py').read_bytes()
    return snapshot


def centered_rows(receipt, csv_bytes):
    spectrum = receipt['spectrum']
    if (receipt.get('all_pass') is not True
            or spectrum.get('positive_visible_lines_per_local_spin_charge_branch') != 228
            or spectrum.get('all_visible_weights_strictly_positive') is not True
            or spectrum.get('finite_energies_identified_as_relativistic_rest_masses') is not False
            or spectrum.get('physical_particle_species_identified') is not False):
        raise ValueError('centered finite-spectrum scope changed')
    rows = spectrum['rows']
    csv_rows = list(csv.DictReader(io.StringIO(csv_bytes.decode())))
    if len(rows) != 228 or len(csv_rows) != 228:
        raise ValueError('all228 centered spectral rows required')
    first = interval(rows[0]['excitation_energy_interval'], True)
    result = []
    previous = F()
    for index, (row, csv_row) in enumerate(zip(rows, csv_rows, strict=True)):
        energy = interval(row['excitation_energy_interval'], True)
        weight = interval(row['local_branch_weight_interval'], True)
        if type(row['index']) is not int or row['index'] != index or energy[0] <= previous:
            raise ValueError('ordered disjoint centered lines required')
        if (csv_row['index'] != str(index) or
                [csv_row['energy_lower'], csv_row['energy_upper']] != row['excitation_energy_interval'] or
                [csv_row['weight_lower'], csv_row['weight_upper']] != row['local_branch_weight_interval']):
            raise ValueError('CSV and exact spectral receipt disagree')
        previous = energy[1]
        ratio = (F(1), F(1)) if index == 0 else (energy[0]/first[1], energy[1]/first[0])
        result.append({'index': index, 'energy_interval': row['excitation_energy_interval'],
                       'relative_to_first_interval': decimal_interval(*ratio),
                       'local_weight_interval': row['local_branch_weight_interval']})
    return result


def assemble(snapshot):
    data = lambda name: loads(snapshot[name])
    mass = data('model/native_mass_channels/results.json')
    cat = mass['catalogue']
    if (mass.get('all_pass') is not True or mass.get('integrated_gate_accepted') is not False
            or type(mass.get('physical_comparison_count_added')) is not int
            or mass['physical_comparison_count_added'] != 0
            or cat.get('physical_species_identified') is not False
            or cat.get('spontaneous_order_proved') is not False
            or cat.get('interacting_mass_shell_proved') is not False
            or cat['exhaustion']['solution_dimension'] != 4
            or cat['clifford']['folded_dimension'] != 16
            or cat['clifford']['generic_signed_level_multiplicities'] != [4, 4, 4, 4]):
        raise ValueError('supplied quadratic catalogue scope changed')
    free = data('model/ballistic_free_boundary/results.json')['theorem']
    if (free.get('rest_masses_in_this_free_class') != 'zero'
            or free.get('full_field_limit') != 'free eight-valley Weyl vacuum hierarchy'
            or free.get('physical_species_identified') is not False):
        raise ValueError('free continuum theorem changed')
    spin = data('model/spin_canonical_energy/results.json')
    theorem = spin['theorem']
    if (theorem.get('actual_continuum_dynamics_nonstatic_and_gapless') is not True
            or theorem.get('sharp_mass_shell_proved') is not False
            or theorem.get('physical_species_or_particle_masses_identified') is not False):
        raise ValueError('spin field energy cannot be promoted to a particle mass')
    charged = data('model/centered_particle_spectrum/results.json')
    charged_lines = centered_rows(charged, snapshot['model/centered_particle_spectrum/lines.csv'])
    ordered = data('model/completion/ordered_spin_energy_readout.json')
    if (ordered.get('total_certified_states') != 256 or ordered.get('cluster_count') != 68
            or ordered.get('relativistic_masses_identified') is not False
            or ordered.get('physical_species_identified') is not False
            or ordered.get('zero_field_source_ground_spectrum') is not False
            or len(ordered['all_ordered_levels']) != 256
            or len(ordered['isolated_clusters']) != 68):
        raise ValueError('finite ordered-field energies changed scope')
    for name, digest in ordered['source_bindings'].items():
        if name not in snapshot or hashlib.sha256(snapshot[name]).hexdigest() != digest:
            raise ValueError('ordered energy source binding changed: '+name)
    if ordered['tool_sha256'] != hashlib.sha256(snapshot['tools/read_ordered_spin_energies.py']).hexdigest():
        raise ValueError('ordered energy producer binding changed')
    classes = [
        {'sector': 'canonical free continuum', 'rest_energy': '0', 'mass': '0',
         'energy_momentum': 'E=v|p|', 'multiplicity': 'eight Weyl valleys',
         'source_condition': free['sufficient_coupling_class'],
         'status': 'proved free scaling class; no physical species identification',
         'calibration_to_eV': None, 'physical_species': None, 'physical_prediction': False},
        {'sector': 'supplied anticommuting quadratic masses',
         'rest_energy': '|m0-|m|| and |m0+|m||', 'mass': 'E_rest/v^2 if v is the limiting speed',
         'energy_momentum': 'E=+/-sqrt(c^2 sum sin(p_j)^2+(m0+/-|m|)^2)',
         'multiplicity': 'four positive and four negative levels per generic branch; coincident levels add',
         'source_condition': 'quadratic h+B with four freely supplied coefficients; not the U=1 interacting centered-cube spectrum',
         'status': 'exact quadratic design catalogue, not a prediction of its coefficients',
         'calibration_to_eV': None, 'physical_species': None, 'physical_prediction': False},
        {'sector': 'actual canonical ordered-spin continuum', 'rest_energy': None, 'mass': None,
         'energy_momentum': 'positive energy windows proportional to momentum; no sharp dispersion proved',
         'multiplicity': 'one canonical pair; no sharp particle count established',
         'source_condition': 'same selected ordered ground; L=lim K/(Jp)',
         'status': 'nonstatic and gapless continuum fields; sharp rest mass undetermined',
         'calibration_to_eV': None, 'physical_species': None, 'physical_prediction': False},
    ]
    return {'schema': 'oph.particle_mass_table.v1', 'physical_prediction_count': 0,
            'new_empirical_comparison_count': 0, 'all_physical_particle_masses_derived': False,
            'classes': classes,
            'coefficient_designs': [designed_masses((0, 0, 0, 1)),
                                    designed_masses((2, 0, 0, 1)), arbitrary_ratio_design(F(7, 2))],
            'arbitrary_ratio_identity': 'm0=(q+1)/2, m1=(q-1)/2, m2=m3=0 gives rest energies1,q for every q>=1',
            'arbitrary_pair_identity': 'm0=(A+B)/2, |m|=(A-B)/2 gives any supplied pair A>=B>=0',
            'maximum_distinct_positive_rest_energies_in_uniform_family': 2,
            'three_distinct_lepton_masses_from_one_uniform_family': False,
            'spin_order_identified_with_mass_coefficients': False,
            'spin_positive_windows': spin['spectral_windows'],
            'centered_finite_CAR': {
                'reference': charged['spectrum']['reference'],
                'energy_unit': charged['spectrum']['energy_unit'], 'visible_lines_per_charge_spin': 228,
                'all_lines': charged_lines, 'source_full_table': 'model/centered_particle_spectrum/lines.csv',
                'finite_gaps_are_relativistic_masses': False,
                'single_scale_charged_lepton_dictionary_excluded': charged['comparison']['all_three_distinct_frequency_assignments_excluded'],
                'conditional_dictionary_only': True},
            'ordered_finite_spin': {
                'Hamiltonian': ordered['Hamiltonian'], 'J': ordered['J'], 'B_over_J': ordered['B_over_J'],
                'states': 256, 'enclosure_groups': 68,
                'first_six_positive_groups': ordered['isolated_clusters'][1:7],
                'source_full_table': 'model/completion/ordered_spin_energy_readout.json',
                'finite_gaps_are_relativistic_masses': False, 'zero_field_spectrum': False},
            'custody': {'kind': 'Reviewed proof outputs and exact rational table arithmetic',
                        'inherited_proof_producers_rerun_here': False,
                        'calibration_or_species_fit_performed': False},
            'source_bindings': {name: hashlib.sha256(value).hexdigest() for name, value in sorted(snapshot.items())}}


def build(root=ROOT):
    snapshot = checked_snapshot(root)
    result = assemble(snapshot)
    if any((Path(root)/name).read_bytes() != raw for name, raw in snapshot.items()):
        raise ValueError('mass-table source changed during build')
    return result


def markdown(table):
    lines = ['# OPH mass and energy readout', '',
             'The table separates proved continuum results, supplied quadratic designs and finite excitation gaps. '
             '**No measured particle species or calibrated physical mass has been derived.**', '',
             '| Sector | Result | Status |', '| --- | --- | --- |',
             '| [Canonical free continuum](../ballistic_free_boundary/README.md) | Eight Weyl valleys, E=v·norm(p), rest energy 0 | Proved when U=o(c^(2/7)); no species assignment |',
             '| [Supplied quadratic masses](../native_mass_channels/README.md) | E_rest=abs(m0−norm(m)) and abs(m0+norm(m)) | Four positive and four negative levels per generic branch; coefficients are inputs |',
             '| [Actual ordered-spin continuum](../spin_canonical_energy/README.md) | Nonzero positive energy bands at every small momentum | Sharp rest mass and dispersion remain undetermined |', '',
             'For a relativistic branch with limiting speed v, mass would be E_rest/v². '
             'Finite lattice gaps and a zero spectral threshold alone do not establish that relation. '
             'No conversion to eV or MeV is supplied.', '',
             '## Exact coefficient designs', '',
             'These are explicitly chosen quadratic Hamiltonians h+B, not predictions or fits. '
             'Their formula is not the spectrum of the U=1 interacting centered cube below.', '',
             '| (m0,m1,m2,m3) | Rest energies | Larger/smaller |', '| --- | --- | --- |']
    for row in table['coefficient_designs']:
        lines.append('| ('+', '.join(row['coefficients'])+') | '+', '.join(row['rest_energy_branches'])+' | '+str(row['heavy_over_light_ratio'])+' |')
    lines += ['', 'Choosing m0=(A+B)/2 and norm(m)=(A−B)/2 produces any chosen pair A≥B≥0. '
              'In particular, m0=(q+1)/2 and m1=(q−1)/2 produce energies 1 and q for any q≥1. '
              'The free coefficients therefore do not predict a mass ratio. This single uniform family '
              'has at most two distinct positive rest energies, so a direct assignment cannot reproduce '
              'three distinct charged-lepton masses. This statement does not exclude other OPH sectors. '
              'The ordered spin expectation has not been identified with these coefficients by a self-consistency proof.', '',
              '## Centered finite charged spectrum', '',
              'The first six of 228 visible lines per local spin/charge branch follow. '
              'Reference: unique centered N8 ground, K=H−N/2, c=1/8, U=1, open 2×2×2 cube. '
              'Intervals use the dimensionless source clock; all nonzero residues remain in the '
              '[full 228-line table](../centered_particle_spectrum/lines.csv).', '',
              '| Line | Excitation energy interval | Relative to first |', '| --- | --- | --- |']
    for row in table['centered_finite_CAR']['all_lines'][:6]:
        lines.append(f"| {row['index']} | [{', '.join(row['energy_interval'])}] | [{', '.join(row['relative_to_first_interval'])}] |")
    lines += ['', 'Treating these finite gaps as charged-lepton rest energies under one positive scale '
              'is already excluded by the retained ratio certificate. This excludes that dictionary, '
              'not all continuum or composite sectors.', '',
              '## Finite ordered-spin reference', '',
              'Here J=1 and B/J=1/4 on the open eight-site cube. The '
              '[complete certified 256-state readout](ordered_spin_energy_readout.json) has 68 '
              'separated enclosure groups; grouping does not assert exact degeneracy. '
              'This is a positive-field finite reference, not the infinite zero-field ground.', '',
              '| Positive group | Excitation energy interval | Relative to first |', '| --- | --- | --- |']
    for row in table['ordered_finite_spin']['first_six_positive_groups']:
        lines.append(f"| {row['cluster']} | [{', '.join(row['gap_interval'])}] | [{', '.join(row['gap_over_first_gap_interval'])}] |")
    lines += ['', 'The [filtered spin-packet readout](../spin_packet_readout/README.md) additionally '
              'gives certified conditional excitation energies and success probabilities for this same '
              'finite reference, including all rejection outcomes. Its continuum packet construction '
              'has an explicit finite resource bound; it does not determine a relativistic rest mass.', '',
              '## Reproduce', '',
              '```sh', 'python -W error tools/build_particle_mass_table.py --check',
              'python -W error -m pytest -q tools/test_particle_mass_table.py', '```', '',
              'The [machine-readable table](particle_mass_table.json) retains all 228 charged lines, '
              'the exact coefficient examples and source hashes. Its builder checks reviewed proof '
              'outputs and performs exact rational arithmetic; it does not rerun every inherited proof producer.', '']
    return '\n'.join(lines)


def validate(value, root=ROOT):
    if canonical(value) != canonical(build(root)):
        raise ValueError('particle mass table or scientific scope changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', action='store_true')
    action.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = build()
    if args.write:
        (ROOT/JSON_OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
        (ROOT/MARKDOWN_OUTPUT).write_text(markdown(result))
    else:
        validate(loads((ROOT/JSON_OUTPUT).read_bytes()))
        if (ROOT/MARKDOWN_OUTPUT).read_text() != markdown(result):
            raise ValueError('particle mass table Markdown changed')
    print('PASS particle mass table; physical predictions0')


if __name__ == '__main__':
    main()
