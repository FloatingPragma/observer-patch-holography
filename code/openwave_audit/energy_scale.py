"""Conditional energy-range test, independent of a momentum/Bloch dictionary.

For unit-weight graphs of maximum degree d, ||L|| <= 2d and
||H_m|| <= sqrt(m^2 + 2d), H_m = [[mI,B],[B*, -mI]].
The anchored one-load energy diameter is <= 2 J sqrt(m^2 + 2d).
This does not constrain a different many-body or interacting source.
"""
import hashlib
import json
import math
from pathlib import Path


def run():
    electron_mev = .51099895
    lep_mev = 104500.
    degree = 6
    rows = []
    for label, mass, anchor in [
        ('lowest_flat_level', 0., math.sqrt(5 - math.sqrt(5))),
        ('supplied_rest_0.23', .23, .23),
        ('supplied_rest_1e-6', 1e-6, 1e-6),
    ]:
        j = electron_mev / anchor
        absolute = j * math.sqrt(mass * mass + 2 * degree)
        rows.append({'anchor': label, 'mass_native': mass, 'anchor_native': anchor,
                     'J_MeV': j, 'absolute_energy_bound_MeV': absolute,
                     'energy_diameter_bound_MeV': 2 * absolute,
                     'cannot_span_104500_MeV_in_this_one_load_source': 2 * absolute < lep_mev})
    hbarc_gev_m = .1973269804e-15
    k_corner = 4 * math.pi / (3 * math.sqrt(3))
    spacing = k_corner * hbarc_gev_m / 104.5
    return {
        'schema': 'oph.openwave.conditional_energy_range.v1',
        'hypotheses': ['unit-weight port graph of maximum degree at most 6',
                       'one-load Hm source, no added field or unbounded oscillator',
                       'physical energy differences equal J times source energy differences',
                       'one physical electron identified with this same one-load source at beam energies',
                       'declared anchor E0 J = 0.51099895 MeV'],
        'rows': rows,
        'bloch_momentum_comparison': {
            'additional_hypothesis': 'physical momentum equals hbar times crystal wavevector restricted to the first Brillouin zone',
            'beam_momentum_GeV_c': 104.5,
            'conditional_spacing_upper_bound_m': spacing,
            'independent_experimental_radius_limit': False,
        },
        'sources': {
            'beam_energy': 'https://cds.cern.ch/record/507409/files/mpph123.pdf',
            'beam_energy_description': 'Brown et al., Ultimate Performance of the LEP RF System, PAC 2001: beam energy up to 104.5 GeV.',
            'candidate': 'Fable L07 Hm source and its stated electron anchor; comparison data are inputs.'},
        'checks': {'flat_anchor_energy_range_insufficient': rows[0]['cannot_span_104500_MeV_in_this_one_load_source'],
                   'small_mass_case_not_excluded_by_this_bound': not rows[-1]['cannot_span_104500_MeV_in_this_one_load_source'],
                   'reported_4p6e18_scale_is_conditional_BZ_arithmetic': abs(spacing / 4.6e-18 - 1) < .01},
        'scope': 'Conditional obstruction for the named one-load energy identification. Passing this necessary bound does not validate a state, dispersion, electron, or row. No experimental compositeness limit or general OPH exclusion is asserted.',
        'auditor_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'official_rows_changed': [],
    }


if __name__ == '__main__':
    result = run()
    if not all(result['checks'].values()):
        raise AssertionError(result['checks'])
    Path(__file__).with_name('energy_scale_results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
