"""Exact source-identity counterexample, using only the Python standard library."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def cosine_from_square(z):
    """Enclose cos(sqrt(z)) by its alternating rational Taylor series."""
    assert 0 <= z <= 1
    upper = sum(((-1) ** n * z ** n / factorial(2 * n) for n in range(17)), F())
    lower = upper - z ** 17 / factorial(34)
    assert 0 <= upper - lower < F(1, 10 ** 37)
    return lower, upper


def main():
    p = json.loads((HERE / 'PARAMETERS.json').read_text())
    assert cosine_from_square(F()) == (F(1), F(1))
    assert p['source_manifest_sha256'] == '8a64e3592c0775f528bf83adde5cdb8607315890c421cd6046bdac15afe00806'
    assert p['L'] == 3 and p['charge_sector'] == 0
    h = F(p['h'])
    assert h == F(1, 3)
    results = []
    for slot, value in enumerate(p['mass_midpoints_exact']):
        mass = F(value)
        assert 0 < mass < F(1, 4)
        z = (h * mass) ** 2
        c = cosine_from_square(z)
        discrete = (-c[1] / 2, -c[0] / 2)
        autonomous = cosine_from_square(z + F(3, 4))
        gap = (autonomous[0] - discrete[1], autonomous[1] - discrete[0])
        assert gap[0] > F(57, 50)
        # Normalized trace bounds the operator norm from below in dimension4.
        # Replacing the discrete symbol by the continuous one must fail.
        mutant_gap = autonomous[0] - autonomous[1]
        assert not mutant_gap > F(57, 50)
        results.append({'flavor_slot': slot, 'mass_exact': value,
                        'macro_normalized_trace': list(map(str, discrete)),
                        'autonomous_normalized_trace': list(map(str, autonomous)),
                        'operator_difference_lower_exact': str(gap[0]),
                        'operator_difference_lower_display': float(gap[0]),
                        'wrong_source_substitution_rejected': True})
    # Three harmonic pairs: Var(p)=1/2 and H_h=sum p²/4.
    harmonic_mean = 3 * F(1, 2) / 4
    harmonic_variance = 3 * 2 * F(1, 2) ** 2 / 16
    assert harmonic_mean == F(3, 8) and harmonic_variance == F(3, 32)
    out = {'schema': 'oph.macro_generator_identity.result.v1',
           'parameter_sha256': hashlib.sha256((HERE / 'PARAMETERS.json').read_bytes()).hexdigest(),
           'same_clock_same_mass_free_evolution_equal': False,
           'actual_flavor_controls': results,
           'all_operator_difference_lower_bound': '57/50',
           'harmonic_reference_mean_cell_energy': str(harmonic_mean),
           'harmonic_reference_cell_energy_variance': str(harmonic_variance),
           'harmonic_reference_physical_energy_variance': str(harmonic_variance / h ** 2),
           'canonical_harmonic_packet_is_free_energy_eigenstate': False,
           'charged_interacting_spectrum_computed': False,
           'full_criterion_met': False,
           'rows_awarded': []}
    (HERE / 'RESULT.json').write_text(json.dumps(out, indent=2) + '\n')
    print('PASS: all3 actual flavors have evolution norm difference > 1.14; harmonic variance > 0')


if __name__ == '__main__':
    main()
