"""Independent source, analytic-enclosure and continuum-domain controls."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import subprocess
import sys

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import expm

from ..interacting_field.lattice import symbol
from . import sector as s, packet


def m(value):
    value = F(value)
    return mp.mpf(value.numerator)/value.denominator


def test_every_original_source_Fourier_column_and_valley_frame():
    result = s.valley_certificate()
    assert result['exact_Fourier_columns_checked'] == 128
    assert result['positive_chirality_count'] == result['negative_chirality_count'] == 4
    assert [x['jacobian_chirality'] for x in result['valleys']] == [-1, 1, 1, -1, 1, -1, -1, 1]
    assert all(x['continuum_mass'] == '0' for x in result['valleys'])


def test_relativistic_Clifford_algebra_and_scoped_mass_obstruction():
    result = s.relativistic_algebra_certificate()
    assert result['opposite_chirality_Dirac_pair_Clifford_checks'] == 16
    assert result['single_two_component_constant_anticommuting_mass_dimension'] == 0
    assert not result['massive_term_added']


def test_full_outgoing_contact_control_and_negative_vacuum_energy():
    result = s.finite_contact_check()
    assert result['all_N2_states_retained'] == 8128
    assert result['contact_norm_squared'] == '1/8'
    assert result['outside_band_contact_norm_squared'] == '7/64'
    assert .04 < result['full_N2_evolution_error_diagnostic'] < float(F(result['rigorous_vector_error_upper']))
    assert result['normalization_error_diagnostic'] < 1e-12
    assert F(result['one_particle_original_energy_below_vacuum']) < 0


@pytest.mark.parametrize('x', [F(-1), F(-1, 64), F(), F(1, 128), F(1, 64), F(1)])
def test_sine_enclosures_against_high_precision_transcendental(x):
    lo, hi = s.sine_interval(x)
    with mp.workdps(100):
        assert m(lo) <= mp.sin(m(x)) <= m(hi)


@pytest.mark.parametrize('x', [F(), F(1, 64), F(2), F(3), F(1001, 17)])
def test_square_root_enclosure_is_outward(x):
    lo, hi = s.sqrt_interval(x)
    assert lo*lo <= x <= hi*hi
    assert hi-lo <= F(1, 2**80)


def test_exact_energy_intervals_and_mass_shell_defects():
    with mp.workdps(100):
        for row in s.energy_table():
            p, a = list(map(F, row['momentum'])), F(row['spacing'])
            energy = mp.sqrt(sum(mp.sin(m(a*x))**2 for x in p))/m(a)
            continuum = mp.sqrt(sum(m(x)**2 for x in p))
            lo, hi = map(F, row['lattice_positive_energy_interval'])
            assert m(lo) <= energy <= m(hi)
            assert abs(energy-continuum) <= m(row['absolute_dispersion_error_upper'])
            lower, upper = map(F, row['lattice_mass_shell_defect_interval'])
            assert m(lower) <= energy**2-continuum**2 <= m(upper)+mp.mpf('1e-98')
            assert row['lattice_negative_energy_interval'] == [str(-hi), str(-lo)]
            assert row['continuum_mass'] == '0'


@pytest.mark.parametrize('a', [F(1, 4), F(1, 16), F(1, 64)])
def test_projectors_and_propagators_against_original_sine_symbol(a):
    pauli = [np.asarray(x, complex) for x in s.PAULI]
    for bits in product((0, 1), repeat=3):
        for p in [np.array([1., .3, -.4]), np.array([.2, -.7, .5]), np.ones(3)]:
            k = [np.pi*b+float(a)*x for b, x in zip(bits, p, strict=True)]
            lattice = symbol(k, speed=1/float(a))
            continuum = -sum(((-1)**b*x*sigma for b, x, sigma in zip(bits, p, pauli, strict=True)), np.zeros((2, 2), complex))
            E = np.linalg.norm(np.sin(float(a)*p))/float(a)
            Ec = np.linalg.norm(p)
            for sign in [-1, 1]:
                P = (np.eye(2)+sign*lattice/E)/2
                Pc = (np.eye(2)+sign*continuum/Ec)/2
                assert np.linalg.norm(P-Pc, 2) <= float(a*a/6)+1e-13
            assert np.linalg.norm(expm(-.7j*lattice)-expm(-.7j*continuum), 2) <= float(F(7, 10)*a*a/2)+1e-13


def test_infinite_band_fraction_and_complete_error_are_outward():
    b = s.continuum_budget(F(1, 64), 1, 2, 1)
    with mp.workdps(100):
        fraction = 8/(64*mp.pi)**3
        assert m(b['band_fraction_lower']) <= fraction <= m(b['band_fraction_upper'])
        assert m(b['contact_on_band_norm_upper']) >= mp.sqrt(fraction)
    assert b['complete_vector_error_upper'] < F(1237, 10**6)
    assert b['one_particle_projector_error_upper'] == F(1, 24576)
    assert b['complete_vector_error_upper'] == b['contact_on_band_norm_upper']+b['free_generator_error_upper']


def test_number_and_signed_time_domains_with_nonnegative_source_interaction():
    for n in [0, 1]:
        b = s.continuum_budget(F(1, 64), 1, n, -1, interaction=10**100)
        assert b['contact_on_band_norm_upper'] == 0
    for n in [2, 3, 4]:
        positive = s.continuum_budget(F(1, 64), 1, n, 1, interaction=3)
        negative = s.continuum_budget(F(1, 64), 1, n, -1, interaction=3)
        assert positive['complete_vector_error_upper'] == negative['complete_vector_error_upper']
    assert s.continuum_budget(1, 1, 100, 100)['complete_vector_error_upper'] == 2
    assert s.continuum_budget(1, 1, 100, 100)['half_trace_distance_upper'] == 1
    assert s.continuum_budget(F(1, 64), 1, 4, 0)['complete_vector_error_upper'] == 0
    with pytest.raises(ValueError):
        s.continuum_budget(F(1, 64), 1, 2, 1, interaction=-1)


def test_fixed_contact_converges_but_critical_pair_bound_does_not():
    sizes = [16, 32, 64, 128]
    fixed = [s.continuum_budget(F(1, n), 1, 2, 1)['complete_vector_error_upper'] for n in sizes]
    assert all(left > 2*right for left, right in zip(fixed, fixed[1:]))
    critical = [s.continuum_budget(F(1, n), 1, 2, 1, interaction=n*n)['contact_on_band_norm_upper'] for n in sizes]
    assert all(left < right for left, right in zip(critical, critical[1:]))


@pytest.mark.parametrize('args', [(0, 1, 2, 1), (1, 0, 2, 1), (2, 1, 2, 1),
                                  (F(1, 2), 1, -1, 1), (F(1, 2), 1, True, 1),
                                  (.5, 1, 2, 1), (F(1, 2), 1, 2, float('nan'))])
def test_invalid_continuum_domain_rejected(args):
    with pytest.raises(ValueError):
        s.continuum_budget(*args)


def test_invalid_energy_and_interval_domains_rejected():
    for call in [lambda: s.energy_row((1, 2)), lambda: s.energy_row((1, 2, 3), 1),
                 lambda: s.energy_row((1, 2, .5)), lambda: s.sine_interval(F(3, 2)),
                 lambda: s.sqrt_interval(-1), lambda: s.sqrt_interval(1, bits=True)]:
        with pytest.raises(ValueError):
            call()


def test_packet_bindings_and_scope():
    result = packet.run()
    assert result['all_pass']
    assert result['physical_boundary']['physical_species_assignments'] == 0
    assert not result['physical_boundary']['positive_energy_interacting_vacuum_constructed']
    assert set(result['source_sha256']) == set(packet.SOURCE_PATHS)
    assert all(len(value) == 64 for value in result['source_sha256'].values())


def test_replay_rejects_bad_supplied_receipt(tmp_path):
    malformed = tmp_path/'bad.json'
    malformed.write_text('{"schema":"wrong"}')
    result = subprocess.run([sys.executable, '-W', 'error', '-m', 'model.relativistic_sector_audit', '--check', str(malformed)],
                            cwd=packet.ROOT, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert 'Receipt keys mismatch' in result.stderr
