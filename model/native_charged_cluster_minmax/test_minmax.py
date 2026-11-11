import json
from fractions import Fraction as F
from pathlib import Path

from .audit import compute


def test_exact_minmax_replay_and_source_completion():
    value = compute()
    assert value == json.loads(Path(__file__).with_name('results.json').read_text())
    assert value['frozen_bridge_receipt'] == '75590a8e62090717d5b91ada6ac472498952c3b4ace19b15de815419ecdadc6c'
    gaps = value['minmax']['full_charged_rank2_gap_by_completion']
    assert F(gaps['21/2']) > F(17175185, 10**8)
    assert F(gaps['23/2']) >= F(gaps['21/2'])
    lower, upper = map(F, value['conservative_actual_energy']['charged_minus_actual_neutral_minus4_over_m'])
    assert F(111134603, 10**8) < lower < upper < F(112278560, 10**8)
    assert value['corner']['actual_free_lowest_projector_rank'] == 2
    assert value['corner']['corner_mirror_times_T_square'] == '+I'
    assert not value['corner']['Kramers_or_exact_actual_doublet_degeneracy_claimed']


def test_perturbative_receipt_preserves_uncertified_transfer():
    value = json.loads(Path(__file__).with_name('self_energy').joinpath('contractions.json').read_text())
    from .self_energy.contractions import calculate
    assert calculate() == value
    status = value['certificate_status']
    assert status['algebraic_coefficient_intervals']
    assert status['exact_candidate_contraction_arithmetic']
    assert not status['compressed_matrix_transfer_allowance_proved']
    assert not status['interacting_low_energy_remainder_certified']
    assert not status['physical_mass_certified']
