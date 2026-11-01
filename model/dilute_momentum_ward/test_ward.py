"""Exact local CAR identities, contact angles and outward error controls."""
from fractions import Fraction as F

import pytest
import sympy as s

from ..coarse_energy_ward import ward as w
from ..periodic_momentum import momentum as p
from . import ward as d, packet


def test_complete_periodic_Ward_source_and_coarse_partitions():
    geometry = p.PeriodicGeometry()
    data = d.columns(geometry)
    assert data['H'] == p.action(geometry)[0]
    assert [len(c['stress']) for c in data['columns']] == [81]*3
    assert len(data['oriented_edges']) == 81
    for column in data['columns']:
        assert all(f == f.adjoint() for f in column['force'])
        assert not w.total(column['force']).is_zero()


@pytest.mark.parametrize('geometry', [p.PeriodicGeometry((3, 1, 1)), p.PeriodicGeometry((1, 3, 1)),
                                    p.PeriodicGeometry((1, 1, 3)), p.PeriodicGeometry((1, 1, 4))])
def test_wrapping_force_identity_on_each_axis(geometry):
    data = d.columns(geometry)
    for column in data['columns']:
        divergence = w.divergence(column['stress'], len(geometry.sites))
        for density, flux, force in zip(column['local'], divergence, column['force'], strict=True):
            assert w.commutator(data['H'], density).scale((0, 1))+flux == force


def test_sharp_N2_contact_force_and_every_total_momentum_sector():
    result = packet.finite_two_particle_certificate()
    assert result['force_norm_squared_in_U_squared_units'] == [0, 1, 2, 1]
    assert result['contact_principal_angles_squared'] == ['1/4', '1/2', '1/4', '0']
    assert result['escape_squared_per_K'] == ['3/16', '1/4', '3/16', '0']


def test_finite_time_cell_retains_nonzero_force():
    result = packet.cell_diagnostic()
    assert result['complete_forced_cell_identity_residual'] < 1e-12
    assert result['unforced_cell_defect_norm'] > 1e-5


def test_N2_occupation_matrix_matches_independent_full_Fock_restriction():
    geometry = p.PeriodicGeometry((1, 1, 3))
    data = d.columns(geometry)
    for operator in [data['H'], data['V'], *data['columns'][2]['force']]:
        matrix, masks = d.occupation_matrix(operator, 6)
        assert matrix == p.clifford_matrix(operator, 6).extract(masks, masks)


def test_contact_overlap_convolution_counts_all_pairs():
    result = d.band_counts(64, 1, 3)
    assert result['fraction'] == F(27, 262144)
    assert len(result['contact_eigenvalues']) == 125
    assert sum(result['contact_eigenvalues'].values()) == F(27**2, 64**3)
    assert result['initial_escape_norm_squared'] == result['fraction']*(1-result['fraction'])


@pytest.mark.parametrize('side,radius,dimension', [(4, 0, 1), (5, 1, 1), (8, 1, 2), (6, 2, 3)])
def test_contact_angle_maximum_is_symmetric_band_fraction(side, radius, dimension):
    result = d.band_counts(side, radius, dimension)
    assert result['contact_eigenvalues'][(0,)*dimension] == result['fraction']
    assert max(result['contact_eigenvalues'].values()) == result['fraction']
    assert result['initial_escape_norm_squared'] <= F(1, 4)


@pytest.mark.parametrize('value', [F(), F(1, 4), F(1, 8), F(27, 262144), F(1001, 17)])
def test_sqrt_majorant_is_outward_and_one_dyadic_step_sharp(value):
    bound = d.sqrt_upper(value)
    assert bound*bound >= value
    if bound:
        assert (bound-F(1, 2**48))**2 < value


def test_bound_improves_previous_crude_leakage_and_keeps_coupling():
    b = d.budget(F(27, 64**3), 2, F(1, 8), cutoff=F(11, 112), block_sites=8)
    assert b['evolved_band_leakage_upper'] < b['previous_leakage_upper']/90
    assert b['local_Ward_cell_defect_upper'] < F(132, 10**6)
    assert b['coarse_block_defect_upper'] == 8*b['local_Ward_cell_defect_upper']
    assert b['rescaled_by_cutoff_momentum_error_upper'] < F(83, 10000)
    high = d.budget(F(27, 64**3), 2, F(1, 8), 64**2)
    assert high['evolved_band_leakage_upper'] == 1
    assert high['local_Ward_cell_defect_upper'] > 1  # no critical-coupling conclusion


def test_zero_contact_and_one_particle_limits_are_exact():
    for b in (d.budget(F(1, 2), 1, 100), d.budget(F(1, 2), 20, 100, 0)):
        assert b['interaction_on_initial_band_norm_upper'] == 0
        assert b['local_Ward_cell_defect_upper'] == 0
        assert b['evolved_band_leakage_upper'] == 0


def test_coarse_error_budget_does_not_assume_cancellation_of_forces():
    a = d.budget(F(1, 64), 4, F(1, 4), block_sites=1)
    b = d.budget(F(1, 64), 4, F(1, 4), block_sites=27)
    assert b['coarse_block_defect_upper'] == 27*a['local_Ward_cell_defect_upper']


def test_later_cells_keep_the_original_preparation_time():
    first = d.budget(F(1, 1024), 2, F(1, 4))
    second = d.budget(F(1, 1024), 2, F(1, 4), start_time=F(1, 4))
    merged = d.budget(F(1, 1024), 2, F(1, 2))
    assert second['local_Ward_cell_defect_upper'] > first['local_Ward_cell_defect_upper']
    assert first['local_Ward_cell_defect_upper']+second['local_Ward_cell_defect_upper'] == merged['local_Ward_cell_defect_upper']
    assert second['evolved_band_leakage_upper'] == merged['evolved_band_leakage_upper']


@pytest.mark.parametrize('args', [(0, 2, 1), (F(3, 2), 2, 1), (0.1, 2, 1),
                                  (F(1, 4), 0, 1), (F(1, 4), True, 1), (F(1, 4), 2, -1)])
def test_nonexact_or_outside_budget_domain_fails_closed(args):
    with pytest.raises(ValueError):
        d.budget(*args)


def test_band_and_smooth_cutoff_domains_are_explicit():
    for args in [(2, 0, 1), (4, 2, 1), (4, 1, 4), (64, 9, 3)]:
        with pytest.raises(ValueError):
            d.band_counts(*args)
    with pytest.raises(ValueError):
        d.budget(F(1, 4), 2, 1, cutoff=F(3, 2))
    with pytest.raises(ValueError):
        d.budget(F(1, 4), 3, 1, cutoff=F(1, 4))


def test_supplied_receipt_path_is_validated(monkeypatch, tmp_path):
    file = tmp_path/'bad.json'
    file.write_text('{"all_pass":false}')
    monkeypatch.setattr(packet, 'run', lambda: {'all_pass': True})
    monkeypatch.setattr('sys.argv', ['dilute_momentum_ward', '--check', str(file)])
    with pytest.raises((AssertionError, ValueError, SystemExit)):
        packet.main()
