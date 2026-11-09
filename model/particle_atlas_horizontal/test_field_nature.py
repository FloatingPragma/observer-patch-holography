"""Checks of fixed-action particle meanings, including countercontrols."""
import json
from pathlib import Path

import sympy as s

from model.particle_atlas_horizontal import field_nature as f
from model.particle_candidate_atlas import conditional_standard_model as owner


def test_exact_owner_hessians_and_neutral_null_direction():
    data = f.electroweak_algebra()
    q = data['owner']
    g, gp, v, lam = q['symbols']
    expected = v*v*s.Matrix([[g*g, 0, 0, 0], [0, g*g, 0, 0],
                             [0, 0, g*g, -g*gp], [0, 0, -g*gp, gp*gp]])/4
    assert q['gauge_mass_squared'] == expected
    assert q['scalar_mass_squared'] == s.diag(0, 0, 0, 2*lam*v*v)
    assert data['neutral_mass_squared'].rank() == 1
    assert data['neutral_mass_squared'].subs(v, 0) == s.zeros(2)
    assert f._clean(q['gauge_mass_squared']*data['photon_normalized']) == s.zeros(4, 1)
    assert data['charged_mass_squared'] == g*g*v*v*s.eye(2)/4


def test_declared_couplings_change_masses_without_changing_photon_nullity():
    q = f.electroweak_algebra()['owner']
    g, gp, v, lam = q['symbols']
    # Gauge symmetry and field counts do not determine numeric masses.
    for inputs in ({g: 1, gp: 1, v: 1, lam: 1},
                   {g: 2, gp: 3, v: 4, lam: 5}):
        m = q['gauge_mass_squared'].subs(inputs)
        assert m.rank() == 3
        assert m*q['photon_vector_unnormalized'].subs(inputs) == s.zeros(4, 1)
    assert q['m_W_squared'].subs({g: 1, v: 1}) == s.Rational(1, 4)
    assert q['m_W_squared'].subs({g: 2, v: 4}) == 16
    assert q['m_H_squared'].subs({lam: 5, v: 4}) == 160


def test_color_result_uses_native_noncommuting_currents():
    color = f.source_color_algebra()
    assert color['independent_color_currents'] == 8
    assert color['noncommuting_unordered_pairs'] > 0
    assert color['color_action_on_Higgs_zero']
    assert color['Higgs_generated_color_mass_squared'] == s.zeros(8)


def test_mass_means_rest_energy_and_chirality_coupling():
    d = f.dirac_algebra()
    assert d['owner']['massive_Dirac_spacetime_Clifford_checks'] == 16
    assert f._clean(d['hamiltonian']**2-d['energy_squared']*s.eye(4)) == s.zeros(4)
    mu = d['rest_energy']
    zero_p = {x: 0 for x in d['energy_squared'].free_symbols if str(x).startswith('p_')}
    assert d['energy_squared'].subs(zero_p) == mu*mu
    massless = d['hamiltonian'].subs(mu, 0)
    upper = s.diag(1, 1, 0, 0)
    assert massless*upper == upper*massless
    assert d['hamiltonian']*upper != upper*d['hamiltonian']


def test_light_speed_statement_excludes_lattice_counterexample():
    speed = f.speed_algebra()
    c = next(x for x in speed['continuum_group_speed'].free_symbols if str(x) == 'c')
    assert speed['continuum_group_speed'] == c
    assert speed['small_spacing_speed_limit'] == c
    substitutions = {x: 1 for x in speed['lattice_ray_group_speed'].free_symbols}
    substitutions[next(x for x in speed['lattice_ray_group_speed'].free_symbols if str(x) == 'p')] = s.pi/3
    assert s.simplify(speed['lattice_ray_group_speed'].subs(substitutions)) == s.Rational(1, 2)


def test_all_targets_and_antiparticles_keep_owner_semantics():
    rows = f.family_rows()
    originals = {r['target']: r for r in owner.benchmark_rows()}
    assert len(rows) == 18
    assert {r['target'] for r in rows} == set(originals)
    for row in rows:
        for key in ('spin', 'charge', 'color', 'antiparticle', 'tree_mass'):
            assert row[key] == originals[row['target']][key]
        assert row['source_predicted_physical_mass'] is None
        assert row['numerical_mass_input'] is None
        assert row['links'] and row['antiparticle_meaning']
    neutrinos = [r for r in rows if r['family'] == 'neutrinos']
    assert len(neutrinos) == 3
    assert all(r['tree_mass'] is None for r in neutrinos)


def test_certificate_and_frozen_content_do_not_equate_distinct_higgs_domains():
    receipt = f.certificate()
    assert receipt['empirical_mass_inputs'] == []
    assert receipt['physical_particle_assignments'] == 0
    assert receipt['bosonic_degrees'] == 28
    assert receipt['finite_Higgs_vs_continuum']['identified_with_measured_Higgs'] is False
    saved = json.loads(Path(f.__file__).with_suffix('.json').read_text())
    assert saved == receipt


def test_all_local_content_links_resolve():
    root = Path(__file__).resolve().parents[2]
    receipt = f.certificate()
    for path in receipt['source_inputs']:
        assert (root/path).is_file(), path
    for row in receipt['families']:
        for reference in row['links']:
            path = root/reference['path']
            assert path.is_file(), reference
            if reference['anchor'].startswith('L') and reference['anchor'][1:].isdigit():
                assert len(path.read_text().splitlines()) >= int(reference['anchor'][1:])
