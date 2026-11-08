"""Independent symmetry, contrast, custody and complete-replay controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import sympy as s

from model.centered_source_ground import dimer
from model.defect_field_bridge import full_domain
from model.defect_localization_readout import algebra, bounds, packet


@pytest.fixture(scope='module')
def current():
    return packet.build()


def test_actual_complete_program_and_energy(current):
    run = current['complete_source_execution']
    assert (run['reference_program_pulses'], run['created_program_pulses']) == (1090, 1092)
    assert run['all_failure_labels_and_zero_pointer_slots_retained']
    assert not run['success_renormalized']
    assert F(run['created_complete_isometry_error']) < F(1, 10**8)
    assert F(run['analytic_complete_contrast_lower']) > F(669, 1000)
    row = current['actual_cube_excitation']
    assert F(row['local_X_pointer_probability_contrast'][0]) > F(67, 100)
    assert F(row['same_H_B_excitation_energy'][0]) > F(row['same_H_B_odd_sector_gap_lower']) > 1
    assert not row['zero_field_input_claimed_stationary']
    assert current['scope']['physical_mass_table'] == []


def test_fermionic_double_sign_is_necessary_for_actual_Hamiltonian():
    correct, wrong = s.zeros(16), s.zeros(16)
    for mask in range(16):
        out, sign = algebra.rotation(mask, 2)
        correct[out, mask] = sign
        wrong[out, mask] = 1
    h = dimer.data()['K']
    assert correct*h == h*correct
    assert wrong*h != h*wrong


def test_independent_local_donor_completeness_and_parity():
    # A local spin-odd branch and its charge-even failure retain every input.
    r = s.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, -1]])
    z = s.diag(0, 1, -1, 0)
    failure = s.diag(1, 0, 0, 1)
    assert z*z+failure*failure == s.eye(4)
    assert r*z == -z*r and r*failure == failure*r
    # A charge-preserving unitary cannot flip R on the one-dimensional vacuum
    # block. The complete local spin-odd instrument retains its failure there.
    assert s.trace(r) == 0
    vacuum = s.Matrix([1, 0, 0, 0])
    assert z*vacuum == s.zeros(4, 1) and failure*vacuum == vacuum


def test_all_cube_parities_and_charge_are_distinct(current):
    for word in range(256):
        mask = algebra.spin_mask(word, 8)
        other, sign = algebra.rotation(mask, 8)
        assert other == algebra.spin_mask(word ^ 255, 8) and sign == 1
        assert mask.bit_count()-8 == 0
    assert current['exact_branch_symmetry']['all_seam_full_Fock_columns'] == 6976
    assert current['exact_branch_symmetry']['protected_cube_columns'] == 12288


def test_local_unitary_keeps_arbitrary_entangled_exterior():
    # Three spin sites plus an untouched reference, exact integer amplitudes.
    v = np.array([1, -2, 3, 7, 5, -11, 13, 17, 19, 23, -29, 31, 37, 41, 43, -47], dtype=object)
    z = np.array([1-2*(i & 1) for i in range(16)], dtype=object)
    rho = np.outer(v, v)
    after = np.outer(z*v, z*v)
    # Bit zero is the local spin. The other two spins and reference stay joint.
    reduced = lambda a: np.array([[a[2*i, 2*j]+a[2*i+1, 2*j+1]
                                  for j in range(8)] for i in range(8)], dtype=object)
    assert np.array_equal(reduced(rho), reduced(after))
    assert not np.array_equal(rho, after)


def test_exact_bond_commutator_norm_and_creation_energy_bound():
    x = s.Matrix([[0, 1], [1, 0]])
    y = s.Matrix([[0, -s.I], [s.I, 0]])
    z = s.diag(1, -1)
    h = (sum((s.kronecker_product(a, a) for a in (x, y, z)), s.zeros(4))-s.eye(4))/4
    xx = s.kronecker_product(x, s.eye(2))
    zz = s.kronecker_product(z, s.eye(2))
    comm = h*xx-xx*h
    assert set((comm.conjugate().T*comm).eigenvals()) == {0, 1}
    assert set((zz*h*zz-h).eigenvals()) == {-1, 0, 1}


def test_lifetime_formula_retains_both_protocol_errors():
    actual = bounds.trajectory_example()
    direct = F(67, 100)-3*F(1, 16)-F(8*11, 12**3)-6*F(1, 16)/12-F(2, 10**8)
    assert F(actual['complete_X_pointer_contrast_lower']) == direct > F(2, 5)
    assert bounds.contrast_lower(F(1), F(1), F(1), 3, 1, F(0), F(0)) == 0
    assert bounds.contrast_lower(F(1), F(0), F(1), 0, 1, F(0), F(0), F(1, 10)) == F(4, 5)


@pytest.mark.parametrize('change', ({'initial': True}, {'time': 0.2},
    {'coupling': F(0)}, {'degree': True}, {'free_intervals': 0},
    {'tau_over_U': F(1, 3)}, {'per_protocol_half_error': F(-1)}))
def test_bound_domains(change):
    args = dict(initial=F(1), time=F(1, 16), coupling=F(1), degree=3,
                free_intervals=1, tau_over_U=F(1, 100), tau3_over_U2=F(1, 100),
                per_protocol_half_error=F(0))
    args.update(change)
    with pytest.raises(ValueError):
        bounds.contrast_lower(**args)


def test_actual_branch_mutation_is_rejected(monkeypatch):
    original = full_domain.branch_action
    def altered(mask, law, coin, **kwargs):
        raw, out, sign, current = original(mask, law, coin, **kwargs)
        if kwargs['sites'] == 2 and mask == 1:
            raw = (raw+1) % 16
        return raw, out, sign, current
    monkeypatch.setattr(full_domain, 'branch_action', altered)
    with pytest.raises(ArithmeticError, match='parity'):
        algebra.branch_certificate()


@pytest.mark.parametrize('section,key,value', (
    ('scope', 'topological_defect_or_infinite_volume_superselection_identified', True),
    ('scope', 'physical_mass_table', [{'mass': 1}]),
    ('complete_source_execution', 'ten_repair_lifetime_tape_executed', True),
    ('full_record_lifetime_bound', 'complete_X_pointer_contrast_lower', '1'),
    ('source_bindings', 'model/defect_localization_readout/algebra.py', '0'*64)))
def test_complete_receipt_not_only_source_hashes(current, monkeypatch, tmp_path, section, key, value):
    changed = deepcopy(current)
    changed[section][key] = value
    target = tmp_path/'requested.json'
    target.write_text(packet.canonical(changed))
    monkeypatch.setattr(packet, 'build', lambda: current)
    with pytest.raises(ValueError):
        packet.verify(target)


def test_check_uses_requested_path(current, monkeypatch, tmp_path):
    good = tmp_path/'alternate.json'
    good.write_text(packet.canonical(current))
    monkeypatch.setattr(packet, 'build', lambda: current)
    packet.verify(good)
    good.write_text('{}')
    with pytest.raises(ValueError):
        packet.verify(good)


def test_runtime_import_custody_in_fresh_full_producer():
    code = r'''
import sys
from pathlib import Path
from model.defect_localization_readout import packet
value=packet.build()
root=packet.ROOT.resolve()
observed=set()
for module in tuple(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if path:
        p=Path(path).resolve()
        if p.is_relative_to(root) and p.suffix=='.py':
            observed.add(p.relative_to(root).as_posix())
assert observed-set(value['source_bindings'])==set(), sorted(observed-set(value['source_bindings']))
assert 'model/defect_field_bridge/full_domain.py' in observed
assert not any(p.startswith(('docs/','flagship/')) for p in value['source_bindings'])
'''
    run = subprocess.run([sys.executable, '-c', code], cwd=packet.ROOT,
                         text=True, capture_output=True)
    assert run.returncode == 0, run.stdout+run.stderr
