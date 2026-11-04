"""Independent full-Fock matrix controls for the local gauge construction."""
import json

import pytest
import sympy as s

from . import algebra, packet


def test_full_fock_nonabelian_ward_identity():
    # Four physical fermion modes, two sites, a five-state link. These
    # matrices are built independently of the sparse one-particle producer.
    d = 2
    annihilators = []
    for k in range(2*d):
        c = s.zeros(16)
        for mask in range(16):
            if mask >> k & 1:
                c[mask ^ (1 << k), mask] = (-1)**((mask & ((1 << k)-1)).bit_count())
        annihilators.append(c)
    link = []
    for a in range(d):
        for b in range(d):
            u = s.zeros(5)
            u[1+d*a+b, 0] = 1
            link.append(u)
    h = s.zeros(80)
    for a in range(d):
        for b in range(d):
            term = s.kronecker_product(annihilators[a].H*annihilators[d+b], link[d*a+b])
            h += term+term.H
    h += s.kronecker_product(s.eye(16), s.diag(0, 1, 1, 1, 1))
    currents = (s.I*s.eye(2), s.Matrix([[0, 1], [-1, 0]]),
                s.I*s.Matrix([[0, 1], [1, 0]]), s.I*s.diag(1, -1))
    for g in currents:
        for site in range(2):
            charge = sum((g[a, b]*annihilators[d*site+a].H*annihilators[d*site+b]
                          for a in range(d) for b in range(d)), s.zeros(16))
            glink = s.zeros(5)
            glink[1:, 1:] = s.kronecker_product(-g.T, s.eye(d)) if site == 0 else s.kronecker_product(s.eye(d), g)
            total = s.kronecker_product(charge, s.eye(5))+s.kronecker_product(s.eye(16), glink)
            assert total*h-h*total == s.zeros(80)


def test_all_actual_native_generators_and_controls():
    result = algebra.exact_certificate()
    assert result['native_current_rank'] == 12
    assert result['matter_dimension'] == 15
    assert result['literal_Gauss_checks'] == 240
    assert result['link_dimension'] == 60
    assert sorted(x['dimension'] for x in result['native_blocks']) == [1, 2, 3, 3, 6]
    assert result['noncommuting_basis_pairs'] > 0
    assert result['failures_when_link_omitted'] > 0
    assert result['failures_when_link_orientation_reversed'] > 0
    assert result['failures_when_static_reference_omitted'] > 0


def test_complete_records_and_detector_work():
    result = algebra.complete_history()
    assert len(result['outcomes']) == 16
    assert sum(s.sympify(x['probability']) for x in result['outcomes']) == 1
    assert result['squared_continuous_gap'] == '25'
    final = s.sympify(result['energy_accounting'][-1]['energy_after_readout'])
    assert final == s.sympify(result['total_detector_work'])
    assert final > 0
    zero = algebra.complete_history(hopping=0, rounds=2)
    assert zero['total_detector_work'] == '0'
    assert sum(s.sympify(x['probability']) for x in zero['outcomes'] if 'L' in x['word']) == 0


@pytest.mark.parametrize('kwargs', [{'duration': .1}, {'link_energy': -1}, {'rounds': 0}, {'rounds': True}, {'d': 0}])
def test_input_domain(kwargs):
    with pytest.raises(ValueError):
        algebra.complete_history(**kwargs)


def test_receipt_rejects_mutated_claim(tmp_path):
    value = packet.produce()
    path = tmp_path/'mutated.json'
    value['scope']['same_compact_defect_band_gauge_dressing'] = True
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='differs'):
        packet.verify(path)
