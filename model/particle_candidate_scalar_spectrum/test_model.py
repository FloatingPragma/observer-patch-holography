"""Independent small-matrix and complete-clock controls for the radial bound."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import model as m


def test_plateau_support_and_exact_transition():
    for cutoff in (F(1, 40), F(2, 3), F(5)):
        for x in (F(), cutoff/2, cutoff):
            assert m.plateau(x, cutoff) == m.plateau(-x, cutoff) == 1
        assert m.plateau(3*cutoff/2, cutoff) == F(1, 2)
        assert m.plateau(2*cutoff, cutoff) == m.plateau(5*cutoff, cutoff) == 0
        assert all(0 <= m.plateau(cutoff*F(k, 37), cutoff) <= 1 for k in range(-100, 101))
    with pytest.raises(ValueError):
        m.plateau(1, 0)


def test_kernel_inverse_transform_and_derivative_bounds():
    # Integrate the independently specified triangle in frequency space.
    from scipy.integrate import quad
    for time in (.1, .7, 2, 9):
        inverse = quad(lambda w: (1-abs(w))*np.cos(w*time)/(2*np.pi), -1, 1)[0]
        assert abs(inverse-(1-np.cos(time))/(np.pi*time*time)) < 2e-15
    def derivative(time):
        if abs(time) < 1e-5:
            return -time/(12*np.pi)
        return np.sin(time)/(np.pi*time*time)-2*(1-np.cos(time))/(np.pi*time**3)
    center = 2*quad(lambda t: abs(derivative(t)), 0, 1)[0]
    finite_tail = 2*quad(lambda t: abs(derivative(t)), 1, 50, limit=200)[0]
    assert center < 1/9 and finite_tail < 2
    cert = m.derivative_certificate()
    assert cert['unit_kernel_derivative_L1_upper'] == F(19,9)
    assert cert['plateau_derivative_L1_coefficient_upper'] == 15
    assert cert['auxiliary_filter_is_executed'] is False


def test_frequency_cut_trace_bound_and_passivity_on_independent_finite_system():
    # This is a separate four-level theorem control, never a replacement source.
    E = np.array([0., .03, 1., 4.])
    H = np.diag(E)
    v = np.array([1., 2., 3., 4.]); v /= np.linalg.norm(v)
    A = .2*np.eye(4)+.7*np.outer(v, v)
    weights = (5-E)**3; weights /= weights.sum()
    rho = np.diag(weights); c = weights.min()
    comm = H@A-A@H; d0 = np.linalg.norm(comm, 2)
    lam = d0/40; D = np.ptp(E)
    frequency = E[:,None]-E[None,:]
    multiplier = np.array([[float(m.plateau(F(str(w)), F(str(lam)))) for w in row]
                           for row in frequency])
    T = multiplier*A; B = A-T
    assert np.linalg.norm(H@T-T@H, 2) <= 15*lam+1e-14
    assert np.linalg.norm(B, 2) >= (d0-15*lam)/D-1e-14
    mu = abs(A)**2*weights[None,:]
    high = abs(frequency) > lam
    positive = frequency > lam
    negative = frequency < -lam
    weighted = np.trace(rho@B@B).real
    assert mu[high].sum() >= weighted-1e-14
    assert weighted >= c*np.linalg.norm(B, 2)**2-1e-14
    assert mu[positive].sum() >= mu[negative].sum()-1e-14
    assert mu[positive].sum() >= c*d0*d0/(8*D*D)-1e-14
    elastic = mu[frequency == 0].sum()
    assert elastic >= np.trace(rho@A).real**2-1e-14
    # Reversing the population invalidates the positive-half premise.
    reversed_mu = abs(A)**2*weights[::-1][None,:]
    assert reversed_mu[positive].sum() < reversed_mu[negative].sum()


@pytest.mark.parametrize('side', [2, 3, 5])
def test_actual_open_source_degree_count_against_explicit_adjacency(side):
    observed = {}
    edges = set()
    for cell in product(range(side), repeat=3):
        neighbors = []
        for axis in range(3):
            for step in (-1, 1):
                nxt = list(cell); nxt[axis] += step
                if 0 <= nxt[axis] < side:
                    nxt = tuple(nxt); neighbors.append(nxt)
                    edges.add(tuple(sorted((cell, nxt))))
        observed[len(neighbors)] = observed.get(len(neighbors), 0)+1
    assert {r['degree']:r['cells'] for r in m.degree_census(side)} == observed
    assert sum(degree*count for degree,count in observed.items()) == 2*len(edges)


def test_reverse_native_Wilson_hop_has_same_beta_changing_norm():
    # The original source uses Wilson spin matrices; conjugation exchanges
    # incoming/outgoing entries without changing their beta-changing norms.
    for axis in range(3):
        matrix = np.array(m.radial.transport.spin_blocks()[axis], complex)
        assert np.linalg.norm(matrix[:2,2:]) == pytest.approx(np.linalg.norm(matrix.conj().T[:2,2:]))


def test_exact_dyadic_ceiling_and_small_clock_resource_inequalities():
    for numerator in range(1, 24):
        for denominator in range(1, 24):
            value = F(numerator, denominator)
            k = m.ceil_log2(value)
            power = F(2)**k
            assert power/2 < value <= power
    R, lam, q, N = F(7), F(1,40), F(1,23), 4
    plan = m.clock_plan(R,lam,q,N)
    J = 2**plan['clock_bins']['exponent']
    r = 2**plan['Cayley_substeps_per_base_tick']['exponent']
    d = plan['turn_margin']; beta = q/(16*2**N)
    tail = min(F(1), 1/(2*(J*d)**2)+1/(2*J*d))
    assert tail <= beta
    assert F(J-1,768*r*r) <= beta
    assert plan['total_two_clock_event_loss_upper'] == m.dyadic(q/4,N)


def test_huge_resource_exponents_do_not_allocate_a_huge_integer():
    N = 10**100
    plan = m.clock_plan(17,F(1,40),F(1,10000),N)
    assert 0 < plan['clock_qubits_each']-N < 100
    assert 0 < plan['Cayley_substeps_per_base_tick']['exponent']-N < 100
    assert not plan['huge_register_or_histogram_allocated']
    for args in ((0,1,F(1,3),2),(1,0,F(1,3),2),(1,1,1,2),(1,1,F(1,3),0)):
        with pytest.raises(ValueError):
            m.clock_plan(*args)


@pytest.fixture(scope='module')
def source_bounds():
    return {n:m.bind_member(m.contract.member(n)) for n in (0,1,2,8)}


def test_actual_primary_source_band_complete_records_and_scope(source_bounds):
    for n, bound in source_bounds.items():
        result = bound['result']; member = m.contract.member(n)
        assert result['member_id'] == member['member_id']
        assert result['frequency_cut']['cutoff'] == '1/40'
        assert result['reported_positive_transfer_interval'][0] == '3/160'
        floor = result['source_density_floor']['negative_exponent']
        assert result['original_ideal_positive_transfer_weight_lower']['negative_exponent'] == floor
        assert result['original_reported_positive_transfer_weight_lower']['negative_exponent'] == floor
        q = F(result['original_ideal_positive_transfer_weight_lower']['coefficient'])
        assert F(result['original_reported_positive_transfer_weight_lower']['coefficient']) == 3*q/4
        assert result['complete_clock_plan']['clock_qubits_each'] > floor
        assert result['complete_probe']['fine_outcomes']['coefficient'] == member['source']['side']**3+1
        assert result['same_source_state_reference_clock']
        assert result['selected_slope_is_positive']
        assert F(result['original_elastic_weight_lower']) >= F(3481,160000)
        assert result['physical_particle_map'] == result['physical_mass_table_GeV'] == []
        assert not result['isolated_scalar_pole_or_mass_gap_proved']
        assert not result['uniform_compact_positive_energy_weight_proved']
        assert not result['complete_probe']['initial_reference_joint_state_claimed_unchanged']
    assert F(source_bounds[0]['result']['original_elastic_weight_lower']) > F(7,10)


def test_replayed_consumer_rejects_band_clock_or_species_promotion(source_bounds, monkeypatch):
    member = m.contract.member(0); original = source_bounds[0]
    # Replay behavior is separate from the source arithmetic checked above.
    monkeypatch.setattr(m, 'bind_member', lambda actual,cutoff='1/40': original)
    m.verify_bound(member,original)
    for path,value in ((('reported_positive_transfer_interval',0),'0'),
                       (('physical_Higgs_field_identified',),True),
                       (('complete_clock_plan','source_clock_changed'),True),
                       (('original_reported_positive_transfer_weight_lower','coefficient'),'1')):
        changed = deepcopy(original); node=changed['result']
        for key in path[:-1]: node=node[key]
        node[path[-1]]=value
        changed['result_sha256']=m.contract.digest(changed['result'])
        with pytest.raises(ValueError):
            m.verify_bound(member,changed)


def test_inadmissible_frequency_cutoff_fails():
    with pytest.raises(ValueError):
        m.bind_member(m.contract.member(0),'1/20')
