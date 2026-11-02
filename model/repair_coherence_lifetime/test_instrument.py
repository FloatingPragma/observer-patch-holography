"""Independent dense instruments verify the complete-history norm contract."""
from fractions import Fraction as F

import numpy as np
import pytest

from model.defect_field_bridge.full_domain import branch_action


def instrument(law, edge, weight):
    for coin in (0, 1):
        for raw in range(4 if law == 'C' else 8):
            K = np.zeros((8, 8), complex)
            for mask in range(8):
                label, target, phase, _ = branch_action(mask, law, coin, sites=3, flavors=1, seam=edge)
                if label == raw:
                    K[target, mask] = float(weight/2)**0.5*phase
            yield (edge, coin, raw), np.kron(K, np.eye(2))


def coherence(rho):
    return sum(np.linalg.svd(rho[2*i:2*i+2, 2*j:2*j+2], compute_uv=False).sum()
               for i in range(8) for j in range(8) if i != j)


def initial():
    psi = np.zeros(16, complex)
    psi[2], psi[5], psi[8] = 1, 1j, -1
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


@pytest.mark.parametrize('law', ['C', 'E'])
def test_adaptive_all_record_bound_with_entangled_reference(law):
    histories = [((), initial())]
    previous = sum(coherence(rho) for _, rho in histories)
    for step in range(3):
        following = []
        for history, rho in histories:
            # A genuine outcome-dependent selector, with minimum cut1/3.
            left = F(1, 3) if not history or history[-1][-1] % 2 else F(2, 3)
            for edge, weight in [((0, 1), left), ((1, 2), 1-left)]:
                for label, K in instrument(law, edge, weight):
                    branch = K@rho@K.conj().T
                    if np.trace(branch).real:
                        following.append((history+(label,), branch))
        assert sum(np.trace(rho).real for _, rho in following) == pytest.approx(1, abs=1e-13)
        current = sum(coherence(rho) for _, rho in following)
        assert current <= previous*2/3+1e-12
        histories, previous = following, current


def test_postselection_cannot_drop_the_probability_cost():
    psi = np.zeros(16, complex)
    psi[2] = psi[8] = 2**-0.5  # source cat between sites0 and2
    rho = np.outer(psi, psi.conj())
    K = dict(instrument('C', (0, 1), F(1, 2)))[((0, 1), 1, 0)]
    branch = K@rho@K.conj().T
    probability = np.trace(branch).real
    assert probability == pytest.approx(1/4)
    assert coherence(branch/probability) == pytest.approx(coherence(rho))
    assert coherence(branch) == pytest.approx(coherence(rho)/4)


def test_full_nonselective_source_reference_converges():
    rho = initial()
    target_source = np.diag([0., 1/3, 1/3, 0., 1/3, 0., 0., 0.])
    reference = np.trace(rho.reshape(8, 2, 8, 2), axis1=0, axis2=2)
    target = np.kron(target_source, reference)
    operations = [K for edge in ((0, 1), (1, 2)) for _, K in instrument('C', edge, F(1, 2))]
    assert np.linalg.norm(sum(K.conj().T@K for K in operations)-np.eye(16)) < 1e-13
    for tick in range(1, 21):
        rho = sum(K@rho@K.conj().T for K in operations)
        distance = np.linalg.svd(rho-target, compute_uv=False).sum()
        # Path3 random-walk gap1/4 and fixed-N cut floor1/2.
        assert distance <= 2**0.5*(3/4)**tick+2*(1/2)**tick+1e-12
