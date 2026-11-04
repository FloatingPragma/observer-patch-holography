"""Complete record reversibility and independent entropy/work controls."""
import hashlib
import math

import numpy as np

from ..matter_geometry_response import source
from ..observer import native


def entropy(state):
    values = np.linalg.eigvalsh(state)
    if values.min() < -1e-12 or abs(values.sum()-1) > 2e-12:
        raise ValueError('a normalized positive density matrix is required')
    positive = values[values > 0]
    return float(-np.dot(positive, np.log(positive)))


def relative_entropy(state, reference):
    values, vectors = np.linalg.eigh(reference)
    if values.min() <= 0:
        raise ValueError('relative-entropy reference must be faithful')
    log_reference = (vectors*np.log(values))@vectors.conj().T
    return -entropy(state)-float(np.trace(state@log_reference).real)


def binary_entropy(p):
    return sum(-x*math.log(x) for x in (p, 1-p) if x > 0)


def native_record_certificate():
    states = native.states()
    digest = hashlib.sha256()
    count = 0
    for bit in (0, 1):
        for word in range(256):
            edge = source.selected_seam(bit, word)
            for coin in (0, 1):
                images = set()
                for state in states:
                    if edge is None:
                        output, current = state, 0
                    else:
                        output, (_, current) = native.native(state, edge, coin)
                    recovered = list(output)
                    if edge is not None:
                        recovered[edge[0]] += current
                        recovered[edge[1]] -= current
                    if tuple(recovered) != state:
                        raise ArithmeticError('complete current record does not recover its native input')
                    key = (output, current)
                    if key in images:
                        raise ArithmeticError('complete native record map is not injective')
                    images.add(key)
                    digest.update((repr((bit, word, coin, state, output, current))+'\n').encode())
                    count += 1
    return {'native_input_states': len(states), 'complete_histories_checked': count,
            'record_map_sha256': digest.hexdigest(),
            'fixed_coin_and_pointer_native_map_injective': True,
            'all_scheduler_words_coin_bits_idle_and_zero_currents_retained': True,
            'independent_native_input_law_required': True,
            'global_complete_entropy_increment': 'H(B)-I(R:B)+9*log(2)',
            'symmetric_regional_complete_entropy_increment': 'h(1/50)+9*log(2)',
            'extra_current_entropy_added': False,
            'rounded_parity_instrument_failure_outputs_deleted': False,
            'formula_domain': 'ideal proper-code parity instrument; all classical records committed'}


def kraus(parity):
    identity = np.eye(len(parity))
    return (7*identity-6*parity)/math.sqrt(50), 1j*(identity+6*parity)/math.sqrt(50)


def quantum_record_control():
    amplitude = np.array([[2, 1j, 0, 1], [1, 2, 1j, 0],
                          [0, 1, 2, 1j], [1j, 0, 1, 3]], complex)
    amplitude /= np.linalg.norm(amplitude)
    rho = amplitude@amplitude.conj().T
    reference = amplitude.conj().T@amplitude
    D = np.diag([1., 0., 0., 1.])
    probs, branch_entropies, reference_entropies = [], [], []
    for operator in kraus(D):
        branch = operator@amplitude
        probability = float(np.vdot(branch, branch).real)
        probs.append(probability)
        branch_entropies.append(entropy(branch@branch.conj().T/probability))
        reference_entropies.append(entropy(branch.conj().T@branch/probability))
    output_entropy = binary_entropy(probs[1])+np.dot(probs, branch_entropies)
    information = entropy(reference)-np.dot(probs, reference_entropies)
    residual = abs(output_entropy-entropy(rho)-binary_entropy(probs[1])+information)
    phi = sum(operator@rho@operator.conj().T for operator in kraus(D))
    offdiagonal_control = np.linalg.norm(phi-rho)
    if residual > 2e-12 or offdiagonal_control <= 0.1:
        raise ArithmeticError('complete-record identity or arbitrary-state coherence control failed')
    local_residuals = []
    for q in (0.25, 0.5, 0.75):
        local = np.diag([q, 1-q, 1-q, q])/2
        branches = [operator@local@operator.conj().T for operator in kraus(D)]
        probabilities = [float(np.trace(branch).real) for branch in branches]
        output = binary_entropy(probabilities[1])+sum(
            p*entropy(branch/p) for p, branch in zip(probabilities, branches, strict=True))
        local_residuals.append(abs(output-entropy(local)-binary_entropy(1/50)))
    if max(local_residuals) > 2e-12:
        raise ArithmeticError('regional complete pointer entropy is incorrect')
    return {'full_field_reference_dimension': [4, 4],
            'correlated_input_record_identity_residual': float(residual),
            'arbitrary_input_unread_coherence_change': float(offdiagonal_control),
            'regional_record_entropy_residual': float(max(local_residuals)),
            'local_entropy_identity_substituted_for_full_source': False}


def _annihilator(modes, mode):
    value = np.zeros((2**modes, 2**modes))
    for mask in range(2**modes):
        if mask & (1 << mode):
            value[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
    return value


def chain_control():
    """Independent complete 64-state source control, without a core truncation."""
    X = np.array([[0, 1], [1, 0]], complex)
    h = np.zeros((6, 6), complex)
    h[:2, :2], h[4:, 4:] = 64*X, -64*X
    for left in (0, 2):
        h[left:left+2, left+2:left+4] = 0.5j*X
        h[left+2:left+4, left:left+2] = -0.5j*X
    annihilators = [_annihilator(6, j) for j in range(6)]
    H, star = np.zeros((64, 64), complex), np.zeros((64, 64), complex)
    for i in range(6):
        for j in range(6):
            term = h[i, j]*annihilators[i].T@annihilators[j]
            H += term
            if ((i in (2, 3)) != (j in (2, 3))):
                star += term
    I = np.eye(64)
    occupations = [a.T@a for a in annihilators]
    for left in (0, 2, 4):
        H += (occupations[left]-I/2)@(occupations[left+1]-I/2)/65536
    D = (occupations[2]+occupations[3]-I)@(occupations[2]+occupations[3]-I)
    Z = I-2*D
    operators = kraus(D)
    adjoint = sum(operator.conj().T@H@operator for operator in operators)
    operator_residual = float(np.linalg.norm(adjoint-H+(18/25)*star, 2))
    if operator_residual > 2e-12 or np.linalg.norm(Z@star@Z+star) > 1e-12:
        raise ArithmeticError('the complete-Fock source work identity failed')
    energies, vectors = np.linalg.eigh(H)
    beta = 0.01
    probabilities = np.exp(-beta*(energies-energies.min()))
    probabilities /= probabilities.sum()
    rho = (vectors*probabilities)@vectors.conj().T
    after = sum(operator@rho@operator.conj().T for operator in operators)
    work = float(np.trace(H@(after-rho)).real)
    change = entropy(after)-entropy(rho)
    divergence = relative_entropy(after, rho)
    residual = abs(beta*work-change-divergence)
    q_before, q_after = float(np.trace(D@rho).real), float(np.trace(D@after).real)
    if (residual > 2e-12 or work <= 0 or change <= 0 or divergence < -1e-12
            or abs(q_before-q_after) > 2e-12):
        raise ArithmeticError('same-H Gibbs work/entropy inequality failed')
    return {'control_sites': 3, 'full_Fock_rows': 64,
            'mass': 64, 'interaction': '1/65536', 'inverse_temperature': '1/100',
            'complete_Fock_work_operator_residual': operator_residual,
            'Gibbs_work_entropy_identity_residual': residual,
            'readout_work': work, 'unread_entropy_increase': change,
            'relative_entropy_to_same_H_Gibbs': divergence,
            'parity_before': q_before, 'parity_after': q_after,
            '27_site_theorem_replaced_by_chain_control': False,
            'work_identified_as_heat': False}


def certificate():
    return {'native_records': native_record_certificate(),
            'quantum_records': quantum_record_control(),
            'full_Fock_source_control': chain_control(),
            'exact_unread_channel': '(16/25)*rho+(9/25)*(I-2D)*rho*(I-2D)',
            'exact_readout_work': '-(18/25)*Tr(H_star*rho)',
            'full_Gibbs_balance': 'beta*W_read-DeltaS=D(Phi(gamma_beta)||gamma_beta)',
            'stationary_affine_family_assumed_Gibbs': False,
            'apparatus_setup_reset_energy_included': False}
