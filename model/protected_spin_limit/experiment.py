"""Independent small full-Fock checks; these numerical runs are not proofs."""
import math
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh

from model.defect_field_bridge.full_domain import branch_action, RAW_SLOTS
from model.strong_coupling_refinement.algebra import hop, exchange_four
from .algebra import spin_mask
from .bounds import history


def full_hamiltonian(sites, edges, c, U):
    dimension = 1 << (2*sites)
    H = np.zeros((dimension, dimension), dtype=float)
    for mask in range(dimension):
        pairs = sum(((mask >> (2*x)) & 3) == 3 for x in range(sites))
        H[mask, mask] = U*(pairs-mask.bit_count()/2+sites/2)
        for x, y, sign in edges:
            for spin in (0, 1):
                for source, target in ((2*x+spin, 2*y+spin), (2*y+spin, 2*x+spin)):
                    moved = hop(mask, source, target)
                    if moved:
                        other, phase = moved
                        H[other, mask] += c*sign*phase/2
    return H


def unitary(H, duration):
    energies, vectors = eigh(H)
    return (vectors*np.exp(-1j*duration*energies))@vectors.conj().T


def instrument(vector, law, seam, sites):
    output = {(coin, raw): np.zeros_like(vector, dtype=complex)
              for coin in (0, 1) for raw in range(RAW_SLOTS)}
    for coin in (0, 1):
        for mask in range(1 << (2*sites)):
            raw, other, sign, _ = branch_action(mask, law, coin, sites=sites, flavors=2, seam=seam)
            output[coin, raw][other] += sign*vector[mask]/math.sqrt(2)
    return output


def difference_of_pure_states(u, v):
    a, b = float(np.vdot(u, u).real), float(np.vdot(v, v).real)
    overlap = abs(np.vdot(u, v))**2
    return math.sqrt(max(0., (a+b)**2-4*overlap))


def recorded_history(*, c=512, U=512**2, duration=.5, repairs=2, law='C', adaptive=False):
    sites, edges = 3, ((0, 1, 1), (1, 2, -1))
    P = np.zeros((64, 8), complex)
    for word in range(8):
        P[spin_mask(word, sites), word] = 1
    source = unitary(full_hamiltonian(sites, edges, c, U), duration/(repairs+1))
    target = unitary(exchange_four(sites, edges).toarray()*c*c/(4*U), duration/(repairs+1))
    initial = P/math.sqrt(8)  # Complete normalized Choi vector; reference dimension eight.
    actual, ideal = {(): initial}, {(): initial}
    total_slots = zero_slots = 0
    for tick in range(repairs+1):
        actual = {h: source@v for h, v in actual.items()}
        ideal = {h: P@target@(P.conj().T@v) for h, v in ideal.items()}
        if tick == repairs:
            break
        updated, target_updated = {}, {}
        for h in actual.keys() | ideal.keys():
            last = h[-1][0] if h else 0
            seam = (1, 2) if adaptive and last else (0, 1)
            selected = ('E' if (last or tick % 2) else 'C') if adaptive else law
            v = actual.get(h, np.zeros_like(initial))
            w = ideal.get(h, np.zeros_like(initial))
            for record, output in instrument(v, selected, seam, sites).items():
                total_slots += 1
                target_output = w/math.sqrt(2) if record[1] == 0 else np.zeros_like(w)
                if np.any(output) or np.any(target_output):
                    updated[h+(record,)] = output
                    target_updated[h+(record,)] = target_output
                else:
                    # An exact zero child has only zero descendants under this linear kernel.
                    zero_slots += 1
        actual, ideal = updated, target_updated
    distance = sum(difference_of_pure_states(v, ideal[h]) for h, v in actual.items())
    probability = sum(float(np.vdot(v, v).real) for v in actual.values())
    bad_raw = sum(float(np.vdot(v, v).real) for h, v in actual.items() if any(r[1] for r in h))
    return {'trace_norm_difference': distance, 'complete_probability': probability,
            'nonzero_raw_probability': bad_raw, 'retained_nonzero_or_ideal_histories': len(actual),
            'explicit_outcome_slots': total_slots, 'zero_children': zero_slots,
            'analytic_trace_bound': float(history(c, U, len(edges), F(duration), repairs+1)),
            'source_native_pulses_executed': False}
