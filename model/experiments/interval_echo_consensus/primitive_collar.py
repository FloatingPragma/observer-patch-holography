"""Generate the homogeneous material coefficient from simultaneous native vacuum-collar fringes."""
from collections import Counter
from fractions import Fraction
from math import pi, sqrt
import json

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import cell, federation
from . import record_action


def responding(board, edge, words):
    state = np.zeros(board.n, dtype=int)
    state[edge[0]] = 1
    addresses = record_action.addresses(board, edge)
    changes = record_action.changes(addresses, edge, state, words)
    return tuple(a for a, change in zip(addresses, changes, strict=True) if change)


def vacuum_block(h):
    # Both coin outcomes of the prepared native write are kept coherently.
    cols = np.zeros((law.DIM, 2))
    for j, sign in enumerate((1,-1)):
        for coin in (0,1):
            label = 1+law.LABELS.index((coin,sign))
            cols[label*law.N_SRC+law.SRC.index((0,0)), j] = 1/sqrt(2)
    np.testing.assert_allclose(law.echo() @ cols, cols, rtol=0, atol=3e-15)
    blank = np.eye(law.DIM)[:,law.SRC.index((0,0))]
    np.testing.assert_allclose(h @ blank, 0., rtol=0, atol=2e-13)
    block = cols.T @ h @ cols
    np.testing.assert_allclose(h @ cols, cols @ block, rtol=0, atol=2e-13)
    return block


def fringe(block, addresses):
    """Separate fresh preparations for the two terminal quadratures at each exposure."""
    size = 32
    if len(addresses) != len(set(addresses)) or not 0 < len(addresses) < size//2:
        raise ValueError("Need a nonempty set of distinct physical collars in the unwrapped range")
    chars = np.array([[1.,1.],[1.,-1.]])/sqrt(2)
    eigen = chars.T @ block @ chars
    np.testing.assert_allclose(eigen, np.diag(np.diag(eigen)), atol=2e-13, rtol=0)
    span = (eigen[1,1]-eigen[0,0])*len(addresses)
    if not -size/2 < span < size/2:
        raise ValueError("The sampled frequency can wrap")
    reads = []
    local_product = []
    for tick in range(size):
        u = chars.T @ expm(-2j*pi*tick/size*block) @ chars
        local_even, local_odd = np.diag(u)
        # A sharp sign splits equally into even and odd parity. Local terminal
        # quadratures reconstruct its coherence without a multi-collar cat state.
        single = [(abs(local_even+z*local_odd)**2-abs(local_even-z*local_odd)**2)/4
                  for z in (1.,1j)]
        local_product.append((single[0]-1j*single[1])**len(addresses))
        even, odd = np.diag(u)**len(addresses)
        # The two product characters are exact eigenvectors. Their complementary
        # outcome has zero norm; both interference ports are retained.
        reads.append([[abs(even+phase*odd)**2/4, abs(even-phase*odd)**2/4, 0.]
                      for phase in (1.,1j)])
    reads = np.array(reads)
    np.testing.assert_allclose(reads.sum(axis=-1), 1., rtol=0, atol=2e-12)
    differences = reads[:,:,0]-reads[:,:,1]
    collective = differences[:,0]-1j*differences[:,1]
    np.testing.assert_allclose(collective, local_product, atol=2e-12, rtol=0)
    fourier = np.fft.ifft(local_product)
    peak = int(np.argmax(abs(fourier)))
    ideal = np.zeros(size)
    ideal[peak] = 1
    np.testing.assert_allclose(fourier, ideal, atol=2e-12, rtol=0)
    return peak if peak < size//2 else peak-size, reads


def material_frequency_pencil():
    """Free material Ramsey frequency relative to a native echo of the same total load."""
    sharp = np.eye(law.DIM)[:,law.SRC.index((1,0))]
    other = np.eye(law.DIM)[:,law.SRC.index((0,1))]
    material = (sharp+other)/sqrt(2)
    reference = law.echo() @ sharp
    np.testing.assert_allclose(reference @ material, 0, rtol=0, atol=0)
    rates = []
    for h in law.local_term_parts():
        rate = float(material @ h @ material)
        np.testing.assert_allclose(h @ material, rate*material, rtol=0, atol=3e-15)
        np.testing.assert_allclose(h @ reference, 0., rtol=0, atol=3e-15)
        rates.append(rate)
    # Compute the same two Rayleigh quotients exactly in the native blank reduction.
    swap, square = law._load_ops()
    word = np.zeros(law.N_SRC, dtype=int)
    word[[law.SRC.index((1,0)),law.SRC.index((0,1))]] = 1
    norm = int(word @ word)
    offset = Fraction(int(word @ (word-swap@word)), 2*norm)
    slope = Fraction(int(word @ square @ word), norm)
    np.testing.assert_allclose(rates, [float(offset),float(slope)], rtol=0, atol=3e-15)
    if not slope:
        raise ArithmeticError("This material measurement cannot identify the coefficient")
    return offset, slope


def calibrated_strength():
    """Solve the material/primitive frequency equality without a candidate value or release receipt."""
    h0, h1 = law.local_term_parts()
    block = vacuum_block(h0)
    np.testing.assert_allclose(vacuum_block(h1), np.eye(2), atol=3e-15, rtol=0)
    board = cell()
    words = record_action.signed_table()
    widths = {len(responding(board, e, words)) for e in board.seams}
    if len(widths) != 1:
        raise ArithmeticError("The elementary contact calibration depends on the chosen seam")
    frequency, _ = fringe(block, responding(board, board.seams[0], words))
    offset, slope = material_frequency_pencil()
    return float((frequency-offset)/slope)


def material_term():
    """The source definition specializes the native family on the whole 99-state seam."""
    h0, h1 = law.local_term_parts()
    return h0 + calibrated_strength()*h1


def controls():
    h0, h1 = law.local_term_parts()
    board = federation(0)
    words = record_action.signed_table()
    counts, settings = Counter(), {}
    for edge in board.seams:
        sensors = responding(board, edge, words)
        counts[len(sensors)] += 1
        settings.setdefault(len(sensors), sensors)
    measured = {}
    block = vacuum_block(h0)
    for size, sensors in sorted(settings.items()):
        frequency, reads = fringe(block, sensors)
        measured[size] = dict(seams=counts[size], frequency=frequency, quadratures=reads.tolist())
    reference = responding(cell(), cell().seams[0], words)
    unit, raw = fringe(block, reference)
    donor_frequencies = []
    offset, slope = material_frequency_pencil()
    for k in (0,.5,1,2,4,8,16,64,128,-8):
        f, observed = fringe(vacuum_block(h0+k*h1), reference)
        np.testing.assert_allclose(observed, raw, atol=2e-12, rtol=0)
        donor_frequencies.append(dict(strength=k, vacuum_frequency=f, material_frequency=float(offset+k*slope),
                                      matches_standard=(offset+k*slope == f)))
    half, _ = fringe(vacuum_block(h0/2), reference)
    if half*2 != unit:
        raise ArithmeticError("Halving the relative collar rate must change the calibration")
    _, single_reads = fringe(block, reference[:1])
    return dict(primitive=unit, donors=donor_frequencies, single_collar_quadratures=single_reads.tolist(),
                half_parity_rate=half, assembled=measured)


if __name__ == "__main__":
    print(json.dumps(controls(), indent=2, allow_nan=False))
