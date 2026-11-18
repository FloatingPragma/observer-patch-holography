"""Exact capacity and current-domain checks for the native scalar construction.

The small matrices are independent witnesses, not a simulation of 432 bits.
The inherited current basis describes vacuum + color triplet + weak triplet.
It does not select an intrinsic gauge-response menu for the scalar process.
"""
from math import comb

import sympy as sp


LAW = "A-native-scalar-overlap-v3"


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(name + " must be a positive integer")


def fock_capacity(modes):
    """Full Fock space and its two parity blocks, with no occupancy truncation."""
    _positive_integer(modes, "modes")
    block = 1 << (modes - 1)
    return {"modes": modes, "full_dimension": 2 * block,
            "parity_block_dimension": block,
            "full_even_algebra_dimension": 2 * block * block,
            "one_parity_image_dimension": block * block,
            "full_fock_qubits": modes, "one_parity_qubits": modes - 1}


def compiler_capacity():
    """Dimension consequences of the export's connected graph and cycle rank.

    The graph has 174 vertices and 432 edges; its independently checked
    260 local cycle operators have rank 259. This count does not establish
    connectivity or independence for a different graph.
    """
    vertices, edges = 150 + 12 + 12, 240 + 12 * 16
    cycles = edges - vertices + 1
    return {"source_modes": 150, "fuel_modes": 12, "reference_modes": 12,
            "total_modes": vertices, "edge_qubits": edges,
            "local_penalty_cycles": 240 + 20, "independent_cycles": cycles,
            "code_qubits": edges - cycles, "code_dimension": 1 << (edges - cycles),
            "fixed_total_parity": "even", "faithful_full_even_CAR": False}


def singlet_capacity(owners):
    """Dimension of collective singlets in (spin0 + spin1)**owners."""
    _positive_integer(owners, "owners")
    singlets = comb(2 * owners, owners) // (owners + 1)
    return {"owners": owners, "raw_dimension": 7 ** owners,
            "neutral_tensor_dimension": 4 ** owners,
            "collective_singlet_dimension": singlets,
            "requested_binary_code_dimension": 1 << owners,
            "binary_code_fits": (1 << owners) <= singlets}


def _unit(size, row, column):
    out = sp.zeros(size)
    out[row, column] = 1
    return out


def owner_currents():
    """Hermitian coordinates for the declared u(3) + so(3) owner action.

    Indices are vacuum0, color1:4 and Cartesian weak4:7. These coordinates
    span the inherited twelve-current representation; they are not twelve
    additional operations admitted by the scalar experiment.
    """
    color = []
    for i, j in ((1, 2), (1, 3), (2, 3)):
        a, b = _unit(7, i, j), _unit(7, j, i)
        color.extend((a + b, -sp.I * (a - b)))
    color.extend((sp.diag(0, 1, -1, 0, 0, 0, 0),
                  sp.diag(0, 1, 1, -2, 0, 0, 0)))
    weak = tuple(-sp.I * (_unit(7, i, j) - _unit(7, j, i))
                 for i, j in ((5, 6), (6, 4), (4, 5)))
    center = sp.diag(0, 1, 1, 1, 0, 0, 0)
    return tuple(color) + weak + (center,)


def span_rank(matrices):
    matrices = tuple(matrices)
    if not matrices:
        return 0
    return sp.Matrix.hstack(*(a.reshape(a.rows * a.cols, 1) for a in matrices)).rank()


def current_algebra_dimensions():
    currents = owner_currents()
    brackets = tuple(sp.I * (a * b - b * a)
                     for i, a in enumerate(currents) for b in currents[:i])
    rank = span_rank(currents)
    derived = span_rank(brackets)
    # For a compact matrix Lie algebra, its center is the orthogonal
    # complement of the derived algebra under the trace inner product.
    return {"current_rank": rank, "bracket_rank": derived,
            "center_rank": rank - derived,
            "closed_span_rank": span_rank((*currents, *brackets))}


def two_owner_witness():
    """Full 49-state frame with a two-dimensional collective singlet code.

    This is an exact small witness of the domain distinctions. The exported
    ten-owner and fifteen-owner isometries use their specified fusion paths.
    """
    vacuum = sp.eye(49)[:, 0]
    singlet = sum((sp.eye(49)[:, 7 * i + i] for i in (4, 5, 6)), sp.zeros(49, 1)) / sp.sqrt(3)
    code = sp.Matrix.hstack(vacuum, singlet)
    projector = code * code.H
    currents = tuple(sp.kronecker_product(q, sp.eye(7))
                     + sp.kronecker_product(sp.eye(7), q) for q in owner_currents())
    scalar_number = singlet * singlet.H
    charged = sp.eye(49)[:, 7]
    individual_weak = sp.kronecker_product(owner_currents()[10], sp.eye(7))
    changed = individual_weak * singlet
    logical_x = sp.Matrix([[0, 1], [1, 0]])
    raw_x = code * logical_x * code.H + sp.eye(49) - projector
    return {"code": code, "projector": projector, "currents": currents,
            "scalar_number": scalar_number, "charged": charged, "raw_x": raw_x,
            "independent_owner_charge_squared": sp.simplify((changed.H * changed)[0])}


def fock_operators(modes):
    """Exact dense CAR matrices for bounded diagnostic frames only."""
    _positive_integer(modes, "modes")
    if modes > 6:
        raise ValueError("dense diagnostic is limited to six modes")
    size = 1 << modes
    operators = []
    for mode in range(modes):
        c = sp.zeros(size)
        for word in range(size):
            if word & (1 << mode):
                sign = (-1) ** (word & ((1 << mode) - 1)).bit_count()
                c[word ^ (1 << mode), word] = sign
        operators.append(c)
    return tuple(operators)


def parity_kernel_witness(modes=3):
    """A nonzero central even-CAR projection killed by the even code."""
    _positive_integer(modes, "modes")
    if modes > 6:
        raise ValueError("dense diagnostic is limited to six modes")
    size = 1 << modes
    parity = sp.diag(*((-1) ** word.bit_count() for word in range(size)))
    odd_projection = (sp.eye(size) - parity) / 2
    even_words = [word for word in range(size) if word.bit_count() % 2 == 0]
    embedding = sp.eye(size)[:, even_words]
    return {"parity": parity, "odd_projection": odd_projection,
            "even_embedding": embedding,
            "restricted_odd_projection": embedding.H * odd_projection * embedding}


def phase_response_witness():
    """Complete finite records distinguish relative and global number phases.

    Source, fuel and spectator reference start in |0,1,1>, in even total
    parity. A source/fuel splitter and its inverse read a relative phase.
    Global total-number phase has no effect on this fixed-number input.
    """
    source, fuel, reference = fock_operators(3)
    numbers = tuple(c.H * c for c in (source, fuel, reference))
    total = sum(numbers, sp.zeros(8))
    hop = source.H * fuel + fuel.H * source
    splitter = sp.eye(8) + (1 / sp.sqrt(2) - 1) * hop ** 2 - sp.I * hop / sp.sqrt(2)
    initial = sp.eye(8)[:, 6]
    phases = {"none": sp.eye(8), "source_pi": sp.eye(8) - 2 * numbers[0],
              "global_pi_over_two": sp.diag(*(sp.I ** word.bit_count() for word in range(8)))}
    records = {}
    for name, phase in phases.items():
        final = sp.simplify(splitter.H * phase * splitter * initial)
        law = tuple(sp.simplify(sp.conjugate(z) * z) for z in final)
        records[name] = {"complete_basis_law": law,
                         "fuel_one": sum(p for word, p in enumerate(law) if word & 2),
                         "reference_one": sum(p for word, p in enumerate(law) if word & 4)}
    return {"splitter": splitter, "hop": hop, "numbers": numbers,
            "total_number": total, "initial": initial, "records": records}
