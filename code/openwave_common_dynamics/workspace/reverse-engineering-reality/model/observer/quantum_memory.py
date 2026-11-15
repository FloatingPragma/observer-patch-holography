"""Exact memory consequences of the uniform native, forgotten-record channel.

The binary-shell stopping count is sufficient for entanglement breaking,
not its exact index.  The full total-two carrier is never entanglement
breaking at a finite attempt count.  Neither statement includes the
discarded record tape, adaptive control, or additional protected source
registers. Both statements account for an untouched external reference.

All certificate arithmetic uses integers and rational numbers.  The
spectral input is the binary mixing theorem in QUANTUM_MIXING.md; the
separability input is Gurvits--Barnum, quant-ph/0204159, Corollary 3.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import isqrt

from . import quantum_mixing as mixing


DIMENSION = 66
TRACELESS_DIAGONALS = 65
OFFDIAGONALS = 4290
COHERENCE_RATES = {"C": F(5, 6), "E": F(7, 10)}


def _candidate(value):
    if type(value) is not str or value not in COHERENCE_RATES:
        raise ValueError("candidate C or E required")
    return value


def _attempts(value):
    if type(value) is not int or value < 0:
        raise ValueError("a nonnegative integer attempt count required")
    return value


def contraction_interval():
    """Rational enclosure of (55 + sqrt(5))/60, checked by squaring."""
    scale = 1000
    root = isqrt(5 * scale**2)
    lower, upper = F(root, scale), F(root + 1, scale)
    if not 0 < lower**2 < 5 < upper**2:
        raise ArithmeticError("strict square-root enclosure failed")
    return (55 + lower) / 60, (55 + upper) / 60


def choi_purity_interval(candidate, attempts):
    """Enclose the spectral *upper bound*, not the actual Choi purity.

    The lower endpoint bounds that same expression from below; it does
    not bound the actual purity from below.  The normalized Choi operator
    lives in dimension 66**2, so its separability threshold is 1/4355.
    """
    candidate, attempts = _candidate(candidate), _attempts(attempts)
    coherence = COHERENCE_RATES[candidate] ** (2 * attempts)
    return tuple((1 + TRACELESS_DIAGONALS * rate ** (2 * attempts)
                  + OFFDIAGONALS * coherence) / DIMENSION**2
                 for rate in contraction_interval())


def stopping_bound(candidate):
    """First count certified by the upper enclosure; later counts also pass."""
    candidate = _candidate(candidate)
    threshold = F(1, DIMENSION**2 - 1)
    attempts = 0
    while choi_purity_interval(candidate, attempts)[1] > threshold:
        attempts += 1
    # Establish minimality for this spectral bound without relying on
    # rounding.  Failure here cannot be interpreted as retained entanglement.
    if attempts and choi_purity_interval(candidate, attempts - 1)[0] <= threshold:
        raise ArithmeticError("the root enclosure cannot resolve the stopping count")
    return attempts


def _rounded(value, *, upper, scale=10**12):
    scaled = value * scale
    quotient = scaled.numerator // scaled.denominator
    if upper and scaled.denominator != 1:
        quotient += 1
    return F(quotient, scale)


def _sector(difference, *, binary=False):
    return next(sector for sector in mixing.sectors(binary=binary)
                if sector.difference == difference)


def _states_on_first_seam():
    u, v = mixing.seams()[0]
    pile, binary = [0] * 12, [0] * 12
    pile[u], binary[u], binary[v] = 2, 1, 1
    return (u, v), tuple(pile), tuple(binary)


def full_carrier_witness(candidate):
    """Check native rows proving a negative partial-transpose minor for all n.

    Input: (|2u>|0> + |u+v>|1>)/sqrt(2), with adjacent u,v.  After
    partial transpose on the reference, the principal minor on
    |2u,1>, |u+v,0> has first diagonal zero and offdiagonal a**n/2.
    Its other diagonal does not enter its determinant -a**(2*n)/4.
    """
    candidate = _candidate(candidate)
    seam, pile, binary = _states_on_first_seam()
    i, j = mixing.INDEX[pile], mixing.INDEX[binary]
    difference = tuple(a-b for a, b in zip(pile, binary, strict=True))
    sector = _sector(difference)
    position = sector.pairs.index((i, j))
    block = mixing.native_block(candidate, difference)
    row = block[position]
    survival = row[position]
    if survival <= 0 or any(x for k, x in enumerate(row) if k != position):
        raise ArithmeticError("the witness coherence row is not isolated and positive")
    expected = F(23, 30) if candidate == "C" else F(7, 10)
    if survival != expected:
        raise ArithmeticError("the native witness survival coefficient changed")

    diagonal = _sector(mixing.ZERO)
    populations = mixing.native_block(candidate, mixing.ZERO)
    # Invariance of the entire binary population face, not only a one-step
    # zero at the selected input, proves that the zero persists for all n.
    for output, (a, _) in enumerate(diagonal.pairs):
        if not mixing.BINARY[a] and any(populations[output][column]
                                       for column, (b, _) in enumerate(diagonal.pairs)
                                       if mixing.BINARY[b]):
            raise ArithmeticError("binary populations can reach a collision state")
    return {"seam": seam, "pile": pile, "binary": binary,
            "coherence_survival": survival,
            "isolated_coherence_row_size": len(row),
            "binary_population_columns_checked": sum(mixing.BINARY)}


def full_carrier_npt_minor(candidate, attempts):
    """Exact negative determinant, including n=0, from the checked native row."""
    attempts = _attempts(attempts)
    survival = full_carrier_witness(candidate)["coherence_survival"]
    return -survival ** (2 * attempts) / 4


def phase_pair_signal(attempts):
    """C minus E plus-readout probability for the uniform phase-pair mixture."""
    attempts = _attempts(attempts)
    return (COHERENCE_RATES["C"] ** attempts
            - COHERENCE_RATES["E"] ** attempts) / 2


def phase_pair_witness():
    """Verify the common-spectator coherence eigenoperator from native blocks.

    A = sum_{j != u,v} |u+j><v+j|.  Prepare the equal mixture of the ten
    normalized plus states and read the binary effect (I+A+A*)/2.  The
    effect is valid because A+A* has eigenvalues -1, 0, 1.  This is an
    analytic apparatus-level witness, not a new executed source pulse word.
    """
    (u, v), _, _ = _states_on_first_seam()
    difference = tuple(int(j == u) - int(j == v) for j in range(12))
    sector = _sector(difference, binary=True)
    expected_pairs = []
    for spectator in range(12):
        if spectator in (u, v):
            continue
        left, right = [0] * 12, [0] * 12
        left[u] = left[spectator] = right[v] = right[spectator] = 1
        expected_pairs.append((mixing.INDEX[tuple(left)], mixing.INDEX[tuple(right)]))
    if set(sector.pairs) != set(expected_pairs) or len(expected_pairs) != 10:
        raise ArithmeticError("phase-pair support differs from the native coherence sector")
    if len({i for pair in expected_pairs for i in pair}) != 20:
        raise ArithmeticError("the phase-pair supports are not disjoint")
    for candidate, rate in COHERENCE_RATES.items():
        block = mixing.native_block(candidate, difference, binary=True)
        if any(sum(row, F(0)) != rate for row in block):
            raise ArithmeticError("the proposed phase operator is not a native eigenoperator")
    # Delta(n+1)-Delta(n) changes sign exactly when (25/21)**n crosses9/5.
    if not F(25, 21)**3 < F(9, 5) < F(25, 21)**4:
        raise ArithmeticError("the phase-pair optimum is not resolved exactly")
    return {"seam": [u, v], "mixture_components": len(expected_pairs),
            "coherence_rates": {key: str(value) for key, value in COHERENCE_RATES.items()},
            "preparation": "uniform mixture of (|u+j>+|v+j>)/sqrt(2), j outside the seam",
            "readout_effect": "(I+A+A*)/2; A=sum_j |u+j><v+j|",
            "probability_difference_formula": "((5/6)^n-(7/10)^n)/2",
            "maximizing_attempts": 4,
            "maximum_probability_difference": str(phase_pair_signal(4)),
            "new_source_execution": False}


def audit():
    """JSON-friendly rational certificates; theorem imports remain explicit."""
    dimension = sum(mixing.BINARY)
    diagonals = len(_sector(mixing.ZERO, binary=True).pairs)
    units = sum(len(s.pairs) for s in mixing.sectors(binary=True))
    if (dimension, diagonals - 1, units - diagonals) != (
            DIMENSION, TRACELESS_DIAGONALS, OFFDIAGONALS):
        raise ArithmeticError("binary operator-space multiplicities changed")
    lower, upper = contraction_interval()
    threshold = F(1, dimension**2 - 1)
    crude = 0
    while upper**crude > threshold:
        crude += 1
    if lower**(crude - 1) <= threshold:
        raise ArithmeticError("the enclosure cannot resolve the crude stopping bound")
    binary_reports, full_reports = {}, {}
    for candidate in COHERENCE_RATES:
        attempts = stopping_bound(candidate)
        before = choi_purity_interval(candidate, attempts - 1)[0]
        after = choi_purity_interval(candidate, attempts)[1]
        # The scaled excess has threshold one.  Conservative decimal
        # rationals keep receipts readable without omitting the inequality.
        scale = dimension**2 * (dimension**2 - 1)
        before_excess = _rounded((before - F(1, dimension**2)) * scale, upper=False)
        after_excess = _rounded((after - F(1, dimension**2)) * scale, upper=True)
        if not after_excess <= 1 < before_excess:
            raise ArithmeticError("rounded rational evidence failed the stopping inequalities")
        binary_reports[candidate] = {
            "sufficient_entanglement_breaking_attempts": attempts,
            "previous_attempt_spectral_bound_scaled_excess_lower": str(before_excess),
            "stopping_spectral_bound_scaled_excess_upper": str(after_excess),
            "scaled_excess_threshold": "1",
            "exact_entanglement_breaking_index_claimed": False,
        }
        witness = full_carrier_witness(candidate)
        rate = witness["coherence_survival"]
        full_reports[candidate] = {
            **{key: value for key, value in witness.items() if key != "coherence_survival"},
            "coherence_survival": str(rate),
            "partial_transpose_principal_minor": f"-({rate})^(2*n)/4",
            "minor_at_one_attempt": str(-rate**2 / 4),
            "negative_for_every_nonnegative_integer_n": True,
            "finite_time_entanglement_breaking": False,
        }
    return {
        "schema": "oph.observer.quantum_memory.v1",
        "scope": "independent uniform seam and fair coin; repair records forgotten",
        "binary_dimension": dimension,
        "traceless_diagonal_dimensions": diagonals - 1,
        "offdiagonal_dimensions": units - diagonals,
        "contraction_interval": [str(lower), str(upper)],
        "normalized_choi_purity_threshold": str(threshold),
        "purity_bound_formula": "[1+65*((55+sqrt(5))/60)^(2*n)+4290*r_coherence^(2*n)]/4356",
        "crude_sufficient_attempts": crude,
        "binary": binary_reports,
        "full_carrier_entangled_input": "(|2u>|0>+|u+v>|1>)/sqrt(2)",
        "partial_transpose_principal_basis": ["|2u,1>", "|u+v,0>"],
        "full_carrier": full_reports,
        "phase_pair": phase_pair_witness(),
        "separability_reference": "https://arxiv.org/pdf/quant-ph/0204159, Corollary 3",
        "retained_records_or_private_memory_covered": False,
    }
