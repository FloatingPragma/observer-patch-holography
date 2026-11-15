"""Exact checkpoint record algebra and controlled-instrument adapter.

This executes the postcommit construction in model/observer/DERIVATION.md.
Channel completeness and Kraus intertwiners imply the dual fixed-observable
criterion of Lean/Dynamics/ProtectedCharge.lean. Writing a new record is a
separate event: raw record projectors are not fixed on the blank sector.

Dense matrices use immutable tuples of Gaussian rationals (real, imaginary),
compatible with the common-source gauge/matter matrix interface. Sparse native
operators retain the actual 3354-dimensional source/tape basis. Arbitrary
operators, including off-diagonal matrix units, are accepted; these routines
never normalize a branch or infer positivity from a numerical tolerance.
"""
from __future__ import annotations

from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import factorial

from model.observer import native

ZERO = (Q(0), Q(0))
ONE = (Q(1), Q(0))


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def _multiply(a, b):
    return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])


def _conjugate(a):
    return (a[0], -a[1])


def _scalar(value):
    if (type(value) is not tuple or len(value) != 2
            or any(type(part) is not Q for part in value)):
        raise ValueError("Expected a pair of exact Fraction components")
    return value


def _matrix(value):
    if (type(value) is not tuple or not value
            or any(type(row) is not tuple or len(row) != len(value)
                   for row in value)):
        raise ValueError("Expected a nonempty immutable square matrix")
    for row in value:
        for entry in row:
            _scalar(entry)
    return value


def _identity(dimension):
    return tuple(tuple(ONE if row == column else ZERO
                       for column in range(dimension))
                 for row in range(dimension))


def _dagger(matrix):
    return tuple(tuple(_conjugate(matrix[column][row])
                       for column in range(len(matrix)))
                 for row in range(len(matrix)))


def _matmul(left, right):
    size = len(left)
    return tuple(tuple(_sum(_multiply(left[row][k], right[k][column])
                            for k in range(size))
                       for column in range(size))
                 for row in range(size))


def _sum(values):
    result = ZERO
    for value in values:
        result = _add(result, value)
    return result


def _shape(operator, record_count):
    operator = _matrix(operator)
    if (type(record_count) is not int or record_count < 1
            or len(operator) % record_count):
        raise ValueError("Positive record count must divide the matrix size")
    return operator, len(operator) // record_count


def block_readback(operator, record_count):
    """Pinch old record labels, retaining EVERY operator within each block.

    This conditional expectation fixes the generally noncommutative algebra
    direct_sum_r B(work). It is distinct from a classical public readback.
    """
    operator, dimension = _shape(operator, record_count)
    return tuple(tuple(value if row//dimension == column//dimension else ZERO
                       for column, value in enumerate(entries))
                 for row, entries in enumerate(operator))


def public_record_readback(operator, record_count):
    """Conditional expectation onto sum_r C Q_r using normalized block trace."""
    operator, dimension = _shape(operator, record_count)
    means = [_sum(operator[r*dimension+k][r*dimension+k]
                  for k in range(dimension)) for r in range(record_count)]
    means = [(a/dimension, b/dimension) for a, b in means]
    return tuple(tuple(means[row//dimension] if row == column else ZERO
                       for column in range(len(operator)))
                 for row in range(len(operator)))


def _controlled_kraus(record_kraus):
    """Build block diagonal Kraus maps; every record retains every outcome."""
    if (type(record_kraus) is not tuple or not record_kraus
            or any(type(branches) is not tuple or not branches
                   for branches in record_kraus)):
        raise ValueError("Expected record-indexed tuples of Kraus matrices")
    dimension = len(_matrix(record_kraus[0][0]))
    outcomes = len(record_kraus[0])
    for branches in record_kraus:
        if len(branches) != outcomes:
            raise ValueError("Every record must retain the complete outcome menu")
        for matrix in branches:
            if len(_matrix(matrix)) != dimension:
                raise ValueError("Every branch must act on the same work space")
        products = [_matmul(_dagger(matrix), matrix) for matrix in branches]
        gram = tuple(tuple(_sum(matrix[a][b] for matrix in products)
                           for b in range(dimension))
                     for a in range(dimension))
        if gram != _identity(dimension):
            raise ValueError("Each record's Kraus family must be exactly complete")
    size = len(record_kraus)*dimension
    return tuple(tuple(tuple(
        record_kraus[a//dimension][outcome][a % dimension][b % dimension]
        if a//dimension == b//dimension else ZERO
        for b in range(size)) for a in range(size)) for outcome in range(outcomes))


def controlled_instrument(operator, record_kraus):
    """Apply the complete old-record-controlled instrument without readback.

    record_kraus[r][y] is K_(r,y). Each returned unnormalized operator is the
    retained outcome y, including zero-probability outcomes. Old-record
    off-diagonal inputs are acted on as K_(r,y) X_rs K_(s,y)^dagger.
    """
    operator = _matrix(operator)
    kraus = _controlled_kraus(record_kraus)
    if len(operator) != len(kraus[0]):
        raise ValueError("Operator and controlled instrument dimensions differ")
    return tuple(_matmul(_matmul(matrix, operator), _dagger(matrix))
                 for matrix in kraus)


def controlled_certificate(record_kraus):
    """Exact Gram and record-projector intertwiners on the WHOLE carrier."""
    kraus = _controlled_kraus(record_kraus)
    count = len(record_kraus)
    size = len(kraus[0])
    dimension = size//count
    products = [_matmul(_dagger(matrix), matrix) for matrix in kraus]
    gram = tuple(tuple(_sum(matrix[a][b] for matrix in products)
                       for b in range(size)) for a in range(size))
    if gram != _identity(size):
        raise AssertionError("Controlled Gram identity failed")
    for matrix in kraus:
        if any(matrix[a][b] != ZERO for a in range(size) for b in range(size)
               if a//dimension != b//dimension):
            raise AssertionError("An old record is changed")
    return {"record_count": count, "work_dimension": dimension,
            "retained_outcomes": len(kraus), "gram_identity_exact": True,
            "old_record_kraus_intertwiners_exact": True,
            "dual_fixes_all_old_record_projectors": True,
            "arbitrary_operator_and_reference_scope": True,
            "whole_work_block_algebra_preservation_claimed": False}


def _state(state):
    if (type(state) is not native.State or type(state.loads) is not tuple
            or len(state.loads) != 12 or any(type(v) is not int or v < 0
                                           for v in state.loads)
            or sum(state.loads) != 2 or type(state.prefix) is not tuple
            or len(state.prefix) > 2
            or any(type(r) is not tuple or len(r) != 2
                   or any(type(v) is not int for v in r) or r not in native.LETTERS
                   for r in state.prefix)):
        raise ValueError("Expected a state of the actual native source/tape")
    return state


def _sparse(operator):
    if type(operator) is not dict:
        raise ValueError("Expected a sparse native operator dictionary")
    for pair, coefficient in operator.items():
        if type(pair) is not tuple or len(pair) != 2:
            raise ValueError("Expected two source/tape basis states")
        for state in pair:
            _state(state)
        _scalar(coefficient)
    return operator


def _accumulate(output, pair, coefficient):
    updated = _add(output.get(pair, ZERO), coefficient)
    if updated == ZERO:
        output.pop(pair, None)
    else:
        output[pair] = updated


def native_tick_operator(operator):
    """The existing coherent T=A0(.)A0*+A1(.)A1*+P2(.)P2, exactly.

    Equal-amplitude native write columns make every matrix-unit coefficient
    rational: nonfull columns have amplitude 1/sqrt(2), giving product 1/2.
    This is NOT an incoherent coin or current measurement. Different phase
    sectors are separate Kraus branches, as in the original construction.
    """
    output = {}
    for (left, right), coefficient in _sparse(operator).items():
        if left.phase != right.phase:
            continue
        if left.phase == 2:
            _accumulate(output, (left, right), coefficient)
            continue
        for weight_left, after_left in native.tick(left):
            for weight_right, after_right in native.tick(right):
                if weight_left != Q(1, 2) or weight_right != Q(1, 2):
                    raise AssertionError("Native write amplitudes changed")
                _accumulate(output, (after_left, after_right),
                            (coefficient[0]/2, coefficient[1]/2))
    return output


def native_block_readback(operator):
    """Pinch the first committed letter, with blank as a separate block."""
    return {pair: coefficient for pair, coefficient in _sparse(operator).items()
            if pair[0].prefix[:1] == pair[1].prefix[:1]}


def native_first_record_matrix(operator):
    """Trace work/phase/second-slot data on the already committed carrier."""
    result = [[ZERO for _ in native.LETTERS] for _ in native.LETTERS]
    index = {record: i for i, record in enumerate(native.LETTERS)}
    for (left, right), coefficient in _sparse(operator).items():
        if not left.prefix or not right.prefix:
            raise ValueError("First record has not been committed")
        if (left.phase == right.phase and left.loads == right.loads
                and left.prefix[1:] == right.prefix[1:]):
            a, b = index[left.prefix[0]], index[right.prefix[0]]
            result[a][b] = _add(result[a][b], coefficient)
    return tuple(tuple(row) for row in result)


def native_record_weight(operator, record):
    """Linear expectation of a first-letter projector, zero on blank tape."""
    if record not in native.LETTERS:
        raise ValueError("Unknown native record")
    return _sum(coefficient for (left, right), coefficient in _sparse(operator).items()
                if left == right and left.prefix[:1] == (record,))


@lru_cache(maxsize=None)
def checkpoint_basis():
    """The actual postcommit carrier, including unreachable tape/load pairs."""
    return tuple(native.State(loads, prefix) for length in (1, 2)
                 for prefix in product(native.LETTERS, repeat=length)
                 for loads in native.states())


def _checkpoint_operator(operator):
    operator = _sparse(operator)
    if any(not state.prefix for pair in operator for state in pair):
        raise ValueError("The persistence observer excludes the blank sector")
    return operator


def checkpoint_record_element(record):
    """A committed first-letter projector on the complete checkpoint space."""
    if record not in native.LETTERS:
        raise ValueError("Unknown native record")
    return {(state, state): ONE for state in checkpoint_basis()
            if state.prefix[0] == record}


def checkpoint_public_readback(operator):
    """Normalized block trace on the actual six equal 546-dimensional blocks."""
    operator = _checkpoint_operator(operator)
    dimension = len(checkpoint_basis()) // len(native.LETTERS)
    means = {record: tuple(part/dimension for part in native_record_weight(operator, record))
             for record in native.LETTERS}
    return {(state, state): means[state.prefix[0]] for state in checkpoint_basis()
            if means[state.prefix[0]] != ZERO}


@lru_cache(maxsize=None)
def _checkpoint_write_preimages(seam):
    # None is the actual record-dependent native policy. E/F are its two
    # already available interventions, selected from a retained letter.
    if seam not in (None, native.E, native.F):
        raise ValueError("Only the native policy and its E/F controls are admitted")
    inverse = {}
    for before in checkpoint_basis():
        if before.phase != 1:
            continue
        if seam is None:
            outputs = [state for _, state in native.tick(before)]
        else:
            outputs = [native.State(loads, before.prefix+(record,))
                       for coin in (0, 1)
                       for loads, record in (native.native(before.loads, seam, coin),)]
        for after in outputs:
            if after in inverse:
                raise AssertionError("Native write columns lost orthogonality")
            inverse[after] = before
    return inverse


def checkpoint_dual(operator, seam=None):
    """Exact Heisenberg continuation A1* X A1 + P2 X P2 on K.

    The same source columns define the persistence generator (seam=None)
    and the E/F interventions; no unitary or autonomous-time substitution.
    """
    inverse, output = _checkpoint_write_preimages(seam), {}
    for (left, right), coefficient in _checkpoint_operator(operator).items():
        if left.phase == right.phase == 2:
            _accumulate(output, (left, right), coefficient)
            if left in inverse and right in inverse:
                _accumulate(output, (inverse[left], inverse[right]),
                            tuple(part/2 for part in coefficient))
    return output


def checkpoint_generator(operator):
    """T_K*−I, an algebraic persistence generator in continuation-step units."""
    output = checkpoint_dual(operator)
    for pair, coefficient in _checkpoint_operator(operator).items():
        _accumulate(output, pair, tuple(-part for part in coefficient))
    return output


def _evidence_label(label):
    if (type(label) is not tuple or len(label) != 2
            or type(label[0]) is not int or label[0] not in (0, 1)
            or type(label[1]) is not tuple or len(label[1]) != 2
            or any(type(value) is not int for value in label[1])
            or label[1] not in native.LETTERS):
        raise ValueError("An issued/retained native evidence label is required")
    return label


def checkpoint_precedes(first, second):
    """Order two reads of the same evidence across one native continuation."""
    a, r = _evidence_label(first)
    b, s = _evidence_label(second)
    return r == s and a < b


def checkpoint_predict(label):
    """Advance issuance to retention of the SAME evidence, not predict a new bit."""
    _, record = _evidence_label(label)
    return (1, record)


def checkpoint_control(operator, label):
    """The retained current selects an actual native E/F continuation map."""
    _, record = _evidence_label(label)
    return checkpoint_dual(operator, native.F if record[1] > 0 else native.E)


def checkpoint_preparation():
    """Actual first-event state with its classical letter committed, no refit."""
    prior = {(native.State(loads), native.State(loads)): (Q(1, 66), Q(0))
             for loads in native.states(binary=True)}
    return native_block_readback(native_tick_operator(prior))


def checkpoint_buffer_lift(operator, bits, *, density=False):
    """Declared exact buffers: X tensor I or rho tensor |0><0|.

    Basis keys are (native checkpoint, buffer address). These are architecture
    refinements, not occupation cutoffs or spatial refinements.
    """
    if type(bits) is not int or not 0 <= bits <= 3 or type(density) is not bool:
        raise ValueError("The executable buffer witness admits zero to three bits")
    return {((left, address), (right, address)): value
            for (left, right), value in _checkpoint_operator(operator).items()
            for address in (range(1) if density else range(2**bits))}


def checkpoint_observer_certificate():
    """All fields of the finite postcommit OperationalObserver, exactly.

    DERIVATION.md supplies the all-operator proofs and arbitrary buffer
    refinement. This finite certificate checks the actual source columns,
    state, every evidence/control label, and a nonzero generator witness.
    It is a mathematical interface realization, not a new Lean declaration.
    """
    labels = tuple(product((0, 1), native.LETTERS))
    state = checkpoint_preparation()
    if any(left != right or value[1] or value[0] < 0
           for (left, right), value in state.items()):
        raise AssertionError("The committed native preparation is not a diagonal state")
    if _sum(state.values()) != ONE:
        raise AssertionError("The native checkpoint preparation is not normalized")
    weights = {record: native_record_weight(state, record)[0]
               for record in native.LETTERS}
    for label in labels:
        stage, record = label
        element = checkpoint_record_element(record)
        if checkpoint_public_readback(element) != element or checkpoint_generator(element):
            raise AssertionError("Native self-readback or durable evidence failed")
        if any(checkpoint_control(element, control) != element for control in labels):
            raise AssertionError("A native control changed committed evidence")
        prediction = checkpoint_predict(label)
        if prediction[0] < stage or prediction[1] != record:
            raise AssertionError("Evidence retention prediction regressed or changed value")
    witness = native.State(native.states()[0], ((0, 0),))
    unit = {(witness, witness): ONE}
    if checkpoint_generator(unit) != {(witness, witness): (Q(-1), Q(0))}:
        raise AssertionError("The persistence generator witness changed")
    # Nonconstant record-conditioned controls on the actual K, rather than
    # copies of one abstract identity intervention.
    before = native.State((1, 1)+(0,)*10, ((0, 0),))
    loads, letter = native.native(before.loads, native.E, 0)
    after = native.State(loads, before.prefix+(letter,))
    effect = {(after, after): ONE}
    if checkpoint_control(effect, (0, (0, 0))) == checkpoint_control(effect, (0, (0, 1))):
        raise AssertionError("The two native evidence-conditioned controls coincided")
    # Actual conditional forecasts are an additional function of the record.
    # Equal record PROBABILITIES do not imply equal forecasts.
    a, b = (0, 0), (1, 0)
    if weights[a] != Q(14, 33) or weights[a] != weights[b] or native.FORECAST[a] == native.FORECAST[b]:
        raise AssertionError("The event-probability/forecast congruence control changed")
    return {
        "private_dimension": len(checkpoint_basis()), "record_block_dimension": 546,
        "evidence_labels": len(labels), "accessible_algebra": "span of first-letter projectors",
        "state": "actual first native event from uniform binary prior, first letter committed",
        "generator": "T_K* - I", "generator_nonzero": True,
        "control": "current-selected native E/F dual continuation with all outcomes",
        "record_conditioned_controls_are_nonconstant": True,
        "prediction": "issuance to retention of the same committed evidence",
        "all_OperationalObserver_clauses_realized_on_declared_postcommit_tower": True,
        "Lean_inhabitant_compiled": False,
        "buffer_refinement_architecture_supplied": True,
        "buffer_refinement": "X -> X tensor I; rho -> rho tensor |0><0|",
        "record_creation_in_persistence_domain": False,
        "native_repair_global_termination_or_physical_clock_inferred": False,
        "forecast_is_not_the_typed_evidence_retention_prediction": True,
        "forecast_congruence_counterexample": {
            "first_records": [list(a), list(b)], "equal_readout_weight": str(weights[a]),
            "different_conditional_forecasts": [str(native.FORECAST[a]), str(native.FORECAST[b])]},
        "first_record_probabilities": [str(weights[record]) for record in native.LETTERS]}


def product_access_generators(coupling=Q(1)):
    """Actual product source terms and the second-mode coordinate, exactly.

    Rational Bargmann matrices use the SAME least-significant first-mode
    ordering as ProductBand. Conjugation by diag(sqrt(n0! n1!)) gives the
    orthonormal work-code matrices; algebra dimensions are unchanged.
    """
    from model.interacting_band_extension.construction import ProductBand
    from model.interacting_band_extension import algebra
    lattice = ProductBand(coupling=coupling)
    dimension = lattice.dimension

    def embed(sites, factors, coefficient=Q(1)):
        local = dict(zip(sites, factors, strict=True))
        return tuple(tuple(coefficient *
            product_entry(row, column, local) for column in range(dimension))
            for row in range(dimension))

    def product_entry(row, column, factors):
        value = Q(1)
        for mode in range(2):
            left, right = (row//4**mode) % 4, (column//4**mode) % 4
            value *= factors[mode][left][right] if mode in factors else int(left == right)
        return value

    coordinate = algebra.polynomial(4, 1)
    return {"coordinate": embed((1,), (coordinate,)),
            "excluded_first_coordinate": embed((0,), (coordinate,)),
            "terms": tuple(embed(term.sites, term.factors, term.coefficient)
                           for term in lattice.terms),
            "term_names": tuple(term.label for term in lattice.terms),
            "first_mode_parity": tuple((-1)**(index % 4) for index in range(dimension)),
            "bargmann_metric": tuple(factorial(index % 4)*factorial(index//4)
                                     for index in range(dimension))}


def product_access_rank(coupling=Q(1), *, first_mode_pointer=False):
    """Certified rational lower rank of the readout/control closure.

    Work modulo the prime 65537. Every source coefficient denominator is
    invertible there. A nonzero reduced minor is a nonzero rational minor;
    a separate exact parity bound supplies equality in the declared case.
    No numerical-rank tolerance or unobserved Hamiltonian read is used.
    """
    import numpy as np
    if type(first_mode_pointer) is not bool:
        raise ValueError("The extra first-mode pointer selector must be boolean")
    data = product_access_generators(coupling)
    modulus, dimension = 65537, 16

    def reduced(matrix):
        return np.array([[value.numerator*pow(value.denominator, -1, modulus) % modulus
                          for value in row] for row in matrix], dtype=np.int64)

    coordinates = [reduced(data["coordinate"])]
    if first_mode_pointer:
        coordinates.append(reduced(data["excluded_first_coordinate"]))
    terms = [reduced(term) for term in data["terms"]]
    basis, queue = {}, []

    def append(matrix):
        vector = matrix.reshape(-1) % modulus
        for pivot, row in basis.items():
            if vector[pivot]:
                vector = (vector-vector[pivot]*row) % modulus
        nonzero = np.flatnonzero(vector)
        if len(nonzero):
            pivot = int(nonzero[0])
            vector = vector*pow(int(vector[pivot]), -1, modulus) % modulus
            basis[pivot] = vector
            queue.append(vector.reshape(dimension, dimension))

    append(np.eye(dimension, dtype=np.int64))
    for matrix in queue:
        for coordinate in coordinates:
            append(coordinate@matrix % modulus)
        for term in terms:
            append((term@matrix-matrix@term) % modulus)
    # Matrix products sum at most sixteen products of residues, far below
    # int64 overflow. The closure uses only multiplication and derivations.
    return len(basis)


def product_access_certificate():
    """Readout-derived noncommutative access for the declared product menu."""
    data = product_access_generators()
    parity, metric = data["first_mode_parity"], data["bargmann_metric"]
    for matrix in (data["coordinate"], *data["terms"]):
        for row, entries in enumerate(matrix):
            for column, value in enumerate(entries):
                if metric[row]*value != metric[column]*matrix[column][row]:
                    raise AssertionError("A product generator lost its physical adjoint")
                if value and parity[row] != parity[column]:
                    raise AssertionError("An admitted product service changes first-mode parity")
    ranks = (product_access_rank(), product_access_rank(Q(0)),
             product_access_rank(first_mode_pointer=True))
    if ranks != (128, 16, 256):
        raise AssertionError("The product readout/control algebra or its controls changed")
    return {"work_dimension": 16, "measured_mode": 1, "unmeasured_parity_mode": 0,
            "source_term_names": list(data["term_names"]),
            "menu": "second-mode pointer and finite words of the three local source services with selectable durations",
            "local_source_service_menu_supplied": True,
            "fixed_two_read_program_algebra_identified": False,
            "rational_closure_lower_rank": ranks[0], "rank_certificate_prime": 65537,
            "exact_first_mode_parity_upper_dimension": 128,
            "accessible_work_algebra": "B(C8) direct_sum B(C8)",
            "disconnected_coupling_access_dimension": ranks[1],
            "extra_first_mode_pointer_access_dimension": ranks[2],
            "full_B_C16_access_from_declared_menu": False,
            "prewrite_effect_is_postwrite_record_pullback": True,
            "old_records_and_arbitrary_references_preserved": True,
            "new_record_durable_before_its_write": False,
            "physical_source_or_continuum_algebra_identification": False}


PRODUCT_RECORD_BITS = ((), (0,), (1,), (0, 0), (0, 1), (1, 0), (1, 1))
PRODUCT_CENTRAL_EVENT = ((0,), (0, 0), (1, 1))


def _rational_product(left, right):
    """Gaussian-rational multiplication with structural zero entries skipped."""
    size = len(left)
    output = [[ZERO for _ in range(size)] for _ in range(size)]
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if value != ZERO:
                for j, other in enumerate(right[k]):
                    if other != ZERO:
                        output[i][j] = _add(output[i][j], _multiply(value, other))
    return tuple(tuple(row) for row in output)


@lru_cache(maxsize=4, typed=True)
def product_central_record_bounds(interval=Q(1, 16), degree=32):
    """Exact enclosures of the existing complete product record law.

    H, beta=1/4, the second-mode pointer and both +Q/4 kicks are fixed.
    The interval is a declared model duration, not an inferred clock unit.
    Taylor remainders are operator bounds in the physical weighted metric.
    No floating eigensystem or numerical rank certifies these intervals.
    """
    if type(interval) is not Q or not 0 < interval <= Q(1, 8):
        raise ValueError("Expected an exact model interval in (0, 1/8]")
    if type(degree) is not int or not 24 <= degree <= 48:
        raise ValueError("Expected a Taylor degree between 24 and 48")
    data = product_access_generators()
    h = tuple(tuple(sum(term[i][j] for term in data['terms'])
                    for j in range(16)) for i in range(16))
    q, metric = data['coordinate'], data['bargmann_metric']
    h_bound = max(sum(abs(x) for x in row) for row in h)
    q_bound = max(sum(abs(x) for x in row) for row in q)

    def taylor(matrix, scalar):
        base = tuple(tuple(_multiply((x, Q(0)), scalar) for x in row)
                     for row in matrix)
        total = term = _identity(16)
        for k in range(1, degree+1):
            term = tuple(tuple((a/k, b/k) for a, b in row)
                         for row in _rational_product(base, term))
            total = tuple(tuple(_add(x, y) for x, y in zip(a, b, strict=True))
                          for a, b in zip(total, term, strict=True))
        return total

    gibbs = taylor(h, (-Q(1, 4), Q(0)))
    evolve = taylor(h, (Q(0), -interval))
    kick = taylor(q, (Q(0), Q(1, 4)))
    phase = taylor(q, (Q(0), Q(1, 2)))
    pointers = tuple(tuple(tuple((a+sign*b, Q(0)) for a, b in row)
                           for row in phase) for sign in (-1, 1))
    factorial_tail = factorial(degree+1)
    eg = (h_bound/4)**(degree+1)/factorial_tail
    eu = (interval*h_bound)**(degree+1)/factorial_tail
    ef = (q_bound/4)**(degree+1)/factorial_tail
    ek = 2*(q_bound/2)**(degree+1)/factorial_tail

    def effect(word, normalization):
        adjoint = tuple(tuple((a*Q(metric[j], metric[i]), -b*Q(metric[j], metric[i]))
            for j, (a, b) in enumerate(column)) for i, column in enumerate(zip(*word)))
        return tuple(tuple((a/normalization, b/normalization) for a, b in row)
                     for row in _rational_product(adjoint, word))

    first = tuple(_rational_product(pointer, _rational_product(evolve, kick))
                  for pointer in pointers)
    words = {(bit,): word for bit, word in enumerate(first)}
    for bit, word in enumerate(first):
        continuation = _rational_product(kick, word) if bit else word
        advanced = _rational_product(evolve, continuation)
        for last, pointer in enumerate(pointers):
            words[bit, last] = _rational_product(pointer, advanced)
    effects = {bits: effect(word, 2**len(bits)) for bits, word in words.items()}
    prefix = {'even': {(): (Q(25, 36), Q(25, 36))},
              'odd': {(): (Q(25, 36), Q(25, 36))}}
    partitions = {}
    scale = 10**15

    def outward(center, error):
        return (max(Q(0), Q(((center-error)*scale).__floor__(), scale)),
                min(Q(1), Q(((center+error)*scale).__ceil__(), scale)))

    for sector, sign in (('even', 1), ('odd', -1)):
        indexes = tuple(i for i, parity in enumerate(data['first_mode_parity']) if parity == sign)
        z = _sum(gibbs[i][i] for i in indexes)
        if z[1] or z[0] <= 8*eg:
            raise ArithmeticError("The exact sector partition lower bound is not positive")
        partitions[sector] = (z[0]-8*eg, z[0]+8*eg)
        for bits, observable in effects.items():
            n = len(bits)
            kicks = 1+int(n == 2 and bits[0] == 1)
            eta = (1+ek)**n*(1+eu)**n*(1+ef)**kicks-1
            effect_error = eta*(2+eta)
            numerator = _sum(_multiply(gibbs[i][j], observable[j][i])
                             for i in indexes for j in range(16))
            if numerator[1]:
                raise ArithmeticError("A weighted self-adjoint effect has a complex expectation")
            error = (8*(eg*(1+effect_error)+effect_error)+8*eg)/(z[0]-8*eg)
            weight = Q(5, 18) if n == 1 else Q(1, 36)
            prefix[sector][bits] = outward(weight*numerator[0]/z[0], weight*error)
    event = {sector: tuple(sum(prefix[sector][bits][i] for bits in PRODUCT_CENTRAL_EVENT)
                           for i in (0, 1)) for sector in prefix}
    if event['even'][0] <= event['odd'][1]:
        raise ArithmeticError("The exact event intervals do not separate the parity sectors")
    return {'sector_prefix_probability_intervals': prefix,
            'sector_partition_intervals': partitions,
            'central_event_probability_intervals': event,
            'central_event_gap_lower': event['even'][0]-event['odd'][1],
            'hamiltonian_norm_upper': h_bound, 'coordinate_norm_upper': q_bound,
            'model_interval': interval, 'inverse_temperature': Q(1, 4), 'Taylor_degree': degree}


def product_central_source_error(interval=Q(1, 16), steps=64):
    """Conditional service/splitting allowance; excludes unknown state error."""
    from model.interacting_band_extension.construction import ProductBand, splitting_bound
    if type(interval) is not Q or not 0 < interval <= Q(1, 8):
        raise ValueError("Expected an exact model interval in (0, 1/8]")
    if type(steps) is not int or steps < 1:
        raise ValueError("Expected a positive integer splitting count")
    lattice = ProductBand()
    split = splitting_bound(lattice, interval, steps)
    gain = (4*steps*len(lattice.terms)+4)*Q(251, 10**12)
    if gain >= 1:
        raise ValueError("The declared source error accumulation is not contractive")
    vector = gain/(1-gain)
    # One incident has half the two-interval bound; all zero-incident outputs
    # are outside the event. The complete two-read event needs only one bound.
    bound = (Q(5, 18)*min(Q(1), split['one_interval_operator_error_upper'])
            + Q(1, 36)*split['two_interval_complete_output_trace_distance_upper']
            + vector*(2+vector)/2)
    return Q((bound*10**15).__ceil__(), 10**15)


def product_central_weight_interval(counts, statistical_radius, probability_error):
    """Certified inverse interval from ALL seven record-pattern counts.

    For independent identically prepared trials, coverage is at least
    1-2 exp(-2 N radius^2), conditional on the supplied probability-error
    allowance. An empty interval rejects the declared calibration contract.
    """
    if (type(counts) is not dict or set(counts) != set(PRODUCT_RECORD_BITS)
            or any(type(v) is not int or v < 0 for v in counts.values()) or not sum(counts.values())):
        raise ValueError("Every complete record pattern needs a nonnegative integer count")
    if any(type(x) is not Q or not 0 <= x <= 1 for x in (statistical_radius, probability_error)):
        raise ValueError("Expected exact probability-error radii in [0,1]")
    event = product_central_record_bounds()['central_event_probability_intervals']
    odd, even = event['odd'], event['even']
    count = sum(counts[bits] for bits in PRODUCT_CENTRAL_EVENT)
    frequency = Q(count, sum(counts.values()))
    low = max(Q(0), (frequency-statistical_radius-probability_error-odd[1])/(even[1]-odd[1]))
    high = min(Q(1), (frequency+statistical_radius+probability_error-odd[0])/(even[0]-odd[0]))
    return None if low > high else (low, high)


def product_central_sample_budget(weight_radius, probability_error):
    """Sufficient 95% IID trial budget, or None when the error floor dominates."""
    if type(weight_radius) is not Q or not 0 < weight_radius < 1:
        raise ValueError("Expected an exact weight radius in (0,1)")
    if type(probability_error) is not Q or not 0 <= probability_error <= 1:
        raise ValueError("Expected an exact probability-error allowance")
    event = product_central_record_bounds()['central_event_probability_intervals']
    centers = {s: sum(interval)/2 for s, interval in event.items()}
    calibration = max((b-a)/2 for a, b in event.values())
    margin = weight_radius*(centers['even']-centers['odd'])-probability_error-calibration
    # N >= 2/r^2 gives 2 exp(-2 N r^2) <= 2 exp(-4) < 1/20.
    return None if margin <= 0 else (2/margin**2).__ceil__()


def product_central_weight_certificate():
    """Known-family identification from one actual event, not state selection."""
    bounds = product_central_record_bounds()
    source_error = product_central_source_error()
    return {'fixed_coupling': '1', 'fixed_inverse_temperature': '1/4',
            'fixed_model_interval': '1/16', 'split_steps': 64,
            'known_family': 'w rho_beta_even direct_sum (1-w) rho_beta_odd',
            'event_bits': [list(bits) for bits in PRODUCT_CENTRAL_EVENT],
            'all_native_tapes_and_timeouts_retained': True,
            'sector_event_probability_intervals': {
                sector: [str(x) for x in interval] for sector, interval in
                bounds['central_event_probability_intervals'].items()},
            'exact_positive_event_gap_lower': str(bounds['central_event_gap_lower']),
            'exact_Taylor_degree': bounds['Taylor_degree'],
            'declared_source_probability_error_upper': str(source_error),
            'unknown_preparation_error_included': False,
            'zero_incident_timeout_probability_for_every_weight': '25/36',
            'fewer_than_two_incidents_timeout_probability_for_every_weight': '35/36',
            'sufficient_IID_trials_95_percent_weight_radius_one_percent_ideal':
                product_central_sample_budget(Q(1, 100), Q(0)),
            'sufficient_IID_trials_95_percent_weight_radius_one_percent_with_source':
                product_central_sample_budget(Q(1, 100), source_error),
            'sufficient_IID_trials_95_percent_weight_radius_ten_percent_with_source':
                product_central_sample_budget(Q(1, 10), source_error),
            'central_weight_selected_by_KMS_or_repair': False,
            'state_tomography_or_joint_parameter_identification': False,
            'finite_samples_or_laboratory_calibration_executed': False}


def native_checkpoint_certificate():
    """Check all write columns and their checkpoint intertwiners exactly.

    Orthogonal column supports and norms prove the full Gram identity, not
    just classical conservation. The phase Kraus split and retained first
    letter then imply the dual projector identity on every coherent input.
    """
    counts = {}
    for phase in (0, 1, 2):
        seen = set()
        columns = 0
        for prefix in product(native.LETTERS, repeat=phase):
            for loads in native.states():
                before = native.State(loads, prefix)
                images = native.tick(before)
                if sum(weight for weight, _ in images) != 1:
                    raise AssertionError("Native column norm differs from one")
                support = {after for _, after in images}
                if len(support) != len(images) or seen.intersection(support):
                    raise AssertionError("Native columns are not orthogonal")
                seen.update(support)
                if phase and any(after.prefix[:1] != prefix[:1]
                                 for _, after in images):
                    raise AssertionError("A committed first record changed")
                columns += 1
        counts[str(phase)] = columns
    blank = native.State((1, 1) + (0,)*10)
    operator = {(blank, blank): ONE}
    record = (0, 0)
    before = native_record_weight(operator, record)
    after = native_record_weight(native_tick_operator(operator), record)
    if before != ZERO or after != (Q(1, 2), Q(0)):
        raise AssertionError("Record-creation boundary changed")
    return {"source_tape_dimension": sum(counts.values()),
            "all_phase_column_counts": counts,
            "postcommit_dimension": counts["1"] + counts["2"],
            "all_phase_gram_identities_exact": True,
            "postcommit_first_record_intertwiners_exact": True,
            "postcommit_dual_projectors_fixed": True,
            "arbitrary_operator_and_reference_scope": True,
            "global_blank_projector_counterexample": {
                "first_record": list(record), "before": "0", "after": "1/2"},
            "full_operational_observer_inhabitant_claimed": False,
            "postcommit_operational_observer": checkpoint_observer_certificate()}
