"""A local realization and transport law for the existing exchange charge.

The additional A/B labels are coherent apparatus, not new bare native species.
On the one-A/one-B sector the two-site, two-flavor source stencil is exactly
the existing ordered cover.  Tensoring that *local* stencil with identities
is essential: a global CAR hop would introduce spectator signs.  A localized
perturbation of the tracial reference carries a conserved charge expectation;
it is not a localized sharp charge sector or a physical particle assignment.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
import json
import math

import numpy as np

from ..defect_field_bridge import full_domain
from ..defect_protection import protection
from ..exchange_sector import statistics as exchange
from ..observer.source_discrimination import circuit


PORTS = 12
PAIRS = exchange.PAIRS
SEAMS = exchange.SEAMS
CORE = (0, 1)
PAIR_DIMENSION = len(PAIRS)
DIMENSION = 2 * PAIR_DIMENSION


def _integer(value, low, high, name):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"{name} must be an exact integer in [{low}, {high}]")
    return value


def _law(value):
    if type(value) is not str or value not in ("C", "E"):
        raise ValueError("law must be C or E")
    return value


def _region(region):
    if (type(region) is not tuple or any(type(x) is not int for x in region)
            or tuple(sorted(set(region))) != region
            or any(not 0 <= x < PORTS for x in region)):
        raise ValueError("region must be a sorted tuple of distinct native ports")
    return region


def _seam(seam):
    _region(seam)
    if seam not in SEAMS:
        raise ValueError("seam must belong to the unchanged native menu")
    return seam


def _rational(value, name, low=F(0), high=F(1)):
    if type(value) is F:
        # Fraction normally enforces these invariants, but its private slots
        # can be altered or bypassed by deserializers/object construction.
        try:
            numerator, denominator = value.numerator, value.denominator
        except AttributeError as error:
            raise ValueError(f"{name} has malformed Fraction storage") from error
        if (type(numerator) is not int or type(denominator) is not int
                or denominator <= 0 or math.gcd(numerator, denominator) != 1):
            raise ValueError(f"{name} must have canonical Fraction storage")
    if type(value) not in (int, F) or not low <= value <= high:
        raise ValueError(f"{name} must be an exact rational in [{low}, {high}]")
    return F(value)


def _coefficients(values):
    if type(values) is not tuple or len(values) != PAIR_DIMENSION:
        raise ValueError("one exact nonnegative coefficient per unordered pair required")
    return tuple(_rational(x, "coefficient") for x in values)


def cover_configuration(index):
    """Actual local qutrit tensor coordinates, with 1=A and 2=B."""
    index = _integer(index, 0, DIMENSION - 1, "cover index")
    orientation, pair = divmod(index, PAIR_DIMENSION)
    a, b = PAIRS[pair]
    result = [0] * PORTS
    result[a], result[b] = ((1, 2) if not orientation else (2, 1))
    return tuple(result)


def cover_index(configuration):
    if (type(configuration) is not tuple or len(configuration) != PORTS
            or any(type(x) is not int or x not in (0, 1, 2) for x in configuration)
            or configuration.count(1) != 1 or configuration.count(2) != 1):
        raise ValueError("twelve qutrits containing exactly one A and one B required")
    pair = tuple(i for i, value in enumerate(configuration) if value)
    orientation = int(configuration[pair[0]] == 2)
    return orientation * PAIR_DIMENSION + exchange.INDEX[pair]


def local_branch_action(mask, law, coin):
    """Unchanged full sixteen-column stencil, including dual-rail collisions."""
    _integer(mask, 0, 15, "local dual-rail mask")
    _law(law)
    _integer(coin, 0, 1, "native coin")
    result = full_domain.branch_action(mask, law, coin, sites=2, flavors=2, seam=(0, 1))
    if (type(result) is not tuple or len(result) != 4
            or any(type(x) is not int for x in result)):
        raise ArithmeticError("local source returned malformed branch metadata")
    raw, output, sign, current = result
    if not (0 <= raw < 16 and 0 <= output < 16 and sign == 1 and current in (-1, 0, 1)):
        raise ArithmeticError("local source stencil changed its sign or output domain")
    if current != (mask & 3).bit_count() - (output & 3).bit_count():
        raise ArithmeticError("local source current disagrees with its actual endpoints")
    return result


def local_completion(law):
    """All 32 coin/raw slots, retaining zero slots and the collision outcome."""
    law = _law(law)
    rows = {(coin, raw): [] for coin in (0, 1) for raw in range(16)}
    for coin in (0, 1):
        for mask in range(16):
            raw, output, sign, _ = local_branch_action(mask, law, coin)
            rows[coin, raw].append((mask, output, sign))
    return {label: tuple(row) for label, row in rows.items()}


def cover_action(index, law, seam, coin):
    """Tensor the local stencil with literal spectator identities."""
    configuration = list(cover_configuration(index))
    law, (a, b) = _law(law), _seam(seam)
    _integer(coin, 0, 1, "native coin")
    raw, output, sign, current = local_branch_action(
        configuration[a] + 4 * configuration[b], law, coin)
    configuration[a], configuration[b] = output & 3, output >> 2
    return raw, cover_index(tuple(configuration)), sign, current


def ordered_cover_branches(law):
    """Actual local maps for the original complete native raw menu."""
    law = _law(law)
    width = 1 << circuit.RECORD_BITS[law]
    rows = {(seam, coin, raw): [] for seam in SEAMS
            for coin in (0, 1) for raw in range(width)}
    for seam in SEAMS:
        for coin in (0, 1):
            for index in range(DIMENSION):
                raw, output, sign, _ = cover_action(index, law, seam, coin)
                if raw >= width or sign != 1:
                    raise ArithmeticError("local tensor stencil does not realize the ordered cover")
                rows[seam, coin, raw].append((index, output))
    return {label: tuple(row) for label, row in rows.items()}


def intertwiner_certificate():
    """Compare fresh generated local maps with every predecessor raw branch."""
    columns, slots, zeros, raw_payload, local_rows = 0, 0, 0, [], []
    for law in ("C", "E"):
        completion = local_completion(law)
        incoming, outgoing = [0] * 16, [0] * 16
        collisions = 0
        for (coin, raw), row in completion.items():
            if len({target for _, target, _ in row}) != len(row):
                raise ArithmeticError("local raw branch is not a partial injection")
            lookup = {origin: target for origin, target, _ in row}
            flip = lambda x: ((x & 5) << 1) | ((x & 10) >> 1)
            for origin, target, sign in row:
                incoming[origin] += 1
                outgoing[target] += 1
                if origin.bit_count() != target.bit_count() or sign != 1:
                    raise ArithmeticError("local completion changed number or phase")
                if lookup.get(flip(origin)) != flip(target):
                    raise ArithmeticError("local completion resolves the transported labels")
                collision = (origin & 3) == 3 or (origin >> 2) == 3
                if collision != (raw == full_domain.COLLISION_RAW) or (collision and origin != target):
                    raise ArithmeticError("local collision complement was changed")
                collisions += collision
        if incoming != [2] * 16 or outgoing != [2] * 16:
            raise ArithmeticError("whole sixteen-column source completion is not bistochastic")
        local_rows.append({"law": law, "raw_slots": len(completion),
                           "zero_slots": sum(not row for row in completion.values()),
                           "input_columns_per_coin": len(incoming),
                           "collision_columns": collisions})

        actual = ordered_cover_branches(law)
        predecessor = protection.branches(law)
        for row in predecessor:
            if row.law != law:
                raise ArithmeticError("predecessor branch law metadata changed")
            _seam(row.seam)
            _integer(row.coin, 0, 1, "predecessor coin")
            _integer(row.raw, 0, (1 << circuit.RECORD_BITS[law]) - 1, "predecessor raw")
        expected_labels = set(actual)
        labels = [(row.seam, row.coin, row.raw) for row in predecessor]
        if len(labels) != len(expected_labels) or set(labels) != expected_labels:
            raise ArithmeticError("predecessor complete raw branch menu is missing or duplicated")
        for row in predecessor:
            label = (row.seam, row.coin, row.raw)
            pairs = actual[label]
            if set(pairs) != set(protection.lifted_pairs(row)):
                raise ArithmeticError("local tensor and ordered-cover branch intertwiner failed")
            lookup = dict(pairs)
            if len(lookup) != len(pairs) or len(set(lookup.values())) != len(pairs):
                raise ArithmeticError("ordered-cover branch lost partial injectivity")
            flip = lambda i: (i + PAIR_DIMENSION) % DIMENSION
            if any(lookup.get(flip(i)) != flip(j) for i, j in pairs):
                raise ArithmeticError("complete native branch fails Q commutation")
            slots += 1
            columns += len(pairs)
            zeros += not pairs
            raw_payload.append([law, row.seam, row.coin, row.raw, pairs])
    return {"native_raw_slots": slots, "supported_cover_columns": columns,
            "input_columns_checked_per_slot": DIMENSION, "zero_native_slots": zeros,
            "local_completions": local_rows,
            "local_maps_sha256": hashlib.sha256(json.dumps(
                raw_payload, separators=(",", ":")).encode()).hexdigest()}


def initial_coefficients(core=CORE):
    core = set(_region(core))
    return tuple(F(set(pair) <= core) for pair in PAIRS)


def branch_transport(coefficients, law, seam, coin, raw):
    """Unnormalized contrast after a complete retained flag; no conditioning."""
    coefficients = _coefficients(coefficients)
    law, seam = _law(law), _seam(seam)
    _integer(coin, 0, 1, "native coin")
    _integer(raw, 0, (1 << circuit.RECORD_BITS[law]) - 1, "native raw slot")
    output = [F(0)] * PAIR_DIMENSION
    for i, value in enumerate(coefficients):
        record, target, sign, _ = cover_action(i, law, seam, coin)
        if sign != 1:
            raise ArithmeticError("contrast channel acquired an unrequested phase")
        if record == raw:
            output[target % PAIR_DIMENSION] += value / 2
    return tuple(output)


def _history(history):
    if type(history) is not tuple:
        raise ValueError("a tuple of complete native history labels required")
    for event in history:
        if type(event) is not tuple or len(event) != 4:
            raise ValueError("history event must be (law, seam, coin, raw)")
        law, seam, coin, raw = event
        _law(law)
        _seam(seam)
        _integer(coin, 0, 1, "native coin")
        _integer(raw, 0, (1 << circuit.RECORD_BITS[law]) - 1, "native raw slot")
    return history


def history_transport(coefficients, history):
    """Finite branch composition, including valid zero and abort maps.

    An abort has zero weight on this exact code. Its zero map remains zero
    under a formal suffix; no nonzero physical continuation is inferred.
    """
    coefficients = _coefficients(coefficients)
    for event in _history(history):
        coefficients = branch_transport(coefficients, *event)
    return coefficients


def pair_kernel():
    """Exact row-stochastic kernel of the actual uniform-seam C/E process."""
    kernels = []
    for law in ("C", "E"):
        matrix = [[F(0)] * PAIR_DIMENSION for _ in PAIRS]
        for i in range(PAIR_DIMENSION):
            for seam in SEAMS:
                for coin in (0, 1):
                    _, target, sign, _ = cover_action(i, law, seam, coin)
                    if sign != 1:
                        raise ArithmeticError("pair kernel encountered a changed local sign")
                    matrix[i][target % PAIR_DIMENSION] += F(1, 2 * len(SEAMS))
        kernels.append(tuple(tuple(row) for row in matrix))
    if kernels[0] != kernels[1]:
        raise ArithmeticError("C and E no longer have the same pair population kernel")
    return kernels[0]


def advance(coefficients, steps):
    coefficients = _coefficients(coefficients)
    if type(steps) is not int or steps < 0:
        raise ValueError("steps must be a nonnegative exact integer")
    kernel = pair_kernel()
    edges = [(i, j, weight) for i, row in enumerate(kernel)
             for j, weight in enumerate(row) if weight]
    for _ in range(steps):
        output = [F(0)] * PAIR_DIMENSION
        for i, j, weight in edges:
            output[j] += coefficients[i] * weight
        coefficients = tuple(output)
    return coefficients


def pair_currents(coefficients):
    """Oriented net contrast flux J_ij=(w_i-w_j)T_ij, with i<j.

    These are currents on the pair configuration graph. Every edge moves
    one endpoint on one native seam; no additive one-port charge is inferred.
    Zero currents remain in the complete movement-edge inventory.
    """
    coefficients = _coefficients(coefficients)
    kernel = pair_kernel()
    if any(kernel[i][j] != kernel[j][i] for i in range(PAIR_DIMENSION)
           for j in range(PAIR_DIMENSION)):
        raise ArithmeticError("the net-current formula requires the actual symmetric kernel")
    currents = []
    for i in range(PAIR_DIMENSION):
        for j in range(i + 1, PAIR_DIMENSION):
            if kernel[i][j]:
                seam = tuple(sorted(set(PAIRS[i]).symmetric_difference(PAIRS[j])))
                if seam not in SEAMS:
                    raise ArithmeticError("pair current does not move along one native seam")
                currents.append((i, j, (coefficients[i] - coefficients[j]) * kernel[i][j]))
    return tuple(currents)


def region_flux(coefficients, region):
    """One-step regional balance, including arrivals and returns.

    The boundary separates pairs wholly inside the region from all other
    pairs. It is not an absorbing first-exit convention. Positive outflow and
    inflow count actual directed transitions; their difference is net flux.
    """
    coefficients = _coefficients(coefficients)
    region = set(_region(region))
    inside = tuple(set(pair) <= region for pair in PAIRS)
    kernel = pair_kernel()
    outgoing, incoming = F(0), F(0)
    for i in range(PAIR_DIMENSION):
        for j in range(PAIR_DIMENSION):
            weight = coefficients[i] * kernel[i][j]
            if inside[i] and not inside[j]:
                outgoing += weight
            if not inside[i] and inside[j]:
                incoming += weight
    net = sum((int(inside[i]) - int(inside[j])) * current
              for i, j, current in pair_currents(coefficients))
    after = advance(coefficients, 1)
    before_inside = sum((value for value, admitted in zip(coefficients, inside, strict=True)
                         if admitted), F(0))
    after_inside = sum((value for value, admitted in zip(after, inside, strict=True)
                        if admitted), F(0))
    if net != outgoing - incoming or after_inside - before_inside != -net:
        raise ArithmeticError("regional pair-current balance failed")
    return {"inside_before": before_inside, "inside_after": after_inside,
            "gross_outflow": outgoing, "gross_inflow": incoming,
            "net_outflow": net, "inside_change": after_inside - before_inside}


def current_certificate(coefficients, region):
    """Check every configuration's divergence, then the regional boundary sum."""
    coefficients = _coefficients(coefficients)
    region = _region(region)
    currents, divergence = pair_currents(coefficients), [F(0)] * PAIR_DIMENSION
    for i, j, current in currents:
        divergence[i] += current
        divergence[j] -= current
    after = advance(coefficients, 1)
    changes = tuple(b - a for a, b in zip(coefficients, after, strict=True))
    if any(flow != -change for flow, change in zip(divergence, changes, strict=True)):
        raise ArithmeticError("pair-current divergence disagrees with the complete channel")
    return {"configuration_vertices": len(divergence), "movement_edges": len(currents),
            "divergence": [str(value) for value in divergence],
            "one_attempt_change": [str(value) for value in changes],
            "nonzero_currents": [{"from_pair": PAIRS[i], "to_pair": PAIRS[j],
                                  "native_seam": tuple(sorted(set(PAIRS[i]) ^ set(PAIRS[j]))),
                                  "net_flux": str(current)}
                                 for i, j, current in currents if current],
            "region": region,
            "regional_balance": {key: str(value) for key, value in region_flux(coefficients, region).items()}}


def region_report(coefficients, region, epsilon=F(1)):
    """Exact regional charge and distance, including split-pair correlations."""
    coefficients = _coefficients(coefficients)
    region = set(_region(region))
    epsilon = _rational(epsilon, "epsilon", F(-1), F(1))
    weights = {"inside": F(0), "outside": F(0), "split": F(0)}
    for pair, value in zip(PAIRS, coefficients, strict=True):
        count = len(region.intersection(pair))
        weights[("outside", "split", "inside")[count]] += value
    total = sum(weights.values(), F(0))
    return {**weights, "total": total,
            "charge_expectation": epsilon * total / PAIR_DIMENSION,
            "regional_charge_expectation": epsilon * weights["inside"] / PAIR_DIMENSION,
            "global_trace_distance": abs(epsilon) * total / DIMENSION,
            "regional_trace_distance": abs(epsilon) * weights["inside"] / DIMENSION,
            "exterior_trace_distance": abs(epsilon) * weights["outside"] / DIMENSION}


def density_matrix(coefficients, epsilon=F(1)):
    """Floating evaluation of the exact rational state (I+epsilon sum f_x P_xQ)/132."""
    coefficients = _coefficients(coefficients)
    epsilon = _rational(epsilon, "epsilon", F(-1), F(1))
    result = np.eye(DIMENSION, dtype=complex) / DIMENSION
    for i, value in enumerate(coefficients):
        result[i, i + PAIR_DIMENSION] = result[i + PAIR_DIMENSION, i] = float(epsilon * value / DIMENSION)
    return result


def _matrix(matrix):
    # These helpers evaluate ordinary binary64 matrices, not arbitrary array
    # coercion protocols or extended-precision numerical certificates.
    if type(matrix) is not np.ndarray:
        raise ValueError("a plain numpy array is required; masked/subclass arrays are unsupported")
    value = matrix
    maximum_size = {"i": 8, "u": 8, "f": 8, "c": 16}
    if (value.shape != (DIMENSION, DIMENSION) or value.dtype.kind not in maximum_size
            or value.dtype.itemsize > maximum_size[value.dtype.kind]
            or not np.isfinite(value).all()):
        raise ValueError("a finite full 132-by-132 standard-precision numeric matrix required")
    if value.dtype.kind in "iu" and (np.any(value > 2**53)
                                     or (value.dtype.kind == "i" and np.any(value < -(2**53)))):
        raise ValueError("integer matrix entries must be exactly representable within the binary64 range")
    converted = np.asarray(value, dtype=np.complex128)
    if not np.isfinite(converted).all():
        raise ValueError("matrix conversion must remain finite")
    return converted


def region_density(matrix, region):
    """Actual qutrit tensor partial trace, in its complete supported local basis.

    The returned basis embeds into (C^3)^|region|. Omitted ambient basis
    states have identically zero rows and columns for every input in H132.
    No conditioning or trace normalization occurs.
    """
    value, region = _matrix(matrix), _region(region)
    outside = tuple(i for i in range(PORTS) if i not in region)
    configurations = tuple(cover_configuration(i) for i in range(DIMENSION))
    local = tuple(tuple(state[i] for i in region) for state in configurations)
    exterior = tuple(tuple(state[i] for i in outside) for state in configurations)
    basis = tuple(sorted(set(local)))
    indices = {state: i for i, state in enumerate(basis)}
    result = np.zeros((len(basis), len(basis)), complex)
    for i, left in enumerate(local):
        for j, right in enumerate(local):
            if exterior[i] == exterior[j]:
                result[indices[left], indices[right]] += value[i, j]
    return basis, result


def label_phase(matrix, port, *, dephase=False):
    """Added local A-sensitive service: a phase, or its retained/ignored read."""
    value = _matrix(matrix)
    port = _integer(port, 0, PORTS - 1, "phase port")
    if type(dephase) is not bool:
        raise ValueError("dephase must be a boolean")
    phases = np.array([-1 if cover_configuration(i)[port] == 1 else 1
                       for i in range(DIMENSION)])
    phased = value * phases[:, None] * phases[None, :]
    return (value + phased) / 2 if dephase else phased


def mixing_certificate():
    """Finite connectedness, holding and exact bistochasticity imply convergence."""
    kernel = pair_kernel()
    if any(sum(row) != 1 for row in kernel):
        raise ArithmeticError("pair kernel is not stochastic")
    if any(sum(row[j] for row in kernel) != 1 for j in range(PAIR_DIMENSION)):
        raise ArithmeticError("uniform pair reference is not stationary")
    if any(kernel[i][j] != kernel[j][i] for i in range(PAIR_DIMENSION) for j in range(PAIR_DIMENSION)):
        raise ArithmeticError("actual pair kernel lost reversal symmetry")
    reached, tree = {0}, []
    queue = [0]
    for i in queue:
        for j, weight in enumerate(kernel[i]):
            if weight and j not in reached:
                reached.add(j)
                queue.append(j)
                tree.append((i, j))
    holding = tuple(kernel[i][i] for i in range(PAIR_DIMENSION))
    if len(reached) != PAIR_DIMENSION or min(holding) <= 0:
        raise ArithmeticError("pair kernel lacks connected positive holding")
    return {"states": len(reached), "spanning_tree": tree,
            "minimum_holding_probability": str(min(holding)),
            "maximum_holding_probability": str(max(holding)),
            "stationary_pair_probability": str(F(1, PAIR_DIMENSION)),
            "limiting_coefficient": "sum(initial_coefficients)/66",
            "proof": "finite connected stochastic matrix with positive diagonal and uniform stationary law"}


def certificate():
    initial = initial_coefficients()
    witness = [{"attempts": step, **{key: str(value) for key, value in
                region_report(advance(initial, step), CORE).items()}}
               for step in (0, 1, 2, 3, 8, 16)]
    uniform = tuple(F(sum(initial), PAIR_DIMENSION) for _ in PAIRS)
    return {"local_intertwiner": intertwiner_certificate(),
            "mixing": mixing_certificate(), "core": CORE, "witness": witness,
            "asymptotic_region": {key: str(value) for key, value in region_report(uniform, CORE).items()},
            "pair_current_witness": current_certificate(initial, CORE),
            "return_flux_witness": current_certificate(advance(initial, 1), CORE),
            "scope": {"bare_binary_carrier": False, "additional_coherent_labels": True,
                      "localized_sharp_charge_sector": False,
                      "localized_conserved_charge_perturbation": True,
                      "physical_particle_identification": False,
                      "arbitrary_label_sensitive_errors_protected": False}}
