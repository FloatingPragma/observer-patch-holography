"""Finite projective radar with supplied preparation, cadence and turnaround services.

The source propagates throughout each half-unit interval. A target occupation
test follows every interval, up to twelve tests. A hit prepares the rank-one
target state, which is the immediate return preparation. The timeout branch
retains the final no-hit state and has no return operation.
"""

import math

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.scalar_channel.carrier import require


ENDPOINTS = (0, 5, 11, 16, 22, 29)
MARK = 2
TICKS = 12
DT = .5
LATTICE = 10**12
DOMAIN_RTOL = 1e-8
CONDITION_FLOOR = 1e-10
HELD_PAIR_RTOL = .05
CONFIDENCE_ALPHA = 1e-6
CONFIDENCE_DYADS = 30
SPECTRAL_NUMERICAL_GUARD = 1e-10
RADAR_SEED_NAMESPACE = 826401


def _settings(protocol):
    require(tuple(protocol['endpoints']) == ENDPOINTS, 'radar endpoints differ from the frozen menu')
    require(protocol['mark'] == MARK, 'radar mark differs from the frozen menu')
    require(protocol['dt'] == DT and protocol['ticks'] == TICKS, 'radar cadence differs from the frozen menu')
    for key in ('shots', 'seed'):
        value = protocol[key]
        require(type(value) is int and value >= (1 if key == 'shots' else 0), f'invalid radar {key}')


def _unitary(hamiltonian, dt):
    h = np.asarray(hamiltonian)
    require(h.ndim == 2 and h.shape[0] == h.shape[1] and np.isfinite(h).all(), 'invalid radar Hamiltonian')
    require(np.max(abs(h - h.conj().T)) < 1e-12, 'radar Hamiltonian is not Hermitian')
    values, modes = eigh(h)
    return (modes * np.exp(-1j * dt * values)) @ modes.conj().T


def first_hit(unitary, endpoints, ticks):
    """Sequential projectors, with the norm of every surviving branch retained.

    Axes are source, target, outcome. Outcome k < ticks means first hit at
    tick k+1. The last outcome is timeout. Diagonal menus are unused zeros.
    """
    unitary = np.asarray(unitary, dtype=complex)
    endpoints = tuple(endpoints)
    laws = np.zeros((len(endpoints), len(endpoints), ticks + 1))
    for target, site in enumerate(endpoints):
        survivors = np.eye(len(unitary), dtype=complex)[:, endpoints]
        for tick in range(ticks):
            survivors = unitary @ survivors
            laws[:, target, tick] = abs(survivors[site])**2
            survivors[site] = 0
        laws[:, target, -1] = np.sum(abs(survivors)**2, axis=0)
    laws[np.arange(len(endpoints)), np.arange(len(endpoints))] = 0
    mask = ~np.eye(len(endpoints), dtype=bool)
    require(np.max(abs(laws.sum(axis=-1)[mask] - 1)) < 1e-10, 'first-hit instrument loses probability')
    return laws


def renewal_prediction(unitary, endpoints, ticks):
    """Unmeasured powers and the renewal amplitude equation, without projections."""
    endpoints = tuple(endpoints)
    size = len(endpoints)
    amplitudes = np.empty((ticks + 1, size, size), complex)
    power = np.eye(len(unitary), dtype=complex)
    for tick in range(ticks + 1):
        amplitudes[tick] = power[np.ix_(endpoints, endpoints)]
        if tick < ticks:
            power = unitary @ power
    laws = np.zeros((size, size, ticks + 1))
    for source in range(size):
        for target in range(size):
            if source == target:
                continue
            phi = np.zeros(ticks + 1, complex)
            for tick in range(1, ticks + 1):
                phi[tick] = amplitudes[tick, target, source] - sum(
                    amplitudes[tick - earlier, target, target] * phi[earlier]
                    for earlier in range(1, tick))
            laws[source, target, :-1] = abs(phi[1:])**2
            timeout = 1 - np.sum(abs(phi[1:])**2)
            require(timeout >= -1e-10, 'renewal law has negative timeout weight')
            laws[source, target, -1] = max(0., timeout)
    return laws


def _frozen_unitary(protocol, coupling):
    from .joint_source import radar_source

    _settings(protocol)
    require(type(coupling) is int and coupling in (0, 1), 'mark intervention must be 0 or 1')
    h = np.asarray(radar_source(protocol, coupling))
    require(h.shape == (150, 150) and np.isrealobj(h), 'radar source must be the real 150-mode restriction')
    return _unitary(h, protocol['dt'])


def leg_probabilities(protocol, coupling):
    """Six by six by thirteen complete leg laws; the six diagonal menus are unused."""
    return first_hit(_frozen_unitary(protocol, coupling), ENDPOINTS, TICKS)


def round_trip_law(legs):
    """Joint first-hit records with one explicit outbound-timeout/no-return bin.

    Both outcome axes have a terminal timeout symbol. A terminal outbound
    symbol and a return hit is an impossible record with zero weight.
    """
    size, _, outcomes = legs.shape
    joint = np.zeros((size, size, outcomes, outcomes))
    for source in range(size):
        for target in range(size):
            if source != target:
                joint[source, target, :-1] = np.outer(legs[source, target, :-1], legs[target, source])
                joint[source, target, -1, -1] = legs[source, target, -1]
    return joint


def acquisition_weights(probabilities):
    """Integer probability lattice; the largest bin absorbs the rounding remainder."""
    p = np.asarray(probabilities, float)
    require(np.isfinite(p).all() and np.all(p >= 0), 'invalid acquisition probabilities')
    require(abs(p.sum() - 1) < 1e-10, 'incomplete acquisition probabilities')
    weights = np.floor(p * LATTICE + .5).astype(np.int64)
    flat = weights.ravel()
    flat[int(np.argmax(p))] += LATTICE - int(flat.sum())
    require(np.all(weights >= 0) and int(weights.sum()) == LATTICE, 'invalid acquisition lattice')
    return weights


def replay_counts(weights, shots, seed, coupling, source, target):
    """Reproduce counts from retained integer weights, independent of floating point."""
    weights = np.asarray(weights)
    require(weights.dtype.kind in 'iu' and np.all(weights >= 0) and int(weights.sum()) == LATTICE,
            'invalid retained acquisition weights')
    generator = np.random.Generator(np.random.PCG64(np.random.SeedSequence(
        [seed, RADAR_SEED_NAMESPACE, coupling, source, target])))
    cumulative = np.cumsum(weights.ravel(), dtype=np.int64)
    counts = np.zeros(weights.size, dtype=np.int64)
    for start in range(0, shots, 65536):
        draws = generator.integers(0, LATTICE, size=min(65536, shots - start), dtype=np.int64)
        bins = np.searchsorted(cumulative, draws, side='right')
        counts += np.bincount(bins, minlength=weights.size)
    return counts.reshape(weights.shape)


def distance_statistics(joint, dt):
    """Symmetric half-round-trip means from successful retained records alone."""
    size = joint.shape[0]
    ticks = joint.shape[-1] - 1
    distances = np.full((size, size), np.nan)
    np.fill_diagonal(distances, 0.)
    successes = np.zeros((size, size))
    half_times = dt * (np.arange(1, ticks + 1)[:, None] + np.arange(1, ticks + 1)[None, :]) / 2
    for left in range(size):
        for right in range(left + 1, size):
            retained = joint[left, right, :-1, :-1] + joint[right, left, :-1, :-1]
            total = float(retained.sum())
            successes[left, right] = successes[right, left] = total
            if total > 0:
                distances[left, right] = distances[right, left] = float(np.sum(retained * half_times) / total)
    return distances, successes


def predict_withheld_pair(distances, left, right):
    """Three-dimensional trilateration using fourteen pair distances, never the fifteenth.

    Four anchors define a Gram matrix. The differences of squared distances
    locate each withheld endpoint in that frame. The distance to the first
    anchor is an independent domain check for each endpoint.
    """
    distances = np.asarray(distances, float).copy()
    require(distances.shape == (6, 6) and 0 <= left < right < 6, 'invalid six-endpoint holdout')
    distances[left, right] = distances[right, left] = np.nan
    anchors = [i for i in range(6) if i not in (left, right)]
    result = dict(pair=[left, right], anchors=anchors, predicted_distance=None,
                  domain_valid=False, reason=None, anchor_eigenvalues=None,
                  radial_relative_errors=None)
    training = [distances[i, j] for i in range(6) for j in range(i + 1, 6) if (i, j) != (left, right)]
    if not np.isfinite(training).all() or np.any(np.asarray(training) <= 0):
        result['reason'] = 'missing_or_nonpositive_training_distance'
        return result
    origin, *basis = anchors
    radii = distances[origin, basis]**2
    gram = (radii[:, None] + radii[None, :] - distances[np.ix_(basis, basis)]**2) / 2
    eigenvalues = np.linalg.eigvalsh(gram)
    result['anchor_eigenvalues'] = eigenvalues.tolist()
    scale = max(float(np.max(abs(gram))), np.finfo(float).tiny)
    if eigenvalues[0] < -DOMAIN_RTOL * scale:
        result['reason'] = 'non_euclidean_anchor_distances'
        return result
    if eigenvalues[0] <= CONDITION_FLOOR * max(float(eigenvalues[-1]), scale):
        result['reason'] = 'singular_or_ill_conditioned_anchor_frame'
        return result
    rhs = np.array([(distances[p, origin]**2 + radii - distances[p, basis]**2) / 2
                    for p in (left, right)])
    coefficients = np.linalg.solve(gram, rhs.T).T
    difference = coefficients[0] - coefficients[1]
    predicted_squared = float(difference @ gram @ difference)
    result['predicted_distance'] = math.sqrt(max(0., predicted_squared))
    radial_errors = [abs(float(rhs[k] @ coefficients[k]) - distances[p, origin]**2)
                     / max(distances[p, origin]**2, scale) for k, p in enumerate((left, right))]
    result['radial_relative_errors'] = radial_errors
    if max(radial_errors) > DOMAIN_RTOL:
        result['reason'] = 'endpoint_distances_incompatible_with_anchor_frame'
        return result
    result['domain_valid'] = True
    result['reason'] = 'inside_trilateration_domain'
    return result


def geometry_test(distances):
    """Evaluate every predeclared held pair, retaining each domain failure."""
    tests = []
    for left in range(6):
        for right in range(left + 1, 6):
            result = predict_withheld_pair(distances, left, right)
            observed = float(distances[left, right])
            result['observed_distance'] = observed if math.isfinite(observed) else None
            result['relative_error'] = None
            result['accepted'] = False
            if result['domain_valid'] and math.isfinite(observed) and observed > 0:
                result['relative_error'] = abs(result['predicted_distance'] - observed) / observed
                result['accepted'] = result['relative_error'] <= HELD_PAIR_RTOL
            tests.append(result)
    return dict(dimension=3, domain_relative_tolerance=DOMAIN_RTOL,
                anchor_condition_floor=CONDITION_FLOOR, held_pair_relative_tolerance=HELD_PAIR_RTOL,
                valid_domains=sum(t['domain_valid'] for t in tests),
                accepted_pairs=sum(t['accepted'] for t in tests),
                all_pairs_accepted=all(t['accepted'] for t in tests), pairs=tests)


def euclidean_error_bound(distances, lower=None, upper=None):
    """Necessary relative fit error in dimensions one through five.

    For S_ij=D_ij^2 and J=I-11'/6, B=-JSJ/2. If every relative
    distance error is at most delta, the Gram perturbation norm is at most
    (2*delta+delta^2)*max_i sum_j S_ij/2. Weyl bounds a rank-d positive
    semidefinite approximation by both lambda_(d+1)(B) and -lambda_min(B).
    Confidence intervals give a separate norm bound on the unknown Gram
    matrix. The numerical guard is deducted from every spectral obstruction.
    """
    distances = np.asarray(distances, float)
    require(distances.shape == (6, 6), 'spectral bound requires six endpoints')
    uncertain = lower is not None or upper is not None
    if uncertain:
        require(lower is not None and upper is not None, 'both distance envelopes are required')
        lower, upper = np.asarray(lower, float), np.asarray(upper, float)
        require(lower.shape == upper.shape == (6, 6) and np.isfinite(lower).all()
                and np.isfinite(upper).all() and np.all(lower >= 0) and np.all(lower <= upper),
                'invalid distance confidence envelopes')
        center = np.where(np.isfinite(distances), distances, (lower + upper) / 2)
        square_error = np.maximum(abs(upper**2 - center**2), abs(center**2 - lower**2))
        uncertainty_norm = float(np.max(square_error.sum(axis=1)) / 2)
        relative_error_scale = float(np.max((upper**2).sum(axis=1)) / 2)
    elif not np.isfinite(distances).all():
        return dict(defined=False, reason='missing_conditional_distance', dimensions=[])
    else:
        center = distances
        uncertainty_norm = 0.
        relative_error_scale = float(np.max((distances**2).sum(axis=1)) / 2)
    require(np.all(center >= 0) and np.max(abs(center - center.T)) < 1e-12,
            'spectral distance matrix must be symmetric and nonnegative')
    centering = np.eye(6) - np.ones((6, 6)) / 6
    gram = -.5 * centering @ center**2 @ centering
    eigenvalues, modes = np.linalg.eigh(gram)
    numerical_guard = SPECTRAL_NUMERICAL_GUARD * max(1., float(np.linalg.norm(gram, ord='fro')))
    reconstruction_error = float(np.linalg.norm(gram - (modes * eigenvalues) @ modes.T, ord='fro'))
    require(reconstruction_error < numerical_guard / 10, 'Gram eigensolver exceeds the numerical guard')
    dimensions = []
    for dimension in range(1, 6):
        positive_tail = max(0., float(eigenvalues[-dimension - 1]))
        negative_part = max(0., -float(eigenvalues[0]))
        obstruction = max(0., max(positive_tail, negative_part) - uncertainty_norm - numerical_guard)
        # This form avoids cancellation for small lower bounds.
        ratio = obstruction / relative_error_scale if relative_error_scale > 0 else 0.
        bound = ratio / (math.sqrt(1 + ratio) + 1)
        dimensions.append(dict(dimension=dimension, positive_tail_eigenvalue=positive_tail,
                               negative_eigenvalue_obstruction=negative_part,
                               surviving_gram_obstruction=obstruction,
                               minimum_possible_maximum_relative_error=bound,
                               excludes_five_percent_fit=bound > HELD_PAIR_RTOL))
    return dict(defined=True, gram_eigenvalues=eigenvalues.tolist(),
                distance_uncertainty_operator_bound=uncertainty_norm,
                relative_error_operator_coefficient=relative_error_scale,
                numerical_guard=numerical_guard, eigensolver_reconstruction_error=reconstruction_error,
                statement='necessary error lower bound for every Euclidean embedding in the stated dimension',
                dimensions=dimensions)


def _lattice_conditional_bias(weights, counts, dt, ticks):
    """Bound quantization bias without reading unrounded source probabilities."""
    bias = np.zeros((6, 6))
    # Rounding a complete outcome law and assigning the remainder to one bin
    # costs at most outcomes/(2*lattice) in total variation. The additive term
    # covers the instrument's allowed probability-sum roundoff.
    unconditional_tv = (ticks + 1)**2 / (2 * LATTICE) + 1e-10
    width = (ticks - 1) * dt
    for left in range(6):
        for right in range(left + 1, 6):
            total, weighted_bias = 0, 0.
            for source, target in ((left, right), (right, left)):
                successes = int(counts[source, target, :-1, :-1].sum())
                success_probability = float(weights[source, target, :-1, :-1].sum() / LATTICE)
                if successes:
                    require(success_probability > 0, 'successful count has zero acquisition weight')
                    weighted_bias += successes * width * min(1., 2 * unconditional_tv / success_probability)
                total += successes
            if total:
                bias[left, right] = bias[right, left] = weighted_bias / total
    return bias


def confidence_envelopes(distances, successful_counts, dt, ticks, lattice_bias=None):
    """Hoeffding intervals, simultaneous over both preparations and all fifteen dyads.

    Conditional on the actual success counts, successful half-trip times are
    independent and lie in [dt,ticks*dt]. Reversing the two legs exchanges
    factors of the same convolution. Integer-lattice sampling contributes the
    supplied additive bias bound. A pair without successes retains the full
    allowed interval and supplies no estimated mean.
    """
    distances = np.asarray(distances, float)
    successful_counts = np.asarray(successful_counts)
    if lattice_bias is None:
        lattice_bias = np.zeros_like(distances)
    lower, upper = np.zeros((6, 6)), np.zeros((6, 6))
    statistical_radii = np.zeros((6, 6))
    for left in range(6):
        for right in range(left + 1, 6):
            count = successful_counts[left, right]
            if count > 0 and math.isfinite(float(distances[left, right])):
                radius = (ticks - 1) * dt * math.sqrt(math.log(2 * CONFIDENCE_DYADS / CONFIDENCE_ALPHA) / (2 * count))
                center = float(distances[left, right])
                total_radius = radius + lattice_bias[left, right]
                lo, hi = max(dt, center - total_radius), min(ticks * dt, center + total_radius)
                statistical_radii[left, right] = statistical_radii[right, left] = radius
            else:
                lo, hi = dt, ticks * dt
            lower[left, right] = lower[right, left] = lo
            upper[left, right] = upper[right, left] = hi
    summary = dict(family_failure_probability=CONFIDENCE_ALPHA, simultaneous_dyads=CONFIDENCE_DYADS,
                   conditional_sample_range=[dt, ticks * dt],
                   lower=lower.tolist(), upper=upper.tolist(),
                   statistical_half_width=statistical_radii.tolist(),
                   acquisition_lattice_bias_bound=np.asarray(lattice_bias).tolist(),
                   no_success_rule='full allowed interval with undefined sample mean',
                   conditioning='actual successful counts; reversed duration laws agree by convolution')
    return summary, lower, upper


def common_multiplier_intersection(base_lower, base_upper, marked_lower, marked_upper):
    """Necessary positive common scale, tested by intersecting all dyad ratio intervals."""
    pairs = []
    for left in range(6):
        for right in range(left + 1, 6):
            entries = [base_lower[left, right], base_upper[left, right],
                       marked_lower[left, right], marked_upper[left, right]]
            if not np.isfinite(entries).all() or min(entries) <= 0:
                return dict(defined=False, reason='missing_or_nonpositive_conditional_distance',
                            interval=None, common_multiplier_excluded=False, pairs=pairs)
            low = float(marked_lower[left, right] / base_upper[left, right])
            high = float(marked_upper[left, right] / base_lower[left, right])
            pairs.append(dict(pair=[left, right], lower=low, upper=high))
    low, high = max(p['lower'] for p in pairs), min(p['upper'] for p in pairs)
    return dict(defined=True, interval=[low, high] if low <= high else None,
                intersection_lower=low, intersection_upper=high,
                common_multiplier_excluded=low > high, pairs=pairs,
                interpretation='empty intersection excludes uniform scaling of these conditional mean records')


def _json_matrix(matrix):
    return [[float(x) if math.isfinite(float(x)) else None for x in row] for row in matrix]


def _pair_summaries(joint, counts, dt):
    size = len(joint)
    ticks = joint.shape[-1] - 1
    durations = np.r_[np.arange(1, ticks + 1), ticks]
    elapsed = durations[:, None] + np.r_[np.ones(ticks), 0][:, None] * durations[None, :]
    rows = []
    for source in range(size):
        for target in range(size):
            if source == target:
                continue
            p, c = joint[source, target], counts[source, target]
            successful = c[:-1, :-1]
            successful_total = int(successful.sum())
            round_trip_mean = None if successful_total == 0 else float(
                dt * np.sum(successful * elapsed[:-1, :-1]) / successful_total)
            rows.append(dict(source=source, target=target,
                             outbound_timeout_probability=float(p[-1, -1]),
                             return_timeout_probability=float(p[:-1, -1].sum()),
                             successful_round_trip_probability=float(p[:-1, :-1].sum()),
                             outbound_timeout_count=int(c[-1, -1]),
                             return_timeout_count=int(c[:-1, -1].sum()),
                             successful_round_trip_count=successful_total,
                             conditional_round_trip_mean=round_trip_mean,
                             acquired_projector_tests=int(np.sum(c * elapsed)),
                             expected_projector_tests_per_trial=float(np.sum(p * elapsed)),
                             acquired_source_action_time=float(dt * np.sum(c * elapsed))))
    return rows


def run(protocol):
    """Return a JSON-ready summary and numerical tables for a compressed receipt.

    This function executes the scientific acquisition. Synthetic tests use the
    matrix-level functions and do not call the contact source.
    """
    _settings(protocol)
    direct, renewal, joint_laws, acquired, lattices = [], [], [], [], []
    settings, exact_distance_matrices, lower_envelopes, upper_envelopes = [], [], [], []
    for coupling in (0, 1):
        unitary = _frozen_unitary(protocol, coupling)
        legs = first_hit(unitary, ENDPOINTS, TICKS)
        predicted = renewal_prediction(unitary, ENDPOINTS, TICKS)
        residual = float(np.max(abs(legs - predicted)))
        require(residual < 1e-10, 'renewal prediction differs from sequential projectors')
        joint = round_trip_law(legs)
        counts = np.zeros(joint.shape, dtype=np.int64)
        weights = np.zeros(joint.shape, dtype=np.int64)
        for source in range(6):
            for target in range(6):
                if source != target:
                    weights[source, target] = acquisition_weights(joint[source, target])
                    counts[source, target] = replay_counts(weights[source, target], protocol['shots'],
                                                          protocol['seed'], coupling, source, target)
        exact_distances, _ = distance_statistics(joint, DT)
        observed_distances, successful_counts = distance_statistics(counts, DT)
        confidence, lower, upper = confidence_envelopes(
            observed_distances, successful_counts, DT, TICKS,
            _lattice_conditional_bias(weights, counts, DT, TICKS))
        pairs = _pair_summaries(joint, counts, DT)
        settings.append(dict(mark_written=coupling, renewal_maximum_absolute_error=residual,
                             maximum_acquisition_lattice_total_variation=float(np.max(
                                 .5 * abs(weights / LATTICE - joint).sum(axis=(-2, -1)))),
                             distances_from_born_laws=_json_matrix(exact_distances),
                             distances_from_counts=_json_matrix(observed_distances),
                             successful_counts_by_unordered_pair=successful_counts.astype(np.int64).tolist(),
                             geometry_from_born_laws=geometry_test(exact_distances),
                             geometry_from_counts=geometry_test(observed_distances),
                             euclidean_error_from_born_laws=euclidean_error_bound(exact_distances),
                             distance_confidence=confidence,
                             euclidean_error_with_simultaneous_confidence=euclidean_error_bound(
                                 observed_distances, lower, upper),
                             pairs=pairs,
                             source_emissions=30 * protocol['shots'],
                             return_turnarounds=sum(protocol['shots'] - p['outbound_timeout_count'] for p in pairs),
                             acquired_projector_tests=sum(p['acquired_projector_tests'] for p in pairs),
                             acquired_source_action_time=sum(p['acquired_source_action_time'] for p in pairs)))
        direct.append(legs)
        renewal.append(predicted)
        joint_laws.append(joint)
        acquired.append(counts)
        lattices.append(weights)
        exact_distance_matrices.append(exact_distances)
        lower_envelopes.append(lower)
        upper_envelopes.append(upper)
    arrays = dict(radar_leg_probabilities=np.asarray(direct),
                  radar_renewal_probabilities=np.asarray(renewal),
                  radar_joint_probabilities=np.asarray(joint_laws),
                  radar_joint_counts=np.asarray(acquired, dtype=np.int64),
                  radar_acquisition_weights=np.asarray(lattices, dtype=np.int64),
                  radar_distance_confidence_lower=np.asarray(lower_envelopes),
                  radar_distance_confidence_upper=np.asarray(upper_envelopes))
    summary = dict(
        endpoints=list(ENDPOINTS), mark=MARK, contact_coupling=1, baseline_endpoint_potential=1,
        ticks_per_leg=TICKS, interval=DT, shots_per_ordered_pair=protocol['shots'],
        intervention_values=[0, 1], diagonal_menu='unused_zero_sentinel',
        leg_outcomes='indices 0 through 11 mean first hit at tick index+1; index 12 means timeout',
        joint_axes=['mark_written', 'source_endpoint', 'target_endpoint', 'outbound', 'return'],
        no_return_record='outbound timeout occupies (12,12); cells (12,0..11) are impossible',
        rng=dict(bit_generator='PCG64', namespace=RADAR_SEED_NAMESPACE, integer_lattice=LATTICE,
                 entropy_order=['seed', 'namespace', 'mark_written', 'source_endpoint', 'target_endpoint'],
                 receipt_replay='regenerate from retained integer weights; compare Born laws with numerical tolerance'),
        supplied_services=[
            'exact scalar site preparation with a fixed parity companion',
            'projective local occupation tests and classical first-hit recording',
            'external cadence of half an exchange-time unit',
            'instantaneous local return after a rank-one hit with no added delay',
            'trial reset and retention of all timeout and no-return records'],
        interpretation='simulation-generated first-hit records; conditional half-round-trip means in supplied cadence units',
        limits=[
            'no implementation of these ideal services by autonomous endpoint matter is asserted',
            'successful-trip conditioning is explicit and timeout weights are retained',
            'trilateration tests three-dimensional Euclidean consistency of retained means',
            'the raw count trilateration is descriptive; confidence claims use the Gram-spectrum bounds',
            'a spectral lower bound below five percent does not establish an embedding',
            'no physical length, light speed, causal cone, spacetime metric or cosmological expansion is assigned'],
        settings=settings,
        common_multiplier_from_born_laws=common_multiplier_intersection(
            exact_distance_matrices[0] - 1e-10, exact_distance_matrices[0] + 1e-10,
            exact_distance_matrices[1] - 1e-10, exact_distance_matrices[1] + 1e-10),
        common_multiplier_born_distance_roundoff=1e-10,
        common_multiplier_with_simultaneous_confidence=common_multiplier_intersection(
            lower_envelopes[0], upper_envelopes[0], lower_envelopes[1], upper_envelopes[1]),
        mark_joint_total_variation=(.5 * abs(arrays['radar_joint_probabilities'][1]
                                           - arrays['radar_joint_probabilities'][0]).sum(axis=(-2, -1))).tolist())
    return summary, arrays
