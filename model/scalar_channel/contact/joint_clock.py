"""Complete joint records from one collar clock and the shared scalar source.

The initial state is blank collar/fuel-present times one declared scalar
mode. Two finite pi/4 pulses retain the source throughout. Between them the
source acts for the declared wait, followed by a phase on the blank arm.
The terminal service copies the orthogonal arm and scalar-mode records.
Its 300 outcomes include both fuel results and every seam and wedge.

``run`` returns a JSON-ready summary and numeric archive arrays. ``analyze``
reads those arrays without propagating a state. Source-on sparse action
executes the records; a separate spectral path sum forecasts them without
fitting probabilities or clock rates. No continuum, laboratory-time or
particle-mass identification is supplied by this finite experiment.
"""

import hashlib
import math

import numpy as np
from scipy import sparse as sp
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from ..carrier import require
from .joint_source import ENDPOINTS, MARK, MODES, arm_sources, radar_source, source_certificate, validate_source

PHASES = (0.0, math.pi / 2, math.pi, 3 * math.pi / 2)
COUPLINGS = (0, 1, -1)
SIGNS = (1, -1)
SEEDS = ENDPOINTS + (MARK,)
SAMPLING_LATTICE = 10 ** 12
FAMILIES = (('joint', 2), ('arm_path', 3), ('pure_arm', 2), ('dephased', 2),
            ('independent_clock', 1), ('independent_scalar', 1))
# Digest of the original six integer arrays at 44e455330.
ACQUISITION_DIGEST = '546d1450118a6d66397e7db5ad25914cf5174da3615252d518ea1557a2a64c8d'
SOURCE_BIN_TOLERANCE = .5 / SAMPLING_LATTICE
SOURCE_TV_TOLERANCE = 1 / SAMPLING_LATTICE
MODELS = ('native_positive_mark', 'marked_contact_suppressed',
          'marked_contact_reversed', 'independent_records_product', 'dephased_before_read')


def validate_protocol(protocol):
    validate_source(protocol)
    require(tuple(protocol['couplings']) == COUPLINGS, 'joint contact menu differs')
    require(tuple(protocol['joint_signs']) == SIGNS, 'joint collar sign menu differs')
    require(np.asarray(protocol['joint_phases']).shape == (4,)
            and np.max(abs(np.asarray(protocol['joint_phases']) - PHASES)) < 1e-14,
            'joint reader requires all four declared phases')
    require(tuple(protocol['joint_times']) == (0.5, 2, 6), 'joint wait menu differs')
    require(protocol['joint_drive_strength'] == 64, 'joint finite-pulse strength differs')
    require(type(protocol['shots']) is int and protocol['shots'] > 0, 'invalid shot count')
    require(type(protocol['seed']) is int and protocol['seed'] >= 0, 'invalid sampling seed')
    require(0 < protocol['joint_familywise_alpha'] < 1, 'invalid joint familywise alpha')
    require(0 < protocol['exact_forecast_tolerance'] <= 1e-8, 'invalid exact forecast tolerance')


def _step(H, state, duration):
    return expm_multiply(-1j * duration * H, state,
                         traceA=-1j * duration * H.diagonal().sum())


def _source(protocol, coupling, sign):
    H0, H1 = arm_sources(protocol, coupling, sign)
    free = sp.block_diag((sp.csr_matrix(H0), sp.csr_matrix(H1)), format='csr')
    I = sp.eye(MODES, format='csr')
    pulse = free + protocol['joint_drive_strength'] * sp.bmat([[None, I], [I, None]], format='csr')
    return free, pulse


def execute(protocol, coupling, sign):
    """Sparse source-on history, including retained intermediate arm records.

    Arrays use seed, time, phase, terminal arm, scalar mode. The arm-path
    control instead uses seed, time, intermediate arm, terminal arm, mode.
    Pure-arm controls use seed, time, prepared arm, terminal arm, mode.
    """
    free, driven = _source(protocol, coupling, sign)
    tau = math.pi / (4 * protocol['joint_drive_strength'])
    seed = np.eye(2 * MODES, dtype=complex)[:, SEEDS]
    prepared = _step(driven, seed, tau)
    shape = (len(SEEDS), len(protocol['joint_times']))
    law = np.zeros((*shape, 4, 2, MODES))
    paths = np.zeros((*shape, 2, 2, MODES))
    pure = np.zeros_like(paths)
    independent_clock = np.zeros((*shape, 4, 2))
    independent_scalar = np.zeros((*shape, 4, MODES))
    # Separate fresh ensembles implement one pointer at a time. They do not
    # turn two independently sampled marginals into coincident observations.
    for ti, wait in enumerate(protocol['joint_times']):
        waiting = _step(free, prepared, wait)
        for pi, phase in enumerate(PHASES):
            phased = waiting.copy()
            phased[:MODES] *= np.exp(-1j * phase)
            final = _step(driven, phased, tau).reshape(2, MODES, len(SEEDS))
            p = abs(final) ** 2
            law[:, ti, pi] = p.transpose(2, 0, 1)
            independent_clock[:, ti, pi] = p.sum(axis=1).T
            independent_scalar[:, ti, pi] = p.sum(axis=0).T
        for arm in (0, 1):
            retained = waiting.copy()
            retained[(1 - arm) * MODES:(2 - arm) * MODES] = 0
            final = _step(driven, retained, tau).reshape(2, MODES, len(SEEDS))
            paths[:, ti, arm] = (abs(final) ** 2).transpose(2, 0, 1)
            start = np.zeros_like(seed)
            start[arm * MODES + np.array(SEEDS), np.arange(len(SEEDS))] = 1
            final = _step(free, start, wait).reshape(2, MODES, len(SEEDS))
            pure[:, ti, arm] = (abs(final) ** 2).transpose(2, 0, 1)
    return law, paths, pure, independent_clock, independent_scalar


def _spectral_unitary(H, duration):
    values, modes = eigh(H)
    return (modes * np.exp(-1j * duration * values)) @ modes.conj().T


def forecast(protocol, coupling, sign):
    """Independent finite-pulse branch composition on the same frozen source.

    For each intermediate arm a, A[b,a] = P[b,a] exp(-it H[a]) P[a,0]
    acts on the scalar seed. The joint law is |exp(-i theta) A[b,0] +
    A[b,1]|^2. Summing |A[b,a]|^2 gives the intervened arm measurement.
    """
    # The forecast constructs its branches from the frozen scalar source
    # and the native empty-load eigenvalues 0 and kappa+(1-sign)/2. The
    # source certificate compares these values with the full native law.
    H0 = radar_source(protocol, 0)
    H1 = radar_source(protocol, coupling) + (protocol['kappa'] + (1 - sign) / 2) * np.eye(MODES)
    source = np.block([[H0, protocol['joint_drive_strength'] * np.eye(MODES)],
                       [protocol['joint_drive_strength'] * np.eye(MODES), H1]])
    pulse = _spectral_unitary(source, math.pi / (4 * protocol['joint_drive_strength']))
    eig = [eigh(H0), eigh(H1)]
    amplitudes = np.zeros((len(SEEDS), len(protocol['joint_times']), 2, 2, MODES), complex)
    for ti, wait in enumerate(protocol['joint_times']):
        for arm in (0, 1):
            values, modes = eig[arm]
            prepared = pulse[arm * MODES:(arm + 1) * MODES, SEEDS]
            waiting = modes @ (np.exp(-1j * wait * values)[:, None] * (modes.conj().T @ prepared))
            final = pulse[:, arm * MODES:(arm + 1) * MODES] @ waiting
            amplitudes[:, ti, arm] = final.reshape(2, MODES, len(SEEDS)).transpose(2, 0, 1)
    law = np.stack([abs(np.exp(-1j * phase) * amplitudes[:, :, 0] + amplitudes[:, :, 1]) ** 2
                    for phase in PHASES], axis=2)
    return law, amplitudes


def canonical_sampling_weights(probabilities, outcome_axes=2):
    """Integer masses sum to 10^12; the first largest rounded bin takes the remainder.

    For K categories the total variation change is at most K/(2*10^12).
    The original full probability law is retained separately in the archive.
    """
    p = np.asarray(probabilities)
    categories = int(np.prod(p.shape[-outcome_axes:]))
    rows = p.reshape(-1, categories)
    require(np.min(rows) >= -1e-14 and np.max(abs(rows.sum(axis=1) - 1)) < 1e-10,
            'cannot sample an incomplete record law')
    normalized = np.maximum(rows, 0) / np.maximum(rows, 0).sum(axis=1, keepdims=True)
    weights = np.rint(normalized * SAMPLING_LATTICE).astype(np.int64)
    largest = np.argmax(weights, axis=1)
    weights[np.arange(len(rows)), largest] += SAMPLING_LATTICE - weights.sum(axis=1)
    require(np.min(weights) >= 0 and np.all(weights.sum(axis=1) == SAMPLING_LATTICE),
            'canonical sampling masses are invalid')
    return weights.reshape(p.shape)


def _sample(probabilities, shots, rng, outcome_axes=2):
    weights = canonical_sampling_weights(probabilities, outcome_axes)
    categories = int(np.prod(weights.shape[-outcome_axes:]))
    result = np.array([rng.multinomial(shots, row / SAMPLING_LATTICE)
                       for row in weights.reshape(-1, categories)], dtype=np.int64)
    return result.reshape(weights.shape)


def authenticated_weights(archive):
    names = sorted(name + '_sampling_weights' for name, _ in FAMILIES)
    require(sorted(archive) == names, 'missing, duplicate or extra clock weight family')
    digest = hashlib.sha256()
    weights = {}
    for name in names:
        value = archive[name]
        require(value.dtype.kind == 'i' and value.dtype.itemsize == 8, 'clock weights require signed int64')
        digest.update(name.encode() + b'\0' + value.astype('<i8', copy=False).tobytes())
        weights[name] = value
    require(digest.hexdigest() == ACQUISITION_DIGEST, 'original integer clock weights differ')
    return weights


def check_acquisition_law(probability, weights, axes):
    """Compare each original lattice row with its separately propagated probability law."""
    require(weights.shape == probability.shape and np.all((weights >= 0) & (weights <= SAMPLING_LATTICE)),
            'clock weight shape or range differs')
    bins = math.prod(weights.shape[-axes:])
    rows = weights.reshape(-1, bins)
    require(np.all(rows.sum(axis=1) == SAMPLING_LATTICE), 'clock weight row has the wrong total')
    p = probability.reshape(-1, bins)
    p = p / p.sum(axis=1, keepdims=True)
    error = abs(rows / SAMPLING_LATTICE - p)
    ordinary = .5 / SAMPLING_LATTICE + SOURCE_BIN_TOLERANCE
    # The first largest rounded bin receives the remainder. Two bins can tie
    # after rounding when their source probabilities differ by one grid unit.
    near_maximum = p >= p.max(axis=1, keepdims=True) - (1 / SAMPLING_LATTICE + 2 * SOURCE_BIN_TOLERANCE)
    exceptional = error > ordinary
    require(np.all(exceptional.sum(axis=1) <= 1) and not np.any(exceptional & ~near_maximum),
            'ordinary clock weights exceed the source bin tolerance')
    require(error.max() <= bins / (2 * SAMPLING_LATTICE) + SOURCE_BIN_TOLERANCE,
            'clock remainder exceeds the source bin tolerance')
    tv = .5 * error.sum(axis=1)
    require(tv.max() <= bins / (2 * SAMPLING_LATTICE) + SOURCE_TV_TOLERANCE,
            'clock weights exceed the source total-variation tolerance')
    return dict(rows=len(rows), bins_per_row=bins, maximum_bin_error=float(error.max()),
                maximum_total_variation=float(tv.max()),
                locally_requantized_bins_different=int(np.count_nonzero(
                    weights != canonical_sampling_weights(probability, axes))))


def coarse_grain(values):
    """Seven named seam outcomes and the sum of every other scalar mode."""
    array = np.asarray(values)
    rest = np.ones(MODES, dtype=bool)
    require(array.shape[-1] == MODES, 'coarse reader requires all scalar modes')
    rest[list(SEEDS)] = False
    return np.concatenate((array[..., list(SEEDS)], array[..., rest].sum(axis=-1, keepdims=True)), axis=-1)


def _tv(left, right, axes=(-2, -1)):
    return 0.5 * np.abs(left - right).sum(axis=axes)


def _family_bound(shots, categories, settings, alpha):
    """Union of finite multinomial L1 concentration bounds, divided by two."""
    return math.sqrt((categories * math.log(2) + math.log(settings / alpha)) / (2 * shots))


def _hypotheses(archive):
    p = np.asarray(archive['forecast_probabilities'])
    native = p[COUPLINGS.index(1)]
    return (native, p[COUPLINGS.index(0)], p[COUPLINGS.index(-1)],
            native.sum(axis=-1, keepdims=True) * native.sum(axis=-2, keepdims=True),
            np.asarray(archive['dephased_probabilities'])[COUPLINGS.index(1)])


def _scores(protocol, archive, tv_bound):
    """Selection uses only five stations at two waits; every other fold is sealed."""
    counts = np.asarray(archive['joint_counts'])[COUPLINGS.index(1)]
    forecasts = _hypotheses(archive)
    # Flatten only the menu. Arm and scalar mode remain complete outcomes.
    groups = dict(selection=[(e, t) for e in range(5) for t in (0, 1)],
                  sixth_station=[(5, t) for t in range(3)],
                  late_first_five=[(e, 2) for e in range(5)],
                  marked_reference=[(6, t) for t in range(3)])
    scores = {}
    for group, pairs in groups.items():
        actual = np.stack([counts[:, e, t] for e, t in pairs])
        rows = {}
        for name, forecast_law in zip(MODELS, forecasts):
            predicted = np.stack([forecast_law[:, e, t] for e, t in pairs])
            nll = float(-(actual * np.log(np.maximum(predicted, 1e-300))).sum() / actual.sum())
            tv = float(np.max(_tv(actual / protocol['shots'], predicted)))
            rows[name] = dict(mean_negative_log_likelihood=nll, maximum_total_variation=tv,
                              passes_simultaneous_probability_test=bool(tv <= tv_bound))
        best = min(row['mean_negative_log_likelihood'] for row in rows.values())
        for row in rows.values():
            row['excess_negative_log_likelihood'] = row['mean_negative_log_likelihood'] - best
        scores[group] = dict(settings=int(np.prod(actual.shape[:-2])), models=rows)
    selection = scores['selection']['models']
    best = min(x['mean_negative_log_likelihood'] for x in selection.values())
    winners = [name for name in MODELS if abs(selection[name]['mean_negative_log_likelihood'] - best) <= 1e-12]
    return dict(selection_rule='minimum complete-record negative log likelihood; numeric tie tolerance 1e-12 nats',
                likelihood_floor=1e-300, selected_models=winners,
                selection_accepted_models=[name for name in MODELS
                                           if selection[name]['passes_simultaneous_probability_test']],
                folds=scores)


def analyze(protocol, archive, *, acquisition_weights=None):
    """Validate complete saved records and score the frozen comparisons."""
    validate_protocol(protocol)
    archived = acquisition_weights is not None
    weight_source = authenticated_weights(acquisition_weights) if archived else archive
    shape = (3, 2, 7, 3, 4, 2, MODES)
    p = np.asarray(archive['joint_probabilities'])
    f = np.asarray(archive['forecast_probabilities'])
    n = np.asarray(archive['joint_counts'])
    require(p.shape == f.shape == n.shape == shape, 'joint archive shape differs')
    require(np.all(np.isfinite(p)) and np.min(p) >= 0
            and np.max(abs(p.sum(axis=(-2, -1)) - 1)) < 1e-10, 'joint law is incomplete')
    require(np.issubdtype(n.dtype, np.integer) and np.min(n) >= 0
            and np.all(n.sum(axis=(-2, -1)) == protocol['shots']), 'joint counts are incomplete')
    discrepancy = float(np.max(abs(p - f)))
    require(discrepancy <= protocol['exact_forecast_tolerance'], 'independent forecast disagrees')
    settings = sum(int(np.prod(np.asarray(archive[name + '_counts']).shape[:-axes]))
                   for name, axes in FAMILIES)
    # The largest family has 2 intermediate arms x 2 terminal arms x 150 modes.
    bound = _family_bound(protocol['shots'], 4 * MODES, settings,
                          protocol['joint_familywise_alpha']) + 4 * MODES / (2 * SAMPLING_LATTICE)
    sampling = {}
    acquisition_checks = {}
    for name, axes in FAMILIES:
        probability = np.asarray(archive[name + '_probabilities'])
        count = np.asarray(archive[name + '_counts'])
        summed = tuple(range(-axes, 0))
        require(probability.shape == count.shape and np.all(np.isfinite(probability))
                and np.min(probability) >= 0
                and np.max(abs(probability.sum(axis=summed) - 1)) < 1e-10,
                f'{name} record probabilities are incomplete')
        require(np.issubdtype(count.dtype, np.integer) and np.min(count) >= 0
                and np.all(count.sum(axis=summed) == protocol['shots']),
                f'{name} record counts are incomplete')
        weights = np.asarray(weight_source[name + '_sampling_weights'])
        if archived:
            check = check_acquisition_law(probability, weights, axes)
            check['maximum_count_total_variation_from_weights'] = float(np.max(
                _tv(count / protocol['shots'], weights / SAMPLING_LATTICE, summed)))
            require(check['maximum_count_total_variation_from_weights'] <= _family_bound(
                protocol['shots'], 4 * MODES, settings, protocol['joint_familywise_alpha']),
                f'{name} archived counts fail the simultaneous acquisition-law test')
            acquisition_checks[name] = check
        else:
            require(np.array_equal(weights, canonical_sampling_weights(probability, axes)),
                    f'{name} canonical sampling law differs')
        sampling[name] = float(np.max(_tv(count / protocol['shots'], probability, summed)))
    require(np.array_equal(coarse_grain(n), archive['coarse_counts']), 'coarse counts dropped outcomes')
    require(np.max(abs(coarse_grain(p) - archive['coarse_probabilities'])) < 1e-13,
            'coarse probabilities dropped outcomes')
    control_metrics = {}
    for ci, coupling in enumerate(COUPLINGS):
        actual = p[ci]
        controls = dict(marked_contact_suppressed=p[0], marked_contact_reversed=p[COUPLINGS.index(-coupling)],
                        independent_records_product=actual.sum(axis=-1, keepdims=True)
                        * actual.sum(axis=-2, keepdims=True),
                        dephased_before_read=np.asarray(archive['dephased_probabilities'])[ci])
        control_metrics[str(coupling)] = {
            name: dict(maximum_exact_total_variation=float(np.max(_tv(actual, alternative))),
                       maximum_sampled_total_variation=float(np.max(_tv(n[ci] / protocol['shots'], alternative))),
                       rejected_by_simultaneous_bound=bool(np.max(_tv(n[ci] / protocol['shots'], alternative)) > bound))
            for name, alternative in controls.items()}
    path = np.asarray(archive['path_amplitudes_real']) + 1j * np.asarray(archive['path_amplitudes_imag'])
    path_law = abs(path) ** 2
    require(path.shape == (3, 2, 7, 3, 2, 2, MODES), 'spectral path shape differs')
    path_error = float(np.max(abs(path_law - archive['arm_path_probabilities'])))
    require(path_error <= protocol['exact_forecast_tolerance'], 'retained arm paths disagree')
    dephased = np.repeat(np.asarray(archive['arm_path_probabilities']).sum(axis=-3)[..., None, :, :],
                         4, axis=-3)
    require(np.max(abs(dephased - archive['dephased_probabilities'])) <= protocol['exact_forecast_tolerance'],
            'the dephased law differs from its retained intermediate arm records')
    require(np.max(abs(p.sum(axis=-1) - archive['independent_clock_probabilities']))
            <= protocol['exact_forecast_tolerance'], 'independent clock read differs from its marginal')
    require(np.max(abs(p.sum(axis=-2) - archive['independent_scalar_probabilities']))
            <= protocol['exact_forecast_tolerance'], 'independent scalar read differs from its marginal')
    reconstructed = (p[..., 0, :, :] - p[..., 2, :, :]
                     + 1j * (p[..., 1, :, :] - p[..., 3, :, :])) / 4
    expected = path[..., 0, :, :] * path[..., 1, :, :].conj()
    fringe_error = float(np.max(abs(reconstructed - expected)))
    require(fringe_error <= protocol['exact_forecast_tolerance'], 'four-phase fringe fails its path forecast')
    scores = _scores(protocol, archive, bound)
    sampled_pass = max(sampling.values()) <= bound
    selected = scores['selected_models']
    controls_separated = scores['selection_accepted_models'] == ['native_positive_mark']
    held_pass = all(scores['folds'][group]['models']['native_positive_mark']['passes_simultaneous_probability_test']
                    for group in ('sixth_station', 'late_first_five', 'marked_reference'))
    result = dict(
        source=source_certificate(protocol), joint_record_shape=list(shape),
        axis_order=['coupling', 'collar_sign', 'scalar_seed', 'wait', 'phase', 'terminal_arm', 'scalar_mode'],
        couplings=list(COUPLINGS), collar_signs=list(SIGNS), scalar_seeds=list(SEEDS),
        times=list(protocol['joint_times']), phases=list(PHASES), terminal_arms=['blank_fuel_present', 'written_fuel_absent'],
        coarse_modes=[*SEEDS, 'all_other_modes'], total_joint_settings=int(np.prod(shape[:-2])),
        shots_per_setting=protocol['shots'], normalization_max_error=float(np.max(abs(p.sum(axis=(-2, -1)) - 1))),
        independent_forecast_max_error=discrepancy, retained_arm_path_max_error=path_error,
        four_phase_path_max_error=fringe_error, simultaneous_sampling_total_variation_bound=bound,
        familywise_alpha=protocol['joint_familywise_alpha'], sampled_family_count=settings,
        sampling_probability_lattice=SAMPLING_LATTICE,
        maximum_sampling_canonicalization_total_variation=4 * MODES / (2 * SAMPLING_LATTICE),
        sampling_max_total_variation=sampling, sampling_within_bound=sampled_pass,
        marked_contact_controls=control_metrics, held_out_model_selection=scores,
        result=('finite_joint_record_closure' if sampled_pass and held_pass
                and controls_separated and selected == ['native_positive_mark']
                else 'finite_joint_record_obstruction'),
        classification_scope='finite preparation, interaction and pointer records on the declared 150-mode source',
        arm_control='intermediate arm outcome retained; no heralded path is discarded',
        independent_read_control='fresh ensembles measure only the collar or only the scalar mode',
        supplied_services=['one scalar excitation at a named seam', 'blank coherent collar and finite strength-64 pulses',
                           'phase settings and wait durations', 'complete orthogonal scalar-mode and collar pointer service'],
        empirical_calibration_inputs=[],
    )
    if archived:
        result['original_acquisition'] = dict(weight_digest=ACQUISITION_DIGEST,
                                              source_bin_tolerance=SOURCE_BIN_TOLERANCE,
                                              source_total_variation_tolerance=SOURCE_TV_TOLERANCE,
                                              families=acquisition_checks)
    return result


def run(protocol, *, acquisition=None):
    """Return (summary, arrays); callers archive arrays without JSON expansion."""
    validate_protocol(protocol)
    rows = [[execute(protocol, coupling, sign) for sign in SIGNS] for coupling in COUPLINGS]
    predicted = [[forecast(protocol, coupling, sign) for sign in SIGNS] for coupling in COUPLINGS]
    arrays = {}
    for index, name in enumerate(('joint', 'arm_path', 'pure_arm', 'independent_clock', 'independent_scalar')):
        arrays[name + '_probabilities'] = np.array([[row[index] for row in group] for group in rows])
    arrays['forecast_probabilities'] = np.array([[row[0] for row in group] for group in predicted])
    amplitude = np.array([[row[1] for row in group] for group in predicted])
    arrays['path_amplitudes_real'] = amplitude.real
    arrays['path_amplitudes_imag'] = amplitude.imag
    dephased = arrays['arm_path_probabilities'].sum(axis=-3)
    arrays['dephased_probabilities'] = np.repeat(dephased[..., None, :, :], 4, axis=-3)
    # Independent deterministic streams isolate each measurement family.
    weights = None
    if acquisition is None:
        streams = np.random.SeedSequence([protocol['seed'], 150, 64]).spawn(len(FAMILIES))
        for (name, axes), stream in zip(FAMILIES, streams):
            arrays[name + '_sampling_weights'] = canonical_sampling_weights(arrays[name + '_probabilities'], axes)
            arrays[name + '_counts'] = _sample(arrays[name + '_probabilities'], protocol['shots'],
                                              np.random.default_rng(stream), axes)
    else:
        measured, weights = acquisition
        for name, _ in FAMILIES:
            count = measured[name + '_counts']
            require(count.dtype == np.int64, f'{name} acquired counts require int64')
            arrays[name + '_counts'] = count
    arrays['coarse_probabilities'] = coarse_grain(arrays['joint_probabilities'])
    arrays['coarse_counts'] = coarse_grain(arrays['joint_counts'])
    report = analyze(protocol, arrays, acquisition_weights=weights)
    if weights is not None:
        arrays.update(weights)
    return report, arrays
