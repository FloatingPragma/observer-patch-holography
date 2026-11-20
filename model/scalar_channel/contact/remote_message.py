"""Repeated private receiver messages through the fixed active/silent echo switch.

Each bit uses a whole finite block. The reader receives two known calibration
click totals and eight payload click totals. The receiver sign, launch bits,
capture records and timed joint histories are absent from that interface.
"""

import math

import numpy as np

from ..carrier import require
from . import radar


TRIALS = 2**42
PAYLOAD_BITS = 8
BLOCKS = 2 + PAYLOAD_BITS
NAMESPACE = 826403
ALPHA = 1e-6
ARCHIVE_ALPHA = 1e-6
ARCHIVE_CATEGORIES = 60 * BLOCKS * 2 * 25 * 25
SAMPLER_PROBABILITY_ROUNDOFF = 1e-12
CALIBRATION_STATEMENTS = 2 * 60
PAYLOAD_STATEMENTS = PAYLOAD_BITS * 60
OUTCOMES = 2 * 25 * 25
LATTICE_TV_BOUND = OUTCOMES / (2 * radar.LATTICE) + 1e-10
SOURCE_PROBABILITY_ROUNDOFF = 1e-10
CONTROL_FOR_BIT = (1, 0)  # The frozen tag array orders active before silent.
ABSTAIN = -1


def sign_measurement_certificate():
    """The supplied written d=+1 state gives a fair native sign readout at A_s=1."""
    from ...polynomial_sv import law as contrast
    from .law import collar_states

    exchange, square = contrast.local_term_parts()
    signs = [collar_states(sign)[:, 1] for sign in (-1, 1)]
    written = np.zeros(contrast.DIM)
    label = contrast.LABELS.index((0, 1))
    written[contrast.N_SRC * (1 + label) + contrast.SRC.index((0, 0))] = 1
    probabilities = [float(abs(np.vdot(arm, written))**2) for arm in signs]
    error = max(float(np.linalg.norm(square @ vector - vector)) for vector in [written] + signs)
    energies = [float(arm @ exchange @ arm) for arm in signs]
    exchange_error = max(float(np.linalg.norm(exchange @ arm - energy * arm))
                         for arm, energy in zip(signs, energies))
    require(max(abs(p - .5) for p in probabilities) < 1e-13 and error < 1e-13
            and exchange_error < 1e-13, 'native receiver sign does not preserve the written contact')
    return dict(preparation='supplied empty-load written memory |0,+1>',
                outcomes=['antisymmetric sign, bit 0, silent', 'symmetric sign, bit 1, active'],
                probabilities=probabilities, contrast_square_eigenvalue=1,
                contact_preservation_residual=error, exchange_eigenvalues=energies,
                exchange_eigenstate_residual=exchange_error,
                source_effect='both fixed sign sectors retain the same +1 endpoint potential; their collar energy is a common phase',
                service='ideal local preparation, sign measurement, and classical register controlling the supplied release switch')


def allowances(trials):
    """Separate simultaneous calibration and future-block errors, including lattice bias."""
    require(type(trials) is int and 0 < trials <= np.iinfo(np.int64).max, 'invalid block trial count')
    calibration = math.sqrt(math.log(2 * CALIBRATION_STATEMENTS / (ALPHA / 2)) / (2 * trials))
    payload = math.sqrt(math.log(2 * PAYLOAD_STATEMENTS / (ALPHA / 2)) / (2 * trials))
    return dict(calibration_sampling=calibration, payload_sampling=payload,
                acquisition_lattice_bias=LATTICE_TV_BOUND,
                calibration=calibration + LATTICE_TV_BOUND,
                payload=payload + LATTICE_TV_BOUND)


def read_message(calibration_A_clicks, payload_A_clicks, trials):
    """The entire A-reader interface; calibration positions have known labels 0 and 1."""
    require(len(calibration_A_clicks) == 2 and len(payload_A_clicks) == PAYLOAD_BITS,
            'reader requires two labelled calibrations and eight payload totals')
    all_counts = list(calibration_A_clicks) + list(payload_A_clicks)
    require(all(type(x) in (int, np.int64) and 0 <= x <= trials for x in all_counts),
            'reader counts must be bounded integers')
    budget = allowances(trials)
    calibration = [int(x) for x in calibration_A_clicks]
    frequencies = [x / trials for x in calibration]
    calibration_intervals = [[max(0., x - budget['calibration']), min(1., x + budget['calibration'])]
                             for x in frequencies]
    future_radius = budget['calibration'] + budget['payload']
    predictive_intervals = [[max(0., x - future_radius), min(1., x + future_radius)] for x in frequencies]
    difference = calibration[1] - calibration[0]
    orientation = 1 if difference > 0 else -1 if difference < 0 else 0
    margin = abs(difference) / trials - 2 * future_radius
    separated = orientation != 0 and margin > 0
    numerator = calibration[0] + calibration[1]
    decoded = []
    for count in payload_A_clicks:
        # Integer arithmetic preserves the midpoint for counts above 2**32.
        signed_side = orientation * (2 * int(count) - numerator)
        decoded.append(ABSTAIN if not separated or signed_side == 0 else int(signed_side > 0))
    return dict(calibration_A_clicks=calibration, payload_A_clicks=[int(x) for x in payload_A_clicks],
                calibration_probability_intervals=calibration_intervals,
                future_frequency_intervals=predictive_intervals,
                confidence_separated=separated, separation_margin=margin,
                bit_one_has_higher_click_rate=orientation > 0 if separated else None,
                threshold_click_numerator=numerator if separated else None,
                threshold_click_denominator=2 if separated else None,
                decoded_bits=decoded, abstentions=decoded.count(ABSTAIN))


def _generator(seed, role, mark, source, receiver, block=0):
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(
        [seed, NAMESPACE, role, mark, source, receiver, block])))


def _private_payload(seed, mark, source, receiver):
    """Reproducible outcomes of eight fresh fair sign reads, isolated from the reader."""
    return _generator(seed, 2, mark, source, receiver).integers(0, 2, size=PAYLOAD_BITS, dtype=np.int64)


def multinomial_block(weights, trials, seed, role, mark, source, receiver, block):
    """A complete large-N histogram; exact zero-weight bins receive no counts."""
    weights = np.asarray(weights)
    require(weights.dtype.kind in 'iu' and np.all(weights >= 0)
            and int(weights.sum()) == radar.LATTICE, 'invalid message acquisition weights')
    support = np.flatnonzero(weights.ravel())
    # Put the largest bin last so NumPy's residual probability cannot put
    # floating-point residue into an impossible outcome.
    largest = int(support[np.argmax(weights.ravel()[support])])
    support = np.r_[support[support != largest], largest]
    probabilities = weights.ravel()[support] / radar.LATTICE
    sampled = _generator(seed, role, mark, source, receiver, block).multinomial(int(trials), probabilities)
    result = np.zeros(weights.size, dtype=np.int64)
    result[support] = sampled.astype(np.int64)
    require(int(result.sum()) == trials, 'large-count message histogram loses trials')
    return result.reshape(weights.shape)


def A_click_total(histogram):
    """Discard the launch bit and every receiver outcome; retain only A click/timeout."""
    histogram = np.asarray(histogram)
    require(histogram.shape == (2, 25, 25), 'invalid tagged block record')
    return int(histogram[:, :, :-1].sum())


def sufficient_block_trials(probability_gap):
    """A fixed-design sufficient budget, without changing the acquired block size.

    A true gap larger than 4*e_cal+2*e_future makes the calibration gate pass
    throughout the simultaneous event. The same event then gives correct
    midpoint decoding for every payload block.
    """
    available = probability_gap - 6 * LATTICE_TV_BOUND
    if available <= 0:
        return None
    calibration_constant = math.sqrt(math.log(2 * CALIBRATION_STATEMENTS / (ALPHA / 2)) / 2)
    payload_constant = math.sqrt(math.log(2 * PAYLOAD_STATEMENTS / (ALPHA / 2)) / 2)
    # Strict separation requires an integer larger than the squared threshold.
    return math.floor(((4 * calibration_constant + 2 * payload_constant) / available)**2) + 1


def _validate_laws(laws):
    laws = np.asarray(laws, float)
    require(laws.shape == (2, 2, 6, 6, 2, 25, 25)
            and np.isfinite(laws).all() and np.all(laws >= 0), 'invalid frozen tagged source laws')
    for mark in (0, 1):
        for source in range(6):
            for receiver in range(6):
                menu = laws[mark, :, source, receiver]
                if source == receiver:
                    require(np.all(menu == 0), 'self-dyad must be an unused zero sentinel')
                else:
                    require(np.max(abs(menu.sum(axis=(-2, -1)) - .5)) < 1e-10,
                            'each source launch bit must retain probability one half')
                    require(np.max(abs(menu[0, :, -1] - menu[1, :, -1])) < 1e-12,
                            'active and silent controls change the missed outbound branch')
                    require(np.max(menu[1, :, :-1, :-1]) < 1e-12,
                            'silent receiver releases a captured scalar')
    return laws


def _validate_saved_records(records):
    records = np.asarray(records)
    require(records.shape == (2, 6, 6, BLOCKS, 2, 25, 25) and records.dtype == np.dtype(np.int64),
            'saved message records require the complete int64 histogram shape')
    require(np.all(records >= 0) and np.all(records <= TRIALS), 'saved message count is outside block bounds')
    for source in range(6):
        require(np.all(records[:, source, source] == 0), 'saved self-dyad must be zero')
        for receiver in range(6):
            if source != receiver:
                require(np.all(records[:, source, receiver].sum(axis=(-3, -2, -1)) == TRIALS),
                        'saved message block has the wrong trial total')
    return records


def _validate_frozen_weights(laws, frozen_weights):
    """Verify retained acquisition masses without requantizing a roundoff-shifted law.

    An ordinary bin incurs at most half a lattice unit of rounding. The bin
    receiving the normalization remainder has the full lattice allowance.
    A roundoff-level tie can change which maximum receives that remainder.
    """
    frozen_weights = np.asarray(frozen_weights)
    require(frozen_weights.shape == laws.shape and frozen_weights.dtype.kind in 'iu',
            'frozen message weights require the complete integer shape')
    require(np.all(frozen_weights >= 0) and np.all(frozen_weights <= radar.LATTICE),
            'frozen message weight is outside lattice bounds')
    retained = frozen_weights.astype(np.int64, copy=True)
    maximum_tv, maximum_bin_error = 0., 0.
    ordinary_bound = .5 / radar.LATTICE + SOURCE_PROBABILITY_ROUNDOFF
    residual_bound = OUTCOMES / (2 * radar.LATTICE) + SOURCE_PROBABILITY_ROUNDOFF
    for source in range(6):
        require(np.all(retained[:, :, source, source] == 0), 'frozen self-dyad weights must be zero')
        for receiver in range(6):
            if source == receiver:
                continue
            for mark in (0, 1):
                for control in (0, 1):
                    weights = retained[mark, control, source, receiver]
                    law = laws[mark, control, source, receiver]
                    require(int(weights.sum()) == radar.LATTICE, 'frozen message weights have the wrong menu total')
                    require(np.all(weights[law == 0] == 0), 'frozen message weights populate an exact physical zero')
                    error = abs(weights / radar.LATTICE - law)
                    tv = float(error.sum() / 2)
                    require(tv <= LATTICE_TV_BOUND, 'frozen message weights exceed the physical-law total-variation bound')
                    near_maximum = law >= float(law.max()) - 2 * SOURCE_PROBABILITY_ROUNDOFF
                    bin_bounds = np.where(near_maximum, residual_bound, ordinary_bound)
                    require(np.all(error <= bin_bounds), 'frozen message weights exceed a physical-law bin bound')
                    maximum_tv = max(maximum_tv, tv)
                    maximum_bin_error = max(maximum_bin_error, float(error.max()))
    return retained, dict(mode='verified frozen acquisition masses',
                          maximum_physical_law_total_variation=maximum_tv,
                          maximum_physical_law_bin_error=maximum_bin_error,
                          total_variation_bound=LATTICE_TV_BOUND,
                          ordinary_bin_bound=ordinary_bound, residual_bin_bound=residual_bound,
                          source_probability_roundoff=SOURCE_PROBABILITY_ROUNDOFF,
                          exact_physical_zero_bins_preserved=True)


def validate_block_record(histogram, weights):
    """Simultaneous Bernstein compatibility with every canonical acquisition category."""
    histogram, weights = np.asarray(histogram), np.asarray(weights)
    require(histogram.shape == weights.shape == (2, 25, 25)
            and histogram.dtype == weights.dtype == np.dtype(np.int64), 'invalid message block dtype or shape')
    require(np.all(histogram >= 0) and np.all(histogram <= TRIALS) and int(histogram.sum()) == TRIALS,
            'invalid message block bounds or trial total')
    require(np.all(weights >= 0) and int(weights.sum()) == radar.LATTICE,
            'invalid canonical message weights')
    require(np.all(histogram[weights == 0] == 0), 'saved message populates a forbidden acquisition bin')
    probability = weights / radar.LATTICE
    logarithm = math.log(2 * ARCHIVE_CATEGORIES / ARCHIVE_ALPHA)
    allowance = (np.sqrt(2 * probability * (1 - probability) * logarithm / TRIALS)
                 + 2 * logarithm / (3 * TRIALS) + SAMPLER_PROBABILITY_ROUNDOFF)
    discrepancy = abs(histogram / TRIALS - probability)
    require(np.all(discrepancy <= allowance), 'saved message fails the simultaneous Bernstein category check')
    return float(np.max(discrepancy / allowance))


def run(spec, tag_laws, saved_records=None, frozen_weights=None):
    """Acquire fixed calibration/payload blocks; grade only after the restricted readout."""
    radar._settings(spec)
    laws = _validate_laws(tag_laws)
    certificate = sign_measurement_certificate()
    if saved_records is not None:
        saved_records = _validate_saved_records(saved_records)
    shape = (2, 6, 6, BLOCKS, 2, 25, 25)
    records = np.zeros(shape, dtype=np.int64)
    if frozen_weights is None:
        weights = np.zeros(laws.shape, dtype=np.int64)
        weight_validation = dict(mode='quantized recomputed physical law')
    else:
        weights, weight_validation = _validate_frozen_weights(laws, frozen_weights)
    private_bits = np.full((2, 6, 6, PAYLOAD_BITS), ABSTAIN, dtype=np.int64)
    decoded_bits = np.full_like(private_bits, ABSTAIN)
    reader_counts = np.zeros((2, 6, 6, BLOCKS), dtype=np.int64)
    dyads = []
    largest_archive_fraction = 0.
    for mark in (0, 1):
        for source in range(6):
            for receiver in range(6):
                if source == receiver:
                    continue
                secret = _private_payload(spec['seed'], mark, source, receiver)
                private_bits[mark, source, receiver] = secret
                if frozen_weights is None:
                    for control in (0, 1):
                        weights[mark, control, source, receiver] = radar.acquisition_weights(
                            laws[mark, control, source, receiver])
                labels = np.r_[np.array([0, 1], dtype=np.int64), secret]
                for block, label in enumerate(labels):
                    control = CONTROL_FOR_BIT[int(label)]
                    role = 0 if block < 2 else 1
                    if saved_records is None:
                        records[mark, source, receiver, block] = multinomial_block(
                            weights[mark, control, source, receiver], TRIALS, spec['seed'], role,
                            mark, source, receiver, block if block < 2 else block - 2)
                    else:
                        retained = saved_records[mark, source, receiver, block]
                        largest_archive_fraction = max(largest_archive_fraction, validate_block_record(
                            retained, weights[mark, control, source, receiver]))
                        records[mark, source, receiver, block] = retained
                    reader_counts[mark, source, receiver, block] = A_click_total(
                        records[mark, source, receiver, block])
                local = reader_counts[mark, source, receiver]
                reading = read_message(local[:2].tolist(), local[2:].tolist(), TRIALS)
                decoded = np.asarray(reading['decoded_bits'], dtype=np.int64)
                decoded_bits[mark, source, receiver] = decoded
                # Secret comparison occurs after the reader returns its immutable decisions.
                errors = int(np.count_nonzero((decoded != ABSTAIN) & (decoded != secret)))
                correct = int(np.count_nonzero(decoded == secret))
                probabilities = [float(laws[mark, control, source, receiver, :, :, :-1].sum())
                                 for control in CONTROL_FOR_BIT]
                gap = abs(probabilities[1] - probabilities[0])
                dyads.append(dict(mark=mark, source=source, receiver=receiver,
                                  source_A_click_probabilities_by_bit=probabilities,
                                  source_probability_gap=gap,
                                  sufficient_trials_per_block=sufficient_block_trials(gap),
                                  reader=reading, private_payload_for_grading=secret.tolist(),
                                  correct_bits=correct, incorrect_bits=errors,
                                  abstentions=int(np.count_nonzero(decoded == ABSTAIN))))
    weakest = min(dyads, key=lambda x: x['source_probability_gap'])
    least_margin = min(dyads, key=lambda x: x['reader']['separation_margin'])
    arrays = dict(remote_message_time_counts=records, remote_message_acquisition_weights=weights,
                  remote_message_A_click_counts=reader_counts,
                  remote_message_private_bits=private_bits, remote_message_decoded_bits=decoded_bits)
    summary = dict(
        experiment='fixed-block receiver-private binary messages through the frozen echo release switch',
        classification='high-resource repeated statistical communication',
        trials_per_calibration_or_bit=TRIALS, calibration_labels=[0, 1],
        calibration_controls=['silent', 'active'], fresh_bits_per_dyad=PAYLOAD_BITS,
        source_launch_bit='independent uniform bit each trial; A discards it',
        receiver_bit_source=certificate,
        reader_interface='two labelled A-click calibration totals, eight A-click payload totals, and the common trial count',
        simultaneous_confidence=dict(family_failure_probability=ALPHA,
                                     calibration_failure_budget=ALPHA / 2, payload_failure_budget=ALPHA / 2,
                                     calibration_statements=CALIBRATION_STATEMENTS,
                                     payload_statements=PAYLOAD_STATEMENTS, allowances=allowances(TRIALS),
                                     scope='all stations, both marks, both calibration labels, and all eight held-out bits',
                                     guarantee='with probability at least 1-alpha every non-abstaining bit is correct; insufficient calibration separation causes abstention'),
        decision='midpoint threshold only after calibration intervals enlarged by future-block allowances are disjoint; its orientation follows the known labels',
        no_tuning='block sizes, slot times, controls, confidence allocation and payload length are fixed before any block is acquired',
        record_validation=dict(mode='seeded acquisition' if saved_records is None else 'retained histogram source compatibility',
                               seeded_exact_replay=saved_records is None,
                               archive_family_failure_probability=ARCHIVE_ALPHA,
                               archive_category_count=ARCHIVE_CATEGORIES,
                               archive_probability_roundoff=SAMPLER_PROBABILITY_ROUNDOFF,
                               largest_bernstein_allowance_fraction=largest_archive_fraction if saved_records is not None else None,
                               scope='category compatibility does not establish the acquisition provenance; decoding uses the retained counts directly'),
        acquisition_weight_validation=weight_validation,
        rng=dict(bit_generator='PCG64', namespace=NAMESPACE, roles={'calibration': 0, 'payload': 1, 'private_sign': 2},
                 entropy_order=['seed', 'namespace', 'role', 'mark', 'source', 'receiver', 'block'],
                 aggregate_sampler='multinomial on positive integer-lattice bins, with largest bin last',
                 reader_receives_seed=False,
                 privacy_scope='private means absent from the reader interface; reproducible pseudorandom sampling makes no cryptographic secrecy claim'),
        resources=dict(trials_per_payload_bit=TRIALS, calibration_trials_per_dyad=2 * TRIALS,
                       payload_trials_per_dyad=PAYLOAD_BITS * TRIALS,
                       total_trials_all_dyads=60 * BLOCKS * TRIALS,
                       elapsed_ticks_per_trial=4 * spec['ticks'],
                       maximum_local_occupation_tests_per_bit=4 * spec['ticks'] * TRIALS,
                       source_action_per_payload_bit=4 * spec['ticks'] * spec['dt'] * TRIALS,
                       receiver_sign_reads=60 * PAYLOAD_BITS,
                       held_receiver_bits_during_each_block=1,
                       timed_outcomes_per_block=OUTCOMES, histogram_counter_bits=43,
                       stored_counter_dtype='int64', raw_first_click_record_bits_per_trial=11,
                       A_count_register_bits=43, active_scalar_modes=150, isolated_local_latch_modes=2),
        weakest_source_dyad=dict(mark=weakest['mark'], source=weakest['source'], receiver=weakest['receiver'],
                                 probability_gap=weakest['source_probability_gap'],
                                 sufficient_trials_per_block=weakest['sufficient_trials_per_block'],
                                 calibration_passed=weakest['reader']['confidence_separated'],
                                 all_bits_correct=weakest['correct_bits'] == PAYLOAD_BITS),
        smallest_observed_calibration_margin=dict(mark=least_margin['mark'], source=least_margin['source'],
                                                  receiver=least_margin['receiver'],
                                                  margin=least_margin['reader']['separation_margin']),
        every_dyad_calibration_passed=all(x['reader']['confidence_separated'] for x in dyads),
        all_payload_bits_correct=all(x['correct_bits'] == PAYLOAD_BITS for x in dyads),
        total_correct_bits=sum(x['correct_bits'] for x in dyads),
        total_incorrect_bits=sum(x['incorrect_bits'] for x in dyads),
        total_abstentions=sum(x['abstentions'] for x in dyads),
        record_axes=['mark', 'source', 'receiver', 'block', 'launch_bit', 'B_first_tick', 'A_first_tick'],
        block_order='silent calibration, active calibration, eight held-out private bits',
        timeout_index=24, diagonal_menu='unused zero count sentinel and private-bit sentinel -1',
        boundaries=[
            'individual A clicks remain ambiguous between replies and surviving missed waves',
            'the reader sees neither receiver captures nor private sign values nor launch-bit records',
            'joint histories and private values serve the independent grader after the readout',
            'supplied local sign readout, isolated latches, absolute timing and switching are not autonomous apparatus proofs',
            'this bounded statistical channel supplies no spacetime or radar-click certification'],
        dyads=dyads)
    return summary, arrays
