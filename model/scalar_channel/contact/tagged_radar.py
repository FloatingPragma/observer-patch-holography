"""A binary time-bin challenge and local echo with retained missed-wave returns.

The challenge bit controls only the launch time. It is absent from the scalar
carrier and both decoders. The receiver derives its tag from its first click,
parks the scalar, and releases it in that tag's declared return slot. A missed
outbound wave survives into the complete return window in both controls.
"""

import math

import numpy as np

from ..carrier import require
from . import radar


SHOTS = 1048576
NAMESPACE = 826402
ALPHA = 1e-6
ACTIVE, SILENT = 0, 1
TIMEOUT = 2
EVENT_NAMES = (
    'receiver_correct', 'receiver_wrong', 'outbound_timeout', 'receiver_capture', 'b0_capture',
    'delivered_correct', 'delivered_wrong', 'return_tag_correct', 'return_tag_wrong',
    'return_timeout_after_capture', 'false_echo', 'any_A_click', 'A_decodes_source_bit',
    'both_timeouts', 'outbound_late', 'return_late')
# Every decoder cell and declared event in all 120 settings, plus three
# conditional witness frequencies in each of the 60 active settings.
STATEMENTS = 120 * (18 + len(EVENT_NAMES)) + 60 * 3


def _monitor(unitary, state, detector, ticks):
    """Unnormalized first-click branches and the complete no-click survivor."""
    state = np.asarray(state, complex).copy()
    hits = np.zeros(ticks)
    for tick in range(ticks):
        state = unitary @ state
        hits[tick] = abs(state[detector])**2
        state[detector] = 0
    return hits, state


def pair_histories(unitary, source, receiver, ticks):
    """Uniform-input laws with axes control, bit, B-first-tick, A-first-tick.

    Both tick axes have 2*ticks click outcomes followed by timeout. B tick
    indices are relative to global zero, A indices to global 2*ticks. A
    captured scalar is isolated until global 2*ticks+B_decoded*ticks.
    """
    unitary = np.asarray(unitary, complex)
    require(unitary.ndim == 2 and unitary.shape[0] == unitary.shape[1]
            and source != receiver and ticks > 0, 'invalid tagged radar instrument')
    horizon = 2 * ticks
    result = np.zeros((2, 2, horizon + 1, horizon + 1))
    site_source = np.eye(len(unitary), dtype=complex)[:, source]
    site_receiver = np.eye(len(unitary), dtype=complex)[:, receiver]
    returns = []
    for decoded in (0, 1):
        hits, survivor = _monitor(unitary, site_receiver, source, horizon - decoded * ticks)
        law = np.r_[np.zeros(decoded * ticks), hits, np.vdot(survivor, survivor).real]
        returns.append(law)
    for bit in (0, 1):
        hits, missed_wave = _monitor(unitary, site_source, receiver, horizon - bit * ticks)
        outbound = np.r_[np.zeros(bit * ticks), hits]
        for first_tick, weight in enumerate(outbound):
            decoded = first_tick // ticks
            result[ACTIVE, bit, first_tick] = .5 * weight * returns[decoded]
            result[SILENT, bit, first_tick, -1] = .5 * weight
        # No removal, reset, or projection onto |source> occurs on this branch.
        false_returns, final_survivor = _monitor(unitary, missed_wave, source, horizon)
        timeout_law = np.r_[false_returns, np.vdot(final_survivor, final_survivor).real]
        result[:, bit, -1] = .5 * timeout_law
    require(np.max(abs(result.sum(axis=(-2, -1)) - .5)) < 1e-10,
            'tagged histories lose a challenge or timeout branch')
    require(np.array_equal(result[ACTIVE, :, -1], result[SILENT, :, -1]),
            'silent control changes the missed outbound wave')
    return result


def decode_histories(histories, ticks):
    """Both decoders use their first-click slot, independently of the source bit."""
    labels = np.zeros((3, 2 * ticks + 1))
    labels[0, :ticks] = labels[1, ticks:2 * ticks] = labels[2, -1] = 1
    return np.einsum('ij,...bjk,lk->...bil', labels, histories, labels)


def factorized_capture_histories(long_legs, short_legs, source, receiver, ticks):
    """Independent first-hit utility prediction for every captured branch."""
    horizon = 2 * ticks
    result = np.zeros((2, horizon + 1, horizon + 1))
    outgoing = [long_legs[source, receiver, :-1],
                np.r_[np.zeros(ticks), short_legs[source, receiver, :-1]]]
    returning = [long_legs[receiver, source],
                 np.r_[np.zeros(ticks), short_legs[receiver, source]]]
    for bit in (0, 1):
        for tick in range(horizon):
            result[bit, tick] = .5 * outgoing[bit][tick] * returning[tick // ticks]
    return result


def _event_masks():
    b, j, k = np.indices((2, 3, 3))
    captured, clicked = j < 2, k < 2
    correct = captured & clicked & (j == b) & (k == b)
    return dict(receiver_correct=j == b, receiver_wrong=captured & (j != b),
                outbound_timeout=j == TIMEOUT, receiver_capture=captured, b0_capture=(b == 0) & captured,
                delivered_correct=correct, delivered_wrong=captured & clicked & ~correct,
                return_tag_correct=captured & (k == j), return_tag_wrong=captured & clicked & (k != j),
                return_timeout_after_capture=captured & ~clicked, false_echo=~captured & clicked,
                any_A_click=clicked, A_decodes_source_bit=k == b,
                both_timeouts=~captured & ~clicked, outbound_late=(b == 0) & (j == 1),
                return_late=(j == 0) & (k == 1))


def event_values(law):
    return {name: float(np.sum(law[mask])) for name, mask in _event_masks().items()}


def mutual_information(joint):
    """Finite-alphabet information in bits; a zero-mass subensemble has no value."""
    joint = np.asarray(joint, float)
    total = float(joint.sum())
    if total <= 0:
        return None
    joint = joint / total
    product = joint.sum(axis=1)[:, None] * joint.sum(axis=0)[None, :]
    positive = joint > 0
    return max(0., float(np.sum(joint[positive] * np.log2(joint[positive] / product[positive]))))


def information(law):
    """Separate source-bit information from the returned receiver tag beyond that bit."""
    law = np.asarray(law, float)
    law = law / law.sum()
    captured = float(law[:, :2].sum())
    per_trial = 0.
    by_bit = []
    for bit in (0, 1):
        subensemble = law[bit, :2]
        value = mutual_information(subensemble)
        by_bit.append(value)
        if value is not None:
            per_trial += float(subensemble.sum()) * value
    return dict(source_bit_to_B=mutual_information(law.sum(axis=2)),
                source_bit_to_A=mutual_information(law.sum(axis=1)),
                receiver_tag_to_A_given_bit_and_capture=per_trial / captured if captured else None,
                receiver_tag_to_A_given_bit_and_capture_by_bit=by_bit,
                receiver_tag_information_weighted_per_trial=per_trial,
                captured_probability=captured)


def _draw(weights, shots, seed, mark, control, source, receiver):
    generator = np.random.Generator(np.random.PCG64(np.random.SeedSequence(
        [seed, NAMESPACE, mark, control, source, receiver])))
    cumulative = np.cumsum(np.asarray(weights, dtype=np.int64).ravel(), dtype=np.int64)
    require(cumulative[-1] == radar.LATTICE, 'invalid tagged acquisition lattice')
    counts = np.zeros(weights.size, dtype=np.int64)
    for start in range(0, shots, 65536):
        draws = generator.integers(0, radar.LATTICE, size=min(65536, shots - start), dtype=np.int64)
        counts += np.bincount(np.searchsorted(cumulative, draws, side='right'), minlength=weights.size)
    return counts.reshape(weights.shape)


def _unconditional_radius(shots, outcomes):
    return math.sqrt(math.log(2 * STATEMENTS / ALPHA) / (2 * shots)) + outcomes / (2 * radar.LATTICE) + 1e-10


def information_lower_bound(counts, weights, outcomes):
    """A simultaneous Pinsker witness for a tag beyond A's known input bit.

    Within b=0 and B-capture, coarse grain the record to X=(B_decoded=0)
    and Y=(A_decoded=0). If covariance has magnitude at least c, binary
    total variation from independence is at least 2c and I(X;Y) >= 8c^2/ln2.
    The conditioning is a joint-record analysis, without an authentication
    flag supplied to A during the experiment.
    """
    shots = int(counts.sum())
    sub = counts[0, :2]
    retained = int(sub.sum())
    result = dict(condition='b=0 and B clicked', retained_trials=retained,
                  covariance_interval=None, mutual_information_lower_bits=0.,
                  capture_weighted_information_lower_bits_per_trial=0.,
                  information_given_bit_and_capture_lower_bits=0.)
    if retained == 0:
        return result
    sub_mass = float(weights[0, :2].sum() / radar.LATTICE)
    require(sub_mass > 0, 'retained tag witness has zero acquisition weight')
    lattice_tv = outcomes / (2 * radar.LATTICE) + 1e-10
    conditional_bias = min(1., 2 * lattice_tv / sub_mass)
    radius = math.sqrt(math.log(2 * STATEMENTS / ALPHA) / (2 * retained)) + conditional_bias
    joint = float(sub[0, 0] / retained)
    x = float(sub[0].sum() / retained)
    y = float(sub[:, 0].sum() / retained)
    lo = max(0., joint - radius) - min(1., x + radius) * min(1., y + radius)
    hi = min(1., joint + radius) - max(0., x - radius) * max(0., y - radius)
    covariance = max(0., lo, -hi)
    bound = 8 * covariance**2 / math.log(2)
    unconditional_radius = _unconditional_radius(shots, outcomes)
    weight_lower = max(0., retained / shots - unconditional_radius)
    capture_upper = min(1., float(counts[:, :2].sum() / shots) + unconditional_radius)
    result.update(covariance_interval=[lo, hi], conditional_frequency_radius=radius,
                  mutual_information_lower_bits=bound,
                  capture_weighted_information_lower_bits_per_trial=weight_lower * bound,
                  information_given_bit_and_capture_lower_bits=weight_lower * bound / capture_upper
                  if capture_upper else 0.)
    return result


def _event_summary(probabilities, counts, radius):
    exact, observed = event_values(probabilities), event_values(counts)
    shots = int(counts.sum())
    return {name: dict(probability=exact[name], count=int(observed[name]),
                       confidence_interval=[max(0., observed[name] / shots - radius),
                                            min(1., observed[name] / shots + radius)])
            for name in EVENT_NAMES}


def _resources(histories, counts, ticks, dt):
    """Occupation tests, monitored propagation and the fixed total action budget."""
    horizon = 2 * ticks
    b, bt, at = np.indices((2, horizon + 1, horizon + 1))
    b_duration = np.minimum(bt + 1, horizon) - b * ticks
    a_duration = np.minimum(at + 1, horizon)
    # Pointer tests can run on vacuum before release and during latch parking.
    b_tests = np.minimum(bt + 1, horizon)
    a_tests = a_duration
    captured = bt < horizon
    return_source_ticks = np.where(captured, np.maximum(0, a_duration - (bt // ticks) * ticks), a_duration)
    result = []
    for control in (ACTIVE, SILENT):
        source_ticks = b_duration + (return_source_ticks if control == ACTIVE
                                     else np.where(captured, 0, a_duration))
        c = counts[control]
        p = histories[control]
        result.append(dict(acquired_B_occupation_tests=int(np.sum(c * b_tests)),
                           acquired_A_occupation_tests=int(np.sum(c * a_tests)),
                           expected_propagation_ticks_until_A_record=float(np.sum(p * source_ticks)),
                           acquired_propagation_action_until_A_record=float(dt * np.sum(c * source_ticks)),
                           elapsed_scheduler_ticks=int(2 * horizon * c.sum()),
                           elapsed_source_action_time=float(2 * horizon * dt * c.sum()),
                           captured_scalar_count=int(c[:, :-1].sum()),
                           released_echo_count=int(c[:, :-1].sum()) if control == ACTIVE else 0))
    return result


def run(protocol):
    """Execute the frozen acquisition and return a summary plus numerical receipt arrays."""
    radar._settings(protocol)
    require(protocol['shots'] == SHOTS, 'tagged radar requires exactly 1048576 trials per setting')
    ticks, dt = protocol['ticks'], protocol['dt']
    horizon = 2 * ticks
    shape = (2, 2, 6, 6, 2, horizon + 1, horizon + 1)
    laws = np.zeros(shape)
    counts = np.zeros(shape, dtype=np.int64)
    weights = np.zeros(shape, dtype=np.int64)
    dyads = []
    outcomes = 2 * (horizon + 1)**2
    radius = _unconditional_radius(SHOTS, outcomes)
    maximum_factorization_error = 0.
    for mark in (0, 1):
        unitary = radar._frozen_unitary(protocol, mark)
        long_legs = radar.first_hit(unitary, radar.ENDPOINTS, horizon)
        short_legs = radar.first_hit(unitary, radar.ENDPOINTS, ticks)
        for source, source_site in enumerate(radar.ENDPOINTS):
            for receiver, receiver_site in enumerate(radar.ENDPOINTS):
                if source == receiver:
                    continue
                p = pair_histories(unitary, source_site, receiver_site, ticks)
                factorized = factorized_capture_histories(long_legs, short_legs, source, receiver, ticks)
                error = float(np.max(abs(p[ACTIVE, :, :-1] - factorized[:, :-1])))
                require(error < 1e-10, 'captured records differ from independent leg factorization')
                maximum_factorization_error = max(maximum_factorization_error, error)
                laws[mark, :, source, receiver] = p
                for control in (ACTIVE, SILENT):
                    w = radar.acquisition_weights(p[control])
                    weights[mark, control, source, receiver] = w
                    counts[mark, control, source, receiver] = _draw(
                        w, SHOTS, protocol['seed'], mark, control, source, receiver)
                c = counts[mark, :, source, receiver]
                decoder_p = decode_histories(p, ticks)
                decoder_c = decode_histories(c, ticks).astype(np.int64)
                decoder_w = decode_histories(weights[mark, :, source, receiver], ticks).astype(np.int64)
                event_summaries = [_event_summary(decoder_p[k], decoder_c[k], radius) for k in (ACTIVE, SILENT)]
                contrast = {}
                for name in ('any_A_click', 'A_decodes_source_bit', 'false_echo', 'delivered_correct'):
                    point = (event_summaries[ACTIVE][name]['count'] - event_summaries[SILENT][name]['count']) / SHOTS
                    contrast[name] = dict(
                        probability_difference=event_summaries[ACTIVE][name]['probability']
                        - event_summaries[SILENT][name]['probability'],
                        observed_difference=point,
                        confidence_interval=[max(-1., point - 2 * radius), min(1., point + 2 * radius)])
                dyads.append(dict(mark_written=mark, source=source, receiver=receiver,
                                  active_events=event_summaries[ACTIVE], silent_events=event_summaries[SILENT],
                                  active_information=information(decoder_p[ACTIVE]),
                                  silent_information=information(decoder_p[SILENT]),
                                  active_information_from_counts=information(decoder_c[ACTIVE]),
                                  silent_information_from_counts=information(decoder_c[SILENT]),
                                  active_tag_information_lower_bound=information_lower_bound(
                                      decoder_c[ACTIVE], decoder_w[ACTIVE], outcomes),
                                  active_minus_silent=contrast,
                                  resources_by_control=_resources(p, c, ticks, dt)))
    decoded_laws = decode_histories(laws, ticks)
    decoded_counts = decode_histories(counts, ticks).astype(np.int64)
    arrays = dict(tagged_time_probabilities=laws, tagged_time_counts=counts,
                  tagged_time_acquisition_weights=weights,
                  tagged_decoder_probabilities=decoded_laws, tagged_decoder_counts=decoded_counts)
    summary = dict(
        protocol='binary local time-bin challenge and fixed-window echo',
        mark_sectors=[0, 1], controls=['active_echo', 'silent_capture'],
        endpoints=list(radar.ENDPOINTS), ticks_per_slot=ticks, interval=dt, shots_per_control_and_dyad=SHOTS,
        challenge_bit_probability=[.5, .5], launch_global_ticks=[0, ticks],
        B_test_global_ticks=[1, horizon], echo_release_global_ticks=[horizon, horizon + ticks],
        A_test_global_ticks=[horizon + 1, 2 * horizon],
        decoders='first-click slot 0 or 1; timeout label 2; neither decoder reads the source bit',
        scalar_tag='only launch time and local click-triggered release time; no bit-labelled scalar degree of freedom',
        time_axes=['mark', 'control', 'source', 'receiver', 'bit', 'B_first_tick', 'A_first_tick'],
        time_timeout_index=horizon,
        tick_index_rule='B click index+1 is a global tick; A click index+1 follows global 2N',
        decoder_axes=['mark', 'control', 'source', 'receiver', 'bit', 'B_decoder', 'A_decoder'],
        diagonal_menu='unused zero sentinel',
        independent_capture_factorization_maximum_error=maximum_factorization_error,
        simultaneous_confidence=dict(family_failure_probability=ALPHA, statements=STATEMENTS,
                                     unconditional_frequency_radius=radius,
                                     method='Hoeffding and Bonferroni with integer-acquisition bias; conditional tag witnesses use actual retained sample sizes'),
        rng=dict(bit_generator='PCG64', namespace=NAMESPACE, lattice=radar.LATTICE,
                 entropy_order=['seed', 'namespace', 'mark', 'control', 'source', 'receiver']),
        finite_memory=dict(active_scalar_modes=150, isolated_local_latch_modes=2,
                           one_excitation_source_plus_latches_dimension=152,
                           inert_occupied_parity_reference_modes=1,
                           challenge_record_bits=1, first_click_record_outcomes_per_endpoint=horizon + 1,
                           first_click_record_bits_per_endpoint=math.ceil(math.log2(horizon + 1)),
                           total_raw_record_bits=1 + 2 * math.ceil(math.log2(horizon + 1)),
                           local_schedule_counters=2, counter_bits_each=math.ceil(math.log2(2 * horizon + 1)),
                           decoder_bits_are_functions_of_first_click_records=True),
        supplied_services=[
            'uniform classical bit from a local binary preparation and readout at A',
            'coherent site release from the isolated source latch at the bit-selected absolute tick',
            'projective local occupation readers with orthogonal classical first-click records',
            'local capture into an isolated B latch on a hit and bit-slot-controlled release',
            'absolute deadlines, synchronized cadence and instantaneous local switches',
            'trial reset only after the complete A window'],
        source_action='the contact source stays on through the full 4N-tick schedule, including scalar-vacuum latch intervals',
        information_scope=[
            'source-bit to A information includes direct late-wave leakage measured by the silent control',
            'receiver-tag information conditions on the known source bit and a receiver capture',
            'capture conditioning uses shared records after acquisition and is not an authentication flag delivered to A',
            'an active-minus-silent click contrast tests the remote reply intervention while retaining the same missed wave',
            'mutual information differences are not used as an additive delivered-information measure',
            'the protocol has no asserted zero-error code, physical speed or autonomous clock implementation'],
        dyads=dyads)
    return summary, arrays
