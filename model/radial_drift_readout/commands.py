"""Timed integrated commands for the unchanged complete native readout word."""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256

from model.certified_reference import certificate as cert
from model.native_spin_observer import circuits
from model.observer.source_scalar import compiler
from model.radial_geometry_join import source

from .bounds import rational


@lru_cache(maxsize=1)
def round_word():
    quantum, error = source.word_and_bound()
    coin = list(circuits.emit([('X', (0,), cert.PiMultiple(F(1, 4)))]+cert.cnot(0, 1)))
    events = list(quantum)
    for j in range(9):
        events += compiler.remap(coin, [7+2*j, 8+2*j])
    if len(quantum) != 40 or len(coin) != 4 or len(events) != 76:
        raise ArithmeticError('the complete parity and nine-coin word changed')
    if F(error['nine_coin_attempt_trace_distance_upper']) >= F(1, 10**9):
        raise ArithmeticError('complete inherited control precision is insufficient')
    return tuple(events), error


def timing_plan(total_time=F(1, 10**6), rounds=32):
    total_time = rational(total_time, positive=True)
    if type(rounds) is not int or not 1 <= rounds <= 1024:
        raise ValueError('an integer round count from 1 to 1024 is required')
    word, parent = round_word()
    duration = total_time/(2*rounds*len(word))
    interval = total_time/(2*rounds)
    commands, maximum, digest = [], F(), sha256()
    for event in word:
        integrated = tuple(F(float(x)) for x in compiler.integrated_controls(event['kind'], event['angle']))
        rates = tuple(value/duration for value in integrated)
        if any(rate*duration != value for rate, value in zip(rates, integrated, strict=True)):
            raise ArithmeticError('timing changed an integrated source coefficient')
        maximum = max(maximum, *map(abs, rates))
        row = {'kind': event['kind'], 'literal_angle_hex': event['angle'].hex(),
               'local_round_targets': list(event['targets']),
               'duration': str(duration),
               'twelve_exact_integrated_port_coefficients': list(map(str, integrated)),
               'twelve_exact_port_rates': list(map(str, rates))}
        commands.append(row)
        digest.update(repr((event['kind'], event['angle'].hex(), tuple(event['targets']), integrated, rates)).encode())
    if rounds*(len(word)*duration+interval) != total_time:
        raise ArithmeticError('a native pulse or interround time interval was omitted')
    return {
        'rounds': rounds, 'total_elapsed_source_time': str(total_time),
        'positive_pulse_duration': str(duration),
        'positive_nonpulse_time_per_round': str(interval),
        'native_pulses_per_round': len(word), 'native_pulses_total': rounds*len(word),
        'integer_repair_attempts': rounds,
        'fresh_logical_apparatus_bits': rounds*23,
        'complete_round_word_sha256': cert.word_digest(word),
        'timed_command_sha256': digest.hexdigest(),
        'one_round_commands': commands,
        'round_rebase': 'field wires 0,1 fixed; each local apparatus wire j>=2 maps to j+23*round_index',
        'required_absolute_port_rate_upper': str(maximum),
        'one_round_complete_native_half_trace_error_upper': parent['nine_coin_attempt_trace_distance_upper'],
        'native_integrated_word_changed_by_time_assignment': False,
        'same_Hamiltonian_clock_for_pulse_rates_and_source_drift_required': True,
        'nonpulse_budget_must_cover_import_export_pointer_repair_and_controller_latency': True,
        'arbitrary_local_port_drift_extension_admitted': False,
        'local_port_drift_extension_must_preserve_the_stated_norm_and_remote_factor': True,
        'physical_port_bandwidth_or_clock_calibration_derived': False,
        'full_timed_27_site_native_process_executed': False,
    }
