"""Complete finite CHSH score from declared coherent source services.

The known A0 B1 input and coherent load/export interface are explicit.
No spacelike separation or autonomous continuum dynamics is claimed.
"""
from __future__ import annotations

from fractions import Fraction as F
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BELL = ROOT / 'model/native_coherent_occupation_control'
from model.native_coherent_occupation_control import bell, derived_control as beam
from model.native_effect_synthesis import intervals as iv
from model.native_program import interpreter as native


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def run():
    graph, edge, spare = bell.layout(0)
    masks = bell.code_masks(edge, spare)
    theta = F(25, 32)
    carrier = {'kind': 'federation', 'level': 0, 'total': 2}
    targets = tuple((endpoint % 12, other % 12) for endpoint, other in zip(edge, spare, strict=True))
    cells = tuple(endpoint // 12 for endpoint in edge)
    rail_edges = tuple(tuple(sorted((endpoint, other))) for endpoint, other in zip(edge, spare, strict=True))
    require(all(local_edge in graph.seams for local_edge in rail_edges),
            'the readout rails must be literal native seams')
    require(all(endpoint > other for endpoint, other in zip(edge, spare, strict=True)),
            'the actual seam order must measure the declared logical Z sign')

    initial = {(masks[1], 0, masks[1], 0): F(1)}
    first = beam.lift(beam.declare(-theta, targets[0], cells[0]), initial, carrier)
    second = beam.lift(beam.declare(theta, targets[1], cells[1]), first['success'], carrier)
    herald_program = native.program(federation={'level': 0, 'total': 2}, biased=False,
        commands=[{'op': 'repair', 'selector': {'kind': 'fixed', 'seam': list(edge)}, 'law': 'E', 'repeat': 1}],
        capacity=1, overflow='export_oldest')
    herald = native.execute(herald_program, second['success'])
    require(herald.status == 'COMPLETE' and len(herald.branches) == 25,
            'the complete native herald must retain all 25 raw slots')
    accepted = [row for row in herald.branches if row.record.law == 'E' and row.record.raw == 0]
    aborted = [row for row in herald.branches if row not in accepted]
    require(len(accepted) == 2 and len(aborted) == 23,
            'all accepted and aborted herald records must remain distinct')
    abort_mass = sum((native.trace(row.density) for row in aborted), F())
    acceptance_mass = sum((native.trace(row.density) for row in accepted), F())
    require(abort_mass + acceptance_mass == native.trace(second['success']),
            'the complete herald must preserve its original unnormalized input mass')

    # The herald precedes either setting choice and is identical in every
    # branch. Aborted runs output (+1,+1), independently of both settings.
    settings = ((F(), theta / 2), (F(), -theta / 2),
                (theta, theta / 2), (theta, -theta / 2))
    readout_program = native.program(federation={'level': 0, 'total': 2}, biased=False,
        commands=[{'op': 'repair', 'selector': {'kind': 'fixed', 'seam': list(local_edge)},
                   'law': 'E', 'repeat': 1} for local_edge in rail_edges],
        capacity=2, overflow='export_oldest')
    correlations = []
    setting_reports = []
    for left, right in settings:
        score = abort_mass
        probability = abort_mass
        complete_readout_slots = 0
        retained_herald_readouts = []
        for row in accepted:
            left_run = beam.lift(beam.declare(left, targets[0], cells[0]), row.density, carrier)
            right_run = beam.lift(beam.declare(right, targets[1], cells[1]), left_run['success'], carrier)
            history = native.execute(readout_program, right_run['success'])
            require(history.status == 'COMPLETE' and history.attempted_steps == 2,
                    'the admitted two-attempt local readout must complete')
            require(len(history.branches) == 625,
                    'the complete two-attempt readout must retain every raw slot')
            complete_readout_slots += len(history.branches)
            readout_law = []
            for branch in history.branches:
                mass = native.trace(branch.density)
                require(mass >= 0, 'an actual complete probability cannot be negative')
                require(len(branch.memory) == 2, 'both local readout records must remain in the FIFO')
                ra, rb = branch.memory
                # On these one-occupation rails the initial native E contrast
                # is precisely logical Z, regardless of the retained coin.
                if mass:
                    require(ra.law == rb.law == 'E', 'only the admitted E readout has nonzero mass')
                    require(ra.label in (-1, 1) and rb.label in (-1, 1),
                            'the one-occupation rail record must give a complete binary output')
                    score += ra.label * rb.label * mass
                probability += mass
                readout_law.append({'left_record': asdict(ra), 'right_record': asdict(rb),
                                    'mass': str(mass)})
            retained_herald_readouts.append({'herald_record': asdict(row.record),
                                            'complete_native_readout_law': readout_law})
        correlations.append(score)
        setting_reports.append({'angles': [str(left), str(right)], 'correlation': str(score),
                                'complete_probability_mass': str(probability),
                                'retained_native_readout_slots': complete_readout_slots,
                                'accepted_herald_histories': retained_herald_readouts})

    score = correlations[0] + correlations[1] + correlations[2] - correlations[3]
    sine, _ = iv.sincos(2 * theta)
    probability = iv.divide(iv.mul(sine, sine), 2)
    _, cosine = iv.sincos(theta)
    _, triple_cosine = iv.sincos(3 * theta)
    bell_score = iv.add(iv.mul(iv.exact(3), cosine), iv.neg(triple_cosine))
    ideal_score = iv.add(iv.exact(2), iv.mul(probability, iv.add(bell_score, iv.exact(-2))))
    lower, upper = (F(x, iv.Q) for x in ideal_score)
    # Up to four coherent services occur on any complete path. Native
    # instruments and retained classical stopping are exact CP instruments.
    error = (1 + beam.GATE)**4 - 1
    halfdiamond = error * (2 + error) / 2
    require(lower - 8 * halfdiamond <= score <= upper + 8 * halfdiamond,
            'the stored complete CHSH score exceeds the source implementation error bound')
    require(lower - 8 * halfdiamond > F(12, 5),
            'the directed source bound does not certify the declared complete score')
    require(score > F(12, 5), 'the actual stored complete score must pass the declared bound')
    for setting in setting_reports:
        require(abs(F(setting['complete_probability_mass']) - 1) <= 2 * halfdiamond,
                'the complete stored probability mass exceeds its certified gain error')

    paths = [BELL / 'derived_control.py', BELL / 'bell.py', Path(__file__),
             ROOT / 'model/native_program/interpreter.py',
             ROOT / 'code/source_selection_model/response.json',
             ROOT / 'code/m1_source_realization/receipt.json']
    return {'status': 'PASS', 'schema': 'oph.finite_complete_source_chsh.v1',
            'input_interface': 'known native A0 B1 basis state and declared coherent M1 load/export',
            'field_or_spacelike_experiment_claimed': False,
            'setting_independent_herald': True, 'abort_output_pair': [1, 1],
            'conditional_acceptance_normalization_used': False,
            'acceptance_mass': str(acceptance_mass), 'abort_mass': str(abort_mass),
            'complete_CHSH_score': str(score), 'complete_CHSH_score_diagnostic': float(score),
            'ideal_CHSH_interval': [str(lower), str(upper)],
            'complete_CHSH_score_certified_lower': '12/5',
            'one_setting_complete_halfdiamond_upper': str(halfdiamond),
            'herald_raw_slots_retained': 25, 'native_readout_raw_slots_per_setting': 1250,
            'maximum_native_attempts_per_path': 3, 'maximum_coherent_source_pulses_per_path': 4,
            'maximum_retained_native_records_per_path': 3,
            'native_records_explained': 'one herald record retained by controller, two local readout FIFO records',
            'source_complement_outputs_per_accepted_path': 8,
            'herald_native_program': herald_program, 'readout_native_program': readout_program,
            'abort_herald_records': [{'record': asdict(row.record),
                                     'mass': str(native.trace(row.density))} for row in aborted],
            'settings': setting_reports,
            'source_bindings': {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}


if __name__ == '__main__':
    result = run()
    Path(__file__).with_name('chsh-results.json').write_text(json.dumps(result, indent=2) + '\n')
    print('PASS: complete CHSH score', result['complete_CHSH_score_diagnostic'], '> 12/5 with all herald failures')
