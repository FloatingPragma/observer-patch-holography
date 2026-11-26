import copy
import itertools
import json
from math import pi
from pathlib import Path

import numpy as np
import pytest

from ... import apparatus as app
from .readback import check
from .service import Node, program, run


def receipt():
    return json.loads(Path(__file__).with_name('results.json').read_text())


def test_exact_native_joint_pre_reads_and_terminal_against_classical_sampling_law():
    states = list(itertools.product((0, 2), repeat=4))
    marginals = np.array([.14, .83, .57, .39])
    init = np.array([np.prod(np.where(np.asarray(state) == 2, marginals, 1-marginals)) for state in states])
    for angle in (pi/2, pi/3):
        local = app.pulse(angle)
        branches = {(): (np.diag(init).astype(complex), init.copy())}
        for a, b in [(0, 1), (0, 2), (0, 1)]:
            gate = np.zeros((16, 16), complex)
            swap = []
            for i, row in enumerate(states):
                new = list(row)
                new[a], new[b] = new[b], new[a]
                swap.append(states.index(tuple(new)))
                col = app.LAW.SRC.index((row[a], row[b]))
                for j, candidate in enumerate(states):
                    if all(row[k] == candidate[k] for k in set(range(4)) - {a, b}):
                        gate[j, i] = local[app.LAW.SRC.index((candidate[a], candidate[b])), col]
            updated = {}
            p = np.sin(angle/2)**2
            for prefix, (density, weight) in branches.items():
                for value in (0, 2):
                    keep = np.array([row[a] == value for row in states])
                    mat = gate @ (density * keep[:, None] * keep[None, :]) @ gate.conj().T
                    flat = weight * keep
                    predicted = (1-p)*flat + p*flat[swap]
                    np.testing.assert_allclose(mat.diagonal(), predicted, atol=2e-14)
                    updated[prefix + (value,)] = mat, predicted
            branches = updated


def test_controller_replay_and_statistical_controls():
    data = receipt()
    assert check(data)['reached_looks'] == 52
    assert all(v['covered'] for row in data['panels'] for v in row['withheld_evaluation'])
    cases = {(row['family'], row['law']): row['panel'] for row in data['panels']}
    assert [cases['irrational', law]['records'][0]['decision_read']['sign'] for law in ('native', 'third_turn')] == ['high', 'low']
    assert [cases[label, 'native']['refusal'] for label in ('matched_a', 'matched_b', 'near_high',
                                                           'second_zero', 'terminal_zero')] == [
        'decision_unresolved', 'decision_unresolved', 'decision_unresolved', 'decision_unresolved', 'terminal_unresolved']
    changed = copy.deepcopy(data)
    word = changed['panels'][0]['panel']['public']['score_words_at_both_owners'][1][0]
    word[0] = word[0][::-1]
    with pytest.raises(ValueError, match='damaged terminal scores'):
        check(changed)


def test_budget_and_unsupported_program_refuse_without_reading():
    for tree, replicas, message in [(program(), 100, 'unsupported_budget_or_audit_depth'),
                                    (Node((2, 5), 0, 0), 4*32768, 'unsupported_layout_or_read')]:
        panel = run(tree, [0.]*12, np.random.SeedSequence(9), replicas=replicas)
        assert panel['refusal'] == message and panel['measured_copies'] == 0
