"""Apply the protected-obstruction classifier to exact phase-pair loss."""
from fractions import Fraction
import importlib.util

from ..observer import native
from ..observer.quantum_mixing import seams


CLASSIFIER_PATH = 'Lean/ObservableNormalForms/tools/verify_protected_obstruction_models.py'


def classify_phase_loss(root):
    states = native.states(binary=True)
    positions = {state: index for index, state in enumerate(states)}
    edges = seams()
    pairs = set()
    for a, b in edges:
        for spectator in range(12):
            if spectator in (a, b):
                continue
            x = tuple(int(port in (a, spectator)) for port in range(12))
            y = tuple(int(port in (b, spectator)) for port in range(12))
            if x[0] != y[0]:
                pairs.update(((positions[x], positions[y]), (positions[y], positions[x])))
    if len(pairs) != 100:
        raise ArithmeticError('the source one-port phase support changed')
    spec = importlib.util.spec_from_file_location('oph_protected_classifier', root / CLASSIFIER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = {}
    for law in ('C', 'E'):
        counts = set()
        for i, j in pairs:
            survived = 0
            for edge in edges:
                for coin in (0, 1):
                    x, (_, current_x) = native.native(states[i], edge, coin)
                    y, (_, current_y) = native.native(states[j], edge, coin)
                    record_x = (current_x, (states[i][edge[0]] + states[i][edge[1]]) % 2)
                    record_y = (current_y, (states[j][edge[0]] + states[j][edge[1]]) % 2)
                    same = current_x == current_y if law == 'C' else record_x == record_y
                    if same:
                        if (positions[x], positions[y]) not in pairs:
                            raise ArithmeticError('a retained phase pair leaves the declared support')
                        survived += 1
            counts.add(survived)
        if len(counts) != 1:
            raise ArithmeticError('phase loss has no common two-state quotient')
        count = counts.pop()
        kernel = ((count, 60 - count), (0, 60))
        result[law] = {'survival': str(Fraction(count, 60)),
                       'pair_count': len(pairs), 'kernel': [list(row) for row in kernel],
                       'classification': module.classify(kernel, 0, frozenset({1}))}
    result['hold_control'] = module.classify(((60, 0), (0, 60)), 0, frozenset({1}))
    return result
