"""Complete coherent source-clock histories and their CP effective instrument.

The exact evaluator is the retained scientific result. The common observer
runtime independently executes Kraus factorizations; its actual binary64
outputs receive an exact entrywise error bound against that result.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
import math
from pathlib import Path
import random

import numpy as np

from . import quantum_clock as q
from ..observer import native
from ..primitive_observer import runtime
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*q.SOURCE_PATHS,
    'model/observer_correspondence/quantum_experiment.py',
    'model/primitive_observer/runtime.py', 'model/common_source/open_protocol.py',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py'))))
MENU = (*q.MARKS, q.TIMEOUT)
NUMERIC_ALLOWANCE = F(1, 10**10)


def _flags(feedback, reference=False):
    if type(feedback) is not bool or type(reference) is not bool:
        raise ValueError('original Boolean protocol and reference choices required')


def initial_density(reference=False):
    _flags(False, reference)
    r = 2 if reference else 1
    result = np.full((12*r, 12*r), F(0), dtype=object)
    indices = (0, 3) if reference else (0, 1)
    for i, j in product(indices, repeat=2):
        result[i, j] = F(1, 2)
    return result


def _phase(rho, reference):
    signs = tuple(-1 if i//reference == 0 else 1 for i in range(len(rho)))
    return np.array([[signs[i]*value*signs[j] for j, value in enumerate(row)]
                     for i, row in enumerate(rho)], dtype=object)


def _analyzers():
    plus, minus = (np.full((12, 12), F(0), dtype=object) for _ in range(2))
    for i, j in product((0, 1), repeat=2):
        plus[i, j] = F(1, 2)
        minus[i, j] = F(1 if i == j else -1, 2)
    rest = np.diag([F(int(i not in (0, 1))) for i in range(12)])
    return plus, minus, rest


def _analyze(rho, reference):
    # Only the two logical rows/columns mix. Preserve the untouched reference.
    outputs = []
    for sign in (1, -1):
        result = np.full(rho.shape, F(0), dtype=object)
        for a, b in product(range(reference), repeat=2):
            value = (rho[a, b]+sign*rho[a, reference+b]
                     +sign*rho[reference+a, b]+rho[reference+a, reference+b])/4
            for i, j in product((0, 1), repeat=2):
                result[i*reference+a, j*reference+b] = sign**(i+j)*value
        outputs.append(result)
    rest = rho.copy()
    rest[:2*reference, :] = F(0)
    rest[:, :2*reference] = F(0)
    return (*outputs, rest)


def exact_history(law, cap, *, effective=False, feedback=True, reference=False):
    """All 81 complete labels, without successful-branch normalization."""
    _flags(feedback, reference)
    maps = q.capped_maps(law, cap, effective)
    r = 2 if reference else 1
    prepared = initial_density(reference)
    outputs = {}
    for first in MENU:
        state = q.apply(maps[first], prepared, r)
        if first == q.TIMEOUT:
            outputs[(q.TIMEOUT,)] = state
            continue
        if feedback and first[:2] == (0, 1):
            state = _phase(state, r)
        for second in MENU:
            final = q.apply(maps[second], state, r)
            if second == q.TIMEOUT:
                outputs[first, q.TIMEOUT] = final
            else:
                for field, density in enumerate(_analyze(final, r)):
                    outputs[first, second, ('phase', field)] = density
    if len(outputs) != 81 or sum(np.trace(value) for value in outputs.values()) != 1:
        raise ArithmeticError('the complete exact quantum history lost outcomes or mass')
    return outputs


def root_output(density, reference=1):
    """Actual partial trace: local vacuum, then occupations at 0,1,5."""
    if (type(reference) is not int or reference not in (1, 2)
            or type(density) is not np.ndarray
            or density.shape != (12*reference, 12*reference)):
        raise ValueError('one-particle density and its original reference dimension required')
    result = np.zeros((4*reference, 4*reference), dtype=density.dtype)
    for i, j in product(range(3), repeat=2):
        result[(i+1)*reference:(i+2)*reference, (j+1)*reference:(j+2)*reference] = density[
            q.UNION[i]*reference:(q.UNION[i]+1)*reference,
            q.UNION[j]*reference:(q.UNION[j]+1)*reference]
    for port in range(12):
        if port not in q.UNION:
            result[:reference, :reference] += density[
                port*reference:(port+1)*reference, port*reference:(port+1)*reference]
    return result


def _kraus(channel):
    """Factor the small diagonal Choi block plus its rank-one jump blocks.

    Rational Schur elimination certifies the small Gram block before square
    roots enter the execution representation. Exact zero pivots/rows remain
    zero, avoiding spurious numerical support in impossible record branches.
    """
    p = np.array(channel.population, dtype=float)
    gram = np.array(channel.coherence, dtype=object)+np.diag(
        [channel.population[i][i] for i in range(12)])
    groups = []
    for k in range(12):
        pivot = gram[k, k]
        if pivot < 0 or (pivot == 0 and any(gram[k, j] for j in range(k, 12))):
            raise ArithmeticError('native quantum instrument lost its exact CP factorization')
        if not pivot:
            continue
        vector = [F(0)]*k+[gram[i, k]/pivot for i in range(k, 12)]
        groups.append(np.diag([float(x)*math.sqrt(float(pivot)) for x in vector]))
        for i in range(k+1, 12):
            for j in range(i, 12):
                gram[i, j] -= pivot*vector[i]*vector[j]
                gram[j, i] = gram[i, j]
    if p.min() < 0:
        raise ArithmeticError('native population jump lost positivity')
    for i, j in product(range(12), repeat=2):
        if i != j and p[i, j]:
            matrix = np.zeros((12, 12))
            matrix[i, j] = math.sqrt(p[i, j])
            groups.append(matrix)
    return tuple(groups) or (np.zeros((12, 12)),)


def runtime_program(law, cap, *, effective=False, feedback=True):
    _flags(feedback)
    maps = q.capped_maps(law, cap, effective)
    event = runtime.Instrument(kraus=tuple(_kraus(maps[mark]) for mark in MENU), visible=True,
        resources=(('native_attempt_budget', cap), ('interval_requests', 1)),
        binding={'law': law, 'cap': cap, 'symmetric': effective, 'marks': MENU})
    identity = np.eye(12)
    zero = np.zeros((12, 12))
    stopped = runtime.Instrument(kraus=tuple((zero,) for _ in q.MARKS)+((identity,),), visible=True,
        binding={'operation': 'no further interval after timeout'})
    unchanged = runtime.Instrument(kraus=((identity,),), binding={'operation': 'identity'})
    z = identity.copy(); z[0, 0] = -1
    phase = runtime.Instrument(kraus=((z,),), resources=(('local_phase_controls', 1),),
        binding={'operation': 'Z at port0, ideal declared local service'})
    analyzer = runtime.Instrument(kraus=tuple((np.array(matrix, dtype=float),) for matrix in _analyzers()),
        visible=True, resources=(('phase_analyzers', 1),), binding={'outcomes': ['plus', 'minus', 'complement']})
    unread = runtime.Instrument(kraus=((zero,), (zero,), (identity,)), visible=True,
        binding={'operation': 'no analyzer after timeout'})
    events = (
        runtime.Event('first incident interval', instruments=(event,)),
        runtime.Event('record-controlled local phase', control=0, instruments=tuple(
            phase if feedback and mark != q.TIMEOUT and mark[:2] == (0, 1) else unchanged for mark in MENU)),
        runtime.Event('second interval unless stopped', control=0,
                      instruments=tuple(event if mark != q.TIMEOUT else stopped for mark in MENU)),
        runtime.Event('terminal coherent field unless stopped', control=1,
                      instruments=tuple(analyzer if mark != q.TIMEOUT else unread for mark in MENU)))
    return runtime.Program('native quantum incident-clock correspondence', 12,
        np.array(initial_density(), dtype=float), events,
        preparation_resources=(('declared_coherent_preparations', 1),),
        binding={'source': 'native twelve-port one-particle carrier', 'law': law,
                 'cap': cap, 'effective': effective, 'feedback': feedback})


def _decode_history(history):
    outputs = {}
    for branch in history.branches:
        first, second, phase = branch.record
        if first == 5:
            valid, word = second == 5 and phase == 2, (q.TIMEOUT,)
        elif second == 5:
            valid, word = phase == 2, (MENU[first], q.TIMEOUT)
        else:
            valid, word = True, (MENU[first], MENU[second], ('phase', phase))
        if not valid:
            if np.any(branch.density != 0):
                raise ArithmeticError('an impossible padding record carries quantum output')
        else:
            outputs[word] = branch.density
    if len(outputs) != 81:
        raise ArithmeticError('runtime record projection lost complete outcomes')
    return outputs


def binary_output_bound(numerical, exact):
    """Rigorous a-posteriori trace-distance bound for the stored binary64 data.

    Each binary float becomes its exact rational value. The sum of absolute
    real/imaginary entry errors bounds the sum of matrix trace norms.
    """
    if type(numerical) is not dict or type(exact) is not dict or numerical.keys() != exact.keys():
        raise ValueError('complete equal quantum output menus required')
    total = F(0)
    for word, target in exact.items():
        actual = numerical[word]
        if (type(actual) is not np.ndarray or type(target) is not np.ndarray
                or actual.dtype not in (np.dtype('float64'), np.dtype('complex128'))
                or target.dtype != object or actual.ndim != 2
                or actual.shape != target.shape or not np.isfinite(actual).all()
                or any(type(x) not in (int, F) for x in target.flat)):
            raise ValueError('finite matching complete output matrices required')
        for value, wanted in zip(actual.flat, target.flat, strict=True):
            if type(wanted) is F:
                numerator, denominator = wanted.numerator, wanted.denominator
                # Fraction's protected slots are writable. Its own abs() can
                # return a negative value for a corrupted negative denominator,
                # cancelling a genuine matrix error. Validate then detach.
                if (type(numerator) is not int or type(denominator) is not int
                        or denominator <= 0 or math.gcd(numerator, denominator) != 1):
                    raise ValueError('canonical exact rational target entries required')
                wanted = F(numerator, denominator)
            total += abs(F(float(value.real))-wanted)+abs(F(float(value.imag)))
    return total/2


def _draw(rng, weights):
    """Exact rational sampling; no float-dependent branch decisions."""
    if not weights or any(w < 0 for w in weights) or sum(weights) != 1:
        raise ArithmeticError('a complete exact sampling distribution is required')
    denominator = math.lcm(*(F(w).denominator for w in weights))
    pick = rng.randrange(denominator)
    for index, weight in enumerate(weights):
        pick -= int(weight*denominator)
        if pick < 0:
            return index
    raise ArithmeticError('sampling lost probability')


def source_trajectory(law, cap, seed, *, effective=False, feedback=True, reference=False):
    """Run actual native repairs; private audit records never control the observer.

    These samples illustrate realized costs and timeouts. The complete exact
    instruments, not sample frequencies, establish the scientific conclusions.
    """
    _flags(feedback, reference)
    source = q.source(law)
    q.capped_maps(law, cap, effective)  # validate original cap/selector before RNG
    if type(seed) is not int or seed < 0:
        raise ValueError('an original nonnegative integer seed is required')
    rng = random.Random(seed)
    incident = tuple(e for e in source.seams if set(e).intersection(q.UNION))
    remote = tuple(e for e in source.seams if not set(e).intersection(q.UNION))
    r = 2 if reference else 1
    rho = initial_density(reference)
    audit, public = [], []
    controls = 0

    def repair(seam):
        nonlocal rho
        coin = rng.randrange(2)
        groups = {}
        for port in range(12):
            loads = tuple(int(i == port) for i in range(12))
            out, (_, current) = native.native(loads, seam, coin)
            label = current if law == 'C' else loads[seam[0]]-loads[seam[1]]
            groups.setdefault((label, current), []).append((port, out.index(1)))
        alternatives = []
        for key, images in sorted(groups.items()):
            value = np.full(rho.shape, F(0), dtype=object)
            for (i, a), (j, b) in product(images, repeat=2):
                value[a*r:(a+1)*r, b*r:(b+1)*r] += rho[i*r:(i+1)*r, j*r:(j+1)*r]
            alternatives.append((key, value))
        chosen = _draw(rng, [np.trace(value) for _, value in alternatives])
        (label, current), selected = alternatives[chosen]
        rho = selected/np.trace(selected)
        audit.append([list(seam), coin, label, current])
        return (*(int(bool(set(seam).intersection(region))) for region in q.REGIONS),
                current*(int(seam[1] == 0)-int(seam[0] == 0)))

    for interval in range(2):
        if effective:
            n = 0
            while n < cap and rng.randrange(5) >= 2:
                n += 1
            if n == cap:
                for _ in range(cap):
                    repair(remote[rng.randrange(len(remote))])
                public.append(q.TIMEOUT)
                break
            split = _draw(rng, [F(math.comb(2*k, k)*math.comb(2*(n-k), n-k), 4**n)
                                for k in range(n+1)])
            for _ in range(split):
                repair(remote[rng.randrange(len(remote))])
            mark = repair(incident[rng.randrange(len(incident))])
            for _ in range(n-split):
                repair(remote[rng.randrange(len(remote))])
        else:
            mark = q.TIMEOUT
            for _ in range(cap):
                seam = source.seams[rng.randrange(len(source.seams))]
                observed = repair(seam)
                if seam in incident:
                    mark = observed
                    break
        public.append(mark)
        if mark == q.TIMEOUT:
            break
        if interval == 0 and feedback and mark[:2] == (0, 1):
            rho = _phase(rho, r)
            controls += 1
    if public[-1] != q.TIMEOUT:
        values = _analyze(rho, r)
        outcome = _draw(rng, [np.trace(v) for v in values])
        rho = values[outcome]/np.trace(values[outcome])
        public.append(('phase', outcome))
    if len(audit) > 2*cap or np.trace(rho) != 1:
        raise ArithmeticError('native trajectory violated its complete resource/state contract')
    return {'law': law, 'cap': cap, 'seed': seed, 'effective': effective,
            'reference_dimension': r, 'feedback': feedback, 'public_word': public,
            'native_repairs': audit, 'native_calls': len(audit), 'local_phase_controls': controls,
            'terminal_root_reference_sha256': _digest({('terminal',): rho}, r)}


def execute_history(law, cap, *, effective=False, feedback=True, reference=False):
    _flags(feedback, reference)
    program = runtime_program(law, cap, effective=effective, feedback=feedback)
    prepared = np.array(initial_density(True), dtype=float) if reference else None
    prefix = runtime.execute(program, initial=prepared, stop_after=2)
    serialized = runtime.checkpoint_to_json(prefix)
    restored = runtime.checkpoint_from_json(program, serialized)
    resumed = runtime.execute(program, checkpoint=restored)
    full = runtime.execute(program, initial=prepared)
    if any(not np.array_equal(a.density, b.density) or a.resources != b.resources
           or a.record != b.record for a, b in zip(resumed.branches, full.branches, strict=True)):
        raise ArithmeticError('source-clock continuation changed retained quantum output')
    outputs = _decode_history(full)
    exact = exact_history(law, cap, effective=effective, feedback=feedback, reference=reference)
    bound = binary_output_bound(outputs, exact)
    if bound > NUMERIC_ALLOWANCE:
        raise ArithmeticError('the common runtime exceeded its exact output error allowance')
    return {'history': full, 'outputs': outputs, 'exact': exact, 'error_bound': bound,
            'checkpoint_bytes': len(serialized.encode())}


def _digest(outputs, reference):
    rows = [[word, [[str(x) for x in row] for row in root_output(value, reference)]]
            for word, value in outputs.items()]
    return hashlib.sha256(canonical(rows).encode()).hexdigest()


def _summary(outputs, reference):
    field = [F(0)]*3
    timeout = F(0)
    for word, value in outputs.items():
        if word[-1] == q.TIMEOUT:
            timeout += np.trace(value)
        else:
            field[word[-1][1]] += np.trace(value)
    return {'complete_output_labels': len(outputs), 'timeout_probability': str(timeout),
            'phase_outcomes_unconditioned': list(map(str, field)),
            'root_reference_output_sha256': _digest(outputs, reference)}


def _json_exact(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, (list, tuple)):
        return [_json_exact(item) for item in value]
    if type(value) is dict:
        return {key: _json_exact(item) for key, item in value.items()}
    return value


def reference_witness(law, cap, effective=False):
    """A genuine pre-analyzer two-tick branch, with no success conditioning."""
    channel = q.capped_maps(law, cap, effective)[(1, 1, 0)]
    state = q.apply(channel, q.apply(channel, initial_density(True), 2), 2)
    local = root_output(state, 2)
    # Partial transpose on R sends [2,5] to [3,4]. Its 2x2 minor
    # has zero diagonal; a negative determinant proves nonseparability.
    minor = local[3, 3]*local[4, 4]-local[2, 5]*local[5, 2]
    expected = -F(625, 1327104)*(1-F(3, 5)**cap)**4 if law == 'C' else F(0)
    if minor != expected:
        raise ArithmeticError('the exact reference-sensitive branch changed')
    return {'public_prefix': [[1, 1, 0], [1, 1, 0]],
            'stage': 'after second incident, before destructive phase analyzer',
            'unnormalized_partial_transpose_minor': str(minor),
            'branch_mass': str(np.trace(state)),
            'negative_minor_certifies_entanglement': minor < 0}


def run():
    cap = 16
    p, remote = F(2, 5), F(3, 5)
    success = 1-remote**cap
    delta = 1-success**2
    rows = []
    for law, feedback, reference in product(('C', 'E'), (False, True), (False, True)):
        pair, states = [], []
        for effective in (False, True):
            result = execute_history(law, cap, effective=effective, feedback=feedback, reference=reference)
            summary = _summary(result['exact'], 2 if reference else 1)
            if F(summary['timeout_probability']) != delta:
                raise ArithmeticError('the full stopping budget differs from its source derivation')
            contrast = F(summary['phase_outcomes_unconditioned'][0])-F(summary['phase_outcomes_unconditioned'][1])
            prediction = {('C', False): F(49, 144), ('E', False): F(1, 16),
                          ('C', True): F(7, 144), ('E', True): -F(1, 16)}[law, feedback]
            if contrast != (F(0) if reference else success**2*prediction):
                raise ArithmeticError('the coherent field differs from its exact source prediction')
            states.append({word: root_output(value, 2 if reference else 1)
                           for word, value in result['exact'].items()})
            pair.append({'representation': 'effective' if effective else 'native', **summary,
                         'unconditional_plus_minus_contrast': str(contrast),
                         'runtime_complete_output_bound': str(NUMERIC_ALLOWANCE),
                         'serialized_continuation_identical': True})
        discrepancy = sum(sum(abs(x) for x in (states[0][word]-states[1][word]).flat)
                          for word in states[0])/2
        if discrepancy > 2*delta:
            raise ArithmeticError('complete local quantum outputs exceed the cap comparison budget')
        rows.append({'law': law, 'feedback': feedback, 'untouched_reference_dimension': 2 if reference else 1,
                     'exact_local_output_entrywise_distance_bound': str(discrepancy),
                     'experiments': pair})
    trajectories = [source_trajectory(law, h, seed, effective=effective, reference=reference)
                    for law, effective, reference, h, seed in
                    product(('C', 'E'), (False, True), (False, True), (1, cap), range(2))]
    if not any(row['public_word'][-1] == q.TIMEOUT for row in trajectories):
        raise ArithmeticError('the retained raw execution omitted its timeout control')
    return _json_exact({'schema': 'oph.quantum_clock_correspondence.v1',
        'source_sha256': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS}),
        'contract': {'states': 12, 'conserved_total': 1, 'root': list(q.UNION),
            'regions': [list(region) for region in q.REGIONS], 'marks': list(q.MARKS),
            'initial_occupation_disclosure': False, 'initial_remote_state': 'vacuum; fixed by J',
            'preparation': 'root plus or root/reference Bell state on ports0,1',
            'feedback': 'local Z0 after first B-only clock mark', 'cap_per_interval': cap,
            'requested_incident_ticks': 2, 'analyzer': 'plus,minus,complement; only after two completed ticks',
            'comparison_algebra': 'root vacuum/one-particle algebra and untouched reference',
            'ideal_preparation_phase_and_detector_services': True},
        'bounds': {'per_representation_complete_output': str(delta),
            'native_vs_effective_complete_output': str(2*delta),
            'runtime_vs_exact_per_execution': str(NUMERIC_ALLOWANCE),
            'maximum_native_calls': 2*cap, 'expected_native_calls': str(success*(1+success)/p),
            'retained_runtime_branches': 108, 'maximum_work_reference_dimension': 24,
            'maximum_terminal_complex_entries': 108*24**2,
            'maximum_public_record_symbols': 3,
            'expected_feedback_controls_if_enabled': str(success/4),
            'expected_terminal_analyzers': str(success**2),
            'joint_reference_preparation': 'one supplied joint-input service; not a native pulse cost'},
        'one_particle_floors': {law: str(q.source(law).incident_hs_lower) for law in ('C', 'E')},
        'full_binary_triangle': q.triangle_certificate(),
        'total_two_floors': _json_exact(q.total_two_certificate()),
        'record_law_boundary': _json_exact(q.record_family_obstruction()),
        'pre_analyzer_reference_witnesses': [dict(law=law, effective=effective,
            **reference_witness(law, cap, effective)) for law, effective in product(('C', 'E'), (False, True))],
        'native_source_trajectories': trajectories,
        'experiments': rows,
        'scope': {'local_quantum_and_reference_correspondence': True,
            'finite_caps_compared_by_error_bound_not_exact_equality': True,
            'full_source_states_identified': False, 'fractional_CPTP_flow_established': False,
            'physical_time_or_relativistic_QFT_established': False}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', type=Path)
    group.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('--output requires a new destination')
    result = run()
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text(encoding='utf-8')), result)
        print('Complete quantum source-clock correspondence replayed.')
    elif args.output:
        with args.output.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result))


if __name__ == '__main__':
    main()
