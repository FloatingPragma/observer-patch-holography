"""Outward finite-word source calculation with actual mixed states and flags."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import mpmath
from interval_math import (iv, real, complex_, require, packed, upper, lower,
    width, root_of_unity, zeros, eye, adj, mm, add, scale, power, trace,
    sandwich, hermitian_trace, difference_upper, diagonal, conjugate)

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[2]/'row_01_charge_quantization/work/2026-11-16_codex_common_reference_menu'
EDGES = ((0, 1), (1, 2), (0, 2))


def native(q, edge, coin, law='native'):
    a, b = edge
    d = q[a]-q[b]
    j = ((d-(2*coin-1)*(abs(d) % 2))//2 if law == 'native' else d)
    out = list(q)
    out[a] -= j
    out[b] += j
    return out


def mean(edges, law='native'):
    out = zeros(3)
    for edge in edges:
        for coin in (0, 1):
            for i in range(3):
                q = [int(j == i) for j in range(3)]
                y = native(q, edge, coin, law).index(1)
                out[y][i] += real(1)/(2*len(edges))
    return out


def basis():
    s2, s6 = iv.sqrt(real(2)), iv.sqrt(real(6))
    e = [[complex_(), complex_(-2/s6)],
         [complex_(1/s2), complex_(1/s6)],
         [complex_(-1/s2), complex_(1/s6)]]
    b = zeros(4)
    for i in range(2):
        b[i][2*i] = b[i][2*i+1] = complex_(1/s2)
        b[2+i][2*i] = complex_(0, -1/s2)
        b[2+i][2*i+1] = complex_(0, 1/s2)
    return e, b


E, B = basis()
BD = adj(B)


def lift(position):
    rho = zeros(4)
    for i in range(2):
        for j in range(2):
            rho[i][j] = position[i][j]
    return rho


def counter(k):
    """Exact geometric finite word, enclosed at every arithmetic operation."""
    q = 256*2**k
    eps = real(1)/(2**k)
    z = []
    for lam in (4, 6):
        p = 1-eps*eps*lam/10
        t = iv.sqrt(1-p*p)
        z.extend([complex_(p, t), complex_(p, -t)])
    numer = [(1-a**q)/q for a in z]
    root = root_of_unity(q)
    phase = complex_(1)
    gram = zeros(4)
    for r in range(1, q//2):
        phase *= root
        vals = [n/(1-a*phase) for n, a in zip(numer, z)]
        for i in range(4):
            for j in range(i, 4):
                gram[i][j] += vals[i]*conjugate(vals[j])
    for i in range(4):
        for j in range(i):
            gram[i][j] = conjugate(gram[j][i])
    return gram, z


def cp(gram, rho):
    spectral = mm(mm(BD, rho), B)
    out = [[spectral[i][j]*gram[i][j] for j in range(4)] for i in range(4)]
    return mm(mm(B, out), BD)


def preparation(gram):
    p = mean(((0, 1), (0, 1), (0, 2), (0, 2), (1, 2)))
    v = [[complex_(0)], [complex_(1)], [complex_(0)]]
    v = scale(mm(add(eye(3), scale(p, -1)), v), 1/(2*iv.sqrt(real(3))))
    for _ in range(8):
        v = scale(mm(p, v), 1/iv.sqrt(real(3)))
    optical = sandwich(adj(E), mm(v, adj(v)))
    rho = cp(gram, lift(optical))
    # D0 is the actual admitted arc decoder. Do not replace its output by f.
    position = [row[:2] for row in rho[:2]]
    require(hermitian_trace(position) > 0, 'nonpositive preparation probability')
    return position


def slot_operator(slot, m, law):
    v, r = divmod(slot % (3*m), m)
    w = (v+1) % 3
    edges = [e for e in EDGES if v in e]*(m-r)+[e for e in EDGES if w in e]*r
    p = mean(edges, law)
    g = [E[1], E[2], E[0]]
    # Four complete accepted target outcomes each have operator P/sqrt(3).
    return scale(power(mm(mm(adj(g), p), g), 4), real(1)/9)


def case(gram, z, position, k, m, reference, law, turns, direction):
    ar, af = eye(2), eye(2)
    frames = 3*m*turns
    for step in range(1, frames+1):
        r = step % m
        rot = (2*m+direction*step) % (3*m)
        ref = ((2*m+direction*min(r, m-r)) % (3*m) if reference == 'matched' else 2*m)
        ar = mm(slot_operator(rot, m, law), ar)
        af = mm(slot_operator(ref, m, law), af)
    arms = (ar, af)
    out = [[cp(gram, lift(scale(sandwich(arms[a], position, arms[b]), real(1)/2)))
            for b in range(2)] for a in range(2)]
    total = hermitian_trace(add(out[0][0], out[1][1]))
    require(total > 0, 'nonpositive analyzer denominator')
    cross = trace(out[0][1]).real
    visibility = 2*cross/total
    sign = -1 if law == 'native' and turns == 1 else 1
    require(visibility < 0 if sign == -1 else visibility > 0, 'visibility sign')
    tau = scale(add(add(out[0][0], out[1][1]), scale(add(out[0][1], out[1][0]), sign)), real(1)/2)
    success = hermitian_trace(tau)
    require(success > 0, 'nonpositive returned-clock probability')
    spectral = scale(mm(mm(BD, tau), B), 1/success)
    weights = [spectral[i][i].real for i in range(4)]
    single_off = weights[2]+weights[3]
    pair_off = 1-weights[0]**2-weights[1]**2
    require(single_off < real('0.002'), 'single impurity gate')
    require(pair_off < real('0.01'), 'pair impurity gate')
    records = {}
    for stride in (1, 8):
        moment = sum((weights[i]*z[i]**(stride*2**k) for i in range(4)), complex_())**2
        records[str(stride)] = {'x_plus': (1+moment.real)/2, 'x_minus': (1-moment.real)/2,
                                'y_plus': (1+moment.imag)/2, 'y_minus': (1-moment.imag)/2}
    plane = real(1)/(48 if law == 'native' else 12)
    hops = (1-real(1)/256)**(frames*(2 if reference == 'matched' else 1))
    selectors = (1-(real(3)/8)**8)**11*(1-(real(1)/4)**8)**frames
    complete_cost = success*plane**frames*hops*selectors
    require(complete_cost > 0, 'complete source probability')
    return {'k': k, 'M': m, 'reference': reference, 'law': law, 'turns': turns, 'direction': direction,
            'visibility': visibility, 'single_off_band': single_off, 'pair_off_band': pair_off,
            'selected_clock_probability': complete_cost, 'independent_pair_probability': complete_cost**2,
            'records': records}


def run(ks, references):
    started = time.perf_counter()
    rows = []
    max_error = 0.
    source_hashes = {}
    for k in ks:
        gram, z = counter(k)
        position = preparation(gram)
        for reference in references:
            path = SOURCE/f'k{k}_{reference}/RECEIPT.json'
            source_hashes[str(path.relative_to(SOURCE))] = hashlib.sha256(path.read_bytes()).hexdigest()
            saved = json.loads(path.read_text())
            indexed = {(x['M'], x['law'], x['turns'], x['direction']): x for x in saved['cases']}
            for m in (4, 8):
                for law in ('native', 'swap'):
                    for turns in (1, 2):
                        for direction in (1, -1):
                            result = case(gram, z, position, k, m, reference, law, turns, direction)
                            old = indexed[m, law, turns, direction]
                            errors = [difference_upper(result['visibility'], old['visibility_direct']),
                                      difference_upper(result['pair_off_band'], old['pair_off_band'])]
                            for stride, probabilities in result['records'].items():
                                errors += [difference_upper(v, old['pair_records'][stride][name])
                                           for name, v in probabilities.items()]
                            error = max(errors)
                            require(error < 1e-10, 'source-center comparison')
                            max_error = max(max_error, error)
                            row = {a: result[a] for a in ('k', 'M', 'reference', 'law', 'turns', 'direction')}
                            for name in ('visibility', 'single_off_band', 'pair_off_band', 'selected_clock_probability', 'independent_pair_probability'):
                                row[name] = packed(result[name])
                            row['pair_records'] = {s: {name: packed(v) for name, v in probs.items()} for s, probs in result['records'].items()}
                            row['source_center_error_upper'] = error
                            rows.append(row)
        print(json.dumps({'progress_k': k, 'cases': len(rows)}), flush=True)
    return {'schema': 'native-source-interval-certificate-v1', 'passed': True,
            'known_result_verification': True, 'mpmath_version': mpmath.__version__, 'interval_dps': iv.dps,
            'source_files': source_hashes, 'k_values': ks, 'references': references, 'case_count': len(rows),
            'max_source_center_error_upper': max_error, 'cases': rows, 'seconds': time.perf_counter()-started,
            'physical_error_model': 'Exact declared native apparatus; certified arithmetic enclosure. Additional implementation perturbations require the separately stated conditional error theorem.',
            'external_G6': False, 'physical_award': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', nargs='+', type=int, default=[2])
    parser.add_argument('--references', nargs='+', choices=['matched', 'stationary'], default=['matched', 'stationary'])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--record-to', type=Path)
    group.add_argument('--verify', type=Path)
    args = parser.parse_args()
    require(all(0 <= k <= 6 for k in args.k), 'clock range')
    if args.record_to is not None:
        require(not args.record_to.exists(), 'refusing to overwrite evidence')
    result = run(args.k, args.references)
    if args.verify is not None:
        saved = json.loads(args.verify.read_text())
        left, right = dict(result), dict(saved)
        left.pop('seconds'); right.pop('seconds')
        require(left == right, 'saved certificate differs')
    else:
        with args.record_to.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('cases', 'source_files')}), flush=True)


if __name__ == '__main__':
    main()
