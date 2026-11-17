"""Fresh complete native target prefix with the actual retained reservoir spectator."""
import math
import numpy as np
from native_source import native
EDGES = ((0, 1), (0, 1), (0, 2), (0, 2), (1, 2))

def columns(background, g, sign):
    target = g['target']
    fine = []
    rows = []
    R = np.eye(33)
    for e, edge in enumerate(EDGES):
        physical = tuple((target[x] for x in edge))
        for c in (0, 1):
            mats = np.zeros((3, 3, 3))
            for p in range(3):
                q = list(background)
                for v in target:
                    q[v] = 0
                q[target[p]] = sign
                out, rec = native(q, physical, c)
                j = rec[3]
                y = next((i for i, v in enumerate(target) if out[v] == sign))
                expected = q.copy()
                expected[target[p]] = 0
                expected[target[y]] = sign
                require(out == expected, 'spectator disturbed')
                mats[j + 1, y, p] = 1
                cc = c ^ bool(j)
                R[:, 11 * p + 2 * e + c] = 0
                R[11 * y + 2 * e + cc, 11 * p + 2 * e + c] = 1
                rows.append([p, e, c, j, y, int(cc)])
            fine.extend(mats)
    return (np.asarray(fine) / math.sqrt(10), R, rows)

def helmert(branches):
    out = [sum(branches) / math.sqrt(len(branches))]
    prefix = branches[0].copy()
    for k in range(1, len(branches)):
        out.append((prefix - k * branches[k]) / math.sqrt(k * (k + 1)))
        prefix += branches[k]
    return np.asarray(out)

def encode(x):
    return {'real': x.real.tolist(), 'imag': x.imag.tolist()} if np.iscomplexobj(x) else x.tolist()

def prefix(background, g, sign):
    require([background[v] for v in g['target']] == [0, sign, 0], 'scientific check failed')
    fine, R, rows = columns(background, g, sign)
    continuing = helmert(fine)
    blank = np.zeros_like(fine)
    blank_records = []
    for e in range(5):
        for c in (0, 1):
            out, rec = native(background, tuple(g['blank']), c)
            require(out == background and rec[3] == 0, 'scientific check failed')
            blank_records.append(rec)
            blank[6 * e + 3 * c + 1] = np.eye(3) / math.sqrt(10)
    reference = helmert(blank)
    contrast = np.asarray([(b + s * a) / 2 for a, b in zip(continuing, reference) for s in (1, -1)])
    a = np.r_[np.full(10, 1 / math.sqrt(10)), 0.0]
    A = np.kron(np.eye(3), a[:, None])
    W = R @ (2 * A @ A.T - np.eye(33))
    live = A.astype(complex)
    waves = []
    for _ in range(256):
        waves.append(live.copy())
        live = W @ live
    counter = np.fft.fft(np.asarray(waves), axis=0) / 256
    rho = np.diag([0.0, 1.0, 0.0])
    records = []
    probability = 1.0
    for name, ops, accept in [('contrast', contrast, [1])] + [('cooling', continuing, [0])] * 8 + [('counter', counter, list(range(1, 128)))]:
        before = rho.copy()
        outs = np.einsum('kxa,ab,kyb->kxy', ops, rho, ops.conj())
        probs = np.trace(outs, axis1=1, axis2=2).real
        selected = sum((outs[k] for k in accept))
        p = float(np.trace(selected).real)
        defect = float(np.max(abs(np.einsum('kba,kbc->ac', ops.conj(), ops) - np.eye(ops.shape[2]))))
        require(min(probs) > -1e-12 and abs(sum(probs) - 1) < 1e-10 and (defect < 1e-10), 'scientific check failed')
        records.append({'stage': name, 'input_density': encode(before), 'all_probabilities': probs.tolist(), 'accepted': accept, 'probability': p, 'completeness_defect': defect})
        rho = selected / p
        probability *= p
    return {'sign': sign, 'actual_source_before': background, 'source_columns': rows, 'blank_records': blank_records, 'R': R.tolist(), 'W': W.tolist(), 'P': (A.T @ R @ A).tolist(), 'records': records, 'clock_density': encode(rho), 'complete_prefix_probability': probability, 'meaning': 'Executed native cooling and complete finite counter on the actual extracted basis state; no full rotation/pair-family award transferred.'}

def require(condition, message='scientific check failed'):
    if not condition:
        raise AssertionError(message)
