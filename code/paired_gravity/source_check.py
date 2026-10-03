"""Independent primitive Gram reconstruction, exact local replay and barrier."""
from fractions import Fraction as F
import importlib.util
import sys
from pathlib import Path
from .format import need, keys, equal, rational, digest

# This verifier is standalone; never import the producer's arithmetic module.
_spec = importlib.util.spec_from_file_location('_paired_gravity_independent_scalar',
    Path(__file__).resolve().parents[2]/'code/source_scalar_execution/verify_source_scalar_execution.py')
parent = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = parent
_spec.loader.exec_module(parent)


def model():
    mass, original, _, _, addresses = parent.model()
    R = parent.R
    rows = [dict(row) for row in original]
    for i, row in enumerate(rows):
        row[i] -= 1
    # Re-prove the maximum-principle hypotheses on every edge and boundary.
    reached = {0}
    for _ in range(64):
        reached |= {j for i in reached for j in rows[i]}
    need(len(reached) == 64, 'connected source support')
    boundary = False
    barrier = []
    for item in addresses:
        xyz = list(map(parent.parse, item['coordinate_Qphi']))
        b = sum((x*(1-x) for x in xyz), R())/6
        need(b.sign() > 0 and (b-F(1, 8)).sign() <= 0, 'positive bounded barrier')
        barrier.append(b)
    for i, row in enumerate(rows):
        need(row[i].sign() > 0 and mass[i].sign() > 0, 'positive diagonal and action measure')
        boundary |= sum(row.values(), R()).sign() > 0
        need(sum(row.values(), R()).sign() >= 0, 'Dirichlet diagonal dominance')
        for j, a in row.items():
            if j != i:
                need(a.sign() < 0 and mass[i]*a == mass[j]*rows[j][i], 'same symmetric edge action')
        need((sum((a*barrier[j] for j, a in row.items()), R())-1).sign() >= 0, 'barrier dominates unity')
    need(boundary, 'nonempty grounded boundary')
    return rows


def rounded(value, denominator):
    # Integer bisection in Q(sqrt(5)), independent of the producer's float seed.
    R = parent.R
    lo, hi = -denominator, denominator
    need((value+1).sign() >= 0 and (value-1).sign() < 0, 'bounded executed field')
    while hi-lo > 1:
        m = (lo+hi)//2
        if (value-F(m, denominator)).sign() >= 0:
            lo = m
        else:
            hi = m
    result = F(lo, denominator)
    need((value-result).sign() >= 0 and (value-result-F(1, denominator)).sign() < 0, 'rounding enclosure')
    return result


def exponential_minus(x):
    need(0 <= x <= F(1, 100), 'nonlinear field domain')
    value = term = F(1)
    for n in range(1, 29):
        term *= -x/n
        value += term
    return value+term*(-x)/29, value


def nonlinear_rounded(a, off, strength):
    lo, hi = 0, (2**36)//100
    def sign(n):
        left, right = exponential_minus(F(n, 2**36))
        lower = a*F(n, 2**36)+off-strength*right
        upper = a*F(n, 2**36)+off-strength*left
        need(lower.sign() > 0 or upper.sign() < 0, 'resolved nonlinear update enclosure')
        return 1 if lower.sign() > 0 else -1
    need(sign(lo) < 0 < sign(hi), 'nonlinear root bracket')
    while hi-lo > 1:
        m = (lo+hi)//2
        if sign(m) > 0:
            hi = m
        else:
            lo = m
    return F(lo, 2**36)


def verify_case(item, strength, reverse, seed, nonlinear, rows):
    keys(item, 'strength reverse seed nonlinear sweeps denominator u w projection chains events reads error_intervals')
    equal([item['strength'], item['reverse'], item['seed'], item['nonlinear'], item['sweeps'], item['denominator']],
          [str(strength), reverse, seed, nonlinear, 32, 2**36], 'complete source protocol')
    R = parent.R
    order = sorted(range(64), reverse=reverse)
    state = {i: [0, F(0), strength/16 if seed and i == 42 else F(0)] for i in range(64)}
    hashes = {str(g): '0'*64 for g in (-1, 0, 1, 2)}
    projected = '0'*64
    event_count = read_count = 0
    for sweep in range(33):
        for i in order:
            inputs = [] if sweep == 0 else [[j, *[state[j][0]], str(state[j][1]), str(state[j][2])] for j in sorted(rows[i])]
            if sweep == 0:
                nu, nw = state[i][1:]
            else:
                force = [sum((a*state[j][k] for j, a in rows[i].items() if j != i), R()) for k in (1, 2)]
                nu = nonlinear_rounded(rows[i][i], force[0], strength) if nonlinear and i == 21 else rounded((R(strength if i == 21 else 0)-force[0])/rows[i][i], 2**36)
                nw = rounded(-force[1]/rows[i][i], 2**36)
            event = dict(site=i, version=state[i][0]+1, reads=inputs, u=str(nu), w=str(nw))
            projected = digest(dict(parent=projected, **event))
            for g in hashes:
                v_reads = [str(int(g)*F(a[2])+F(a[3])) for a in inputs]
                hashes[g] = digest(dict(parent=hashes[g], gamma=int(g), v=str(int(g)*nu+nw), v_reads=v_reads, **event))
            state[i] = [state[i][0]+1, nu, nw]
            event_count += 1
            read_count += len(inputs)
    equal([item['u'], item['w']], [[str(state[i][k]) for i in range(64)] for k in (1, 2)], 'every final source record')
    equal([item['projection'], item['chains'], item['events'], item['reads']],
          [projected, hashes, event_count, read_count], 'all writers, reads and constitutive branches')
    need(len(set(hashes.values())) == 4, 'physical v records distinguish constitutive laws')
    need(type(item['error_intervals']) is list and len(item['error_intervals']) == 2, 'both field error bounds')
    for k, interval in enumerate(item['error_intervals'], 1):
        need(type(interval) is list and len(interval) == 2, 'residual interval shape')
        lo, hi = map(rational, interval)
        need(0 <= lo < hi < F(1, 10**6) and hi-lo <= F(2, 10**16), 'useful outward residual bound')
        maximum, minimum = R(), R()
        for i, row in enumerate(rows):
            raw = sum((a*state[j][k] for j, a in row.items()), R())
            force = (exponential_minus(state[i][k]) if nonlinear else (F(1), F(1))) if k == 1 and i == 21 else (F(0), F(0))
            a, b = raw-strength*force[1], raw-strength*force[0]
            aa, bb = (-a if a.sign() < 0 else a), (-b if b.sign() < 0 else b)
            low, high = (aa, bb) if (aa-bb).sign() <= 0 else (bb, aa)
            if a.sign() <= 0 <= b.sign():
                low = R()
            if (high-maximum).sign() > 0:
                maximum = high
            if (low-minimum).sign() > 0:
                minimum = low
        need((minimum/8-lo).sign() >= 0 and (hi-maximum/8).sign() >= 0, 'independent exact maximum-principle enclosure')
        need(maximum.sign() > 0 if k == 1 or seed else maximum.sign() == 0,
             'finite residual and invariant zero-sector control')


def verify(items):
    need(type(items) is list and len(items) == 16, 'two strengths, both schedules, both seeds and both source laws')
    rows = model()
    for item, (s, reverse, seed, nonlinear) in zip(items, [(s, r, z, n) for s in (F(1, 32), F(1, 16))
                  for r in (False, True) for z in (False, True) for n in (False, True)]):
        verify_case(item, s, reverse, seed, nonlinear, rows)
