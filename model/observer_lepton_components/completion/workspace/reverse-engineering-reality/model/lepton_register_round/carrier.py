"""Read counts on the carrier with the sixty addresses of the public program (section 6).

A cell is one carrier with its own loads, record and sixty-address register; the walk is the plain law of
`Observer.repair` with the public local completion at the addressed seam. Two cells are independent registers.
"""
from fractions import Fraction as F
from itertools import combinations

from model.coherent_seam_repair.program import Observer
from model.observer.native import native
from . import walk as wk

N_PORTS = 12


def carrier_walk():
    return wk.Walk(Observer().addresses)


def one_load_sector(w):
    """The walk on one load equals Observer.repair and Observer.reflect on all 1440 basis states."""
    obs = Observer()
    D = obs.address_count
    repair_ok = reflect_ok = True
    for port in range(12):
        for parity in (0, 1):
            for k in range(D):
                state = obs.basis_state(port, parity, k)
                key = (wk.config(port), parity, k)
                mine = w.repair({key: F(1)})
                theirs = {(wk.config(i // (2 * D)), (i // D) % 2, i % D): F(x)
                          for i, x in enumerate(obs.repair(state).amplitudes) if x != 0}
                repair_ok &= mine == theirs
                mine = w.reflect({key: F(1)})
                theirs = {(wk.config(i // (2 * D)), (i // D) % 2, i % D): F(x)
                          for i, x in enumerate(obs.reflect(state).amplitudes) if x != 0}
                reflect_ok &= mine == theirs
    return repair_ok, reflect_ok


def double_load(w, port):
    """One tick from a blank double at the port with the uniform address; weights and the native current."""
    start = w.uniform({wk.config(port, port): F(1)})
    ticked = w.step(start)
    out = {'blank': F(0), 'record 2': F(0), 'record 3': F(0)}
    current_ok, positive = True, True
    for (cfg, rec, k), amp in ticked.items():
        label = 'blank' if rec == 0 else f'record {rec}'
        out[label] = out.get(label, F(0)) + amp * amp / w.size
        a, b = w.addresses[k]
        before = wk.config(port, port)
        for coin in (0, 1):
            moved, (_, current) = native(before, (a, b), coin)
            expected = 0 if rec == 0 else 1 if rec == 2 else -1 if rec == 3 else None
            current_ok &= current == expected and moved == cfg
        if rec in (2, 3):
            positive &= amp == 1
    return out, current_ok, positive, ticked


def two_cells(w, port_x, port_y):
    """P(R+) = P(record 2 in X), P(H) = P(record of Y written), and the joint law over the 60 x 60 address pairs."""
    _, _, _, tx = double_load(w, port_x)
    _, _, _, ty = double_load(w, port_y)
    p_r = p_h = joint = total = F(0)
    for (cx, rx, kx), ax in tx.items():
        for (cy, ry, ky), ay in ty.items():
            weight = (ax * ay) ** 2 / (w.size * w.size)
            total += weight
            rp, hold = rx == 2, ry in (2, 3)
            p_r += weight * rp
            p_h += weight * hold
            joint += weight * (rp and hold)
    return dict(p_r=p_r, p_h=p_h, joint=joint, total=total, count=p_r + p_h - 2 * joint)


def native_pair(w, ports):
    """Weights of a cell with one load on each of two ports: untouched, agree, repair with the load first or second."""
    ticked = w.step(w.uniform({wk.config(*ports): F(1)}))
    out = {'untouched': F(0), 'agree': F(0), 'repair, load first': F(0), 'repair, load second': F(0), 'other': F(0)}
    signs = {'repair, load first': set(), 'repair, load second': set()}
    for (cfg, rec, k), amp in ticked.items():
        a, b = w.addresses[k]
        touched = (a in ports) + (b in ports)
        if touched == 0:
            label = 'untouched'
        elif touched == 1 and rec == 0:
            label = 'agree'
        elif touched == 1 and rec == 1:
            label = 'repair, load first' if a in ports else 'repair, load second'
            signs[label].add(1 if amp > 0 else -1)
        else:
            label = 'other'
        out[label] += amp * amp / w.size
    return out, {k: sorted(v) for k, v in signs.items()}


def neighbours():
    obs = Observer()
    return {p: {q for e in obs.seams for q in e if p in e and q != p} for p in range(12)}


def one_register(w):
    """Entries of the completion, branch weights of every definite two-load configuration, and two doubles in one
    register on non-adjacent ports."""
    entries = sorted({str(x) for row in wk.completion() for x in row})
    multiples = True
    for a in range(12):
        for b in range(a, 12):
            ticked = w.step(w.uniform({wk.config(a, b): F(1)}))
            multiples &= all((amp * amp / w.size * 240).denominator == 1 for amp in ticked.values())
    nb = neighbours()
    rows = set()
    for a, b in combinations(range(12), 2):
        if b in nb[a]:
            continue
        ticked = w.step(w.uniform({wk.config(a, a, b, b): F(1)}))
        p_r = sum((amp * amp for (cfg, rec, k), amp in ticked.items() if rec == 2 and cfg[a] == 1), F(0)) / w.size
        p_h = sum((amp * amp for (cfg, rec, k), amp in ticked.items() if rec in (2, 3) and cfg[b] == 1), F(0)) / w.size
        joint = sum((amp * amp for (cfg, rec, k), amp in ticked.items()
                     if rec == 2 and cfg[a] == 1 and rec in (2, 3) and cfg[b] == 1), F(0)) / w.size
        rows.add((p_r, p_h, joint))
    return entries, multiples, rows


def weak_round(w, port):
    """The weak round on two cells with doubles at the port: X read for record 2, Y read for a written record."""
    start = w.uniform({wk.config(port, port): F(1)})
    n2 = F(w.size)
    ws = w.step(start)
    a = w.inverse_step(wk.restrict(ws, lambda k: k[1] == 2))
    b = w.inverse_step(wk.restrict(ws, lambda k: k[1] in (2, 3)))
    s_a, s_b = wk.dot(start, a) / n2, wk.dot(start, b) / n2
    pair = wk.tensor(start, start)
    flagged = {k: v for k, v in wk.tensor(ws, ws).items() if (k[0][1] == 2) != (k[1][1] in (2, 3))}
    X = wk.on_both(flagged, w.inverse_step)
    overlap, norm = wk.dot(pair, X) / (n2 * n2), wk.norm2(X) / (n2 * n2)

    def content(v, s_):
        perp = wk.combine(v, start, -s_)
        off = sum((x * x for (cfg, rec, k), x in perp.items() if rec != 0 or cfg != wk.config(port, port)), F(0)) / n2
        tot = sum(perp.get((wk.config(port, port), 0, k), F(0)) for k in range(w.size))
        return dict(norm2=wk.norm2(perp) / n2, off=off, uniform=tot * tot / w.size / n2)
    return dict(inverse=w.inverse_step(ws) == start, s_a=s_a, s_b=s_b, overlap=overlap, norm=norm,
                A=content(a, s_a), B=content(b, s_b))


def result():
    w = carrier_walk()
    repair_ok, reflect_ok = one_load_sector(w)
    per_port = {p: double_load(w, p)[:3] for p in range(12)}
    weights_ok = all(x[0] == {'blank': F(5, 6), 'record 2': F(1, 12), 'record 3': F(1, 12)} for x in per_port.values())
    current_ok = all(x[1] for x in per_port.values())
    positive_ok = all(x[2] for x in per_port.values())
    cells = {p: two_cells(w, p, p) for p in range(12)}
    c2 = cells[2]
    n = N_PORTS
    nb = neighbours()
    pairs = [(a, b) for a, b in combinations(range(12), 2) if b not in nb[a]]
    np_rows, np_signs = set(), set()
    for ports in pairs:
        wt, sg = native_pair(w, ports)
        np_rows.add(tuple(sorted((k, v) for k, v in wt.items())))
        np_signs.add(tuple((k, tuple(v)) for k, v in sorted(sg.items())))
    np_expected = {'untouched': F(2, 3), 'agree': F(1, 6), 'repair, load first': F(1, 12),
                   'repair, load second': F(1, 12), 'other': F(0)}
    entries, multiples, one_reg = one_register(w)
    rounds = {p: weak_round(w, p) for p in (2, 10)}
    r2 = rounds[2]
    checks = dict(
        walk_equals_observer_repair_and_reflect_on_1440_basis_states=repair_ok and reflect_ok,
        double_load_blank_5_6_record_2_and_3_1_12_at_every_port=weights_ok,
        records_2_and_3_are_the_native_current_signs=current_ok,
        emission_amplitudes_all_positive=positive_ok,
        two_cells_read_1_12_1_6_joint_1_72_at_every_port=all(
            (c['p_r'], c['p_h'], c['joint'], c['total']) == (F(1, 12), F(1, 6), F(1, 72), 1) for c in cells.values()),
        two_cell_count_is_2_9_equal_to_3_over_n_minus_4_over_n2=(
            all(c['count'] == F(2, 9) for c in cells.values()) and F(3, n) - F(4, n * n) == F(2, 9)),
        native_pairs_give_2_3_1_6_1_12_1_12=np_rows == {tuple(sorted(np_expected.items()))} and len(pairs) == 36,
        second_endpoint_repair_has_the_minus_sign=np_signs == {(('repair, load first', (1,)), ('repair, load second', (-1,)))},
        completion_entries_are_0_1_and_plus_minus_one_half=entries == ['-1/2', '0', '1', '1/2'],
        one_register_weights_are_multiples_of_1_240=multiples and (F(2, 9) * 240).denominator != 1,
        one_register_two_doubles_count_1_4=one_reg == {(F(1, 12), F(1, 6), F(0))},
        weak_round_on_sixty_addresses_has_P_2_9=all(
            r['inverse'] and r['s_a'] == F(1, 12) and r['s_b'] == F(1, 6) and r['overlap'] == r['norm'] == F(2, 9)
            for r in rounds.values()),
        sixty_address_non_returns_11_144_and_5_36_need_only_the_address_consensus=all(
            r['A']['norm2'] == F(11, 144) and r['B']['norm2'] == F(5, 36) and r['A']['off'] == r['B']['off'] == 0
            and r['A']['uniform'] == r['B']['uniform'] == 0 for r in rounds.values()))
    data = dict(
        double_load_at_port_2={k: str(v) for k, v in per_port[2][0].items()},
        two_cells_at_port_2={k: str(v) for k, v in c2.items()},
        count_formula='3/n - 4/n^2 at n = 12',
        native_pair_weights={k: str(v) for k, v in np_expected.items()}, native_pairs_checked=len(pairs),
        completion_entries=entries,
        one_register_two_doubles=[{'P(R+)': str(a), 'P(H)': str(b), 'P(R+ and H)': str(c), 'count': str(a + b - 2 * c)}
                                  for a, b, c in sorted(one_reg)],
        two_nine_times_240=str(F(2, 9) * 240),
        weak_round_port_2={'<D|W^-1 Pi_R+ W|D>': str(r2['s_a']), '<D|W^-1 Pi_H W|D>': str(r2['s_b']),
                           '<start|X>': str(r2['overlap']), '<X|X>': str(r2['norm']),
                           'non_return_X': {k: str(v) for k, v in r2['A'].items()},
                           'non_return_Y': {k: str(v) for k, v in r2['B'].items()}})
    return data, checks
