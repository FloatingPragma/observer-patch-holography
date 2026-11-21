"""Certify a finite-bank clock's mean rate and replay its capacity boundaries on the native sphere."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import sqrt
from pathlib import Path

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
TOL = Fraction(1, 2000)
ITINERARY = (8, 16, 31, 32, 2, 1)


def mean_rate(q, bits):
    if q == 1:
        return Fraction(0)
    if q >= 2**bits:
        return Fraction(1)
    pool = 1 << (q-1).bit_length()
    match = (1 - Fraction(pool-q, pool)**clock.ADDRESS_TRIES) / q
    live = {-1: Fraction(1, 2), 1: Fraction(1, 2)}
    returned = Fraction(0)
    for _ in range(1, clock.STEPS):
        following = {}
        for x, weight in live.items():
            branch = weight * match / 2
            for dest in (x-1, x+1):
                if dest == 0:
                    returned += branch
                else:
                    following[dest] = following.get(dest, 0) + branch
        live = following
    return 1 - returned / q


def within_target(q, rate):
    squared_target = 1 - Fraction(1, q*q)
    return max(0, rate-TOL)**2 <= squared_target <= (rate+TOL)**2


def profile(bits):
    overflow = 2**bits
    rows = []
    for q in range(1, overflow+1):
        rate = mean_rate(q, bits)
        float_rate = float(rate)
        if abs(clock.interval_response(q, age_bits=bits)['finite_rate'] - float_rate) > 1e-14:
            raise ArithmeticError('finite-walk certificate and registered kernel disagree')
        rows.append(dict(q=q, mean_rate=float_rate, passes=within_target(q, rate),
                         absolute_error=abs(float_rate - sqrt(1-1/q**2))))
    # Beyond the overflow threshold the computed rate is exactly one; the discrepancy decreases with q.
    rejected = [r['q'] for r in rows[:-1] if not r['passes']]
    q = overflow
    while not within_target(q, Fraction(1)):
        rejected.append(q)
        q += 1
    worst = max(rows, key=lambda r: r['absolute_error'])
    return dict(bits=bits, first_overflow=overflow, first_passing_tail_interval=q,
                outside_tolerance=rejected, worst_interval=worst['q'],
                max_absolute_rate_error=round(worst['absolute_error'], 12),
                first_tail_error=round(rows[-1]['absolute_error'], 12))


def witness():
    boards = []
    for level in (0, 1, 2):
        for seed in (31, 47):
            row = clock.material_replay(federation(level), ITINERARY, age_bits=5, seed=seed)
            events = row['local_events']
            offset = 0
            for q in ITINERARY:
                for step, event in enumerate(events[offset:offset+q], start=1):
                    moved = step == q
                    if event['handoff'] != moved or event['age'] != (0 if moved else min(step, 31)):
                        raise ArithmeticError('the physical age/transfer record differs from the boundary test')
                offset += q
            counts = dict(home=sum(e['home_primary'] for e in events), traveller=sum(e['tick'] for e in events))
            if row['reunion'] != counts or row['home_retained_traveller'] != counts['traveller']:
                raise ArithmeticError('the material seam failed to communicate the two local counts')
            boards.append(dict(cells=federation(level).cells, seed=seed, moving_cells=row['moving_cells'],
                               transfer_cells=row['transfer_cells'], reunion=counts))
    return boards


def tape_controls():
    rows = []
    for bits in (3, 4, 5):
        for q in (8, 16, 31, 32):
            width = (q-1).bit_length()
            tape = iter([1, *map(int, format(q-1, f'0{width}b')), 0])
            accepted = clock.closing_tick(q, 1, tape.__next__, age_bits=bits)
            overflow = q >= 2**bits
            remaining = list(tape)
            if accepted != overflow or len(remaining) != (width+2 if overflow else 0):
                raise ArithmeticError('the closing mark failed the selected first-return/overflow control')
            rows.append(dict(bits=bits, q=q, accepted=accepted))
    return rows


def run():
    profiles = [profile(bits) for bits in (3, 4, 5)]
    if profiles[-1]['outside_tolerance']:
        raise ArithmeticError('the five-bit profile violates the frozen tolerance')
    return dict(protocol_sha256=hashlib.sha256((HERE/'CAPACITY_PROTOCOL.md').read_bytes()).hexdigest(),
                tolerance=str(TOL), itinerary=list(ITINERARY), profiles=profiles,
                closing_tape=tape_controls(), material=witness())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', action='store_true', help='create the capacity receipt on first measurement')
    result = run()
    path = HERE/'capacity.json'
    if parser.parse_args().record:
        with path.open('x', encoding='utf-8') as f:
            f.write(json.dumps(result, indent=2) + '\n')
    elif result != json.loads(path.read_text()):
        raise SystemExit('capacity replay differs from stored observations')
    for row in result['profiles']:
        print(row)
