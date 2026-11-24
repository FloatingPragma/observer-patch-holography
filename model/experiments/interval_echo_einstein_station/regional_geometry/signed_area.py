"""Two native rulers at common fixed volumes, with sharp source-word measurements."""
import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import isqrt

import mpmath as mp
import numpy as np
import sympy as sp

from . import both_volumes, measure

ROOT = measure.HERE
CORNERS = [(0, 0), (0, 1), (1, 0), (1, 1)]
# Free corner velocities from the sum and difference of unit projected gradients on 20 cells.
FREE = {"both_increase": ("-.030310", "-.461984"), "opposite": ("-.521716", ".419664")}
STEP = Fraction(1, 32)
SHOTS = 2**23
BASIS = []
for i, j in CORNERS:
    a = sp.zeros(3)
    a[i, j] = a[2, 2] = 1
    a[2, j] = a[i, 2] = -1
    BASIS.append(a)
REF = sp.Matrix([[4, 6, 2], [2, 4, 6], [6, 2, 4]])/36


def gram(distances, sites):
    first, *rest = sites
    def sq(i, j):
        return 0 if i == j else distances[tuple(sorted((i, j)))]**2
    return sp.Matrix([[(sq(first, i)+sq(first, j)-sq(i, j))/2 for j in rest] for i in rest])


def conditional_ranges(raw, menu):
    result = {}
    for k, pair in enumerate(measure.PAIRS):
        if max(pair) >= 4:
            continue
        result[pair] = tuple(sp.Matrix(3, 3, lambda i,j: sp.Rational(raw[i][j][k][menu][m]))
                             for m in range(2))
    return result


def geometry_gradient(law):
    # All range quotients and all log-volume derivatives are rational on the deposited moments.
    distances = {}
    slopes = [{} for _ in BASIS]
    for pair, (success, hand) in law.items():
        chance = sum(REF.multiply_elementwise(success))
        d = sum(REF.multiply_elementwise(hand))/chance
        distances[pair] = d
        for slope, transfer in zip(slopes, BASIS, strict=True):
            slope[pair] = sum(transfer.multiply_elementwise(hand-d*success))/chance
    parts = []
    for sites in [tuple(range(4)), *combinations(range(4), 3)]:
        g = gram(distances, sites)
        det = g.det()
        if not g.is_positive_definite:
            raise ValueError("The supplied ranges do not define a positive simplex")
        inv = g.inv()
        row = []
        for slope in slopes:
            # The Gram matrix is quadratic in the ranges. Polarization gives its exact derivative.
            plus = {k:d+slope[k] for k,d in distances.items()}
            minus = {k:d-slope[k] for k,d in distances.items()}
            dg = (gram(plus, sites)-gram(minus, sites))/2
            row.append((inv*dg).trace()/2)
        parts.append((det, sp.Matrix([row])))
    return parts[0][1], parts[1:]


def root_interval(x):
    # Adjacent dyadic bounds on sqrt(x); all endpoint and sign operations remain rational.
    scale = 1 << 88
    floor = isqrt((int(x.p)*scale**2)//int(x.q))
    return sp.Rational(floor, scale), sp.Rational(floor+1, scale)


def area_interval(faces, tangent):
    lower = upper = 0
    for det, dlog_area in faces:
        a, b = root_interval(det)
        slope = (dlog_area*tangent)[0]/2
        first, second = a*slope, b*slope
        lower += min(first, second)
        upper += max(first, second)
    if lower <= 0 <= upper:
        raise ValueError("Area derivative sign is unresolved")
    return [float(np.nextafter(float(lower), -np.inf)), float(np.nextafter(float(upper), np.inf))]


def charge_interval(delta):
    # The constant in -log(36 p) vanishes because delta has zero sum.
    lo = hi = 0
    for n in (2, 4, 6):
        t = sp.Rational(n-1, n+1)
        left = 2*sum(t**(2*j+1)/(2*j+1) for j in range(64))
        right = left + 2*t**129/(129*(1-t**2))
        coefficient = -sum(d for d,p in zip(delta,REF,strict=True) if 36*p == n)
        a,b = coefficient*left, coefficient*right
        lo += min(a,b)
        hi += max(a,b)
    if lo <= 0 <= hi:
        raise ValueError("The reference charge sign is unresolved")
    return [float(np.nextafter(float(lo), -np.inf)), float(np.nextafter(float(hi), np.inf))]


def floating(x):
    return mp.mpf(int(x.p))/int(x.q)


def finite_reader(law):
    pairs = sorted(law)
    banks = [[mp.matrix([floating(v) for v in arr]) for arr in law[pair]] for pair in pairs]
    def read(state):
        d = {pair:mp.fdot(state, hand)/mp.fdot(state, success)
             for pair, (success, hand) in zip(pairs, banks, strict=True)}
        def simplex(sites):
            o, *points = sites
            def sq(i, j):
                return 0 if i == j else d[tuple(sorted((i, j)))]**2
            g = mp.matrix([[(sq(o,i)+sq(o,j)-sq(i,j))/2 for j in points] for i in points])
            return mp.sqrt(mp.det(g))/mp.factorial(len(points))
        return simplex(range(4)), sum(simplex(sites) for sites in combinations(range(4), 3))
    return read


def finite_arms(readers, tangent):
    reference = mp.matrix([floating(x) for x in REF])
    basis = [mp.matrix([int(a) for a in x]) for x in BASIS]
    slopes = [floating(x) for x in tangent]
    targets = [read(reference)[0] for read in readers]
    arms = []
    for step in (-mp.mpf(1)/32, mp.mpf(1)/32):
        def state(x, y):
            return reference + x*basis[0]+y*basis[1]+step*(slopes[2]*basis[2]+slopes[3]*basis[3])
        def equations(x, y):
            p = state(x,y)
            return tuple(read(p)[0]/target-1 for read,target in zip(readers, targets, strict=True))
        x,y = mp.findroot(equations, (step*slopes[0],step*slopes[1]), tol=mp.mpf('1e-44'))
        p = state(x,y)
        if min(p) <= 0:
            raise ValueError("finite transfer leaves the preparation family")
        arms.append((p, [read(p) for read in readers]))
    minus, plus = arms
    delta = plus[0]-minus[0]
    k = [-mp.log(9*x) for x in reference]
    secants = [float((hi[1]-lo[1])*16) for hi,lo in zip(plus[1],minus[1], strict=True)]
    return dict(probabilities=[[float(x) for x in side[0]] for side in arms],
                max_fractional_volume_error=float(max(abs(records[i][0]/v-1)
                    for _,records in arms for i,v in enumerate(targets))),
                area_secants=secants, fixed_reference_charge_secant=float(16*mp.fdot(delta,k)))


def measure_charge(arms, rng):
    ref = np.array(REF, dtype=float).ravel()
    reference = rng.multinomial(SHOTS, ref)
    sigma = (reference/reference.sum()).reshape(3,3)
    score = np.log(np.outer(sigma.sum(axis=1),sigma.sum(axis=0))/sigma)
    counts = [rng.multinomial(SHOTS, p) for p in arms]
    a,b = [(n/SHOTS).reshape(3,3) for n in counts]
    slope = 16*(b-a)
    z = -slope/sigma + (slope.sum(axis=1)/sigma.sum(axis=1))[:,None] + slope.sum(axis=0)/sigma.sum(axis=0)
    def variance(p,x):
        return np.sum(p*(x-np.sum(p*x))**2)/SHOTS
    se = np.sqrt(256*(variance(a,score)+variance(b,score))+variance(sigma,z))
    return dict(shots_per_bank=SHOTS, reference_counts=reference.tolist(),
                source_counts=[n.tolist() for n in counts], charge_secant=float(np.sum(score*slope)),
                counting_standard_error=float(se))


def run():
    bank = json.loads((ROOT/'results.json').read_text(), parse_float=Fraction)
    boards, frozen = [], {}
    rng = np.random.default_rng(367024)
    with mp.workdps(60):
        for board in bank['screens']:
            laws = [conditional_ranges(board['word_moments'], i) for i in range(2)]
            gradients = [geometry_gradient(law) for law in laws]
            volumes = gradients[0][0].col_join(gradients[1][0])
            pivot = volumes[:, :2]
            if pivot.det() == 0:
                raise ValueError("volume constraints do not admit the two free coordinates")
            readers = [finite_reader(law) for law in laws]
            ref = np.array(REF, dtype=float)
            floats = np.array(board['word_moments'], dtype=float)
            base = both_volumes.cycle_basis()
            op = np.stack([both_volumes.gradients(floats,ref,d) for d in base], axis=1)
            kernel = np.linalg.svd(op[:2]/np.linalg.norm(op[:2],axis=1)[:,None])[2][2:].T
            a = op[2:] @ kernel
            unit = a/np.linalg.norm(a,axis=1)[:,None]
            rows = {}
            for label, free in FREE.items():
                v = sp.Matrix([sp.Rational(x) for x in free])
                tangent = (pivot.inv()*(-volumes[:,2:]*v)).col_join(v)
                assert volumes*tangent == sp.zeros(2,1)
                delta = sum((x*d for x,d in zip(tangent,BASIS,strict=True)),sp.zeros(3))
                assert delta*sp.ones(3,1) == delta.T*sp.ones(3,1) == sp.zeros(3,1)
                dlog = -sum(floating(x)*mp.log(9*floating(p)) for x,p in zip(delta,REF,strict=True))
                row = dict(area_slope_enclosures=[area_interval(faces,tangent) for _,faces in gradients],
                           charge_slope=float(dlog), charge_slope_enclosure=charge_interval(delta),
                           finite=finite_arms(readers,tangent))
                if label not in frozen:
                    frozen[label] = tangent
                    row['sharp_measurement'] = measure_charge(row['finite']['probabilities'],rng)
                row['transferred_volume_slopes'] = [float(x) for x in volumes*frozen[label]]
                row['transferred_area_enclosures'] = [area_interval(faces,frozen[label]) for _,faces in gradients]
                rows[label] = row
            boards.append(dict(cells=board['cells'], fixed_volume_area_cosine=float(unit[0]@unit[1]), witnesses=rows))
    return dict(geometry='first-return, free-pulse tetrahedra', finite_step=float(STEP), boards=boards)


if __name__ == '__main__':
    from ..collar_pairing.measure import compare
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    result = run()
    dest = ROOT/'signed_area.json'
    if parser.parse_args().write:
        dest.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    else:
        compare(result,json.loads(dest.read_text()))
    for board in result['boards']:
        print(board['cells'], [(k, v['finite']['area_secants'],v['finite']['fixed_reference_charge_secant'])
                               for k,v in board['witnesses'].items()])
