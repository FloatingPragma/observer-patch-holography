"""Uniform finite and continuous weak-detector tests at fixed six-tick deadline.

Split the enrolled Hamiltonian into its within-cell terms and V, the exchange
between different cells. Each missed inspection multiplies the wave by
I+(sqrt(1-p)-1)P at its target. This is a contraction preserving cell blocks.
In continuous monitoring the generator acquires -P/2, also within one cell.
Both full and uncoupled propagators contract. The Duhamel expansion with a
remainder after d insertions of V therefore has remainder norm <=(t||V||)^d/d!.
This works across any fixed set of discrete missed inspections: compose their
contractions between the propagators in the time-ordered integrals. If the
stations are d or more cell hops apart, the terms with fewer than d insertions
vanish between the two stations. No sharp projection on misses is used.

For a continuous relay at s and return after t the corresponding density is
at most gamma² ||V||^(4d) (st)^(2d)/(d!)^4. This holds even on entangled
head/strain states and after conditioning on the first click. Integrating on
s+t<=T gives the two Dirichlet integrals in continuum(). At spacing h the
joint probability at the two ordered appointments has the same bound with
p² in place of gamma². The code uses p=1-exp(-gamma h)<=gamma h.

The six stations below are an additional detector challenge on the same live
80-cell source. Their placement uses cell incidence only. Certifying almost
equal timeout ranges addresses space inferred from these actual readings;
the clock deadline and frozen weak detector remain apparatus choices.
"""

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import ceil, factorial
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.special import jv

from .. import live_radar


T = 6


def station_challenge(board):
    owners = np.asarray(board.external) // 12
    adj = coo_matrix((np.ones(len(owners)), owners.T), shape=(board.cells, board.cells))
    hop = shortest_path(adj+adj.T, directed=False, unweighted=True).astype(int)

    def first_six(chosen, options):
        if len(chosen) == 6:
            return chosen
        for j, cell in enumerate(options):
            if len(chosen) + len(options) - j < 6:
                break
            remaining = [c for c in options[j+1:] if hop[cell, c] >= 5]
            result = first_six([*chosen, cell], remaining)
            if result is not None:
                return result
        return None

    cells = first_six([], list(range(board.cells)))
    if cells is None:
        raise ValueError("No six separated cells on the supplied board")
    return [12*c for c in cells], hop[np.ix_(cells, cells)]


def intercell_norm(board, charge):
    route = live_radar.source.formation.region(board)
    integer_h = live_radar.source.system(board, route, charge, coupled=True)[0].tocoo()
    cell = (np.arange(integer_h.shape[0]) % (4*board.n)) // 48
    which = cell[integer_h.row] != cell[integer_h.col]
    a, b = integer_h.row[which], integer_h.col[which]
    weights = integer_h.data[which]
    # Combining duplicate sparse entries and dropping explicit zeros makes
    # the matching test count actual links. A symmetric matching has exact
    # operator norm equal to its greatest absolute entry.
    v = coo_matrix((weights,(a,b)), shape=integer_h.shape).tocsr()
    v.eliminate_zeros()
    assert (v != v.T).nnz == 0 and np.diff(v.indptr).max() == 1
    return Fraction(int(abs(v.data).max()),40)


def continuum(d, v):
    coefficient = v**(4*d) / factorial(d)**4
    joint = coefficient * T**(4*d+2)*factorial(2*d)**2/factorial(4*d+2)
    # Score rebate=(T-s-t)/2; its exact simplex integral adds T/[2(4d+3)].
    rebate = joint*Fraction(T, 2*(4*d+3))
    return min(Fraction(1), joint), Fraction(T,2)-rebate


def discrete(d, v, reciprocal_step):
    n = T*reciprocal_step
    coefficient = v**(4*d) / (factorial(d)**4 * reciprocal_step**(4*d+2))
    joint = rebate = 0
    for k in range(1, n):
        for m in range(1, n-k+1):
            term = (k*m)**(2*d)
            joint += term
            rebate += (n-k-m)*term
    return min(Fraction(1), coefficient*joint), Fraction(T,2)-coefficient*rebate/(2*reciprocal_step)


def geometry(low, high):
    # Six Euclidean R4 points must have distance aspect squared <=2/3:
    # on 1-perp their centered squared Gram has five eigenvalues bounded
    # below by (6l²-4u²)/4 (row-sum perturbation from equal distances).
    # S³ chords give at least the geodesic aspect by concavity of sin.
    correction = next(Fraction(c,100) for c in range(30,-1,-1)
                      if 3*(low-Fraction(c,100))**2 > 2*(high+Fraction(c,100))**2)
    assert 3*(low-correction)**2 > 2*(high+correction)**2
    # Five points in R3 instead require aspect squared <=3/5. In a
    # strongly convex normal ball with |sec|<=K, x=sqrt(K)R<pi, the
    # geodesic distances relative to tangent-space distances are between
    # sin(x)/x and sinh(x)/x by Rauch comparison. The first factor uses
    # a minimizing geodesic remaining in the ball, the second a straight
    # chord in its tangent coordinates. sinh(x)/sin(x) increases on (0,pi).
    for i in range(87,-1,-1):
        x = Fraction(i,100)
        sin_lower = x-x**3/6+x**5/120-x**7/5040
        sinh_upper = x+x**3/6+x**5/120+x**7/5040 + x**9/362880/(1-x*x/110)
        if 5*low**2*sin_lower**2 > 3*high**2*sinh_upper**2:
            break
    assert 5*low**2*sin_lower**2 > 3*high**2*sinh_upper**2
    return dict(round_three_sphere_max_distance_change_strict_lower=str(correction),
                convex_normal_ball_abs_sectional_curvature_times_radius_squared_strict_lower=str(x*x))


def panel(complete, low):
    return dict(completion_upper=str(complete), range_lower=str(low), range_upper="3",
                geometry=geometry(low, Fraction(3)))


SAMPLE_RATE = 64
DEGREE = 260
COEFF_BITS, MOMENT_BITS, PRODUCT_BITS, OUT_BITS = 42, 48, 23, 19


def checked_coefficients():
    """Propose dyadic Chebyshev coefficients, certify by their polynomial ODE.

    For p on [-1,1], |p(s)-exp(-izs)| <= |p(0)-1|+sup|p'+izp|.
    Differentiation in the Chebyshev basis is exact; its l1 residual bounds
    the sup norm. Bessel library errors consequently cannot pass unchecked.
    z=(128/5)t=(2/5)k on the grid t=k/64.
    """
    scale = 1 << COEFF_BITS
    n = np.arange(DEGREE+1)
    ticks = np.arange(T*SAMPLE_RATE+1)
    proposals = np.rint(jv(n[None,:], (2./5)*ticks[:,None])
                         *np.where(n == 0,1.,2.)*scale).astype(np.int64)
    assert abs(proposals).max() <= 2*scale
    # Exponential has real even and imaginary odd coefficients.
    real = proposals*np.where(n%4 == 0,1,np.where(n%4 == 2,-1,0))
    imag = proposals*np.where(n%4 == 1,-1,np.where(n%4 == 3,1,0))
    worst = 0
    for k,(a,b) in enumerate(zip(real, imag, strict=True)):
        a, b = a.astype(object), b.astype(object)
        def derivative(v):
            out = np.zeros(DEGREE+3, dtype=object)
            for j in range(DEGREE-1,-1,-1):
                out[j] = out[j+2]+2*(j+1)*v[j+1]
            out[0] //= 2          # even: only 2*j*v_j enters the recurrence
            return out[:DEGREE+1]
        residual = (sum(abs(5*derivative(a)-2*k*b)) + sum(abs(5*derivative(b)+2*k*a)))
        value_at_zero = sum(a[::4])-sum(a[2::4])
        worst = max(worst, residual+5*abs(value_at_zero-scale))
    error = Fraction(int(worst),5*scale)
    assert error < Fraction(1,100000)
    return (real >> (COEFF_BITS-PRODUCT_BITS), imag >> (COEFF_BITS-PRODUCT_BITS)), error


def arrival_grid(h, rows, weights, weight_error):
    """Exact integer products and a uniform 8×8 operator error on free arrivals."""
    dim = h.shape[0]
    row = int(abs(h).sum(axis=1).max())
    assert dim <= 88**2 and row <= 1024 and (h != h.T).nnz == 0
    indices = rows.ravel()
    scale = 1 << MOMENT_BITS
    previous = np.zeros((dim,len(indices)), dtype=np.int64)
    previous[indices] = np.eye(len(indices),dtype=np.int64)*scale
    current = (h @ previous)//1024
    shift = MOMENT_BITS-PRODUCT_BITS
    projected = np.empty((DEGREE+1,len(indices),len(indices)), dtype=np.int64)
    projected[0], projected[1] = previous[indices] >> shift, current[indices] >> shift
    for k in range(2,DEGREE+1):
        if int(abs(current).max()) > 2*scale or row*2*scale > (1<<60):
            raise OverflowError("Chebyshev recurrence exceeds its checked integer bound")
        previous, current = current, (h @ current)//512-previous
        projected[k] = current[indices] >> shift
    unit = 1 << PRODUCT_BITS
    # Moment truncation: each column recurrence floor costs <=88/scale;
    # ||U_j(H)||<=j+1 for the Chebyshev polynomials of the second kind.
    assert (abs(projected).max() <= 2*unit and
            (DEGREE+1)*4*unit**2 < (1<<62))
    result = []
    for w in weights:
        packed = (w @ projected.reshape(DEGREE+1,-1)) >> (2*PRODUCT_BITS-OUT_BITS)
        result.append(packed.reshape(-1,len(rows),8,len(rows),8).transpose(0,1,3,2,4))
    moment_error = Fraction(3*88*DEGREE*(DEGREE+1),2*scale) + Fraction(8,unit)
    # Each moment block has norm <=1 before truncation. The coefficients
    # have just one active real/imaginary part each and magnitude <=2.
    error = (weight_error+(DEGREE+1)*(1/Fraction(unit)+3*moment_error)
             + Fraction(16,1 << OUT_BITS))
    return result,error


def gram_upper(real, imag, i, j, weights, divisor):
    """Exact integer quadrature and an eighth-power Schur norm upper bound."""
    a, b = real[:,j,i], imag[:,j,i]
    mx = int(max(abs(a).max(),abs(b).max()))
    assert 16*mx**2*int(np.sum(abs(weights))) < 1<<62
    at, bt = a.swapaxes(-1,-2), b.swapaxes(-1,-2)
    re = np.einsum("t,tij->ij",weights,at@a+bt@b).astype(object)
    im = np.einsum("t,tij->ij",weights,at@b-bt@a).astype(object)
    assert np.array_equal(re,re.T) and np.array_equal(im,-im.T)
    for _ in range(3):
        re,im = re@re-im@im,re@im+im@re
    norm = Fraction(int(max((abs(re)+abs(im)).sum(axis=1))),
                    (divisor*(1 << (2*OUT_BITS)))**8)
    upper = Fraction(ceil(10000*float(norm)**(1/8)),10000)
    if upper**8 < norm:
        upper += Fraction(1,10000)
    assert upper**8 >= norm
    return upper


def local_six():
    """The original six local station ranges, through weak-monitor passivity.

    Write B(t)=P_target U(t)P_start, S(t)=P_target U(t)P_target, and W
    for the no-click arrival; B-W=(gamma/2) int_0^t S(t-s)W(s) ds.
    For every initial vector v, the real L² inner product of Wv with Bv-Wv
    is (gamma/4)||int U(-s)P_target W(s)v ds||² >=0. Thus the weak first-hit
    probability is <=gamma ||int B*B dt||, including every possible head.
    For discrete pre-click W_k the Volterra sum omits the current diagonal.
    The same Gram-square identity gives (1-eta/2)||Wv|| <= ||Bv|| with
    eta=1-exp(-h/2). The hit prefactor p/(1-eta/2)²=4 tanh(h/4)<=h.
    The two one-leg upper bounds multiply even with the conditioned relay
    state, since each is uniform over all states and shorter time windows.
    """
    ports = (0,39,5,7,40,42)
    board = live_radar.source.formation.federation(1)
    route = live_radar.source.formation.region(board)
    rows = live_radar.station_rows(board,ports)
    weights, error = checked_coefficients()
    reports = {}
    for charge in (-1,1):
        h = live_radar.source.system(board,route,charge,coupled=True)[0]
        (re,im), radius = arrival_grid(h,rows,weights,error)
        # Composite Simpson for each unit vector. The fourth derivative of
        # U*PU is bounded by (2||H||)^4 via four commutators.
        hnorm = Fraction(int(abs(h).sum(axis=1).max()),40)
        quad = Fraction(T,180*SAMPLE_RATE**4)*(2*hnorm)**4
        matrix_error = T*(2*radius+radius**2)
        per_kind = {}
        edges = list(combinations(range(6),2))
        for name in ("continuous","spacing_1/16"):
            legs = {}
            for k in range(1,9):
                # Eight windows in units T/8. Scalar operator bounds increase
                # with the window because these are cumulative probabilities.
                count = k*T*SAMPLE_RATE//8
                w = np.array([1]+[4 if t%2 else 2 for t in range(1,count)]+[1])
                for i,j in edges:
                    for start,end in ((i,j),(j,i)):
                        if name == "continuous":
                            norm = gram_upper(re[:count+1],im[:count+1],start,end,w,3*SAMPLE_RATE)
                            norm += quad*k/8
                        else:
                            norm = gram_upper(re[4:count+1:4],im[4:count+1:4],start,end,
                                              np.ones(count//4,dtype=int),16)
                        legs[start,end,k] = min(1,norm+matrix_error*k/8)
            def by_deadline(i,j,k):
                # A return by t has at least one leg <=t/2. Split the
                # domain; bounds apply also to the conditional relay state.
                outward, back = legs[i,j,k], legs[j,i,k]
                return min(outward*back, outward*legs[j,i,k//2]+legs[i,j,k//2]*back)
            worst = max(by_deadline(i,j,8) for i,j in edges)
            # The time rebate averages one half the cumulative return law.
            # The upper Riemann sum works for the continuous and finite
            # detector, including the atom at the deadline.
            rebate = max(Fraction(T,8)*sum(by_deadline(i,j,k) for k in (2,4,6,8)) for i,j in edges)
            record = panel(worst,Fraction(T,2)-rebate)
            record["ordered_one_leg_upper"] = {f"{ports[i]}_to_{ports[j]}":str(v)
                                                for (i,j,k),v in legs.items() if k == 8}
            per_kind[name] = record
        reports[str(charge)] = dict(kernel_op_error=str(radius),simpson_error=str(quad),**per_kind)
    return dict(stations=list(ports),coefficient_polynomial_uniform_error=str(error),by_charge=reports)


def run():
    board = live_radar.source.formation.federation(1)
    stations, hops = station_challenge(board)
    floor = min(hops[i,j] for i,j in combinations(range(6),2)).item()
    norm = {str(q):intercell_norm(board, q) for q in (-1,1)}
    v = max(norm.values())
    return dict(cells=board.cells, stations=stations, cell_distances=hops.tolist(),
                minimum_cell_distance=floor, intercell_operator_norm={k:str(x) for k,x in norm.items()},
                deadline=T, detector_rate=1, continuous=panel(*continuum(floor,v)),
                discrete={str(Fraction(1,h)):panel(*discrete(floor,v,h)) for h in (2,4,8,16,32,64)},
                original_local_six=local_six())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    output = json.dumps(run(), indent=2) + "\n"
    if parser.parse_args().write:
        Path(__file__).with_suffix(".json").write_text(output)
    else:
        print(output, end="")
