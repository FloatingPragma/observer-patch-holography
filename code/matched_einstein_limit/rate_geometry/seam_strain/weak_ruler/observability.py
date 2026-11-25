"""Information in every reciprocal leading contact of graph orders one to four.

The incoming two-head variation is sum_j d_j P_j/4 in lexicographic order
XX,XY,...,ZZ. Both one-head marginals and the 2:1 rail lottery are fixed.
The output observable is the negative log of the retained reference density.

An exact six-row linear system corrects a proposed normalized source direction.
An independent nine-row system expresses it as a signed sum of native rotated
and re-paired source preparations. Positive and negative parts give two physical
lotteries; the remaining probability in each prepares the original reference.
With a common admixture t their reciprocal leading coefficients are identical
on every registered contact at orders one through four. They have different
retained entropy slopes. Only the lottery weights are calibrated; the optical
operators, source channel and entropy are measured independently. Higher-order
time coefficients or full timestamp histories may supply additional readings.
"""

import argparse
from functools import reduce
from itertools import product
import json
from math import gcd
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import qr
from scipy.sparse.csgraph import shortest_path

from . import threshold
from .. import live_modular, live_source


AXES = live_source.qubit_hidden.correlation_axes()
SELECTED = ((3,0,39),(2,1,10),(2,39,46),(4,1,37),(4,2,38),(1,1,38))
NATIVE_ARMS = ("independent_repairing","L2R2","L1R3","L0R2","L3R3","L2R0","L1R1","L1R2","L0R1")
HERE = Path(__file__).parent


def norm_inf(a):
    return max(sum(abs(v) for v in a.row(i)) for i in range(a.rows))


def enumerate_fronts():
    """Remove identity offsets and scalar multiples by exact integer equality."""
    board = threshold.source.formation.federation(1)
    route = threshold.source.formation.region(board)
    hops = shortest_path(board.adjacency(),unweighted=True,directed=False)
    minus, plus = (threshold.hop(board,route,q) for q in (-1,1))
    operators, chosen, counts = {}, {}, []
    for order in threshold.ORDERS:
        a,b = np.where(np.triu(hops == order,1))
        m = threshold.blocks(minus,a,b,order)
        np.testing.assert_array_equal(m.transpose(0,2,1),threshold.blocks(plus,a,b,order))
        distinct, index, inverse = np.unique(m.reshape(-1,64),axis=0,return_index=True,return_inverse=True)
        keys = []
        for x, first in zip(distinct,index,strict=True):
            op, _ = threshold.effect(x.reshape(8,8),order)
            # Subtract the trace using integers even for large order-four powers.
            centered = 8*op-np.trace(op)*np.eye(8,dtype=np.int64).astype(object)
            flat = tuple(map(int,centered.ravel()))
            divisor = reduce(gcd,flat,0)
            if not divisor:
                keys.append(None)
                continue
            orientation = 1 if next(v for v in flat if v)>0 else -1
            key = tuple(orientation*v//divisor for v in flat)
            operators.setdefault(key, [order,int(a[first]),int(b[first])])
            keys.append(key)
        for label in SELECTED:
            n,i,j = label
            if n == order:
                pair, = np.flatnonzero((a==i)&(b==j))
                chosen[label] = keys[inverse[pair]]
        span = sp.Matrix(list(operators)).rank()
        counts.append(dict(order=order,links=len(a),cumulative_classes=len(operators),cumulative_rank=span))
    keys = [chosen[label] for label in SELECTED]
    assert sp.Matrix(keys).rank() == counts[-1]["cumulative_rank"] == 6
    return [sp.Matrix(8,8,key) for key in operators], [sp.Matrix(8,8,key) for key in keys], counts


def input_column(a):
    return sp.Matrix([sp.re(4*sp.trace(x*a)) for x in AXES])


def marginal_generator(formed, ferr, side):
    m = sp.Matrix(2,2,lambda i,j:sum(formed[2*i+k,2*j+k] if side==0
                                     else formed[2*k+i,2*k+j] for k in range(2))) - sp.eye(2)/2
    gap2 = sp.re(2*sp.trace(m*m))
    assert sp.Rational(1,4)+4*ferr < gap2 < 1
    nu = sp.Rational(float(1/np.sqrt(float(gap2))))
    # Weyl and trace contraction under a partial trace, as in live_modular.
    error = 10*ferr + 2*abs(nu**2*gap2-1)
    g = sp.kronecker_product(nu*m,sp.eye(2)) if side==0 else sp.kronecker_product(sp.eye(2),nu*m)
    return g,error


def quarter_turn(rho, g, k):
    if k==0:
        return rho
    if k==2:
        return (4*g*rho*g).applyfunc(sp.expand)
    return ((rho+4*g*rho*g + (1 if k==1 else -1)*2*sp.I*(rho*g-g*rho))/2).applyfunc(sp.expand)


def native_bank(formed, product_head, ferr):
    """Use rotations about the *actual* marginal axes and independent re-pairing."""
    (gl,el),(gr,er) = (marginal_generator(formed,ferr,i) for i in (0,1))
    names, columns = [], []
    for a,b in product(range(4),repeat=2):
        if a==b==0:
            continue
        r = quarter_turn(quarter_turn(formed,gl,a),gr,b)
        names.append(f"L{a}R{b}")
        columns.append(input_column(r-formed))
    names.append("independent_repairing")
    columns.append(input_column(product_head-formed))
    stuff = sp.Matrix.hstack(*columns)
    # The true unitary differs from its approximation by <=2e per head;
    # a product then costs u, and conjugating a norm-one density costs
    # 2u+u². Include both copies of the reference and the re-pairing error.
    u = 2*(el+er)+4*el*er
    error = 3*ferr+2*u+u**2
    idx = [names.index(name) for name in NATIVE_ARMS]
    return [names[i] for i in idx], stuff[:,idx], error


def dyadic_vector(x):
    return sp.Matrix([sp.Rational(round(float(v)*2**48),2**48) for v in x])


def run():
    classes,ops,counts = enumerate_fronts()
    (formed, product_head), ferr, _ = live_source.formation_center()
    _,_,cerr,grams = live_source.pump_center(True)
    channel_error = 4*cerr+4*cerr**2
    outputs = [live_source.rational_read(grams,x,(2,1)) for x in AXES]
    ref = live_source.rational_read(grams,formed,(2,1))
    ref, rround = live_modular.rounded(ref)
    log, logerr, _, _, _ = live_modular.logarithm(ref,ferr+channel_error+rround)
    # On trace-zero inputs the identity offset is immaterial. Row scaling
    # to operator norm <=1 gives the same trace-norm error on every entry.
    ops = [o/norm_inf(o) for o in ops]
    center = sp.Matrix([[sp.re(sp.trace(o*y)) for y in outputs] for o in ops])
    entropy = sp.Matrix([sp.re(sp.trace(log*y)) for y in outputs])
    entropy_error = logerr + (7+logerr)*channel_error
    all_rows = np.array([[float(sp.re(sp.trace(o*y)))/np.linalg.norm(np.array(o,float))
                          for y in outputs] for o in classes])
    _, singular, vh = np.linalg.svd(all_rows,full_matrices=False)
    h = np.array(entropy,float).ravel()
    residual = h-(h@vh[:6].T)@vh[:6]
    v = dyadic_vector(residual/np.linalg.norm(residual))
    columns = qr(np.array(center,float),pivoting=True)[2][:6].tolist()
    invnorm = norm_inf(center[:,columns].inv())
    defect = 6*invnorm*channel_error
    assert defect < sp.Rational(1,1000)
    # Correct exactly six coordinates to kill every true read. The other
    # three coordinates stay fixed; the perturbation estimate is in l∞.
    error = invnorm*(max(abs(z) for z in center*v)+channel_error*sum(abs(z) for z in v))/(1-defect)
    norm_hi = 1+abs(v.dot(v)-1)+6*error
    h_lower = (entropy.dot(v) - entropy_error*sum(abs(z) for z in v)
               -error*sum(abs(entropy[i])+entropy_error for i in columns))/norm_hi
    assert h_lower > sp.Rational(53,100000)
    names, bank, bank_error = native_bank(formed, product_head, ferr)
    inv = bank.inv()
    bank_defect = 9*norm_inf(inv)*bank_error
    assert bank_defect < sp.Rational(1,1000)
    # Exact weights on nine independent native preparations implement the
    # unnormalized exact blind change. Split positive and negative weights
    # into two lotteries and fill each to one with unrotated reference.
    weights = inv*v
    werr = norm_inf(inv)*(error+bank_error*sum(abs(z) for z in weights))/(1-bank_defect)
    scale = int(sp.ceiling(max(sum(max(0,z) for z in weights),sum(max(0,-z) for z in weights))+9*werr))
    native_gain = h_lower*(1-abs(v.dot(v)-1)-6*error)/scale
    assert native_gain > 0 and min(abs(z) for z in weights)>werr
    # Interpolate the two output endpoints at common admixture t. Along the
    # segment their minimum eigenvalue is >=(1-t)m and its distance from the
    # reference in operator norm is <=t. The log resolvent bound then costs
    # <=t/[(1-t)m] per unit trace norm of the difference. The trace norm of
    # sum d_j P_j/4 is <=||d||_2; the retained channel contracts that norm.
    admix = sp.Rational(1,10**7)
    ref_floor = sp.Rational(14,10000)
    difference = admix*(native_gain - admix*norm_hi/((1-admix)*ref_floor*scale))
    assert difference > 0
    left = [float(max(0,z)/scale) for z in weights]
    right = [float(max(0,-z)/scale) for z in weights]
    # A separate, exactly encodable menu. Verify its actual near-equality
    # directly; the inverse bound on the ideal lottery is unnecessary here.
    size = 1 << 40
    plus = [int(sp.floor(max(0,z)*size/scale)) for z in weights]
    minus = [int(sp.floor(max(0,-z)*size/scale)) for z in weights]
    assert max(sum(plus),sum(minus)) <= size
    signed = (sp.Matrix(plus)-sp.Matrix(minus))/size
    direction = bank*signed
    coord_error = bank_error*sum(abs(z) for z in signed)
    out = sum((d*y for d,y in zip(direction,outputs,strict=True)),sp.zeros(8))
    # Every Pauli/4 input has trace norm one. Bound all ten distinct
    # centered output effects after scaling their absolute row-sum norm to one.
    front_error = (channel_error*sum(abs(z) for z in direction)
                   +9*coord_error)
    residuals = [abs(sp.re(sp.trace(o*out)))/norm_inf(o)+front_error for o in classes]
    dyadic_gain = (entropy.dot(direction) - entropy_error*sum(abs(z) for z in direction)
                   -coord_error*sum(abs(z)+entropy_error for z in entropy))
    assert max(residuals) < sp.Rational(1,10**12) and dyadic_gain > sp.Rational(2,10**8)
    return dict(
        fixed_rail_weights=[2,1],basis=[a+b for a,b in product("XYZ",repeat=2)],
        contacts=counts,independent_links=[list(p) for p in SELECTED],
        rank=6,rank_with_retained_modular=7,
        singular_values_unit_HS_per_distinct_centered_operator=singular.tolist(),
        modular_covector_center=h.tolist(),projected_residual_norm=float(np.linalg.norm(residual)),
        blind_direction_center=[float(z) for z in v],correction_columns=columns,
        maximum_coordinate_correction_upper=float(error)*1.001,
        front_minor_inverse_defect_upper=float(defect)*1.001,
        entropy_nats_per_unit_blind_direction_lower=0.00053,
        native_lotteries=dict(
            quarter_turns_about_formed_marginals=names, common_scale=scale,
            positive_center_probabilities=left,negative_center_probabilities=right,
            positive_reference_center=1-sum(left),negative_reference_center=1-sum(right),
            exact_probability_error_upper=float(werr/scale)*1.001,
            bank_inverse_defect_upper=float(bank_defect)*1.001,
            retained_entropy_derivative_difference_lower=float(native_gain)*.999,
            common_admixture=float(admix),finite_entropy_difference_lower=float(difference)*.999,
            finite_weight_control=dict(denominator=size, positive_numerators=plus,
                negative_numerators=minus,positive_reference_numerator=size-sum(plus),
                negative_reference_numerator=size-sum(minus),
                centered_unit_row_sum_front_difference_upper=1e-12,
                retained_entropy_derivative_difference_lower=2e-8),
        ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    output = json.dumps(run(),indent=2,allow_nan=False)+"\n"
    if parser.parse_args().write:
        Path(__file__).with_suffix(".json").write_text(output)
    else:
        print(output,end="")
