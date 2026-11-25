"""Bounded integer replay of the half-tick, cap-six witness on the full 80-cell board."""

from fractions import Fraction
from math import factorial, isqrt
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import qr

from model.experiments.interval_echo_einstein_station.collar_pairing.ranging.certificate import integer_flight
from .measure import basis, components, pairing, permutations

DEGREE = 640
LIMB = 50
BASE = 2**LIMB
UNIT = 2**102
FRAC = 2**100
MASK = BASE-1


def weights():
    # Taylor polynomial of exp(-izx) in the Chebyshev basis, z = 32t.
    # After extracting (-i)^k and dividing by DEGREE!, all summands are integers.
    fact = factorial(DEGREE)
    assert Fraction(192**(DEGREE+1), factorial(DEGREE+1)) < Fraction(1,10**59)
    rows = []
    for tick in range(1,13):
        half_z = 8*tick
        co = []
        for k in range(DEGREE+1):
            term = fact//factorial(k)*half_z**k
            total = term
            for m in range(1, (DEGREE-k)//2+1):
                term = term*half_z**2//(m*(m+k))
                total += (-1 if m%2 else 1)*term
            c = (1 if k==0 else 2)*total*FRAC//fact
            assert abs(c) <= 3*FRAC
            co.append(c)
        rows.append(co)
    return np.array(rows,dtype=object)


def divided_product(h, high, low, shift):
    a, b = h @ high, h @ low
    a += b >> LIMB
    return a >> shift, ((a % (1<<shift)) << (LIMB-shift)) + ((b & MASK) >> shift)


def certified_kernels(board, route, charge, co):
    h = integer_flight(board, route, charge)
    d = h.shape[0]
    row = int(abs(h).sum(axis=1).max())
    assert row <= 41
    # Rounding a component down by one UNIT quantum costs <=sqrt(d)/UNIT in column norm.
    col_error = Fraction(3*(1+isqrt(d))*DEGREE*(DEGREE+1)*(DEGREE+2),6*UNIT)
    col_error += Fraction(3*(DEGREE+1),FRAC) + Fraction(1,10**59)
    assert 4*col_error < Fraction(1,10**19)
    stations = np.array([4*p+i for p in route for i in range(4)])
    left, right = pairing.formation.starting_coins(route,1)
    initial = 2*left+right
    state = np.zeros((2,d),dtype=object)
    output = np.zeros((2,12,16,16),dtype=object)
    high = np.zeros((d,16),np.int64)
    low = high.copy()
    high[stations] = (UNIT//BASE)*np.eye(16,dtype=np.int64)
    prev_hi = prev_lo = None
    for k in range(DEGREE+1):
        # The l2 norm of an exact Chebyshev column is <=1. Check integer overflow separately.
        assert row*(int(np.max(np.abs(high)))+1) < 2**62 and row*BASE < 2**62
        partial = high[stations].astype(object)*BASE+low[stations]
        sign = 1 if k%4 in (0,3) else -1
        output[k%2] += sign*co[:,k,None,None]*partial[None,:,:]
        if charge == 1:
            column = high[:,initial].astype(object)*BASE+low[:,initial]
            state[k%2] += sign*co[7,k]*column
        if k==DEGREE:
            break
        nh, nl = divided_product(h,high,low,6 if k==0 else 5)
        if k:
            nl -= prev_lo
            nh += (nl >> LIMB)-prev_hi
            nl &= MASK
        prev_hi, prev_lo, high, low = high, low, nh, nl
    # Scale every scalar kernel entry to the common denominator FRAC.
    output //= UNIT
    output = output.reshape(2,12,4,4,4,4).transpose(0,1,2,4,3,5)
    return output, state//UNIT, float(col_error)


def multiply(a,b):
    return np.stack(((a[0]@b[0]-a[1]@b[1])//FRAC, (a[0]@b[1]+a[1]@b[0])//FRAC))


def adjoint(a):
    return np.stack((a[0].T, -a[1].T))


def direction_effects(u, cap):
    def arrivals(i,j):
        f = [np.zeros((2,4,4),dtype=object)]
        for n in range(1,cap):
            a = u[:,n-1,j,i].copy()
            for k in range(1,n):
                a -= multiply(u[:,n-k-1,j,j],f[k])
            f.append(a)
        return f

    results = []
    for i,j in permutations(range(u.shape[2]),2):
        out, back = arrivals(i,j), arrivals(j,i)
        total = np.zeros((2,4,4),dtype=object)
        for n in range(2,cap):
            for k in range(1,n):
                a = multiply(back[n-k],out[k])
                total += (n-cap)*multiply(adjoint(a),a)
        results.append(total)
    # Each numerator still needs division by 4 (cap time weight), before charge pooling.
    return results


def error_effect():
    # Block kernel error includes fixed rounding of both complex parts. Each rounded
    # 4-by-4 product costs at most 8/FRAC in operator norm by the Frobenius norm.
    eps = Fraction(1,10**19) + Fraction(8,FRAC)
    rnd = Fraction(8,FRAC)
    f = [Fraction(0)]
    for n in range(1,12):
        # Every true U block and every true first-hit f has operator norm <=1.
        f.append(eps + sum(eps+(1+eps)*f[k]+rnd for k in range(1,n)))
    total = Fraction(0)
    for n in range(2,12):
        for k in range(1,n):
            a = f[k]+f[n-k]+f[k]*f[n-k]+rnd
            total += Fraction(12-n,4)*(2*a+a*a+rnd)
    return total


def sym_matrix(z):
    return sp.Matrix(z.shape[0],z.shape[1],lambda i,j: sp.Rational(float(z[i,j].real))
                     + sp.I*sp.Rational(float(z[i,j].imag)))


def one_norm(a):
    return max(sum(abs(complex_val.as_real_imag()[0])+abs(complex_val.as_real_imag()[1])
                   for complex_val in a[:,j]) for j in range(a.cols))


def log_bound(w, approx):
    # Bound log(w) by two rational atanh series on [1,2], including the log(2) term.
    e, s = 0, w
    while s < 1:
        s *= 2
        e -= 1
    while s > 2:
        s /= 2
        e += 1
    def series(x):
        t = (x-1)/(x+1)
        return 2*sum(t**(2*k+1)/(2*k+1) for k in range(20))
    err = 3*(1+abs(e))*sp.Rational(1,3)**41
    assert abs(series(s)+e*series(sp.Integer(2))-approx)+err < sp.Rational(1,10**12)


def entropy_coordinates(reference, input_error, lower, operators):
    w,V = np.linalg.eigh(np.array(reference,dtype=complex))
    axes = sym_matrix(V)
    diag = sp.diag(*[sp.Rational(float(x)) for x in w])
    logs = sp.diag(*[sp.Rational(float(np.log(x))) for x in w])
    for k in range(4):
        log_bound(diag[k,k],logs[k,k])
    defect = one_norm(axes.H*axes-sp.eye(4))
    residue = one_norm(reference-axes*diag*axes.H)
    assert defect < sp.Rational(1,10**13) and residue < sp.Rational(1,10**13)
    assert max(diag.diagonal())<1 and min(logs.diagonal())>-21
    reference_error = residue+3*defect+input_error
    assert min(diag.diagonal())-reference_error > lower
    elog = reference_error/lower+63*defect+sp.Rational(2,10**12)
    K = -axes*logs*axes.H
    return sp.Matrix([sp.re(sp.trace(K*p)) for p in operators]), elog


def certificate():
    board = pairing.formation.federation(1)
    route = pairing.formation.region(board)
    co = weights()
    data = []
    for charge in (-1,1):
        kernels, vector, bound = certified_kernels(board,route,charge,co)
        data.append(np.array([direction_effects(kernels,cap) for cap in (8,12)],dtype=object))
    # Last charge is the forming excess. Its phase convention cancels on reducing the state.
    a, b = (v.reshape(-1,4) for v in vector)
    rho = (sp.Matrix(a.T@a+b.T@b)+sp.I*sp.Matrix(b.T@a-a.T@b))/FRAC**2
    assert bound < 1e-20
    combined = (data[0]+data[1])
    ops = [sym_matrix(op) for op in basis(False)]
    def coordinates(panel):
        exact_effects = [(sp.Matrix(z[0])+sp.I*sp.Matrix(z[1]))/(8*FRAC) for z in panel]
        return sp.Matrix([[sp.re(sp.trace(m*p)) for p in ops] for m in exact_effects])
    A = coordinates(combined[1])
    eps = error_effect()
    assert eps < Fraction(1,10**13)
    # Each Pauli/2 has trace norm 2. The true 12x15 matrix 1-norm is within 24 eps.
    ea = 24*sp.Rational(eps.numerator,eps.denominator)
    num = np.array(A,dtype=float)
    ent, elog = entropy_coordinates(rho,sp.Rational(1,10**18),sp.Rational(6,10**9),ops)
    B = sym_matrix(np.linalg.pinv(num))
    _,_,vt = np.linalg.svd(num,full_matrices=True)
    e = np.array(ent,dtype=float).ravel()
    tmp = vt[12:].T@(vt[12:]@e)
    d0 = sym_matrix((tmp/np.linalg.norm(tmp))[:,None])
    bnd = one_norm(B)
    mb = one_norm(sp.eye(12)-A*B)+ea*bnd
    assert mb<sp.Rational(1,10)
    correction = bnd*(one_norm(A*d0)+ea*one_norm(d0))/(1-mb)
    # Each true entropy coordinate differs by <=2 elog (trace norm of Pauli/2).
    slope = (ent.T*d0)[0] - 2*elog*one_norm(d0) - (max(abs(v) for v in ent)+2*elog)*correction
    assert slope > sp.Rational(1,10)
    # Certify that a known marginal, or an extra cap on the same tick, removes the ambiguity.
    names = [(i,j) for i in range(4) for j in range(4) if i or j]
    def full_column_rank(mat):
        left = sym_matrix(np.linalg.pinv(np.array(mat,dtype=float)))
        uncertainty = one_norm(sp.eye(mat.cols)-left*mat) + one_norm(left)*2*mat.rows*sp.Rational(eps)
        assert uncertainty < 1
        return float(uncertainty)
    left_known = full_column_rank(A[:,[k for k,(i,j) in enumerate(names) if j]])
    right_known = full_column_rank(A[:,[k for k,(i,j) in enumerate(names) if i]])
    two_menus = full_column_rank(coordinates(combined[0]).col_join(A))

    # A native sharp address read of the forming excess supplies a spanning ensemble.
    # Pick candidate columns using floating QR, then certify their exact independence.
    probs = (a*a+b*b).sum(axis=1)
    elig = np.flatnonzero(np.array(probs,dtype=float)/float(FRAC**2) > 1e-7)
    ops16 = [sp.eye(4)/2]+ops
    complex_branches = (np.array(a[elig],dtype=float)+1j*np.array(b[elig],dtype=float))/float(FRAC)
    approx = np.einsum("xi,xj->xij",complex_branches,complex_branches.conj())
    approx /= np.trace(approx,axis1=1,axis2=2).real[:,None,None]
    frame = components(approx,np.array([np.array(v,dtype=complex) for v in ops16]))
    selected = elig[qr(frame.T, pivoting=True,mode="economic")[2][:16]]
    prep = []
    for x in selected:
        v = sp.Matrix(a[x]+sp.I*b[x])
        conditional = (v*v.H)/int(probs[x])
        prep.append([sp.re(sp.trace(conditional*op)) for op in ops16])
    bank = sp.Matrix(prep).T
    inverse = sym_matrix(np.linalg.inv(np.array(bank,dtype=float)))
    # The exact normalized vector is within 2*1e-20/sqrt(1e-7) of the simulated one.
    # Its rank-one projector has trace-norm error <= twice that distance.
    prep_uncertainty = one_norm(sp.eye(16)-bank*inverse)+16*sp.Rational(2,10**16)*one_norm(inverse)
    assert prep_uncertainty < sp.Rational(1,100)
    pmin = min(sp.Rational(int(probs[x]),FRAC**2) for x in selected)-sp.Rational(1,10**19)
    assert pmin > sp.Rational(1,10**7)
    cinv = one_norm(inverse)/(1-prep_uncertainty)
    amplitude = pmin/(3*cinv*(one_norm(d0)+correction))
    # Test the very same true D against a separately selected, mixed parity source.
    c = sp.Rational("-0.004650970878318474")
    band = sp.Rational(1,10**6)
    even_weight = sp.Rational(3,4)
    s_even = sp.Rational(1,2)+c/even_weight
    marginal_a = sp.Matrix(2,2,lambda i,j:sum(rho[2*i+k,2*j+k] for k in range(2)))
    marginal_b = sp.Matrix(2,2,lambda i,j:sum(rho[2*k+i,2*k+j] for k in range(2)))
    independent = sp.kronecker_product(marginal_a,marginal_b)
    tau = (1-s_even)*independent + s_even*rho
    # Real input may also move within |c-c_center|<=1e-6, with ||rho-pi||<=2.
    mixed_row, e_mixed = entropy_coordinates(tau,sp.Rational(1,10**16)+2*band/even_weight,
                                              sp.Rational(8,1000),ops)
    mixed_slope = ((mixed_row.T*d0)[0] - 2*e_mixed*one_norm(d0)
                   - (max(abs(v) for v in mixed_row)+2*e_mixed)*correction)
    formed_fraction = even_weight*s_even-band
    never_dephased = sp.Rational(1,4)+2*c-2*band
    assert mixed_slope > sp.Rational(1,100)
    receipt = json.loads(Path(__file__).resolve().parents[1].joinpath("full_station.json").read_text(),
                         parse_float=sp.Rational)["independent_ordinary_calibration"]
    u0,s0,v0 = receipt["selected_coordinates"]
    strength = (s0-v0)/(1-u0)
    halo = sp.Rational(1,10**7)
    lo = (s0-v0-2*halo)/(1-u0+halo)
    hi = (s0-v0+2*halo)/(1-u0-halo)
    calibrated_tau = (1-strength)*independent+strength*rho
    hc, he = entropy_coordinates(calibrated_tau,sp.Rational(1,10**16)+2*max(strength-lo,hi-strength),
                                  sp.Rational(8,1000),ops)
    base_signal = (hc.T*d0)[0]
    hc_floor = base_signal-2*he*one_norm(d0)-(max(abs(v) for v in hc)+2*he)*correction
    accessible = s0-2*v0-3*halo
    assert hc_floor > sp.Rational(1,100)
    return dict(chebyshev_degree=DEGREE,fractional_bits=102,column_error_bound=bound,
                range_operator_error_bound=float(eps),entropy_operator_error_bound=float(elog),
                right_inverse_error=float(mb),exact_kernel_tangent_l1_correction=float(correction),
                certified_entropy_slope_floor=float(slope),exact_rank=12,
                known_left_inverse_error=left_known,known_right_inverse_error=right_known,
                two_menu_left_inverse_error=two_menus,span_address_ids=selected.tolist(),
                preparation_inverse_error=float(prep_uncertainty),reachable_step_lower=float(amplitude),
                reachable_entropy_slope_floor=float(amplitude*slope),held_c_band=float(band),
                held_even_entropy_slope_floor=float(mixed_slope),
                held_total_entropy_slope_floor=float(formed_fraction*mixed_slope),
                held_reachable_total_slope_floor=float(amplitude*formed_fraction*mixed_slope),
                held_never_dephased_slope_floor=float(never_dephased*mixed_slope),
                held_reachable_never_dephased_floor=float(amplitude*never_dephased*mixed_slope),
                calibrated_coordinates=[float(v) for v in (u0,s0,v0)], calibrated_box=float(halo),
                calibrated_s_even=float(strength), calibrated_even_center=float(base_signal),
                calibrated_even_floor=float(hc_floor), calibrated_untouched_weight=float(s0-2*v0),
                calibrated_total_center=float((s0-2*v0)*base_signal),
                calibrated_total_floor=float(accessible*hc_floor),
                calibrated_reachable_floor=float(amplitude*accessible*hc_floor))


if __name__ == "__main__":
    print(json.dumps(certificate(),indent=2))
