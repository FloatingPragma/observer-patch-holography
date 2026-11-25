"""Certify the Lie image of the twelve public current strengths on both native sectors."""
import argparse
import json
from pathlib import Path
from fractions import Fraction

import numpy as np
import sympy as sp
from scipy.linalg import expm

from model.candidate_laws.interval_echo import currents
from model.executable_model import cell
from ..__main__ import ordered_word

SCALE = 10**11


def integers(a):
    n = np.rint(a).astype("int64")
    if not np.array_equal(n,a):
        raise ValueError("nonintegral native operator")
    return n.astype(object)


def squares(a):
    return int(np.sum(a*a))


def positive(a):
    m = sp.Matrix(a)
    return all(m[:k,:k].det() > 0 for k in range(1, m.rows+1))


def count(factors, lo, hi):
    return sum(m*p.count_roots(lo,hi) for p,m in factors)


def band_certificate(h, factors, energy, v, source, detector):
    """Isolate the triple, bound its isometry error and certify the source's admission sign."""
    b = np.rint(SCALE*v).astype("int64").astype(object)
    t = round(SCALE*2*energy)
    center = sp.Rational(t, SCALE)
    gap = sp.Rational(2,25)
    assert count(factors,center-gap,center+gap) == 3
    gram = b.T@b - SCALE**2*np.eye(3,dtype=object)
    residual = SCALE*(h@b)-t*b
    assert max(squares(gram),squares(residual))*10**16 < SCALE**4
    projected = b@(b.T@source)
    signed = ((2*detector.astype(int)-1)[:,None]*projected**2).sum(axis=0)
    source_size = (source**2).sum(axis=0)
    scores = [sp.Rational(int(s),int(n)*SCALE**4) for s,n in zip(signed,source_size)]
    assert all(abs(s)>sp.Rational(1,10**5) for s in scores)
    assert all((s>0)==(scores[0]>0) for s in scores)
    return b, bool(scores[0]>0)


def inspect(board, pair):
    sector, source, h, poles, active, _, class_sum, _, exposed, _ = currents._service(board,pair,1.)
    h2 = integers(2*h)
    klass = integers(class_sum)
    assert np.array_equal(klass@h2,h2@klass)
    factors = sp.Matrix(h2).charpoly().as_poly().factor_list()[1]
    sharp = integers(np.rint(source*np.sqrt(10 if pair else 2)))
    assert np.allclose(source,sharp.astype(float)/np.sqrt(10 if pair else 2),atol=1e-14,rtol=0)
    # All sourced eigenvalues belong to one irreducible native characteristic factor.
    first = sp.Rational(str(2*poles[0]["energy"]))
    factor, = [p for p,m in factors if p.count_roots(first-sp.Rational(1,20),first+sp.Rational(1,20))]
    stay = np.zeros_like(sharp)
    for c in factor.all_coeffs():
        stay = h2@stay + int(c)*sharp
    assert not np.any(stay) and factor.degree()==len(poles)
    detector = sector.adjacent if pair else np.ones(board.n, bool)
    eig, vec = np.linalg.eigh(h)
    retained = iter(active)
    bases = []
    for pole in poles:
        energy = pole["energy"]
        band = next(retained) if pole["retained"] else vec[:,np.abs(eig-energy)<1e-8]
        b, admitted = band_certificate(h2,factors,energy,band,sharp,detector)
        assert admitted==pole["retained"]
        if admitted:
            sign = 1 if not bases else -1
            assert positive(sign*b.T@klass@b-SCALE**2*np.eye(3,dtype=object))
            bases.append(b)
    assert len(bases)==2
    antipode, turns = currents._geometry(board)
    if not pair:
        adj = integers(board.adjacency().toarray())
        odd = np.eye(board.n,dtype=int)-np.eye(board.n,dtype=int)[antipode]
        assert np.array_equal(adj@adj@odd,5*odd)
        assert not np.any(np.asarray(adj@sharp,dtype=object)*sharp)
    j = [integers(a.imag) for a in turns]
    if pair:
        j = [integers(currents.port.second_quantized(sector,a.astype(float)).real) for a in j]
    assert not np.any(sum(j))
    flat, norms = [], []
    for p, rotation in enumerate(j):
        load = ((sector.aa==p).astype(int)+(sector.bb==p)) if pair else (np.arange(board.n)==p)
        assert max((np.abs(rotation).sum(axis=1)+load).tolist()) <= 5
        first, second = bases
        real = first.T@(load[:,None]*first)
        imag = first.T@rotation@first
        weak = second.T@rotation@second
        flat.append(np.r_[real.ravel(),imag.ravel(),weak.ravel()])
        if p==0:
            norms = [(sp.Rational(squares(real)+squares(imag),SCALE**4)
                     - sp.Rational(int(real.trace())**2,3*SCALE**4)),
                     sp.Rational(squares(weak),SCALE**4)]
    flat = np.array(flat,dtype=object)
    gram = flat@flat.T
    assert positive(20*gram-SCALE**4*np.eye(12,dtype=object))
    # Each normalized trace-square differs by < .001 from this rational center.
    enclosures = []
    for value in norms:
        k = int(sp.floor(1000*value))
        enclosures.append([(k-1)/1000,(k+2)/1000])
    return exposed, active, h, source[:,currents.SOURCE_PORT], antipode, {
        "loads": 1+int(pair), "sourced_characteristic_factor_of_2H": [int(c) for c in factor.all_coeffs()],
        "positive_Gram_lower": ".05 for approximate columns",
        "true_current_span": 12, "port_zero_traceless_squared_HS": enclosures}


def outcomes(addresses, basis, h, source, antipode, loads):
    u = np.column_stack(basis)
    ideal = np.eye(len(source)) + ((-1)**loads-1)*(basis[0]@basis[0].T)
    first = currents.SOURCE_PORT
    data = {}
    sector = currents.PairSector(cell()) if loads == 2 else None
    def read(psi):
        weights = abs(psi)**2
        if sector is None:
            return weights
        return np.r_[weights[sector.seam_states],weights[~sector.adjacent].sum()]
    # By proper-rotation covariance, each port has the same bounds as port zero.
    # Centre traces N/4 and the certified traceless squared norms give these bounds.
    exposure = Fraction(22,7)*12*(Fraction(7,5) if loads==1 else Fraction(81,100))
    for steps in (512,2048):
        w = ordered_word(addresses,np.ones(12),np.pi,steps)
        whole = np.eye(len(source))+u@(w-np.eye(6))@u.T
        entry = {"operator_error": float(np.linalg.norm(whole-ideal,2))}
        step_bound = exposure**3/(3*steps**2)/(1-exposure/steps)
        entry["rigorous_operator_error_upper"] = float(np.nextafter(float(step_bound),np.inf))
        for time in (0,16):
            free = expm(-1j*time*h)
            sham = read(free@source)
            receipts = read(free@whole@source)
            limit = read(free@ideal@source)
            entry[str(time)] = {
                "sham": sham.tolist(),"word": receipts.tolist(),
                "total_variation": float(abs(receipts-sham).sum()/2),
                "ideal_total_variation": float(abs(limit-sham).sum()/2)}
            if loads==1:
                entry[str(time)]["input_antipodal_pair_probability"] = float(receipts[[first,antipode[first]]].sum())
        if steps==2048:
            control = np.array(addresses)-np.array(addresses).real
            u0 = ordered_word(control,np.ones(12),np.pi,steps)
            off = np.eye(len(source))+u@(u0-np.eye(6))@u.T
            pulsed = read(whole@source)
            rng = np.random.default_rng(104731+loads)
            comparison = {}
            for name,arm in (("pulse",pulsed),("sham",read(source)),("occupation_gain_zero",read(off@source))):
                comparison[name] = rng.multinomial(8192,arm/arm.sum()).tolist()
            entry["zero_flight_shots_per_setting"]=8192
            entry["counts"]=comparison
            entry["zero_gain_operator_deviation"]=float(np.linalg.norm(off-np.eye(len(source)),2))
        data[str(steps)] = entry
    return data


def run():
    board = cell()
    panels=[]
    for paired in (False,True):
        ops, bases, h, source, antipode, item = inspect(board,paired)
        item["common_turn"] = outcomes(ops,bases,h,source,antipode,item["loads"])
        panels.append(item)
    for single, pair in zip(*(p["port_zero_traceless_squared_HS"] for p in panels)):
        assert pair[1] < single[0]
    return {"sectors":panels, "joint_lie_dimension":23, "central_weights":[1,2],
            "simply_connected_semisimple_times_circle_kernel_order":12}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    ans=json.dumps(run(),indent=2,allow_nan=False)+"\n"
    if parser.parse_args().write:
        Path(__file__).with_name("results.json").write_text(ans)
    else:
        print(ans,end="")
