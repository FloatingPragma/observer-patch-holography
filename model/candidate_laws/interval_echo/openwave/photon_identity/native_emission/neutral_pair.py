"""Release one independently prepared excess–vacancy pair in the full blank-memory sector."""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.special import jv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
sys.path.insert(0,str(ROOT))
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from model.candidate_laws.interval_echo.openwave.photon_identity.native_emission.measure import four_h

K = int(candidate.KAPPA)
TIMES = (0, 1, 2, 4, 8)
SHOTS = 65536
SEED = 713401


def sector(board):
    n = board.n
    excess, vacancy = np.where(~np.eye(n, dtype=bool))
    idx = np.full((n,n),-1,dtype=int)
    idx[excess,vacancy] = np.arange(len(excess))
    adjacency = board.adjacency().toarray()
    touch = adjacency[excess,vacancy].astype(bool)
    contrast = 12+2*(board.degree[vacancy]-board.degree[excess])+2*touch
    degree = board.degree[excess]+board.degree[vacancy]-touch
    rows,cols = [],[]
    for a,b in board.seams:
        others = np.flatnonzero((np.arange(n)!=a)&(np.arange(n)!=b))
        first = np.r_[idx[a,others],idx[others,a],idx[a,b]]
        second = np.r_[idx[b,others],idx[others,b],idx[b,a]]
        rows.extend((first,second))
        cols.extend((second,first))
    rr,cc = np.concatenate(rows),np.concatenate(cols)
    hopping = sparse.coo_matrix((-0.5*np.ones(len(rr)),(rr,cc)),shape=(len(excess),len(excess))).tocsr()
    np.testing.assert_array_equal(np.diff(hopping.indptr),degree)
    kinetic = hopping+sparse.diags(degree/2)
    np.testing.assert_array_equal(kinetic@np.ones(len(excess)),0)
    assert np.all(contrast >= 10)
    return idx,excess,vacancy,touch,contrast,kinetic


def full_cell_reference():
    board = cell()
    idx,p,q,_,contrast,kinetic = sector(board)
    term = four_h(K)[:candidate.native_law.N_SRC,:candidate.native_law.N_SRC]/4
    labels = candidate.native_law.SRC
    check = np.zeros(kinetic.shape)
    for i,(excess,vacancy) in enumerate(zip(p,q)):
        loads = np.ones(board.n,dtype=int)
        loads[excess],loads[vacancy] = 2,0
        check[i,i] += K*((6-board.degree)@(loads**2-1))
        for a,b in board.seams:
            local = labels.index((loads[a],loads[b]))
            for j in np.flatnonzero(term[:,local]):
                out = loads.copy()
                out[a],out[b] = labels[j]
                outp,outq = np.flatnonzero(out==2).item(),np.flatnonzero(out==0).item()
                check[idx[outp,outq],i] += term[j,local]
    np.testing.assert_array_equal(check,kinetic.toarray()+K*np.diag(contrast))
    spectrum = np.linalg.eigvalsh(check)
    return {"all_columns_exactly_matching":len(p),"lowest_formation_energy":float(spectrum[0]),
            "highest_formation_energy":float(spectrum[-1])}


def singleton_check(board):
    """A single defect on unit filling and the registered one-load branch share kinetic hops."""
    terms = four_h(K)[:candidate.native_law.N_SRC,:candidate.native_law.N_SRC]/4
    free = candidate.generator(board).toarray()
    for q in (-1,1):
        # Subtract the uniform fixed-material background from each local contact and degree term.
        actual = np.diag(K*(6-board.degree)*(2*q+q*q)).astype(float)
        ix = [candidate.native_law.SRC.index(pair) for pair in ((1+q,1),(1,1+q))]
        local = terms[np.ix_(ix,ix)]
        for edge in board.seams:
            actual[np.ix_(edge,edge)] += local
        np.testing.assert_array_equal(actual,free+np.diag(2*q*K*(6-board.degree)))
    return {"signed_one_defect_columns_checked":2*board.n,
            "reference_energy_of_unit_filling":int(K*sum(6-board.degree))}


def launch_settings(board):
    rows = []
    for a,b in board.external:
        if a//12!=0:
            continue
        local = min(x for x in board.neighbours[a] if x//12==0 and board.degree[x]==5)
        rows.append((a,local,b//12))
    assert len(rows)==3
    return rows


def records(amplitude,p,q,touch,receiver,rng):
    probs = abs(amplitude)**2
    norm_error = float(abs(probs.sum()-1))
    assert norm_error < 3e-11
    a = np.where(p//12==receiver,p%12,12)
    b = np.where(q//12==receiver,q%12,12)
    outcome = a*13+b
    terminal = np.bincount(outcome,weights=probs,minlength=169)
    counts = rng.multinomial(SHOTS,terminal/terminal.sum()).reshape(13,13)
    two = (a<12)&(b<12)
    near = two & touch
    pair_indices = outcome[near]
    unique = np.unique(pair_indices)
    hit_count = counts.ravel()[unique].sum()
    return {"norm_error":norm_error,"joint_detector_counts":counts.tolist(),
            "receiver_probabilities":{"both_adjacent":float(probs[near].sum()),
                "both_nonadjacent":float(probs[two & ~touch].sum()),
                "excess_only":float(probs[(a<12)&(b==12)].sum()),
                "vacancy_only":float(probs[(a==12)&(b<12)].sum()),
                "neither":float(probs[(a==12)&(b==12)].sum())},
            "both_adjacent_count":int(hit_count),
            "global_adjacency_probability":float(probs[touch].sum()),
            "adjacent_in_any_cell_outside_home":float(probs[touch&(p//12==q//12)&(p//12!=0)].sum())}


def chebyshev_flight(h, kinetic_diagonal, sharp, time):
    """Direct late flight from launch with a scaled Chebyshev expansion."""
    low = min(h.diagonal()-kinetic_diagonal)
    high = max(h.diagonal()+kinetic_diagonal)
    center, radius = (high+low)/2, (high-low)/2
    z = (h-sparse.eye(h.shape[0],format="csr")*center)/radius
    order = int(np.ceil(radius*time+25*(radius*time)**(1/3)+40))
    # For real input both Chebyshev recurrences use real sparse matvecs.
    # Repeating with 16 further terms checks stabilisation of the truncated sum.
    n = np.arange(order+17)
    weights = jv(n,radius*time)*(-1j)**n
    f = sharp.real
    g = z@f
    out = weights[0]*f + 2*weights[1]*g
    for j in range(2,order+17):
        f,g = g, 2*(z@g)-f
        out += 2*weights[j]*g
        if j==order:
            first = out.copy()
    tail = float(np.max(abs(out-first)))
    assert tail < 1e-14
    return out*np.exp(-1j*center*time), order, tail


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    board=federation(0)
    idx,p,q,touch,contrast,kinetic=sector(board)
    launches=launch_settings(board)
    sharp=np.zeros((len(p),2*len(launches)),dtype=complex)
    for i,(gateway,interior,receiver) in enumerate(launches):
        sharp[idx[gateway,interior],2*i]=sharp[idx[interior,gateway],2*i+1]=1
    rng=np.random.default_rng(SEED)
    runs, late = {}, {}
    for k in (K,0):
        # Global phase shift removes the constant contrast cost 12 kappa during flight.
        h=kinetic+sparse.diags(k*(contrast-12),dtype=float)
        wave=sharp.copy()
        previous=0
        windows=[]
        for t in TIMES:
            if t>previous:
                wave=expm_multiply(-1j*(t-previous)*h,wave,traceA=-1j*(t-previous)*h.diagonal().sum())
            np.testing.assert_allclose(wave.conj().T@wave,np.eye(wave.shape[1]),atol=7e-12,rtol=0)
            settings=[]
            for i,(_,_,receiver) in enumerate(launches):
                plus,minus=wave[:,2*i:2*i+2].T
                for tag,psi in (("excess_at_gateway",plus),("vacancy_at_gateway",minus),
                                ("even_orientation",(plus+minus)/np.sqrt(2)),
                                ("odd_orientation",(plus-minus)/np.sqrt(2))):
                    settings.append({"station":i,"source":tag,**records(psi,p,q,touch,receiver,rng)})
            windows.append({"exposure":t,"settings":settings})
            previous=t
        runs[str(k)]=windows
        # The early survey chose station zero for independent, later read appointments.
        # Keep the other sharp orientation, the intake law, and fresh count draws.
        late_rng = np.random.default_rng(SEED+1+k)
        calibration, _, _ = chebyshev_flight(h,kinetic.diagonal(),sharp[:, :2],8)
        error = float(np.max(abs(wave[:, :2]-calibration)))
        assert error < 2e-10
        held = []
        for t in (16,32):
            held_wave, order, tail = chebyshev_flight(h,kinetic.diagonal(),sharp[:, :2],t)
            np.testing.assert_allclose(held_wave.conj().T@held_wave,np.eye(2),atol=3e-11,rtol=0)
            held.append({"exposure":t,"station":0,"polynomial_order":order, "max_change_16_extra_terms":tail,
                         "independent_expm_calibration_max_absolute_amplitude_error_at_8":error,"settings":[
                {"source":tag,**records(psi,p,q,touch,launches[0][2],late_rng)}
                for tag,psi in zip(("excess_at_gateway","vacancy_at_gateway"),held_wave.T)]})
            previous = t
        late[str(k)] = held
    print("Validated full cell projection and evolved all 57,360 ordered positions.",flush=True)
    result={"settings":{"cells":board.cells,"full_pair_dimension":len(p),"law":"registered recovered seam",
                "baseline":"one load per port, all seam memories blank","native_kappa":K,
                "receiver_alphabet":"0..11 are local ports; 12 is quiet, excess row and vacancy column",
                "trials_per_source_window":SHOTS,"seed":SEED,
                "launches":[{"gateway":a,"interior":z,"receiver_cell":r} for a,z,r in launches]},
            "one_cell_direct_all_native_terms":full_cell_reference(),
            "signed_one_defect_check_on_same_board":singleton_check(board),
            "energy_certificate":{"minimum_diagonal_contrast_cost_on_20_cells":int(K*min(contrast)),
                "all_boards_lower_bound_above_uniform":10*K,
                "isolated_pair_excess_energy_at_fixed_degrees_on_touching":2*K},
            "runs":runs,"later_selected_station":{"selection":"largest adjacency at exposure 8 in native early survey, with opposite sharp orientation as control",
                "late_count_seed_for_each_kappa":"early seed + 1 + kappa", "runs":late},"code_sha256":sha256(Path(__file__).read_bytes()).hexdigest()}
    with args.output.open("x") as out:
        json.dump(result,out,allow_nan=False,separators=(",",":"))
        out.write("\n")


if __name__=="__main__":
    main()
