"""Neutral dipole-write flight and a mobile-defect receiver control."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
sys.path.insert(0, str(ROOT))
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

LAW = candidate.native_law
K = int(candidate.KAPPA)
NM = 1 + LAW.N_LAB
TIMES = (0., .5, 1., 2., 4., 8.)
SHOTS = 8192
SEED = 18420


def four_h(k):
    """Four times the full native 99-state operator, over the integers."""
    w = LAW._write_pattern()
    r, d2 = LAW._label_ops()
    swap, cost = LAW._load_ops()
    twice_g = np.kron(np.eye(LAW.N_LAB, dtype=int) - r + 2*k*d2,
                      np.eye(LAW.N_SRC, dtype=int))
    twice_q = 2*np.eye(w.shape[0], dtype=int) - w @ w.T
    np.testing.assert_array_equal(twice_q @ twice_g, twice_g @ twice_q)
    full = np.zeros((LAW.DIM, LAW.DIM), dtype=int)
    full[:LAW.N_SRC, :LAW.N_SRC] = 2*(np.eye(LAW.N_SRC, dtype=int)-swap) + 4*k*cost
    full[LAW.N_SRC:, LAW.N_SRC:] = twice_q @ twice_g
    exchange, agreement = LAW.local_term_parts()
    np.testing.assert_allclose(full, 4*(exchange+k*agreement), atol=1e-13, rtol=0)
    return full


def columns(loads):
    return LAW.N_SRC*np.arange(NM) + LAW.SRC.index(tuple(loads))


def certification(h):
    out = {}
    for n in (0, 1, 2):
        ix = columns((n,n))
        leak = np.delete(h[:, ix], ix, axis=0)
        assert np.count_nonzero(leak) == 0
        out[str(n)] = {"memory_columns": len(ix), "off_block_nonzero_entries": 0,
                       "four_times_memory_hamiltonian": h[np.ix_(ix, ix)].tolist()}
    # Magnitude sectors are invariant even when loads are allowed to move.
    labels = np.repeat([-1, *(abs(d) for c,d in LAW.LABELS)], LAW.N_SRC)
    np.testing.assert_array_equal((labels[:, None] - labels[None, :])*h, 0)
    for target in (0, 1, 2):
        number = np.tile([(a == target)+(b == target) for a,b in LAW.SRC], NM)
        np.testing.assert_array_equal((number[:,None]-number[None,:])*h, 0)
    # A homogeneous port load has no current transition on any memory word.
    endpoint_number = np.tile([b for a,b in LAW.SRC],NM)
    current_over_i = (endpoint_number[None,:]-endpoint_number[:,None])*h
    for n in (0,1,2):
        np.testing.assert_array_equal(current_over_i[:,columns((n,n))],0)
    # A spectator prepared by writing a balanced pair has zero action at every load.
    idle = np.zeros(NM, dtype=int)
    idle[1+LAW.LABELS.index((0,0))] = idle[1+LAW.LABELS.index((1,0))] = 1
    np.testing.assert_array_equal(h @ np.kron(idle[:,None], np.eye(LAW.N_SRC, dtype=int)), 0)
    return out


def geography(board):
    a, b = min(e for e in board.external if e[0]//12 == 0)
    c = min(p for p in board.neighbours[b] if p//12 == b//12)
    d = min(p for p in board.neighbours[c] if p//12 == b//12 and p != b and p not in board.neighbours[b])
    source = tuple(sorted((a, min(p for p in board.neighbours[a] if p//12 == a//12))))
    return (a,b,c,d), source


def memory_vector(entries):
    out = np.zeros(NM, dtype=complex)
    for coin, contrast, amplitude in entries:
        out[1 + LAW.LABELS.index((coin, contrast))] = amplitude
    return out/np.linalg.norm(out)


def neutral_read(h):
    unitary = LAW.echo()
    input_pair = np.eye(LAW.DIM)[:, [LAW.SRC.index((2,0)), LAW.SRC.index((0,2))]]
    emitted = unitary @ input_pair
    uniform = columns((1,1))
    np.testing.assert_allclose(np.delete(emitted, uniform, axis=0), 0., atol=1e-14)
    bound = emitted[uniform]
    expected = np.column_stack([memory_vector([(0,d,1),(1,d,1)]) for d in (2,-2)])
    np.testing.assert_allclose(bound, expected, atol=1e-14)
    local = h[np.ix_(uniform, uniform)]/4
    np.testing.assert_allclose(local @ bound, 0., atol=1e-14)
    source_settings = {"positive": (1,0), "negative": (0,1),
                       "linear_superposition": np.array((1,1))/np.sqrt(2),
                       "phase_superposition": np.array((1,1j))/np.sqrt(2)}
    receiver_zero = memory_vector([(0,0,1),(1,0,-1)])
    receiver_excited = memory_vector([(0,2,1),(0,-2,1),(1,2,-1),(1,-2,-1)])
    ramsey = (receiver_zero + receiver_excited)/np.sqrt(2)
    receiver_blank = np.eye(NM)[:,0]
    results = []
    rng = np.random.default_rng(SEED)
    for t in (0., np.pi/(4*K), np.pi/(2*K), 1., 4.):
        u = expm(-1j*t*local)
        us = u @ bound
        np.testing.assert_allclose(us, bound, atol=2e-13)
        unwrite = unitary[:, uniform] @ us
        np.testing.assert_allclose(unwrite, input_pair, atol=3e-13)
        u_pair = np.kron(u,u)
        stats = {}
        for name, amplitudes in source_settings.items():
            s = bound @ amplitudes
            data = {}
            for recv_name, v in (("blank", receiver_blank), ("active_Ramsey", ramsey)):
                psi = (u_pair @ np.kron(s,v)).reshape(NM,NM)
                rho = psi.T @ psi.conj()
                reference = np.outer(u @ v, (u @ v).conj())
                np.testing.assert_allclose(rho, reference, atol=4e-14)
                pointer = v
                yes = np.clip(float(np.vdot(pointer,rho @ pointer).real),0.,1.)
                yes_count = int(rng.binomial(SHOTS,yes))
                data[recv_name] = {"memory_probabilities": rho.diagonal().real.tolist(),
                                   "pointer_yes_probability": yes, "yes_no_counts": [yes_count,SHOTS-yes_count]}
            stats[name] = data
        expected_ramsey = np.cos(2*K*t)**2
        for data in stats.values():
            np.testing.assert_allclose(data["active_Ramsey"]["pointer_yes_probability"], expected_ramsey, atol=1e-13)
            np.testing.assert_allclose(data["blank"]["memory_probabilities"], np.eye(NM)[0], atol=1e-13)
        results.append({"exposure":t,"source_settings":stats})
    direct = abs(expected.conj().T @ bound)**2
    np.testing.assert_allclose(direct,np.eye(2),atol=2e-14)
    return {"endpoint_loads_after_repair":[1,1],"receiver_direct_write_sign_read":direct.tolist(),
            "source_inverse_echo_roundtrip_verified_for_all_four_settings":True,
            "source_mark_energy_in_local_units":0, "active_receiver_gap":4*K,
            "terminal_reads":results}


def active_path(h, board, route, q, park, k):
    """Restrict native terms to all defect positions and all receiving collar states."""
    defects = np.eye(4, dtype=int)*q + 1
    matrix = np.kron(np.diag(k*(2*q + q*q)*(6-board.degree[list(route)])),np.eye(NM))
    for edge_slot, (x,y) in enumerate(zip(range(3),range(1,4))):
        if park and edge_slot == 0:
            continue  # zero term by the exact idle-write identity
        # Orient port order consistently with the microscopic stored edge.
        pair = [x,y] if route[x] < route[y] else [y,x]
        outside = [p for p in range(4) if p not in pair]
        current = [columns(load[pair]) for load in defects]
        for i, ix in enumerate(current):
            destinations = [j for j in range(4) if np.array_equal(defects[j, outside],defects[i,outside])]
            if edge_slot == 1:
                # Test the complete outgoing columns, including unrepresented load configurations.
                allowed = np.concatenate([current[j] for j in destinations])
                assert np.count_nonzero(np.delete(h[:,ix], allowed, axis=0)) == 0
            else:
                expected = [current[j][0] for j in destinations]
                assert np.count_nonzero(np.delete(h[:,ix[0]], expected)) == 0
            for j in destinations:
                if edge_slot == 1:
                    block = h[np.ix_(current[j],ix)]/4
                else:
                    block = h[current[j][0],ix[0]]/4*np.eye(NM)
                matrix[j*NM:(j+1)*NM,i*NM:(i+1)*NM] += block
    np.testing.assert_array_equal(matrix, matrix.T)
    return matrix


def reachable_records(operator, start_coin):
    """Exact connectivity test in the sign-even unit-contrast subspace."""
    coin_basis = np.zeros((NM, 2), dtype=int)
    for coin in (0, 1):
        for sign in (-1, 1):
            coin_basis[1+LAW.LABELS.index((coin, sign)), coin] = 1
    embedding = np.kron(np.eye(4, dtype=int), coin_basis)
    four = (4*operator).astype(int)
    np.testing.assert_array_equal(4*operator, four)
    small = embedding.T @ four @ embedding // 2
    np.testing.assert_array_equal(four @ embedding, embedding @ small)
    seen = {start_coin}
    todo = [start_coin]
    while todo:
        here = todo.pop()
        for other in np.flatnonzero(small[:, here]):
            if int(other) not in seen:
                seen.add(int(other))
                todo.append(int(other))
    return [[i//2, i%2] for i in sorted(seen)]


def material_control(board, route, full, exchange):
    rng = np.random.default_rng(SEED + 1)
    rows = []
    for q in (-1,1):
        for arm, park, k, h in (("free",False,K,full),("parked",True,K,full),
                                ("exchange_only",False,0,exchange)):
            operator = active_path(h,board,route,q,park,k)
            eigenvalues, modes = np.linalg.eigh(operator)
            for coin in (0,1):
                source = np.kron(np.eye(4)[:,0],memory_vector([(coin,-1,1),(coin,1,1)]))
                flipped = [1+j for j,(c,d) in enumerate(LAW.LABELS) if c!=coin]
                reachable = reachable_records(operator, coin)
                if park or k == 0:
                    assert [3, 1-coin] not in reachable
                windows = []
                for t in TIMES:
                    psi = modes @ (np.exp(-1j*t*eigenvalues)*(modes.conj().T @ source))
                    probs = (abs(psi)**2).reshape(4,NM)
                    np.testing.assert_allclose(probs.sum(),1.,atol=1e-13)
                    counts = rng.multinomial(SHOTS, probs.ravel()/probs.sum()).reshape(4,NM)
                    hits = float(probs[3,flipped].sum())
                    if park or k == 0:
                        assert hits < 1e-26
                    windows.append({"exposure":t,"probabilities":probs.tolist(),"counts":counts.tolist(),
                                    "flipped_coin_and_defect_beyond":hits,
                                    "flipped_coin_and_defect_beyond_count":int(counts[3,flipped].sum())})
                rows.append({"charge_change":q,"initial_receiver_coin":coin,"arm":arm,"exact_connected_position_coin_pairs":reachable,
                             "terminal_reads":windows})
        active = [r for r in rows if r["arm"]=="free" and r["charge_change"]==q]
        strengths = [max(w["flipped_coin_and_defect_beyond"] for w in r["terminal_reads"]) for r in active]
        assert min(strengths)<1e-25 and max(strengths)>0.0001
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    board = federation(0)
    route, source = geography(board)
    h, exchange = four_h(K), four_h(0)
    result = {"law":"registered Interval Echo free Hamiltonian after addressed source preparation",
              "cells":board.cells,"number_physical_seams":len(board.seams),"strength":K,
              "uniform_background_certificate":certification(h),"charge_route":list(route),
              "dipole_source_seam":list(source),"independent_receiver_seam":sorted(route[1:3]),
              "receiver_basis":["blank", *[[c,d] for c,d in LAW.LABELS]],
              "all_endpoint_occupancy_counts_commute_on_all_99_states":True,
              "endpoint_load_current_annihilates_all_33_uniform_columns":True,
              "shots_per_setting":SHOTS,"seeds":[SEED, SEED+1],
              "neutral_emission":neutral_read(h),
              "four_port_material_controls":material_control(board, route, h, exchange),
              "source_sha256":{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()
                               for p in (HERE/"DESIGN.md",Path(__file__),Path(LAW.__file__))}}
    with args.output.open("x") as out:
        json.dump(result,out,indent=2,allow_nan=False)
        out.write("\n")


if __name__ == "__main__":
    main()
