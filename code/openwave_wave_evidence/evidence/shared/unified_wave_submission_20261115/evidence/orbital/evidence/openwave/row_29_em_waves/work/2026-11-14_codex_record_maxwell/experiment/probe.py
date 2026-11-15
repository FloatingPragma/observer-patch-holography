#!/usr/bin/env python3
"""Preregistered finite native-current and compressed population experiment."""
import gzip
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import random
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
SHAPE = (3, 3, 3, 8)
CHANNEL_STEPS = 16384
SEED = 2026111429
BLOCKS = (10_000_000, 1_000_000_000)
FACES = ((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),
         (1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
         (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),
         (4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))


def encode(obj):
    return (json.dumps(obj, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode()


def load_native():
    path = HERE / "inputs/native.py"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == \
        "dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d"
    spec = importlib.util.spec_from_file_location("frozen_maxwell_native", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod.native


def main():
    native = load_native()
    cells = list(itertools.product(*(range(n) for n in SHAPE)))
    ci = {z: i for i, z in enumerate(cells)}
    nc = len(cells)
    nv = 12 * nc
    internal = sorted({tuple(sorted((f[j], f[(j+1) % 3])))
                       for f in FACES for j in range(3)})
    assert len(internal) == 30
    edges = []
    for z in cells:
        c = ci[z]
        edges.extend((12*c+a, 12*c+b, -1) for a, b in internal)
        for axis in range(4):
            zz = list(z)
            zz[axis] = (zz[axis]+1) % SHAPE[axis]
            d = ci[tuple(zz)]
            edges.extend((12*c+p, 12*d+p, axis) for p in range(12))
    w = len(edges)
    assert (nc, nv, w) == (216, 2592, 16848)
    assert len({tuple(sorted((a,b))) for a,b,_ in edges}) == w
    branches = []
    displacement = {}
    full_checks = 0
    for eid, (a,b,axis) in enumerate(edges):
        spectator = next(i for i in range(3) if i not in (a,b))
        for slot in range(3):
            loads = [0,0,0]
            loads[slot] = 1
            for coin in (0,1):
                out, (returned_coin, current) = native(loads, (0,1), coin)
                rebuilt = list(loads)
                rebuilt[0] -= current
                rebuilt[1] += current
                assert tuple(rebuilt) == out and returned_coin == coin
                target = out.index(1)
                assert sum(out) == 1
                branches.append([eid, slot, coin, current, target])
                displacement[eid,slot,coin] = current
                if eid < 64:
                    global_loads = [0] * nv
                    mapping = (a,b,spectator)
                    global_loads[mapping[slot]] = 1
                    full_out, rec = native(global_loads, (a,b), coin)
                    assert rec == (coin,current)
                    assert full_out.index(1) == mapping[target]
                    full_checks += 1
    # Actual branch records determine every off-diagonal macro count. The
    # spectator census and untouched-source algebra compress identity branches.
    transfers = np.zeros((12,nc,nc), dtype=np.int64)
    for eid,(a,b,axis) in enumerate(edges):
        if axis < 0:
            continue
        for slot,source in ((0,a),(1,b)):
            for coin in (0,1):
                j = displacement[eid,slot,coin]
                target = (b if j == 1 else a) if slot == 0 else \
                         (a if j == -1 else b)
                if target != source:
                    transfers[source % 12,source//12,target//12] += 1
    for port in range(12):
        transfers[port] += np.diag(2*w-transfers[port].sum(axis=1))
    assert all(np.array_equal(transfers[0],m) for m in transfers)
    counts = transfers[0]
    assert np.array_equal(counts, counts.T)
    assert np.all(counts.sum(axis=1) == 2*w)
    lap = 2*w*np.eye(nc) - counts
    vals, vecs = np.linalg.eigh(lap)
    assert vals[0] > -1e-10 and vals[1] > 0.5
    assert np.max(np.abs(lap@vecs-vecs*vals)) < 1e-10

    # A finite actual tape, separate from the population-compression calculation.
    rng = random.Random(SEED)
    tapes, initial, source_counts = [], [], []
    for mu in range(4):
        pos = rng.randrange(nv)
        initial.append(pos)
        records = []
        observed = pos
        occupancy = np.zeros(nc, dtype=np.int64)
        for _ in range(CHANNEL_STEPS):
            eid = rng.randrange(w)
            coin = rng.randrange(2)
            a,b,_ = edges[eid]
            out,rec = native((int(pos == a), int(pos == b)), (0,1), coin)
            j = rec[1]
            if j == 1:
                assert pos == a and observed == a
                pos,observed = b,b
            elif j == -1:
                assert pos == b and observed == b
                pos,observed = a,a
            assert out == (int(pos == a),int(pos == b))
            records.append([eid,rec[0],j])
            occupancy[observed//12] += 1
        tapes.append(records)
        source_counts.append(occupancy.tolist())
    rebuilt = []
    for start, tape in zip(initial,tapes):
        pos = start
        occ = [0] * nc
        for eid,coin,j in tape:
            a,b,_ = edges[eid]
            if j:
                assert pos == (a if j == 1 else b)
                pos = b if j == 1 else a
            occ[pos//12] += 1
        rebuilt.append(occ)
    assert rebuilt == source_counts
    a_fields = math.sqrt(nc/(4*w*CHANNEL_STEPS)) * \
        (np.asarray(rebuilt)-CHANNEL_STEPS/nc)
    a_fields = a_fields.reshape((4,)+SHAPE)
    directions = (3,0,1,2)
    f_fields = []
    for mu in range(4):
        for nu in range(mu+1,4):
            f_fields.append((np.roll(a_fields[nu],-1,axis=directions[mu])
                             -a_fields[nu]-np.roll(a_fields[mu],-1,
                             axis=directions[nu])+a_fields[mu]).ravel().tolist())

    # This is the exact finite-block population formula evaluated numerically
    # from the transfer census, not a Monte Carlo covariance estimate.
    spectra = []
    for n in BLOCKS:
        g = np.zeros(nc)
        lam = vals[1:]
        z = 1-lam/(2*w)
        correction = z*(-np.expm1(n*np.log(z))) / \
                     (2*w*n*(1-z)**2)
        g[1:] = 1/lam-1/(4*w)-correction
        cov = (vecs*g) @ vecs.T
        for spatial in ((1,0,0),(1,1,0),(1,1,1)):
            p = np.array(spatial)*2*np.pi/3
            k = 2*np.sin(p/2)
            q2 = float(k@k)
            transverse = np.eye(3)-np.outer(k,k)/q2
            wave = np.array([np.exp(1j*sum(p[j]*z[j] for j in range(3)))
                             /math.sqrt(27) for z in cells])
            temporal = []
            v0 = wave * np.array([z[3] == 0 for z in cells])
            for tau in range(4):
                vt = wave * np.array([z[3] == tau for z in cells])
                value = np.vdot(vt,cov@v0)
                assert abs(value.imag) < 2e-10
                temporal.append(float(value.real))
            omega = 2*math.asinh(math.sqrt(q2)/2)
            target = [math.cosh(omega*(4-tau)) /
                      (2*math.sinh(omega)*math.sinh(4*omega))
                      for tau in range(4)]
            t1,t2,t3 = temporal[1:4]
            recurrence = (t1+t3-(2+q2)*t2) / ((2+q2)*t2)
            fitted = math.acosh((t1+t3)/(2*t2))
            bb = q2*t2*transverse
            minus_ee = -(2*t2-t1-t3)*np.eye(3)-np.outer(k,k)*t2
            b_eigs = np.linalg.eigvalsh(bb)
            e_eigs = np.linalg.eigvalsh(minus_ee)
            # IID-position control P=11^T/nc, same normalization and curl.
            control_cov = (np.eye(nc)-np.ones((nc,nc))/nc)/(4*w)
            v2 = wave*np.array([z[3] == 2 for z in cells])
            control = float(np.vdot(v2,control_cov@v0).real)
            if n == BLOCKS[-1]:
                assert max(abs(temporal[t]-target[t])/target[t]
                           for t in (1,2,3)) < 0.01
                assert abs(recurrence) < 0.01
                assert abs(fitted/omega-1) < 0.01
                assert abs(b_eigs[0]) < 1e-10
                assert min(b_eigs[1:]) > 0
                assert max(abs(e_eigs-b_eigs)) / b_eigs[-1] < 0.01
                assert abs(control) < 1e-12 and t2 > 1e-5
            spectra.append(dict(block_events=n,spatial=list(spatial),
                native_slice_covariance=temporal,target_green=target,
                finite_block_recurrence_relative_error=recurrence,
                fitted_slice_energy=fitted,target_slice_energy=omega,
                magnetic_eigenvalues=b_eigs.tolist(),
                reflected_electric_eigenvalues=e_eigs.tolist(),
                changed_law_separated_covariance=control))
    data = dict(shape=SHAPE,edges=edges,branches=branches,
                channel_initial=initial,tapes=tapes,occupation_counts=rebuilt,
                field_strength_readouts=f_fields)
    raw = encode(data)
    (HERE/"RECORDS.json.gz").write_bytes(gzip.compress(raw,mtime=0))
    receipt = dict(schema="oph.record_maxwell.finite.v1",
        graph=dict(cells=nc,vertices=nv,seams=w,shape=SHAPE),
        actual_native_calls=len(branches)+full_checks+4*CHANNEL_STEPS,
        endpoint_branch_calls=len(branches),full_vector_checks=full_checks,
        actual_tape_events=4*CHANNEL_STEPS,
        actual_tape_nonzero_currents=sum(bool(r[2]) for t in tapes for r in t),
        record_only_replay=True,macro_lumping_equal_for_all_ports=True,
        population_compression="finite-block covariance geometric sum from branch-derived transfer",
        compressed_population_block_lengths=BLOCKS,
        population_blocks_are_actual_long_tapes=False,
        raw_uncompressed_sha256=hashlib.sha256(raw).hexdigest(),
        spectra=spectra,checks="all preregistered checks pass",
        physical_time_identified=False,official_row_status="untested")
    (HERE/"RECEIPT.json").write_bytes(encode(receipt))
    print(json.dumps({k:v for k,v in receipt.items() if k != "spectra"},
                     indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
