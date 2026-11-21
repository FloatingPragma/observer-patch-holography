"""Read elapsed intervals in a mature cohort without reading its birth times.

The neutral return counter is an explicitly separate, equipped free channel.
Its own two-cycle triggers an appointment. Only terminal material port counts
reach the elapsed-time decoder.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from . import clock_factory as factory
from . import clock_inference as inference
from . import clock_loader as loader
from . import clock_records, clock_stream, loaded_clock
from .. import formation, law


HERE = Path(__file__).parent
ARCHIVE = HERE / "clock_comparison.npz"
FIXTURES = (((3,3,3,2), (0.75,0.5)), ((3,3,3,4), (0.375,0.875)))
# Appointment labels are furnished only to the scorer of the frequency.
# The elapsed decoder receives two unlabeled count arrays and the even-lapse domain.
RETURNS = (80, 81, 101, 239, 683, 1095, 1915, 2997, 4043)
EARLY_WINDOW = 16
WIRE_PERIOD = 8193
SHOTS = 1 << 36
SEED = 609173
ALPHA = 1e-6
MENUS = ("installed", "passive", "dephased")
PUBLISHED = "8f986c25a4ffc82804f74911c128b0d23c9135d823df06a6997f78c27fc14227"


def optical_counter(shape, station, appointments, *, tag="parity"):
    """Trigger on local private returns. A displaced pulse must give no return.

    For period two the return port is one full hardware site. For period four
    the same pointer reads the two even private positions without retaining
    which-position information. Reference and material occupy declared
    independent channels.
    """
    if shape[-1] not in (2,4) or tag not in ("parity", "displaced"):
        raise ValueError("the reference uses an even private ring of period two or four")
    if not appointments or list(appointments) != sorted(set(appointments)) or min(appointments) < 1:
        raise ValueError("book increasing positive optical returns")
    actual = station if tag == "parity" else ((station[0]+1)%shape[0], *station[1:])
    init = np.zeros((*shape,law.SLOTS), complex)
    init[(*actual[:3], slice(0,None,2), 6)] = 1/np.sqrt(shape[-1])
    init[(*actual[:3], slice(0,None,2), 7)] = -1/np.sqrt(shape[-1])
    one = law.step(init, loader.SOURCE)
    two = law.step(one, loader.SOURCE)
    if np.max(abs(two-init)) > 2e-13:
        raise ArithmeticError("two native free services failed to return the optical state")
    def incidence(state):
        local = state[(*station[:3],slice(0,None,2))]
        return float(np.vdot(local, local).real)
    pair = (incidence(one),incidence(two))
    if max(min(abs(p),abs(1-p)) for p in pair) > 2e-13:
        raise ArithmeticError("the optical pointer would disturb this preparation")
    # This is the reversible pointer addition on two exactly deterministic
    # native parity classes. The pointer starts at zero and cannot wrap here.
    level = 0
    records = {}
    for service in range(1, 2*max(appointments)+1):
        hit = int(round(pair[(service-1)%2]))
        level = (level+hit) % 8192
        if hit and level in appointments:
            records[level] = service
    return records, pair


def stations(shape):
    return ((0,0,0,0), (1,0,0,0), (1,1,0,0), (1,0,1,shape[-1]-1))


def _local_port_kernel(waves):
    return np.stack([np.sum(abs(waves[0]+a*waves[1])**2, axis=-1)/4 for a in loaded_clock.PORTS], axis=-1)


def source():
    rows = np.empty((len(MENUS),len(FIXTURES),4,len(RETURNS),7))
    triggers = np.empty((len(FIXTURES),4,len(RETURNS)), dtype=np.int32)
    enrolled_weights = []
    cursor_controls = []
    for fi, (shape,focus) in enumerate(FIXTURES):
        sites = stations(shape)
        writes, free, err = clock_stream.ordinary_history(shape, sites[0], focus, EARLY_WINDOW, formation.ANGLE)
        q = float(writes.sum())
        if abs(q+free[-1]-1) > 2e-11 or err > 2e-13:
            raise ArithmeticError("native first-write enrollment is inconsistent")
        enrolled_weights.append(q)
        for s,site in enumerate(sites):
            events, _ = optical_counter(shape,site,RETURNS)
            triggers[fi,s] = [events[c] for c in RETURNS]
        none, p = optical_counter(shape,sites[0],RETURNS,tag="displaced")
        cursor_controls.append({"displaced_triggers": len(none), "odd_even_power": list(p)})
        all_ticks = sorted(set(triggers[fi].ravel().tolist()))
        if max(all_ticks) >= WIRE_PERIOD or max(all_ticks) >= 8192 or min(all_ticks)-EARLY_WINDOW < 128:
            raise ValueError("the mature, pre-return clock domain does not hold")
        required = sorted({t-b for t in all_ticks for b in range(EARLY_WINDOW)})
        powers = loaded_clock.mature_powers(shape,max(required)-128)
        for m,menu in enumerate(MENUS[:2]):
            gates = factory.uniform_schedule() if menu == "installed" else {}
            waves = np.zeros((2,*shape,8),complex)
            waves[(slice(None),0,0,0,0)] = law.EQUAL/np.sqrt(2)
            for age in range(factory.END):
                waves = factory.advance(waves,age,gates)
            kernels = {}
            for age, outcome in zip(required, loaded_clock.at_appointments(waves,required,powers,start_age=factory.END)):
                kernels[age] = _local_port_kernel(outcome)
            for j,t in enumerate(all_ticks):
                # The independent natal wires record alternatives, without an
                # apparatus erase. Terminal stations receive a pooled cohort.
                words = sum(clock_stream._convolve(w,kernels[t-b]) for b,w in enumerate(writes))
                for s,site in enumerate(sites):
                    ports = words[site]
                    rows[m,fi,s,j] = np.r_[ports, q-ports.sum(), free[-1], 0.]
        rows[2,fi] = rows[0,fi]
        rows[2,fi,...,:4] = rows[0,fi,...,:4].mean(axis=-1,keepdims=True)
    return rows, triggers, enrolled_weights, cursor_controls


def draw(weights):
    counts = np.empty_like(weights, dtype=np.int64)
    for key in np.ndindex(weights.shape[:-1]):
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED,*key])))
        counts[key] = rng.multinomial(SHOTS, weights[key]/(1<<38))
    return counts


def _contains(bands, value):
    return any(low <= value <= high for low,high in bands)


def compare(words, counts, triggers, enrolled, controls):
    if (counts.shape != words.shape or counts.dtype != np.int64 or (counts < 0).any()
            or (counts > SHOTS).any() or (counts.sum(axis=-1) != SHOTS).any()):
        raise ValueError("invalid comparison counts or resources")
    rate = -loaded_clock.DELTA_MU
    family = int(np.prod(words.shape[:-1]))
    stats = {"enrolled": enrolled, "cursor_controls": controls, "trials": int(counts.sum()),
             "menus": {}, "failures": []}
    fail = stats["failures"]
    if any(c["displaced_triggers"] or any(c["odd_even_power"]) for c in controls):
        fail.append("displaced cursor")
    if abs(words.sum(axis=-1)-1).max() > 2e-11 or words.min() < -2e-11:
        fail.append("original-start closure")
    for m,menu in enumerate(MENUS):
        panel = counts[m].reshape(-1,len(RETURNS),7)
        ticks = triggers[0,0]
        if not np.all(triggers == ticks):
            raise ArithmeticError("station reference records disagree")
        options = dict(alpha=ALPHA,family_size=family,cell_error=2e-11)
        rates = inference.infer_rate(panel,ticks,**options)["per_service_intervals"]
        decoded = []
        for j in range(1,len(RETURNS)):
            data = inference.decode_elapsed(panel[:,0],panel[:,j],rate=rate,**options)
            runs = data["even_runs_step_two"]
            actual = int(ticks[j]-ticks[0])
            decoded.append({"reference_elapsed":actual,"decoded":runs,"number":data["eligible_even_count"]})
            if menu == "installed" and (data["eligible_even_count"] != 1 or not _contains(runs, actual)):
                fail.append(f"blind elapsed {j}")
            if menu == "dephased" and data["eligible_even_count"] != 4096:
                fail.append(f"unresolved dephased {j}")
        if menu == "installed" and not (_contains(rates,rate) and _contains(rates,rate+np.pi)):
            fail.append("source rate or even-sampling alias missing")
        if menu == "passive" and all(r["number"] == 1 and _contains(r["decoded"],r["reference_elapsed"]) for r in decoded):
            fail.append("passive comparison equally predictive")
        raw = words[m,...,:4]
        hand = raw[...,0]-raw[...,2]+1j*(raw[...,1]-raw[...,3])
        expected = hand[...,:1]*np.exp(1j*rate*(ticks-ticks[0]))
        residual = float(abs(hand-expected).max())
        if menu == "installed" and residual > 2e-11:
            fail.append("native two-service hand")
        stats["menus"][menu] = {"rate_intervals":rates,"elapsed":decoded,"max_quadrature_error":residual}
        if menu == "installed":
            individual = [[inference.decode_elapsed(station[:1],station[j:j+1],rate=rate,**options)["eligible_even_count"]
                           for j in range(1,len(ticks))] for station in panel]
            stats["menus"][menu]["individual_receiver_candidate_counts"] = individual
    return stats


def experiment(path=ARCHIVE, *, acquire=False, seed_replay=False):
    path = Path(path)
    blob = None if acquire else clock_records.snapshot(path, PUBLISHED)
    words, ticks, q, controls = source()
    if acquire:
        weights = np.stack([loaded_clock._grid(w) for w in words.reshape(-1,7)]).reshape(words.shape)
        blob = clock_records.snapshot(path, PUBLISHED, arrays={"weights": weights,"counts": draw(weights),"ticks": ticks})
    with np.load(io.BytesIO(blob), allow_pickle=False) as saved:
        if (sorted(saved.files) != ["counts","ticks","weights"] or saved["ticks"].dtype != ticks.dtype
                or not np.array_equal(ticks,saved["ticks"])):
            raise ValueError("incorrect optical comparison record")
        weights,counts = saved["weights"],saved["counts"]
    scale = 1 << 38
    if (weights.dtype != np.int64 or weights.shape != words.shape or (weights<0).any() or (weights>scale).any()
            or (weights.sum(axis=-1)!=scale).any() or abs(weights/scale-words).max()>2e-11):
        raise ValueError("recorded sampling weights do not match the complete quantum experiment")
    if counts.shape != weights.shape or np.any(counts[weights==0] != 0):
        raise ValueError("a forbidden comparison count was recorded")
    if seed_replay and not np.array_equal(counts, draw(weights)):
        raise ValueError("comparison counts disagree with the local seeded acquisition")
    result = compare(words, counts, ticks, q, controls)
    result["archive_sha256"] = hashlib.sha256(blob).hexdigest()
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire",action="store_true")
    parser.add_argument("--seed-replay",action="store_true")
    parser.add_argument("--archive",type=Path,default=ARCHIVE)
    args = parser.parse_args()
    result = experiment(args.archive,acquire=args.acquire,seed_replay=args.seed_replay)
    print(json.dumps(result,indent=2))
    raise SystemExit(int(bool(result["failures"])))
