"""Cross-record counts and angular integrations for the finite two-clock channel."""

import argparse
from fractions import Fraction as F
import json
from math import isclose, log, sqrt
from pathlib import Path

import numpy as np

from .channel import ClockPair, LagGrid

HERE = Path(__file__).resolve().parent
WORDS = ((-1,-1), (-1,1), (1,-1), (1,1))


def conditional_counts(pair, ga, gb, sign):
    """Enumerate the actual two hidden signs and both local attenuation coins."""
    t = pair.tickets
    result = [0]*4
    for first, second in WORDS:
        # Shared table entries supply one fair bit, independent entries supply two.
        multiplicity = 1 if sign == 0 else 2*int(second == sign*first)
        for i, (a,b) in enumerate(WORDS):
            result[i] += multiplicity*(t+a*first*ga)*(t+b*second*gb)
    return result


def local_word_counts(pair, ta, tb, scale_a, scale_b, *, shared_lag=True):
    """Terminal urn from all source-ticket classes, phases, table bits and gate coins."""
    t, grid = pair.tickets, pair.grid
    ga, gb = pair.gate(scale_a), pair.gate(scale_b)
    free = conditional_counts(pair, ga, gb, 0)
    tied = conditional_counts(pair, ga, gb, -1 if pair.weight >= 0 else 1)
    base = conditional_counts(pair, t, t, 1)
    lanes = [(i,pair.class_ticks(i)) for i in range(-grid.max_lag,grid.max_lag+1)]
    independent = t//2-sum(count for _, count in lanes)
    counts = [0]*4
    for z in range(grid.subphases):
        ia = (ta+z)//grid.subphases
        ib = (tb+z)//grid.subphases
        for i in range(4):
            counts[i] += (t//2)*base[i] + independent*free[i]
        for offset, size in lanes:
            frequencies = tied if shared_lag and ib-offset == ia else free
            for i in range(4):
                counts[i] += size*frequencies[i]
    if sum(counts) != grid.subphases*16*t**3:
        raise ArithmeticError("The word counts must include every source class and gate record")
    return counts


def bank_measurements(frozen=None):
    grid = LagGrid(3,6,8)
    source, blank = ClockPair(grid, F(1,256)), ClockPair(grid, F(0))
    cfg = {"vacuum": (blank, True, 1,1), "ordinary": (source, True, 1,1),
           "attenuated": (source, True, F(9,16),F(4,9)), "independent_tables": (source,False,1,1)}
    banks = {}
    rng = np.random.default_rng(26112633)
    # Conservation of every local sign and of every failed match is tested below.
    starts = 1 << 25
    radius = sqrt(log(2*4*len(cfg)/0.01)/(2*starts))
    for name, (material,shared,ga,gb) in cfg.items():
        counts = local_word_counts(material, 11, 45, ga, gb, shared_lag=shared)
        total = sum(counts)
        if counts[0]+counts[1] != total//2 or counts[0]+counts[2] != total//2:
            raise ArithmeticError("One-ended census departed from a fair sign")
        c = F(counts[0]-counts[1]-counts[2]+counts[3], total)
        model = (F(1,2) + material.shift(F(11,64),F(45,64),1,ga,gb)) if shared else F(1,2)
        if c != model:
            raise ArithmeticError("Independent local table enumeration failed the lag law")
        probs = np.array([float(F(n,total)) for n in counts])
        sample = list(map(int,rng.multinomial(starts,probs))) if frozen is None else frozen[name]["sample"]
        if len(sample) != 4 or sum(sample) != starts or min(sample) < 0:
            raise ValueError("The terminal bank must retain all four sign outcomes")
        if max(abs(np.array(sample)/starts-probs)) > radius:
            raise ValueError("An unfiltered sign count lies outside the simultaneous 99 percent intervals")
        banks[name] = {"correlation": float(c), "exact_words": list(map(str,counts)), "sample": sample}
    return {"starts_per_setting": starts, "terminal_words": [list(w) for w in WORDS],
            "frequency_radius": radius, "banks": banks}


def angular_holdouts():
    # Neither the angles nor velocities occur in the specification of the clock source.
    nodes, weights = np.polynomial.legendre.leggauss(32)
    phi = np.arange(64)*(2*np.pi/64)
    rings = np.sqrt(1-nodes*nodes)[:,None]
    unit = np.stack(np.broadcast_arrays(rings*np.cos(phi), rings*np.sin(phi), nodes[:,None]),axis=-1)
    settings = {"rest": (0,0,0), "oblique": (0.31,-0.44,0.17), "transverse": (-0.51,0.28,0.43),
                "fast": (0.69,-0.31,0.29)}
    results = {}
    for name, velocity in settings.items():
        v = np.asarray(velocity)
        x = unit @ v
        gamma = 1/sqrt(1-v@v)
        projected_chord = gamma/np.sqrt(1-x*x)
        levels = []
        for bits in (4,6,8):
            bins = 1 << bits
            weight = F(1, 1 << (bits+9))
            pair = ClockPair(LagGrid(bits,7*bins//8), weight)
            # A common conformal-time origin away from the source clock origin.
            fine = pair.grid.subphases*bins
            ta, tb = np.floor(0.273*fine), np.floor((0.273+x)*fine)
            signal = pair.entropy_from_tick_differences((tb-ta).astype(int))/float(weight)
            answer = float(np.sum(weights*np.mean(signal*projected_chord,axis=1))/2)
            levels.append({"bits": bits, "measured": answer, "absolute_error": abs(answer-gamma**2)})
        flat = float(np.sum(weights*np.mean(projected_chord,axis=1))/2)
        results[name] = {"velocity": list(velocity), "target": gamma**2, "flat_lag_control": flat, "levels": levels}
    return results


def run(frozen=None):
    return {"schema":"oph.clock_pair_lag.v1", "sign_readings":bank_measurements(frozen),
            "tangent_ellipsoid_holdouts":angular_holdouts()}


def verify():
    saved = json.loads((HERE/"results.json").read_text())
    new = run(saved["sign_readings"]["banks"])
    records, previous = new["sign_readings"], saved["sign_readings"]
    if any(records[k] != previous[k] for k in ("banks","starts_per_setting","terminal_words")):
        raise ArithmeticError("Frozen word counts or terminal calibration changed")
    if new["schema"] != saved["schema"] or not isclose(records["frequency_radius"],previous["frequency_radius"],rel_tol=1e-14):
        raise ArithmeticError("Frozen coverage or schema changed")
    for name, row in new["tangent_ellipsoid_holdouts"].items():
        old = saved["tangent_ellipsoid_holdouts"][name]
        if row["velocity"] != old["velocity"] or not isclose(row["target"], old["target"], abs_tol=2e-14):
            raise ArithmeticError("Registered boost changed")
        if abs(row["flat_lag_control"]-old["flat_lag_control"]) > 1e-12:
            raise ArithmeticError("Flat-lag control changed")
        for a,b in zip(row["levels"],old["levels"],strict=True):
            if a["bits"] != b["bits"] or abs(a["measured"]-b["measured"]) > 2e-10:
                raise ArithmeticError("Held boost returned a different endpoint entropy")
    return new


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    if parser.parse_args().write:
        (HERE/"results.json").write_text(json.dumps(run(),indent=2,allow_nan=False)+"\n")
    else:
        verify()
        print("Verified four complete two-clock count banks and four held boost directions.")
