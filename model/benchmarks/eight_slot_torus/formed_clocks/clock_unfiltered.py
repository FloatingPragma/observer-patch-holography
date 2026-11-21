"""Test two-service phase advance when the station admits every first-write history.

The authenticated 193-service stream on the 3**4 chart supplies the original
first wires. Fresh ideal local measurements use the external station service
counts (129,131) and (191,193). Each four-port reader accepts every material
age. A partition by birth history is used only to prove the residual bound.
These words are quantum probabilities. They contain no sampled counts.
"""
import argparse
import io
import json

import numpy as np

from . import clock_records, clock_stream


INTERVALS = ((129, 131), (191, 193))
SHAPE = (3, 3, 3, 3)
MATURITY = 128
DELTA = np.pi/4096
TOL = 2e-11
OUTCOMES = ("+x", "+y", "-x", "-y", "marked at another private position of this public address",
            "marked at another public address", "unmarked")


def hand(ports):
    return ports[...,0]-ports[...,2]+1j*(ports[...,1]-ports[...,3])


def original_words(ports, unmarked):
    """Seven mutually exclusive outcomes for a station at each complete address."""
    local = ports.sum(axis=-1)
    public = local.sum(axis=-1,keepdims=True)
    return np.concatenate([ports, (public-local)[...,None],
                           np.broadcast_to(1-unmarked-public,local.shape)[...,None],
                           np.full((*local.shape,1),unmarked)],axis=-1)


def _local(writes, impulse, tick, begin, end):
    # Distinct (wire, birth time) alternatives carry orthogonal marker records.
    ports = np.zeros((*SHAPE,4))
    for birth in range(begin, min(tick,end)):
        ports += clock_stream._convolve(writes[birth],impulse[tick-birth])
    return ports


def evaluate(writes, unmarked, impulse):
    results = []
    for t,u in INTERVALS:
        phases = np.exp(-1j*(u-t)*DELTA)
        stop = t-MATURITY+1
        all0,all1 = (_local(writes,impulse,a,0,a) for a in (t,u))
        old0,old1 = (_local(writes,impulse,a,0,stop) for a in (t,u))
        young0,young1 = (_local(writes,impulse,a,stop,a) for a in (t,u))
        first0,first1 = (_local(writes,impulse,a,1,2) for a in (t,u))
        residual = hand(all1)-phases*hand(all0)
        young_residual = hand(young1)-phases*hand(young0)
        bound = (young1.sum(axis=-1)+young0.sum(axis=-1))/2
        # Each mature impulse returns its spatial ray after two services,
        # and its signed hand picks up the same known rail phase. The
        # triangle inequality and |hand| <= incidence/2 bound every other
        # first-wire history, including both intervening birth opportunities.
        records = np.stack([original_words(all0,unmarked[t]),original_words(all1,unmarked[u])])
        norm_error = max(abs(all0.sum()+unmarked[t]-1),abs(all1.sum()+unmarked[u]-1))
        earliest = float(np.max(abs(hand(first1)-phases*hand(first0))))
        mature = float(np.max(abs(hand(old1)-phases*hand(old0))))
        decomposition = float(np.max(abs(residual-young_residual)))
        excess = float(np.max(abs(residual)-bound))
        if (records.min() < -TOL or abs(records.sum(axis=-1)-1).max() > TOL or norm_error > TOL
                or max(earliest,mature,decomposition,excess) > TOL):
            raise ArithmeticError("complete records or the mature/young decomposition failed")
        nonzero = (abs(hand(all0)) > TOL) & (abs(hand(all1)) > TOL)
        increments = np.angle(hand(all1)[nonzero]*hand(all0)[nonzero].conj())
        summary = {"station_ticks":[t,u], "source":"native original-start probabilities",
                   "tested_sites":int(np.prod(SHAPE)), "sites_outside_stationary_tolerance":int(np.sum(abs(residual)>TOL)),
                   "max_unfiltered_residual":float(abs(residual).max()), "max_young_power_bound":float(bound.max()),
                   "max_residual_minus_bound":excess, "max_decomposition_error":decomposition,
                   "max_mature_residual":mature, "max_first_birth_residual":earliest,
                   "first_birth_mass":float(writes[1].sum()), "mature_cohort_mass":float(writes[:stop].sum()),
                   "intervening_birth_mass":float(writes[t:u].sum()), "resolved_phases":int(nonzero.sum()),
                   "min_signed_increment":float(increments.min()) if len(increments) else None,
                   "max_signed_increment":float(increments.max()) if len(increments) else None,
                   "mature_signed_increment":float(-(u-t)*DELTA)}
        results.append((summary,records,residual,bound))
    return results


def experiment():
    blob = clock_records.snapshot(clock_stream.ARCHIVE,clock_stream.PUBLISHED)
    with np.load(io.BytesIO(blob),allow_pickle=False) as saved:
        writes, unmarked = saved["writes"], saved["unmarked"]
    if writes.shape != (193,*SHAPE) or unmarked.shape != (194,):
        raise ValueError("the enrolled 193-service first-wire stream is required")
    impulse,_ = clock_stream.impulse(SHAPE,193,True)
    return evaluate(writes,unmarked,impulse)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all-sites", action="store_true", help="include every original-start probability word and bound")
    args = parser.parse_args()
    report = []
    for summary,words,r,bound in experiment():
        if args.all_sites:
            summary.update(outcomes=OUTCOMES, original_start_words=words.tolist(),
                           residual_real=r.real.tolist(), residual_imag=r.imag.tolist(),
                           young_power_bound=bound.tolist())
        report.append(summary)
    print(json.dumps(report,indent=2,allow_nan=False))
