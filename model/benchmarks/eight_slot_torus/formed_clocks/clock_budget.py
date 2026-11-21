"""Shot costs for distinguishing adjacent even lapses and a fixed one-mouth reader.

The two original fixture geometries and the first-write population at their
second launch mouth determine the enrollment. The comparison cohort's
source and results were available when this design was chosen. Its counts
are never read by the acquisition. Every menu books that mouth on each
fixture separately at each of the nine original optical appointments.
"""
import argparse
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
from mpmath.ctx_iv import MPIntervalContext

from . import clock_comparison as primary
from . import clock_inference as reader

HERE = Path(__file__).parent
ARCHIVE = HERE / "clock_budget_counts.npz"
SOURCE_DIGEST = "8f986c25a4ffc82804f74911c128b0d23c9135d823df06a6997f78c27fc14227"
SEED = 8369251
SHOTS = 2**34
FAMILY = 54
ALPHA = 5e-7
ERROR = 2e-11
GRID = 2**38
STATION = 1


def design():
    """Verify the source incidence floor and simultaneous one-mouth precision."""
    iv = MPIntervalContext()
    iv.dps = 50
    theta = sum(iv.pi*(1-iv.cos(iv.pi*(j+1)/128))/4 if j < 128 else iv.pi/2
                for j in range(35,159))
    even = iv.cos(theta)**2
    # Only the first two literal writes, at service one, are needed.
    for _, (fa,fb) in primary.FIXTURES:
        first = (fa*fb*(1-even) + (1-fa)*(1-fb)*even)/4
        if not first > iv.mpf(".065"):
            raise ArithmeticError("first-record incidence floor failed")
    s = iv.mpf(".065")
    delta = iv.pi/4096
    visibility = iv.cos(15*delta/2)
    k = 557787402
    l = iv.log(6*FAMILY/iv.mpf("0.0000005"))
    error = iv.mpf("0.00000000002")
    incidence_lower = iv.mpf(k)/SHOTS - iv.sqrt(l/(2*SHOTS))
    coordinate = iv.sqrt(2*l/k)+4*error/incidence_lower+iv.mpf("0.00000000004")
    offset = iv.sqrt(2)*coordinate
    ratio = offset/(visibility-offset)
    angle = iv.atan2(ratio,iv.sqrt(1-ratio**2))
    spacing = 2*delta
    mean = SHOTS*(s/2-2*error)
    outage = 2*FAMILY*iv.exp(-mean*(1-k/mean)**2/2)
    # Angular and modular bookkeeping in the reader rounds its limits outwards.
    if not (4*angle+iv.mpf("0.00000001") < iv.mpf(".001527") < spacing
            and outage < iv.mpf("1e-110")):
        raise ArithmeticError("one-mouth budget cannot guarantee a unique even lapse")
    return {"incidence_lower":.065, "visibility_lower":float(visibility.a),
            "axis_count_minimum":k, "difference_arc_width_upper":.001527,
            "adjacent_even_phase_spacing":float(spacing.mid), "count_coverage_failure":ALPHA,
            "axis_count_outage_upper":1e-110}


def necessary_equal_shots(incidences, epsilon=1e-6):
    """Necessary independent terminal launches per equal-size group.

    On a miss all available records have the same conditional law under
    both hypotheses. On a hit, even allowing a known phase and a quantum
    receiver, root fidelity is at least cos(delta). Condition on all hit
    flags and apply the pure-state Helstrom error lower bound. Convexity
    of (1-sqrt(1-x))/2 gives product(1-s*sin(delta)^2)^N <= 4e(1-e).
    Visibility zero can make discrimination impossible at every budget.
    """
    s = tuple(float(x) for x in incidences)
    if not s or not all(0 <= x <= 1 for x in s) or not 0 < epsilon < .5:
        raise ValueError("registered incidences and binary error below one half required")
    if not any(s):
        return None
    gap = math.sin(math.pi/4096)**2
    return math.ceil(math.log(4*epsilon*(1-epsilon))/sum(math.log1p(-p*gap) for p in s))


def inputs():
    """Authenticate the source grid and regenerate every original seven-port law."""
    blob = primary.ARCHIVE.read_bytes()
    if hashlib.sha256(blob).hexdigest() != SOURCE_DIGEST:
        raise ValueError("comparison-source archive checksum mismatch")
    words, ticks, _, _ = primary.source()
    with np.load(io.BytesIO(blob),allow_pickle=False) as archive:
        weights = archive["weights"]
        if (weights.shape != words.shape or weights.dtype != np.int64 or
                not np.array_equal(archive["ticks"], ticks) or np.any(weights < 0) or
                np.any(weights > GRID) or not np.all(weights.sum(-1) == GRID) or
                np.max(abs(words-weights/GRID)) > ERROR):
            raise ValueError("seven-port inputs differ from the regenerated source")
    registered = weights[0,...,:4].sum(axis=-1)/GRID
    # Upper source incidence gives a conservative necessary shot count.
    upper = np.minimum(1., registered.max(axis=-1) + 4*ERROR)
    original_cost = {"first_word_launches":8*primary.SHOTS,
                     "each_later_word_launches":8*primary.SHOTS,
                     "installed_panel_launches":72*primary.SHOTS,
                     "all_menu_launches":216*primary.SHOTS,
                     "registered_incidences":registered[:,:,0].tolist(),
                     "necessary_shots_equal_group":necessary_equal_shots(upper.ravel()),
                     "necessary_shots_best_group":necessary_equal_shots([upper.max()]),
                     "necessary_shots_unit_incidence":necessary_equal_shots([1.])}
    return weights[:,:,STATION].copy(), ticks[:,STATION].copy(), original_cost


def draw(weights):
    counts = np.empty_like(weights)
    for key in np.ndindex(weights.shape[:-1]):
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, *key])))
        counts[key] = rng.multinomial(SHOTS,weights[key]/GRID)
    return counts


def score(weights, counts, ticks):
    if (counts.shape != (3,2,9,7) or counts.dtype != np.int64 or (counts < 0).any()
            or (counts > SHOTS).any() or (counts.sum(-1) != SHOTS).any()
            or ((weights == 0) & (counts != 0)).any()):
        raise ValueError("invalid original-start counts")
    records, failures = [], []
    for m, menu in enumerate(primary.MENUS):
        for f in range(2):
            bank = counts[m,f]
            hits = []
            for j in range(1,9):
                candidates = reader.decode_elapsed(bank[:1],bank[j:j+1],rate=-math.pi/4096,
                                                   alpha=ALPHA,family_size=FAMILY,cell_error=ERROR,
                                                   interval="hoeffding")
                n = int(ticks[f,j]-ticks[f,0])
                compatible = candidates["even_runs_step_two"]
                unique = compatible == [[n,n]]
                hits.append(unique)
                records.append({"menu":menu,"fixture":f,"reference_elapsed":n,
                                "decoded_even_runs":compatible, "unique_correct":unique})
                if (menu == "installed" and not unique) or (menu == "dephased"
                        and candidates["eligible_even_count"] != 4096):
                    failures.append([menu,f,j])
            if menu == "passive" and all(hits):
                failures.append([menu,f,"all blind intervals mimic storage"])
    return {"counted_launches":int(counts.sum()),"early_active_launches":2*SHOTS,
            "each_later_active_launches":2*SHOTS,"records":records,"failures":failures}


def experiment(path=ARCHIVE, *, acquire=False, seed_replay=False):
    budget = design()
    weights, ticks, original = inputs()
    path = Path(path)
    if acquire:
        with path.open("xb") as file:
            np.savez_compressed(file,weights=weights,counts=draw(weights),optical_ticks=ticks)
    blob = path.read_bytes()
    digest = hashlib.sha256(blob).hexdigest()
    if not acquire and digest != ARCHIVE.with_suffix(".sha256").read_text(encoding="ascii").split()[0]:
        raise ValueError("reduced acquisition checksum mismatch")
    with np.load(io.BytesIO(blob),allow_pickle=False) as archive:
        if (sorted(archive.files) != ["counts","optical_ticks","weights"]
                or not np.array_equal(archive["weights"],weights)
                or not np.array_equal(archive["optical_ticks"],ticks)):
            raise ValueError("reduced acquisition has different source inputs")
        counts = archive["counts"]
    if seed_replay and not np.array_equal(draw(weights),counts):
        raise ValueError("different draws on this multinomial backend")
    return {"design":budget,"primary_resource_audit":original,"archive_sha256":digest,
            **score(weights,counts,ticks)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--records",default=ARCHIVE,type=Path)
    parser.add_argument("--acquire",action="store_true")
    parser.add_argument("--seed-replay",action="store_true")
    args = parser.parse_args()
    result = experiment(args.records,acquire=args.acquire,seed_replay=args.seed_replay)
    print(json.dumps(result,indent=2))
    raise SystemExit(int(bool(result["failures"])))
