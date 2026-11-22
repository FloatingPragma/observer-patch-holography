"""Row 4 from a signed-record source: contact calibration, spectral sectors and every internal launch."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import exp, log, sqrt
from pathlib import Path

import numpy as np
from scipy import sparse

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector, doubled_hamiltonian
from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier, cell, federation
from . import adjacency_bound, record_action, record_release
from .bounded_clock_flight import DetectorEvolution

TIMES = (1, 5, 7, 17, 31, 37)
MENUS = ((1, 7, 37), (5, 17, 31))
ARMS = ("record", "exchange", "reversed")
BOARD_HASH = "99870398e3d516f899c600e9bf25ba53c7f513142ba8facbe1cedd94b23f915b"


def source_costs():
    board, words = cell(), record_action.signed_table()
    costs = {}
    for metric in ("square", "unary", "word"):
        tables = [record_action.local_flags(board, edge, words, metric) for edge in board.seams]
        profiles = {tuple(sum(table[d]) for d in range(-2, 3)) for table in tables}
        if len(profiles) != 1:
            raise ValueError("The primitive seams do not admit the same instrument")
        profile, = profiles
        if profile[2] != 0 or profile != profile[::-1] or profile[3] <= 0:
            raise ValueError("The source lost its zero or signed symmetry")
        costs[metric] = profile[2:]
    return costs


def cold_contact(costs, h):
    """A single occupied port, with total load one or two. Neither input is an adjacent pair."""
    answer = {}
    for load in (1, 2):
        idx = [law.SRC.index(word) for word in ((load,0),(0,load))]
        columns = h[:,idx]
        remainder = np.delete(columns, idx, axis=0)
        block = columns[idx,:]
        local_cost = costs["square"][load]
        expected = local_cost*np.eye(2) + np.array([[.5,-.5],[-.5,.5]])
        np.testing.assert_allclose(remainder, 0, atol=2e-13, rtol=0)
        np.testing.assert_allclose(block, expected, atol=2e-13, rtol=0)
        answer[load] = [local_cost, local_cost+1]
    return answer


def check_free_blocks(board, costs, pin_strength):
    """Entrywise comparison on both conserved load multisets of the registered carrier."""
    release, _, _, _ = record_release.source("primitive")
    sector = PairSector(board)
    distinct = sector.hamiltonian(pin_strength)
    if (distinct != release).nnz:
        raise ArithmeticError("The complete distinct-load free generators differ")
    _, unit, double = costs["square"]
    # A doubled port changes places with a vacancy on every incident seam.
    # Every missing incidence contributes the same cost against a vacant reference.
    edges = np.array(board.seams)
    a, b = edges.T
    hops = sparse.csr_matrix((-np.ones(2*len(edges)), (np.r_[a,b],np.r_[b,a])),
                             shape=(board.n,board.n))
    twice_exchange = hops + sparse.diags(board.degree, dtype=np.int64)
    cost = np.bincount(edges.ravel(), minlength=board.n)*double + (6-board.degree)*double
    rest = 2*6*unit
    derived_double = twice_exchange/2 + sparse.diags(cost-rest, dtype=np.int64)
    if (derived_double != doubled_hamiltonian(board,pin_strength)).nnz:
        raise ArithmeticError("The complete doubled-load free generators differ")
    return dict(distinct_columns=release.shape[0], doubled_columns=board.n, two_single_load_rest=rest)


def spectral(board, costs):
    blocks = adjacency_bound.certificate(board)
    degree = int(board.degree.max())
    if blocks["minimum_strength"] != costs["square"][1]:
        raise ValueError("The deployed width differs from the source cell")
    # (A-5I)(A+I)(A^2-5I)=0 on the primitive tile; the external matching adds
    # at most one to ||L/2||. Use sqrt(5) < 9/4 to keep all bounds rational.
    a = cell().adjacency().astype(np.int64).toarray()
    eye = np.eye(len(a), dtype=np.int64)
    if np.any((a-5*eye) @ (a+eye) @ (a@a-5*eye)):
        raise ArithmeticError("The primitive adjacency spectrum changed")
    one_upper = Fraction(29,8) + bool(board.external)
    two_upper = 2*one_upper
    # K>=0; on adjacent loads K <= limit/2. Min-max on P and its orthogonal complement
    # counts the entire negative band. For doubled loads the exchange is L/2 in [0,dmax].
    upper_twice = blocks["twice_adjacent_upper"] - 4*costs["square"][1]
    if upper_twice >= 0:
        raise ArithmeticError("Adjacent and separated bands need not be separated")
    doubled = {}
    for metric, (_, unit, double) in costs.items():
        # The pinned degree-six reference includes a vacant comparison at each missing incidence.
        shift = 6*(double-2*unit)
        doubled[metric] = [shift, float(shift+one_upper)]
    lo = doubled["square"][0]
    width = costs["square"][1]
    if degree > 6 or lo <= 0 or -lo+one_upper >= 0 or 2*width <= two_upper:
        raise ArithmeticError("The squared source does not separate doubled occupations by sign")
    cross_squared = Fraction(blocks["four_cross_norm_squared_upper"],4)
    attractive_gap = Fraction(-upper_twice,2)
    reversed_gap = 2*width-two_upper
    return dict(cells=board.cells, distinct_states=board.n*(board.n-1)//2, double_states=board.n,
                attractive_adjacent_band=[-2*width, upper_twice/2],
                reversed_upper_band=[2*width, float(2*width+two_upper)],
                reversed_upper_levels=len(board.seams), exchange_upper=float(two_upper),
                attractive_eigenstate_adjacency_lower=str(attractive_gap**2/(attractive_gap**2+cross_squared)),
                reversed_upper_eigenstate_adjacency_lower=str(reversed_gap**2/(reversed_gap**2+cross_squared)),
                negative_levels_attractive=len(board.seams), negative_levels_exchange=0,
                negative_levels_reversed=board.n, positive_sign_doubled_enclosures=doubled,
                reversed_doubled_enclosure=[-lo,float(-lo+one_upper)])


def home_rows(board, sector, unit):
    for sign in (1, 0, -1):
        h = sector.hamiltonian(sign*unit)
        flight = DetectorEvolution(h, TIMES)
        for home in range(board.cells):
            seams = [(a,b) for a,b in board.seams if a//12 == b//12 == home]
            if len(seams) != 30:
                raise ValueError("Incomplete home detector")
            detector = np.array([sector.index(a,b) for a,b in seams])
            for initial in detector:
                yield flight.read(initial, detector), flight.probability_error


def gates(probabilities, error):
    """Uniform over one chosen sharp launch, separately for the two appointment menus."""
    if probabilities.shape != (3,20,30,len(TIMES),31) or not np.isfinite(probabilities).all():
        raise ValueError("Incomplete census")
    if (not 0 <= error < 1e-5 or probabilities.min() < 0 or probabilities.max() > 1
            or np.max(abs(probabilities.sum(axis=-1)-1)) > 1e-12):
        raise ValueError("Invalid probability law or enclosure")
    home = 1-probabilities[..., -1]
    by_time = {}
    for it,t in enumerate(TIMES):
        by_time[t] = {}
        for a,name in enumerate(ARMS):
            val = home[a,:,:,it]
            low = np.unravel_index(val.argmin(), val.shape)
            high = np.unravel_index(val.argmax(), val.shape)
            by_time[t][name] = dict(minimum=float(val[low]), maximum=float(val[high]),
                                    min_cell_and_seam=list(map(int,low)), max_cell_and_seam=list(map(int,high)))
    shots = 4096
    margin = sqrt(log(800)/(2*shots))
    menus = []
    for menu in MENUS:
        time_slots = [TIMES.index(t) for t in menu]
        lower = home[0][:,:,time_slots].min(axis=2)-error
        upper = home[1,:,:,time_slots[-1]]+error
        fails = (lower <= .80+margin) | (upper >= .65-margin)
        d1 = max(0., float(lower.min())-.80-margin)
        d2 = max(0., .65-margin-float(upper.max()))
        # Bound for fresh independent releases, not an interpretation of a chosen PRNG seed.
        risk = min(1., 4*exp(-int(2*shots*min(d1,d2)**2)))
        menus.append(dict(times=menu, failing_launches=np.argwhere(fails).tolist(),
                          native_lower=float(lower.min()), control_upper=float(upper.max()),
                          fresh_sampling_failure_upper=risk))
    final = home[0,:,:,TIMES.index(37)]
    return dict(arms=ARMS, times=TIMES, terminal_probabilities=probabilities.size,
                complete_distributions=probabilities.size//31, l1_error_upper=error,
                below_94_percent_at_37=np.argwhere(final+error < .94).tolist(),
                unresolved_at_94_percent=np.argwhere(abs(final-.94) <= error).tolist(),
                home_extremes=by_time, menus=menus)


def compare_outcomes(actual, independent):
    if independent.shape != actual.shape or not np.isfinite(independent).all():
        raise ValueError("Incomplete independent detector laws")
    if independent.min() < 0 or independent.max() > 1 or np.max(abs(independent.sum(axis=-1)-1)) > 1e-12:
        raise ValueError("Invalid independent detector law")
    errors = abs(actual-independent).sum(axis=-1)
    worst = np.unravel_index(errors.argmax(), errors.shape)
    if errors[worst] > 1e-9:
        raise ValueError(f"Different terminal law at {worst}: {errors[worst]}")
    return dict(maximum_l1=float(errors[worst]), arm=ARMS[worst[0]], cell=int(worst[1]),
                seam_index=int(worst[2]), exposure=TIMES[worst[3]])


def census(costs):
    path = carrier.HERE / "carriers/FEDERATION_LEVEL0.npz"
    if sha256(path.read_bytes()).hexdigest() != BOARD_HASH:
        raise ValueError("Incorrect launch carrier")
    board = federation(0)
    sector = PairSector(board)
    p = np.empty((3,20,30,len(TIMES),31))
    errors = []
    for dest, (outcomes, error) in zip(p.reshape(-1,len(TIMES),31),
                                       home_rows(board, sector, costs["square"][1]), strict=True):
        dest[:] = outcomes
        errors.append(error)
    return p, gates(p, max(errors))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all-launches", action="store_true")
    parser.add_argument("--outcomes", type=Path, help="write every outcome to a new NumPy file in arm/cell/seam/time order")
    parser.add_argument("--compare", type=Path, help="check independent outcomes in the same order")
    args = parser.parse_args()
    if (args.outcomes or args.compare) and not args.all_launches:
        parser.error("--outcomes and --compare require --all-launches")
    costs = source_costs()
    h0, h1 = law.local_term_parts()
    result = dict(costs=costs, single_occupied_port_lines=cold_contact(costs, h0+candidate.KAPPA*h1),
                  binary_and_doubled=check_free_blocks(federation(0), costs, candidate.KAPPA), spectra=[])
    for level in range(3):
        board = federation(level)
        result["spectra"].append(spectral(board, costs))
    if args.all_launches:
        p, result["census"] = census(costs)
        if args.outcomes:
            with args.outcomes.open("xb") as out:
                np.save(out, p, allow_pickle=False)
        if args.compare:
            result["independent_comparison"] = compare_outcomes(p, np.load(args.compare, allow_pickle=False))
    print(json.dumps(result, indent=2, allow_nan=False))
    if args.all_launches and any(m["failing_launches"] for m in result["census"]["menus"]):
        raise SystemExit("An enrolled launch misses the uniform forecast")


if __name__ == "__main__":
    main()
