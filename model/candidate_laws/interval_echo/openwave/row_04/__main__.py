"""Release the registered two-load excitation, with all external seams open."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply
from scipy.special import jv

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector, theorem_bound
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation


HERE = Path(__file__).resolve().parent
TIMES = ((0, 1, 4, 16, 32, 75, 150), (0, 1, 4, 16, 32))
ARMS = ("sharp", "perturbed", "exchange_only")
SHOTS = 16384


def precision(z):
    """Contour bound for the omitted Bessel coefficients, valid uniformly on [-1,1]."""
    degree = max(1, int(np.ceil(z)))
    while True:
        s = np.arccosh((degree + 1) / z)
        bound = 2 * np.exp(z * np.sinh(s) - (degree + 1) * s) / (-np.expm1(-s))
        if bound < 1e-12:
            return degree, bound
        degree += 1


def evolve(h, kappa, source, times):
    # K is positive, and its largest row sum in absolute value is at most 12.
    # The shifted potential is either -2*kappa or zero on this sector.
    lo, hi = -2 * kappa, 12
    mid, span = (hi + lo) / 2, (hi - lo) / 2
    op = (h - mid * sparse.eye(h.shape[0], format="csr")) / span
    coefficients = []
    for t in times:
        if t == 0:
            coefficients.append(np.ones(1))
        else:
            degree, _ = precision(span * t)
            terms = 2 * jv(np.arange(degree + 1), span * t)
            terms[0] /= 2
            coefficients.append(terms)
    outputs = [np.array(source * weights[0], dtype=complex) for weights in coefficients]
    previous = source
    current = op @ previous
    for k in range(1, max(map(len, coefficients))):
        direction = (-1j)**(k % 4)
        for result, values in zip(outputs, coefficients, strict=True):
            if k < len(values):
                # Each input basis vector is real, so each term updates just one quadrature.
                component = result.real if k % 2 == 0 else result.imag
                component += (direction.real + direction.imag) * values[k] * current
        previous, current = current, 2 * (op @ current) - previous
    return [v * np.exp(-1j * mid * t) for v, t in zip(outputs, times, strict=True)]


def perturbed_amplitudes():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    duration = 2 * np.arcsin(0.2)
    U = expm(-1j * duration * (exchange + candidate.KAPPA * contrast))
    source, dest = [law.SRC.index(x) for x in ((1, 0), (0, 1))]
    a, b = U[[source, dest], source]
    if abs(abs(a)**2 + abs(b)**2 - 1) > 1e-12 or abs(abs(b)**2 - .04) > 1e-12:
        raise ArithmeticError("partial contact branches differ from the enrolled preparation")
    return a, b


def verify_local_reader():
    """Undo the terminal echo after copying its pair flag on each binary input."""
    law = candidate.native_law
    echo = law.echo()
    pair_flag = np.arange(law.DIM) % len(law.SRC) == law.SRC.index((1, 1))
    for a, b in ((0, 0), (0, 1), (1, 0), (1, 1)):
        ket = np.eye(law.DIM)[:, law.SRC.index((a, b))]
        branch = echo @ ket
        event = a == b == 1
        np.testing.assert_allclose(echo.T @ (branch * pair_flag), ket * event, atol=1e-12)
        np.testing.assert_allclose(echo.T @ (branch * ~pair_flag), ket * (not event), atol=1e-12)


def radii(sector):
    car = sector.car
    a, b = np.asarray(car.external).T // 12
    graph = sparse.csr_matrix((np.ones(len(a)), (a, b)), shape=(car.cells, car.cells))
    distance = shortest_path(graph, directed=False, unweighted=True, indices=0).astype(int)
    return np.maximum(distance[sector.aa // 12], distance[sector.bb // 12])


def read_patch(sector, amplitudes, cone, local, origin, seed):
    p = abs(amplitudes)**2
    if abs(p.sum() - 1) > 2e-10:
        raise ArithmeticError("unitary release changed normalization")
    seam = p[local]
    target = float(seam.sum())
    table = np.r_[seam, max(0., 1 - target)]
    shell = np.bincount(cone, weights=p)
    cdf = np.cumsum(shell)
    pr_both = min(1., float(cdf[1]))
    rng = np.random.default_rng(seed)
    return dict(any_seam=float(p[sector.adjacent].sum()), original_seam=float(p[origin]),
                home_internal_seam=target, both_in_home=float(cdf[0]),
                both_within_one_cell_step=pr_both, radial_probabilities=shell.tolist(),
                adjacent_inside_four_cell_patch=float(p[sector.adjacent & (cone <= 1)].sum()),
                half_radius=int(np.searchsorted(cdf, .5)),
                ninety_percent_radius=int(np.searchsorted(cdf, .9)),
                home_seam_or_quiet_counts=rng.multinomial(SHOTS, table/table.sum()).tolist(),
                four_cell_both_count=int(rng.binomial(SHOTS, pr_both)))


def measure(level):
    car = federation(level)
    sector = PairSector(car)
    internal = np.asarray([e for e in car.seams if max(e) < 12], dtype=int)
    a, b = internal[0]
    c = next(n for n in car.neighbours[b] if n < 12 and n != a and n not in car.neighbours[a])
    origin, separate = sector.index(a, b), sector.index(a, c)
    local = np.array([sector.index(*e) for e in internal])
    cone = radii(sector)
    seed = np.zeros(sector.P)
    seed[origin] = 1
    h = sector.hamiltonian(candidate.KAPPA)
    sharp = evolve(h, candidate.KAPPA, seed, TIMES[level])
    source2 = np.zeros_like(seed)
    source2[separate] = 1
    other = evolve(h, candidate.KAPPA, source2, TIMES[level])
    x, y = perturbed_amplitudes()
    perturbed = [x*v + y*w for v, w in zip(sharp, other, strict=True)]
    del other
    control = evolve(sector.hamiltonian(0), 0, seed, TIMES[level])
    panels = {}
    for i, (key, states) in enumerate(zip(ARMS, (sharp, perturbed, control), strict=True)):
        panels[key] = [dict(time=t, **read_patch(sector, state, cone, local, origin, [471305, level, i, j]))
                       for j, (t, state) in enumerate(zip(TIMES[level], states, strict=True))]
    receipt = dict(level=level, cells=car.cells, seams=len(car.seams), external_links=len(car.external),
                   active_binary_pair_dimension=sector.P, source_cell=0, original_seam=[int(a), int(b)],
                   perturbation_contact=[int(b), int(c)], home_internal_addresses=internal.tolist(),
                   arms=panels)
    if level == 0:
        independent = expm_multiply(-1j * h, seed)
        err = float(np.linalg.norm(independent - sharp[1]))
        if err > 1e-10:
            raise ArithmeticError("independent free propagator disagrees")
        starts = [sector.index(i, j) for i, j in combinations(range(12), 2) if j not in car.neighbours[i]]
        # Real symmetric H gives the same transition amplitude with launch and receiver exchanged.
        receipt["herald_from_separated_starts"] = dict(starts=len(starts), time=1,
                 probability=float(np.mean(abs(independent[starts])**2)),
                 independent_expm_error=err)
    receipt["finite_binding_pass"] = (
        min(p["any_seam"] for p in panels["sharp"][1:]) > .9 and
        min(p["any_seam"] for p in panels["perturbed"][1:]) > .8 and
        min(p["any_seam"] for p in panels["exchange_only"][1:]) < .5)
    return receipt


def portable(value):
    if isinstance(value, dict):
        return {k: portable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [portable(v) for v in value]
    if isinstance(value, float):
        return round(value, 10)
    return value


def run():
    verify_local_reader()
    return portable(dict(candidate=candidate.NAME, kappa=candidate.KAPPA, shots_per_setting=SHOTS,
                         free_pair_lower_bound_all_times=1-theorem_bound(candidate.KAPPA),
                         boards=[measure(level) for level in range(len(TIMES))]))


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, help="save a new evidence receipt")
    output = cli.parse_args().output
    result = run()
    if output:
        with output.open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
    elif result != json.loads((HERE / "results.json").read_text()):
        raise SystemExit("Free-release receipt differs from the archive")
    for board in result["boards"]:
        print("cells", board["cells"], "finite free-binding pass", board["finite_binding_pass"])
        for arm in ARMS:
            rows = board["arms"][arm]
            print(arm, "minimum adjacent", min(row["any_seam"] for row in rows),
                  "final home patch", rows[-1]["home_internal_seam"])
