"""Read an arm-permutation singlet selected by ordinary native energy."""

import argparse
import json
from fractions import Fraction
from itertools import combinations, product
from math import factorial
from pathlib import Path

import numpy as np
from scipy import sparse

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.executable_model import carrier
from model.candidate_laws.interval_echo import currents
from .. import measure as original
from ...quark_winding.bound_face.measure import colour_words
from ..frames import EPS
from ..measure import detector
from ..open_release.flight_check import neighbourhood
from ..open_release.measure import dyadic_spectral_audit

HERE = Path(__file__).parent
TARGET = 100.8894585203
ARMS = 2**20
STARTS = 2**30
WAITS = (0, 0.25, 0.5, 1)
STEPS = 64


def apparatus():
    cell = carrier.cell()
    sec = BinarySector(cell, 3)
    gates = set(np.flatnonzero(carrier.federation(0).degree[:12] == 6))
    autos = carrier.automorphisms(cell.seams)
    for face in cell.faces():
        group = [p for p in autos if {p[i] for i in face} == set(face)]
        if all({p[i] for i in gates} == gates for p in group):
            break
    if len(group) != 6:
        raise ValueError("Incomplete group of arm permutations")
    orbits = sorted({tuple(sorted({p[i] for p in group})) for i in range(12)})
    arms = np.array(
        [
            [
                next(x for x in orbit if all(p[x] == x for p in group if p[f] == f))
                for orbit in orbits
            ]
            for f in face
        ]
    )
    sign = np.zeros((220, 220), dtype=int)
    for p in group:
        seq = [p[i] for i in face]
        parity = (-1) ** sum(a > b for a, b in combinations(seq, 2))
        idx = [sec.index[tuple(sorted(p[i] for i in w))] for w in sec.states]
        sign[idx, np.arange(220)] += parity
    gateway_count = np.array([len(gates.intersection(w)) for w in sec.states])
    return sec, gates, arms, sign, gateway_count


def bank(rng, probabilities):
    p = np.round(np.asarray(probabilities), 13)
    if p.min() < -1e-12 or abs(p.sum() - 1) > 1e-10:
        raise ValueError("Incomplete material records")
    p = np.maximum(0, p)
    return rng.multinomial(STARTS, p / p.sum()).tolist()


def filters(energy, eigenvectors, line, launch):
    x = (energy - line) / (64 * np.pi)
    filt = np.exp(-1j * np.pi * (ARMS - 1) * x) * np.sinc(ARMS * x) / np.sinc(x)
    coeff = eigenvectors[launch].astype(complex)
    masses = [1.0]
    for _ in range(4):
        coeff *= filt
        masses.append(float(np.vdot(coeff, coeff).real))
    return eigenvectors @ coeff / np.sqrt(masses[-1]), masses


def source():
    sec, gates, arms, sign, count = apparatus()
    h = sec.hamiltonian(8, False) + np.diag(count / 2)
    if np.any(sign @ h != h @ sign) or np.any(sign @ sign != 6 * sign):
        raise ValueError("Arm sign is not an exact native symmetry")
    en, v = np.linalg.eigh(h)
    target = abs(en - TARGET) < 1e-7
    if target.sum() != 1:
        raise ValueError("Energy does not resolve the declared line")
    # The previous verifier uses a common energy-origin cutoff of 100.5.
    # Shifting every level by two preserves all timings, residuals and gaps.
    cert, delta, proj_err, _ = dyadic_spectral_audit(
        h - 2 * np.eye(220), en - 2, v, target
    )
    pos = int(np.flatnonzero(target)[0])
    ideal = v[:, pos]
    sign_weight = float(ideal @ sign @ ideal / 6)
    if sign_weight - proj_err < 0.99:
        raise ValueError("Could not establish unit sign weight on the simple line")
    launch = sec.index[(1, 8, 9)]
    finite, masses = filters(en, v, TARGET, launch)
    lo_amp = np.sqrt(ideal[launch] ** 2 - proj_err)
    detuning = abs(en[pos] - TARGET) + delta
    gap = min(abs(en[~target] - TARGET)) - delta
    if max(abs(en - TARGET)) + delta >= 32 * np.pi:
        raise ValueError("Timing filter side lobe bound unavailable")
    filter_min = 1 - detuning * (ARMS - 1) / 64
    off = 1 / (ARMS * np.sin(gap / 64))
    lower_amp = lo_amp * filter_min**4
    err = float(2 * off**4 / lower_amp + 2 * proj_err / lo_amp + 2e-10)
    if err >= 1e-7:
        raise ValueError("Unresolved normalized preparation")
    # Exact-line output lies in sign. Bound the actual source using the
    # certified line's vector instead of accumulating total filter duration.
    summary = dict(
        energy=float(en[pos]),
        gap_lower=gap,
        spectral_shift_for_audit=-2,
        certificate=cert,
        sign_weight=sign_weight,
        arm_ports=arms.tolist(),
        native_gateways=sorted(map(int, gates)),
        source_word=[1, 8, 9],
        filter_acceptances=masses,
        original_start_source_lower=lower_amp**2,
        original_start_source_upper=float(
            ideal[launch] ** 2 + proj_err + off**8 + 1e-12
        ),
        finite_to_ideal_vector_upper=err,
        separate_guard_variants={},
    )
    for name, mat in [("guarded", h), ("shuttered", sec.hamiltonian(8, False))]:
        e, vec = np.linalg.eigh(mat)
        summary["separate_guard_variants"][name] = [
            dict(
                energy=float(e[mask].mean()),
                multiplicity=int(mask.sum()),
                sign_rank=int(round(np.sum(vec[:, mask] * (sign @ vec[:, mask])) / 6)),
            )
            for level in np.unique(e.round(8))
            for mask in [abs(e - level) < 6e-9]
        ]
    return sec, arms, en, v, pos, ideal, finite, summary


def colour(sec, arms, state):
    read = []
    for omit in range(4):
        groups = np.delete(arms, omit, axis=1)
        eps, admitted = detector(sec, groups)
        amp = state[[sec.index[tuple(sorted(t))] for t in product(*groups)]].reshape(
            (3, 3, 3)
        )
        admission = float(np.vdot(amp, amp).real)
        singlet = float(abs(np.vdot(eps, state)) ** 2)
        # Fundamental half-turn diag(-1,-1,1), applied to receiver one.
        # Dark output = probability in its -1 eigenspace.
        dark_one = float(np.vdot(amp[:2], amp[:2]).real)
        signs = np.array([-1, -1, 1])
        common = np.einsum("i,j,k->ijk", signs, signs, signs)
        dark_all = float(np.vdot(amp[common == -1], amp[common == -1]).real)
        read.append(
            dict(
                omitted=omit,
                triplet=admission,
                singlet=singlet,
                nonsinglet=float(np.linalg.norm(amp - np.vdot(EPS, amp) * EPS) ** 2),
                single_turn_dark=dark_one,
                common_turn_dark=dark_all,
            )
        )
    return read


def survey(level, gates):
    b = carrier.federation(level)
    edges = b.adjacency().astype(np.int64)
    v = np.asarray(b.degree == 6, dtype=np.int64)
    for step in range(1, 12):
        v = edges @ v
        vals = [int(v[i]) for i in gates]
        if len(set(vals)) > 1:
            return dict(length=step, gate_walk_counts=vals)
    raise ValueError("The finite walk survey did not resolve the gateway asymmetry")


def release(level, kappa, initial, sec, arms, vector_error, centre):
    sector, missing = neighbourhood(level)
    # All native degree services are retained in this compression of the full generator.
    h = sector.hamiltonian(kappa) + sparse.diags(missing / 2)
    h = (h - centre * sparse.eye(sector.dim)).tocsr()
    home = np.flatnonzero(np.all(sector.tri < 12, axis=1))
    if [tuple(w) for w in sector.tri[home]] != sec.states:
        raise ValueError("Changed home-address ordering")
    state = np.zeros(sector.dim, complex)
    state[home] = initial
    # Each exit target has a unique parent and each moving load at most one outward seam.
    boundary = np.sqrt(missing) / 2
    norm_b = np.sqrt(3) / 2
    dt = 1 / STEPS
    fourth = np.linalg.norm(h @ (h @ (h @ (h @ state))))
    integrated = 0.0
    steps = 16
    # The exact generator norm is bounded by 66 for the selected native and sham origins.
    step_error = np.exp(66 * dt) * (66 * dt) ** (steps + 1) / factorial(steps + 1)
    panels = {}
    for tick in range(STEPS + 1):
        t = tick * dt
        if t in WAITS:
            # Two applications of Duhamel. For an initially home terminal
            # vector the backward boundary flux is at most (3/4) u.
            numerical = tick * (step_error + 2e-12)
            home_err = (3 * t * t / 8) * (integrated + 2e-9) + numerical + vector_error
            rows = colour(sec, arms, state[home])
            for r in rows:
                q = max(0, np.sqrt(r["triplet"]) - home_err) ** 2
                bad = (np.sqrt(r["nonsinglet"]) + home_err) ** 2
                r["triplet_lower"] = q
                r["conditional_singlet_lower"] = max(0.0, 1 - bad / q) if q else 0.0
                r["joint_singlet_lower_given_source"] = (
                    max(0, np.sqrt(r["singlet"]) - home_err) ** 2
                )
            panels[str(t)] = dict(home_vector_error_upper=home_err, colours=rows)
        if tick == STEPS:
            break
        term = state
        next_state = state.copy()
        upper = float(np.linalg.norm(boundary * state)) * dt
        for j in range(1, steps + 1):
            term = (-1j * dt / j) * (h @ term)
            next_state += term
            if j <= 3:
                upper += dt * np.linalg.norm(boundary * term) / (j + 1)
        # Roundoff allowance for the centered sparse recurrence and its
        # first three derivatives; row sums and the term density bound it.
        upper += (
            norm_b * fourth * dt**5 / 120 + dt * (tick + 1) * (step_error + 3e-12) * 2
        )
        integrated += upper
        state = next_state
    return panels


def interval_bank(rng, row, err, source_bounds):
    """Inverse-CDF bank on all exact evolutions inside the propagation enclosure."""
    bounds = []
    for key in ("nonsinglet", "triplet"):
        norm = np.sqrt(row[key])
        bounds.append(
            (
                max(0, source_bounds[0] * max(0.0, norm - err) ** 2 - 1e-12),
                min(1, source_bounds[1] * (norm + err) ** 2 + 1e-12),
            )
        )
    bounds.append(tuple(source_bounds))
    cuts = sorted({0.0, 1.0, *(x for band in bounds for x in band)})
    counts = rng.multinomial(STARTS, np.diff(cuts))
    labels = ["other_colour", "singlet", "outside_triplets", "source_refused"]
    bins = dict.fromkeys(labels, 0)
    ambiguities = []
    for low, high, count in zip(cuts[:-1], cuts[1:], counts, strict=True):
        u = (low + high) / 2
        first = sum(u > right for _, right in bounds)
        last = sum(u > left for left, _ in bounds)
        if first == last:
            bins[labels[first]] += int(count)
        else:
            ambiguities.append(
                dict(possible=labels[first : last + 1], count=int(count))
            )
    return dict(
        cumulative_probability_bounds=bounds, certain=bins, unresolved=ambiguities
    )


def native_word(sec, arms, state, source_p, rng):
    groups = np.delete(arms, 3, axis=1)
    settings = original.meter_calibration(groups.tolist())
    _, _, words, _, _, _ = colour_words()
    basis = currents._service(carrier.cell(), False, 1.0)[4][0]
    tensor = state[[sec.index[tuple(sorted(t))] for t in product(*groups)]].reshape(
        3, 3, 3
    )
    p = float(np.vdot(tensor, tensor).real)
    banks = {}
    for name in ("sham", "centre", "half_turn"):
        u = basis.T @ words[name] @ basis
        for collective in (True, False):
            v = u if collective else np.eye(3)
            turned = np.einsum("ia,jb,kc,abc->ijk", u, v, v, tensor)
            # Native unitaries on the complete twelve-port receiver preserve
            # p even when a small component leaves its fundamental triplet.
            visibility = np.vdot(tensor, turned).real
            dark = max(0.0, float((p - visibility) / 2))
            bins = bank(
                rng,
                [
                    1 - source_p,
                    source_p * (1 - p),
                    source_p * dark,
                    source_p * (p - dark),
                ],
            )
            banks[name + ("_all" if collective else "_one")] = dict(
                dark_given_triplets=dark / p, counts=bins
            )
    return dict(receiver_compilation=settings, fringes=banks)


def run():
    sec, arms, en, vec, pos, ideal, finite, meta = source()
    rng = np.random.default_rng(13116019)
    native = colour(sec, arms, finite)
    source_p = meta["filter_acceptances"][-1]
    initial = []
    for r in native:
        d, q = r["singlet"], r["triplet"]
        initial.append(
            bank(
                rng,
                [1 - source_p, source_p * (1 - q), source_p * d, source_p * (q - d)],
            )
        )
    main = native[3]
    sharp = np.eye(220)[:, sec.index[(1, 8, 9)]]
    shutters, sv = np.linalg.eigh(sec.hamiltonian(8, False))
    close = shutters[np.argmin(abs(shutters - 100.390599))]
    shuttered, _ = filters(shutters, sv, close, sec.index[(1, 8, 9)])
    ground, ground_p = filters(en, vec, en[0], sec.index[(1, 8, 9)])
    # Dephasing by fully resolving the sharp port word before extraction.
    eps, adm = detector(sec, np.delete(arms, 3, axis=1))
    dephased = float((abs(finite) ** 2) @ (eps**2))
    controls = dict(
        sharp=colour(sec, arms, sharp),
        shuttered=colour(sec, arms, shuttered),
        symmetric_ground=colour(sec, arms, ground),
        symmetric_ground_source=ground_p[-1],
        dephased_admission=main["triplet"],
        dephased_singlet=dephased,
        common_half_turn_counts=bank(
            rng,
            [
                1 - source_p,
                source_p * (1 - main["triplet"]),
                source_p * main["common_turn_dark"],
                source_p * (main["triplet"] - main["common_turn_dark"]),
            ],
        ),
        single_half_turn_counts=bank(
            rng,
            [
                1 - source_p,
                source_p * (1 - main["triplet"]),
                source_p * main["single_turn_dark"],
                source_p * (main["triplet"] - main["single_turn_dark"]),
            ],
        ),
    )
    # Phase alignment does not affect any terminal read or bound.
    physical = native_word(sec, arms, finite, source_p, rng)
    source_bounds = (
        meta["original_start_source_lower"],
        meta["original_start_source_upper"],
    )
    boards = {}
    for level in range(3):
        boards[str(20 * 4**level)] = dict(
            walk_survey=survey(level, meta["native_gateways"]),
            native=release(
                level,
                8,
                ideal,
                sec,
                arms,
                meta["finite_to_ideal_vector_upper"],
                en[pos] - 144,
            ),
            exchange=release(
                level, 0, ideal, sec, arms, meta["finite_to_ideal_vector_upper"], 0
            ),
        )
        for channel in ("native", "exchange"):
            for panel in boards[str(20 * 4**level)][channel].values():
                panel["full_board_count_bank"] = interval_bank(
                    rng,
                    panel["colours"][3],
                    panel["home_vector_error_upper"],
                    source_bounds,
                )
    # Full compression QHQ >= -16. One external hop has a unique parent,
    # so a home eigenvector has variance <= 3/4 after opening. B=QHP has
    # norm <= sqrt(7) by exact complete-board incidence enumeration.
    cert = json.loads((HERE.parent / "open_release" / "results.json").read_bytes())
    gap = -16 - (
        float(Fraction(meta["certificate"]["selected_energy_upper"])) + 2 - 144
    )
    full_bound = {}
    for board in cert["boards"]:
        square = board["boundary"]["max_integer_Gram_row_sum"] / 4
        fraction = gap**2 / (gap**2 + (np.sqrt(0.75) + np.sqrt(square)) ** 2)
        full_bound[str(board["cells"])] = float(
            fraction - meta["finite_to_ideal_vector_upper"] - 1e-10
        )
    return dict(
        start_count=STARTS,
        source=meta,
        material_half_turn=physical,
        all_time_connected_lower_given_source=full_bound,
        filter_stage_counts=bank(
            rng, [-x for x in np.diff(meta["filter_acceptances"])] + [source_p]
        ),
        preparations=native,
        guarded_counts=initial,
        controls=controls,
        open_boards=boards,
    )


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = run()
    if args.write:
        (HERE / "results.json").write_text(
            json.dumps(data, indent=2, allow_nan=False) + "\n"
        )
    else:
        # Allow machine roundoff only for floating columns. Finite banks and
        # exact source certificates are deterministic on the archived toolchain.
        from ..open_release.measure import compare

        compare(
            json.loads(json.dumps(data)),
            json.loads((HERE / "results.json").read_text()),
        )
    print("Three-arm source and three released boards verified.")


if __name__ == "__main__":
    main()
