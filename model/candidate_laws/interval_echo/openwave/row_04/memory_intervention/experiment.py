"""Retrospective memory erasure, restored memory and permanently cut controls."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.executable_model.carrier import Carrier
from .. import open_pair, open_pair_release
from ..corridor_pair import neutral_column
from ..topology_census import action


HERE = Path(__file__).parent
ARCHIVE = HERE / "results.json"
STRENGTHS, TIMES = (0, .5, 8), (1, 7, 37)
TOLERANCE = 2e-10


def local_certificate():
    """Exact null action on every neutral memory word, with a 99-state cross-check."""
    import sympy as sp

    law = candidate.native_law
    wp = sp.Matrix(law._write_pattern().tolist())
    w = wp * sp.sqrt(2) / 2
    u = sp.BlockMatrix([[sp.zeros(law.N_SRC), w.T],
                       [w, sp.eye(law.N_LAB * law.N_SRC) - w * w.T]]).as_explicit()
    r, d2 = (sp.Matrix(x.tolist()) for x in law._label_ops())
    g0 = sp.diag(sp.zeros(law.N_SRC), sp.kronecker_product((sp.eye(law.N_LAB)-r)/2, sp.eye(law.N_SRC)))
    g1 = sp.diag(sp.zeros(law.N_SRC), sp.kronecker_product(d2, sp.eye(law.N_SRC)))
    ids = [law.N_SRC * (1 + i) + j for i, (_, d) in enumerate(law.LABELS) if d == 0
           for j in range(law.N_SRC)]
    columns = u[:, ids]
    if any(sp.expand(x) != 0 for x in (g0 * columns).col_join(g1 * columns)):
        raise ArithmeticError("Neutral memory projector fails the exact local null-action test")
    h0, h1 = law.local_term_parts()
    native = law.echo()
    error = max(float(np.max(abs(h[:, ids]))) for h in (h0, h1))
    if error > 5e-14 or not np.allclose(native, np.array(u).astype(float), rtol=0, atol=5e-14):
        raise ArithmeticError("Native implementation differs from the exact memory operator")
    # The intervention occurs while this contact has load word 00.
    empty = np.eye(law.DIM)[:, law.SRC.index((0, 0))]
    marked = neutral_column((0, 0))
    if (not np.allclose(native @ marked, empty, rtol=0, atol=5e-14)
            or not np.allclose(native @ native @ marked, marked, rtol=0, atol=5e-14)
            or not np.allclose(native @ empty, marked, rtol=0, atol=5e-14)):
        raise ArithmeticError("One echo fails to reset, or two echoes fail to restore, empty memory")
    return {"neutral_memory_dimension": len(ids), "exact_null_action": True,
            "native_null_action_max_error": error,
            "memory_address_projector_commutes_with_full_free_generator": True,
            "other_terms_argument": "All other seam terms and degree references act trivially on this seam's memory register.",
            "empty_contact_echo": {"native_one_echo_blank": 1, "native_two_echo_neutral": 1,
                                   "fixed_cut_one_echo_neutral": 1}}


def intervention_contact():
    cuts = open_pair.geometry()[0]
    vacant_tip = open_pair.PATH[2]
    contact = min(edge for edge in cuts if vacant_tip in edge)
    sharp = set(open_pair.BASIS[0])
    if sharp.intersection(contact) or contact != (152, 166):
        raise ArithmeticError("Intervention contact is occupied or differs from the declared selection")
    return contact


def generator(written):
    board = federation(0)
    live = Carrier("blank contacts after declared memory intervention", board.n,
                   [edge for edge in board.seams if edge not in written])
    sector = PairSector(live)
    kinetic = sector.exchange()
    costs = live.degree + 6 - board.degree
    contrast = costs[sector.aa] + costs[sector.bb] - 2 * sector.adjacent
    home = np.array([sector.index(*edge) for edge in open_pair.geometry()[1]])
    source = np.array([sector.index(*pair) for pair in open_pair.BASIS])
    return sector, kinetic, contrast, home, source


def flight(parts, density, exhaustion, kappa, exposure):
    sector, kinetic, contrast, home, source = parts
    active = np.flatnonzero(np.any(density != 0, axis=0) | np.any(density != 0, axis=1))
    launch = np.zeros((sector.P, len(active)))
    launch[source[active], np.arange(len(active))] = 1
    h = kinetic + sparse.diags(kappa * (contrast - 3), format="csr", dtype=np.float64)
    columns = expm_multiply((-1j * exposure) * h, launch)
    population = np.sum((columns @ density[np.ix_(active, active)]) * columns.conj(), axis=1).real
    quiet = np.ones(sector.P, dtype=bool)
    quiet[home] = False
    bins = np.r_[population[home], population[quiet].sum() + exhaustion]
    if not np.isfinite(bins).all() or bins.min() < 0 or abs(bins.sum()-1) > TOLERANCE:
        raise ArithmeticError("Intervention loses a source or terminal outcome")
    return bins


def short_exposure_certificate(density):
    """Exact leading quiet probability, retaining the entire accepted source."""
    retry = sum(Fraction(364, 729) ** attempt for attempt in range(8))
    vector = np.array([-13, 0, 14]) / 27
    if not np.allclose(density, float(retry) * np.outer(vector, vector), rtol=0, atol=5e-14):
        raise ArithmeticError("Source differs from the exact accepted three-query ensemble")
    board, home = federation(0), set(open_pair.geometry()[1])
    written = set(open_pair.geometry()[0])
    state = {open_pair.BASIS[0]: -13, open_pair.BASIS[2]: 14}
    coefficients = []
    for cut in (written, written - {intervention_contact()}):
        kinetic, contrast = action(board, cut, state)
        if any(pair not in home for pair in contrast):
            raise ArithmeticError("Contrast changes the leading quiet coefficient")
        numerator = sum(value * value for pair, value in kinetic.items() if pair not in home)
        coefficients.append(Fraction(numerator, 4 * 27**2) * retry)
    difference = coefficients[1] - coefficients[0]
    if difference != Fraction(49, 729) * retry or difference <= 0:
        raise ArithmeticError("Empty-contact reset lacks the exact short-exposure distinction")
    return {"quantity": "Coefficient of exposure squared in quiet probability above source exhaustion",
            "no_echo_or_fixed_cut": str(coefficients[0]), "native_one_echo": str(coefficients[1]),
            "difference": str(difference), "difference_decimal": float(difference),
            "all_real_couplings": True,
            "proof": "The initial accepted vector is wholly detected. Every contrast output remains detected. The squared norm of the complete undetected kinetic derivative gives the leading coefficient."}


def experiment(strengths=STRENGTHS, times=TIMES):
    local = local_certificate()
    contact = intervention_contact()
    written = set(open_pair.geometry()[0])
    baseline, reset = generator(written), generator(written - {contact})
    density, source = open_pair_release.source()
    exhaustion = source["exhaustion"]
    appointments = []
    short = short_exposure_certificate(density)
    for kappa in strengths:
        for exposure in times:
            fixed = flight(baseline, density, exhaustion, kappa, exposure)
            changed = flight(reset, density, exhaustion, kappa, exposure)
            variation = float(abs(changed - fixed).sum()/2)
            appointments.append({"kappa": kappa, "exposure": exposure,
                                 "bins": {"native_no_echo": fixed.tolist(),
                                          "native_one_echo": changed.tolist(),
                                          "native_two_echo": fixed.tolist(),
                                          "permanently_cut_one_echo": fixed.tolist()},
                                 "native_one_echo_home": float(changed[:-1].sum()),
                                 "permanently_cut_one_echo_home": float(fixed[:-1].sum()),
                                 "total_variation": variation})
    primary = next(row for row in appointments if row["kappa"] == 0 and row["exposure"] == 37)
    gate = primary["total_variation"] > .01
    return {"schema": "oph.interval_echo.memory_intervention.v1",
            "design_sha256": sha256((HERE/"DESIGN.md").read_bytes()).hexdigest(),
            "experiment_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "source_class": "Declared eleven-memory material on the frozen Interval Echo free law; supplied local echo and original source apparatus.",
            "scope": "Retrospective supplementary calculation; no public prospective registration, intact composite translation, particle identity or platform award.",
            "local_certificate": local, "contact": list(contact),
            "short_exposure_certificate": short,
            "source_exhaustion": exhaustion,
            "source_protected_weight": source["guaranteed_joint_home"],
            "outcome_order": [list(edge) for edge in open_pair.geometry()[1]] + ["quiet including source exhaustion"],
            "equivalent_arms_reason": "Two native echoes restore the initial neutral record; a permanent external cut keeps its contact term zero under any memory echo. The no-echo, two-echo and permanent-cut operators coincide on this preparation, so their flight is evaluated once.",
            "appointments": appointments,
            "primary_gate": {"kappa": 0, "exposure": 37, "total_variation_strictly_above": .01,
                             "value": primary["total_variation"], "passed": gate}}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--record", type=Path, help="exclusive first-record creation")
    mode.add_argument("--headline", action="store_true", help="recompute only the primary appointment and exact certificates")
    args = cli.parse_args()
    result = experiment((0,), (37,)) if args.headline else experiment()
    if args.record:
        with args.record.open("x", encoding="utf-8") as destination:
            json.dump(result, destination, indent=2, allow_nan=False)
            destination.write("\n")
    else:
        archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
        if args.headline:
            archive["appointments"] = [row for row in archive["appointments"]
                                       if row["kappa"] == 0 and row["exposure"] == 37]
        compare(result, archive)
    if not result["primary_gate"]["passed"]:
        raise ArithmeticError("The registered local distinguishing gate failed")
    print(json.dumps(result["primary_gate"]))


def compare(actual, archived):
    """Replay tolerates numerical propagation error and preserves exact metadata."""
    # Only the propagation bins and their derived metrics vary across solvers.
    lhs, rhs = dict(actual), dict(archived)
    arows, brows = lhs.pop("appointments"), rhs.pop("appointments")
    agate, bgate = lhs.pop("primary_gate"), rhs.pop("primary_gate")
    if lhs != rhs or len(arows) != len(brows):
        raise ArithmeticError("Memory-intervention source or metadata differs from its archive")
    for a, b in zip(arows, brows, strict=True):
        if (a["kappa"], a["exposure"]) != (b["kappa"], b["exposure"]) or a["bins"].keys() != b["bins"].keys():
            raise ArithmeticError("Changed intervention appointment or missing arm")
        for arm in a["bins"]:
            if not np.allclose(a["bins"][arm], b["bins"][arm], rtol=0, atol=TOLERANCE):
                raise ArithmeticError("Memory-intervention terminal bins differ from their archive")
        for key in ("native_one_echo_home", "permanently_cut_one_echo_home", "total_variation"):
            if not np.isfinite(b[key]) or abs(a[key]-b[key]) > TOLERANCE:
                raise ArithmeticError("Memory-intervention metric differs from its archive")
    if agate.keys() != bgate.keys() or any(agate[k] != bgate[k] for k in agate if k != "value") or not np.isfinite(bgate["value"]) or abs(agate["value"]-bgate["value"]) > TOLERANCE:
        raise ArithmeticError("Memory-intervention gate differs from its archive")


if __name__ == "__main__":
    main()
