"""Replicate the one-link assay in site, sharp coin and signed-memory coordinates."""

import argparse
import json

import numpy as np
from scipy.sparse import block_diag, csr_matrix, eye, kron
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from . import measure as m

WAITS = (1, 2, 4, 8, 12)
SHOTS = 2**30
TIMER_SHOTS = 2**26
SEED = 27112511
FAMILY = len(WAITS) * 2 * 2 * 3 * 6
TAIL = .01 / (2 * FAMILY)
EPS = np.sqrt(-np.log(TAIL) / (2 * SHOTS))


def bounds(c, n=SHOTS):
    return [float(beta.ppf(TAIL, c, n-c+1)) if c else 0.,
            float(beta.isf(TAIL, c+1, n-c)) if c < n else 1.]


def device(control):
    law = candidate.native_law
    u = law.echo()
    h0, h1 = law.local_term_parts()
    h = h0 + candidate.KAPPA * h1
    if control:
        # Remove just the dressed reversal generator at d=+/-2.
        select = np.repeat([False] + [abs(d) == 2 for _, d in law.LABELS], law.N_SRC)
        ex, _ = law.generator_parts()
        h -= u[:, select] @ ex[np.ix_(select, select)] @ u[:, select].T
    board = federation(0)
    # Infer all transport blocks from physical seam columns in the sharp basis.
    twin = [law.SRC.index(s) for s in m.READS[2:]]
    pos = np.diag(3 * candidate.KAPPA * (6 - board.degree)).astype(float)
    for edge in board.seams:
        if tuple(edge) != (m.A, m.B):
            pos[np.ix_(edge, edge)] += h[np.ix_(twin, twin)]
    tags = [law.LABELS.index((c, s)) + 1 for c in (0, 1) for s in (2, -2)]
    ends = (m.A, m.B)
    blocks, decoders = [], {}
    for site in range(board.n):
        local = tuple(1 + int(site == v) for v in ends)
        columns = [tag * law.N_SRC + law.SRC.index(local) for tag in tags]
        np.testing.assert_allclose(np.delete(h[:, columns], columns, axis=0), 0, atol=1e-13)
        blocks.append(csr_matrix(h[np.ix_(columns, columns)]))
        # The final detector sums over eleven private memory states. It reads
        # the material port pair (or the collective bin of other outcomes).
        decoders[local] = u[:, columns].reshape(1 + law.N_LAB, law.N_SRC, len(tags))
    generator = block_diag(blocks, format="csr") + kron(csr_matrix(pos), eye(4), format="csr")
    start = np.zeros((board.n, len(tags), 2), complex)
    input_columns = [law.SRC.index(w) for w in m.READS[:2]]
    written_columns = [tag * law.N_SRC + law.SRC.index((1, 1)) for tag in tags]
    start[1] = u[np.ix_(written_columns, input_columns)]
    np.testing.assert_allclose(start[1].T.conj() @ start[1], np.eye(2), atol=1e-13)
    return generator, start, decoders


def material_counts_law(waves, decoders):
    law = candidate.native_law
    sets = {(1, 1): np.setdiff1d(np.arange(waves.shape[0]), [m.A, m.B]),
            (2, 1): [m.A], (1, 2): [m.B]}
    pairs = np.zeros((law.N_SRC, 2))
    for local, sites in sets.items():
        decoded = np.einsum("alc,icb->ialb", decoders[local], waves[sites], optimize=True)
        pairs += (abs(decoded) ** 2).sum(axis=(0, 1))
    sel = [law.SRC.index(r) for r in m.READS]
    p = np.vstack((pairs[sel], np.delete(pairs, sel, axis=0).sum(axis=0)))
    np.testing.assert_allclose(p.sum(0), 1, atol=1e-11, rtol=0)
    return p / p.sum(0)


def run():
    rng = np.random.default_rng(SEED)
    banks, decisions = {}, {}
    timing = {}
    for name, control in (("native", False), ("changed", True)):
        h, start, decoders = device(control)
        mobile = []
        for wait in WAITS:
            waves = expm_multiply((-1j * wait) * h, start.reshape(-1, 2)).reshape(start.shape)
            mobile.append(material_counts_law(waves, decoders))
        for k, sign in enumerate(("plus", "minus")):
            for arm in ("blank", "park", "mobile"):
                key = f"{name}_{sign}_{arm}"
                rows = [p[:, k] if arm == "mobile" else np.eye(5)[k] for p in mobile]
                bank = [m.clock.enroll(rng, p, SHOTS) for p in rows]
                banks[key] = bank
                timer = rng.binomial(TIMER_SHOTS, np.sin(np.asarray(WAITS) / 16)**2)
                banks[key + "_timer"] = timer.tolist()
                timer_intervals = 16 * np.arcsin(np.sqrt([bounds(c, TIMER_SHOTS) for c in timer]))
                timing[key] = timer_intervals.tolist()
                reverse, blank = np.array(bank)[:, 1-k], np.array(banks[f"{name}_{sign}_blank"])[:, 1-k]
                decisions[key] = {
                    "opposite_order": reverse.tolist(),
                    "opposite_intervals": [bounds(c) for c in reverse],
                    "resolved_against_blank": [bounds(c)[0] > bounds(b)[1] for c, b in zip(reverse, blank)],
                    "resolved_hoeffding": ((reverse - blank)/SHOTS > 2 * EPS).tolist(),
                }
    return {"counts": banks, "decisions": decisions, "dimer_duration_intervals": timing,
            "waits": list(WAITS), "shots": SHOTS, "timer_shots": TIMER_SHOTS, "seed": SEED,
            "count_columns": [*map(list, m.READS), "other"], "comparison_marginals": FAMILY,
            "simultaneous_99_probability_radius": float(EPS)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = m.HERE / "replication.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
        compare(result, json.loads(path.read_text()), "sharp-order replication")
    for sign in ("plus", "minus"):
        entry = result["decisions"][f"native_{sign}_mobile"]
        print(sign, entry["opposite_order"], entry["resolved_against_blank"])


if __name__ == "__main__":
    main()
