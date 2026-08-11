"""Conditional-repair thermodynamics certificate.

The finite thermodynamic completion of the observer framework rests on
one repaired distinction and two optimizer instantiations. The
deterministic strict-descent normalizer settles public facts and carries
no entropy inequality; the certificate exhibits the explicit
counterexample. The state-side information projection of the third axiom
selects the exponential family; the transition-side information
projection onto the fibre of the complete repaired visible datum selects
weighted conditional resampling. The resulting kernel is stochastic,
idempotent, reversible, stationary, and fixes every fibre-measurable
observable, and relative entropy to the reference law contracts under
it. The Lean modules `Thermodynamics/FiniteConditionalRepair.lean`,
`Thermodynamics/FirstLawIdentity.lean`, and
`EventAlgebra/PartitionPinchingCP.lean` carry the exact statements; this
certificate replays the exact finite kernel algebra over the rationals,
checks the entropy statements at high precision on randomized instances,
and records the open physical receipts.

Run with --write to refresh the committed receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from fractions import Fraction
from pathlib import Path
from typing import Any

SCHEMA = "oph.conditional_repair_thermodynamics.v1"
STATUS = (
    "FINITE_FOUR_LAW_PACKAGE_CERTIFIED__"
    "TRANSITION_MAXENT_IS_CONDITIONAL_RESAMPLING__"
    "SECOND_LAW_CONTRACTION_AND_GIBBS_PROJECTION_REPLAYED__"
    "FLUCTUATION_AND_CAP_FIRST_LAW_IDENTITIES_ATTAINED__"
    "STRICT_DESCENT_ENTROPY_COUNTEREXAMPLE_RECORDED__"
    "COMMON_REFERENCE_REALIZATION_AND_ENERGY_CLOCK_OPEN"
)
RECEIPT_PATH = Path(__file__).resolve().parent / "runtime" / (
    "conditional_repair_receipt.json"
)


class ThermoError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ThermoError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode()


def tagged_sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------------------
# Strict-descent control: a terminating idempotent normalizer that
# decreases Shannon entropy
# ---------------------------------------------------------------------------


def strict_descent_control() -> dict[str, Any]:
    rep = {"a": "b", "b": "b"}
    require(
        all(rep[rep[x]] == rep[x] for x in rep),
        "normalizer is not idempotent",
    )
    p = {"a": Fraction(1, 2), "b": Fraction(1, 2)}
    push: dict[str, Fraction] = {"a": Fraction(0), "b": Fraction(0)}
    for x, mass in p.items():
        push[rep[x]] += mass
    require(
        push == {"a": Fraction(0), "b": Fraction(1)},
        "pushforward drift",
    )
    # The source law has two positive atoms, so its Shannon entropy is
    # strictly positive; the pushforward is deterministic, so its
    # entropy is exactly zero. The comparison is structural and needs no
    # numerical logarithm.
    support_before = sum(1 for v in p.values() if v > 0)
    support_after = sum(1 for v in push.values() if v > 0)
    require(support_before == 2 and support_after == 1, "support drift")
    return {
        "statement": (
            "the terminating idempotent normalizer a -> b, b -> b maps "
            "the uniform two-point law to a deterministic law, so its "
            "entropy drops from log 2 to zero; strict-descent settlement "
            "alone carries no second law, and the repair entropy "
            "inequality lives on the separate conditional-resampling "
            "kernel"
        ),
        "normalizer": {"a": "b", "b": "b"},
        "entropy_before": "log 2 (two positive atoms)",
        "entropy_after": "0 (deterministic)",
    }


# ---------------------------------------------------------------------------
# Exact conditional-resampling kernel algebra over the rationals
# ---------------------------------------------------------------------------


def random_instance(rng: random.Random, size: int, blocks: int):
    labels = [rng.randrange(blocks) for _ in range(size)]
    # every block nonempty is not required; empty blocks never appear as
    # fibres of a point
    pi = [
        Fraction(rng.randint(1, 30), rng.randint(1, 7))
        for _ in range(size)
    ]
    return labels, pi


def heat_bath(labels, pi):
    size = len(labels)
    mass = {}
    for lab, weight in zip(labels, pi):
        mass[lab] = mass.get(lab, Fraction(0)) + weight
    kernel = [
        [
            (pi[y] / mass[labels[x]]) if labels[y] == labels[x]
            else Fraction(0)
            for y in range(size)
        ]
        for x in range(size)
    ]
    return kernel, mass


def kernel_algebra_certificate() -> dict[str, Any]:
    rng = random.Random(20260803)
    checked = 0
    for _ in range(25):
        size = rng.randint(4, 9)
        blocks = rng.randint(1, 4)
        labels, pi = random_instance(rng, size, blocks)
        kernel, _mass = heat_bath(labels, pi)
        # stochastic rows
        for row in kernel:
            require(sum(row) == 1, "row sum is not one")
        # idempotence
        for x in range(size):
            for y in range(size):
                entry = sum(
                    kernel[x][z] * kernel[z][y] for z in range(size)
                )
                require(entry == kernel[x][y], "kernel is not idempotent")
        # detailed balance and stationarity
        for x in range(size):
            for y in range(size):
                require(
                    pi[x] * kernel[x][y] == pi[y] * kernel[y][x],
                    "detailed balance fails",
                )
        for y in range(size):
            require(
                sum(pi[x] * kernel[x][y] for x in range(size)) == pi[y],
                "reference law is not stationary",
            )
        # fibre-measurable observables are fixed
        values = {lab: Fraction(rng.randint(-9, 9)) for lab in set(labels)}
        observable = [values[lab] for lab in labels]
        for x in range(size):
            image = sum(
                kernel[x][y] * observable[y] for y in range(size)
            )
            require(image == observable[x], "protected observable moves")
        # random-scan mixture of two fibre maps stays stochastic and
        # stationary
        labels2 = [rng.randrange(blocks + 1) for _ in range(size)]
        kernel2, _ = heat_bath(labels2, pi)
        a = Fraction(rng.randint(1, 5), 7)
        mix = [
            [
                a * kernel[x][y] + (1 - a) * kernel2[x][y]
                for y in range(size)
            ]
            for x in range(size)
        ]
        for row in mix:
            require(sum(row) == 1, "mixture row sum drift")
        for y in range(size):
            require(
                sum(pi[x] * mix[x][y] for x in range(size)) == pi[y],
                "mixture stationarity drift",
            )
        checked += 1
    return {
        "statement": (
            "on every random exact instance the fibre kernel is "
            "stochastic, idempotent, reversible, stationary, fixes "
            "fibre-measurable observables, and positive mixtures of "
            "fibre kernels stay stochastic and stationary"
        ),
        "instances": checked,
        "arithmetic": "exact rational",
    }


# ---------------------------------------------------------------------------
# High-precision entropy statements
# ---------------------------------------------------------------------------


def entropy_certificate() -> dict[str, Any]:
    import mpmath

    mp = mpmath.mp
    mp.dps = 60
    rng = random.Random(20260804)

    def kl(q, p):
        total = mp.mpf(0)
        for qi, pi_ in zip(q, p):
            if qi > 0:
                total += qi * mp.log(qi / pi_)
        return total

    contraction_margin = None
    optimizer_margin = None
    pythagorean_residual = mp.mpf(0)
    for _ in range(20):
        size = rng.randint(4, 8)
        blocks = rng.randint(1, 3)
        labels, pi_frac = random_instance(rng, size, blocks)
        total = sum(pi_frac)
        pi = [mp.mpf(w.numerator) / w.denominator / (
            mp.mpf(total.numerator) / total.denominator) for w in pi_frac]
        kernel, _ = heat_bath(labels, pi_frac)
        kernel_mp = [
            [mp.mpf(e.numerator) / e.denominator for e in row]
            for row in kernel
        ]
        raw = [mp.mpf(rng.randint(1, 40)) for _ in range(size)]
        norm = sum(raw)
        p = [r / norm for r in raw]
        pushed = [
            sum(p[x] * kernel_mp[x][y] for x in range(size))
            for y in range(size)
        ]
        before = kl(p, pi)
        after = kl(pushed, pi)
        margin = before - after
        require(
            margin >= -mp.mpf("1e-45"),
            "second-law contraction fails numerically",
        )
        if contraction_margin is None or margin < contraction_margin:
            contraction_margin = margin
        # transition-side optimality: the kernel row beats every random
        # fibre-supported competitor
        x0 = rng.randrange(size)
        fibre = [y for y in range(size) if labels[y] == labels[x0]]
        row = kernel_mp[x0]
        row_score = kl(row, pi)
        for _ in range(6):
            weights = [mp.mpf(rng.randint(1, 20)) for _ in fibre]
            wtot = sum(weights)
            competitor = [mp.mpf(0)] * size
            for pos, y in enumerate(fibre):
                competitor[y] = weights[pos] / wtot
            gap = kl(competitor, pi) - row_score
            require(
                gap >= -mp.mpf("1e-45"),
                "conditional row is not optimal",
            )
            if optimizer_margin is None or gap < optimizer_margin:
                optimizer_margin = gap
        # Gibbs Pythagorean identity on a random exponential family
        sigma_raw = [mp.mpf(rng.randint(1, 30)) for _ in range(size)]
        sigma = [s / sum(sigma_raw) for s in sigma_raw]
        charge = [mp.mpf(rng.randint(-4, 4)) for _ in range(size)]
        lam = mp.mpf(rng.randint(-2, 2)) / 3
        tau_raw = [s * mp.e ** (-lam * c) for s, c in zip(sigma, charge)]
        tau = [t / sum(tau_raw) for t in tau_raw]
        moment = sum(t * c for t, c in zip(tau, charge))
        # draw rho on the same moment surface: mix tau with a zero-mean
        # perturbation along a direction orthogonal to the charge
        rho = list(tau)
        if size >= 3:
            i, j, kdx = rng.sample(range(size), 3)
            ci, cj, ck = charge[i], charge[j], charge[kdx]
            # solve e_i + a e_j + b e_k with 1 + a + b = 0 and
            # ci + a cj + b ck = 0
            det = cj - ck
            if abs(det) > mp.mpf("1e-9"):
                a = (ck - ci) / det
                b = -1 - a
                eps = mp.mpf("1e-3")
                cap = min(tau[i], tau[j], tau[kdx]) / (
                    1 + abs(a) + abs(b)) / 2
                eps = min(eps, cap)
                rho[i] += eps
                rho[j] += eps * a
                rho[kdx] += eps * b
        moment_rho = sum(r * c for r, c in zip(rho, charge))
        require(
            abs(moment_rho - moment) < mp.mpf("1e-40"),
            "constraint surface drift",
        )
        residual = abs(kl(rho, sigma) - kl(rho, tau) - kl(tau, sigma))
        require(
            residual < mp.mpf("1e-40"),
            "Gibbs Pythagorean identity fails",
        )
        if residual > pythagorean_residual:
            pythagorean_residual = residual
    return {
        "statement": (
            "relative entropy to the reference contracts under the "
            "fibre kernel, the kernel row minimizes relative entropy "
            "among fibre-supported competitors, and the Gibbs "
            "information-projection identity holds on the constrained "
            "surface, at sixty-digit precision on randomized instances"
        ),
        "instances": 20,
        "worst_contraction_margin": mpf_str(contraction_margin),
        "worst_optimizer_gap": mpf_str(optimizer_margin),
        "worst_pythagorean_residual": mpf_str(pythagorean_residual),
    }


def mpf_str(x) -> str:
    import mpmath

    return mpmath.nstr(x, 8, strip_zeros=False)


def low_temperature_certificate() -> dict[str, Any]:
    import mpmath

    mp = mpmath.mp
    mp.dps = 50
    rng = random.Random(20260805)
    rows = []
    for _ in range(8):
        d = rng.randint(3, 9)
        g0 = rng.randint(1, d - 1)
        gap = Fraction(rng.randint(1, 5), rng.randint(1, 3))
        spectrum = [Fraction(0)] * g0 + [
            gap + Fraction(rng.randint(0, 6), 2) for _ in range(d - g0)
        ]
        worst = None
        for beta_int in (1, 2, 5, 11, 23):
            beta = mp.mpf(beta_int)
            weights = [
                mp.e ** (-beta * (mp.mpf(s.numerator) / s.denominator))
                for s in spectrum
            ]
            z = sum(weights)
            excited = sum(weights[g0:]) / z
            bound = mp.mpf(d - g0) / g0 * mp.e ** (
                -beta * (mp.mpf(gap.numerator) / gap.denominator)
            )
            require(
                excited <= bound + mp.mpf("1e-45"),
                "low-temperature bound fails",
            )
            # trace distance to the ground-sector uniform state equals
            # the excited mass
            tv = (
                sum(
                    abs(weights[i] / z - mp.mpf(1) / g0)
                    for i in range(g0)
                )
                + sum(abs(w / z) for w in weights[g0:])
            ) / 2
            require(
                abs(tv - excited) < mp.mpf("1e-44"),
                "trace-distance identity fails",
            )
            slack = bound - excited
            if worst is None or slack < worst:
                worst = slack
        rows.append(
            {
                "dimension": d,
                "ground_degeneracy": g0,
                "gap": str(gap),
                "worst_bound_slack": mpf_str(worst),
            }
        )
    return {
        "statement": (
            "the excited Gibbs mass obeys the finite gap bound, and the "
            "trace distance from the Gibbs state to the uniform ground "
            "state equals the excited mass, on randomized exact spectra"
        ),
        "spectra": rows,
    }


def pinching_tower_certificate() -> dict[str, Any]:
    import mpmath

    mp = mpmath.mp
    mp.dps = 40
    rng = random.Random(20260806)

    def entropy(eigs):
        total = mp.mpf(0)
        for e in eigs:
            if e > mp.mpf("1e-30"):
                total -= e * mp.log(e)
        return total

    checked = 0
    for _ in range(6):
        n = rng.randint(3, 5)
        # random density matrix
        raw = mpmath.matrix(
            [[mp.mpf(rng.randint(-5, 5)) + 1j * rng.randint(-5, 5)
              for _ in range(n)] for _ in range(n)]
        )
        rho = raw * raw.H
        rho = rho / mpmath.mp.mpf(sum(rho[i, i].real for i in range(n)))
        # coordinate partition into two blocks
        cut = rng.randint(1, n - 1)
        blocks = [list(range(cut)), list(range(cut, n))]
        pinched = mpmath.matrix(n, n)
        for block in blocks:
            for i in block:
                for j in block:
                    pinched[i, j] = rho[i, j]
        averaged = mpmath.matrix(n, n)
        for block in blocks:
            weight = sum(rho[i, i].real for i in block)
            for i in block:
                averaged[i, i] = weight / len(block)
        # entropies through eigenvalues
        def eigs(mat):
            return [
                complex(v).real
                for v in mpmath.mp.eig(mat, left=False, right=False)
            ]
        s_rho = entropy(sorted(eigs(rho)))
        s_pinched = entropy(sorted(eigs(pinched)))
        s_avg = entropy(sorted(eigs(averaged)))
        require(
            s_rho <= s_pinched + mp.mpf("1e-25")
            and s_pinched <= s_avg + mp.mpf("1e-25"),
            "entropy tower fails",
        )
        checked += 1
    return {
        "statement": (
            "the two-stage record tower holds on randomized density "
            "matrices: pinching removes inter-record coherence without "
            "lowering entropy, and public-record averaging raises it "
            "further; the Kraus completeness of pinching is the Lean "
            "theorem kraus_complete"
        ),
        "instances": checked,
    }


def zeroth_law_certificate() -> dict[str, Any]:
    import mpmath

    mp = mpmath.mp
    mp.dps = 50
    rng = random.Random(20260807)

    def gibbs(spectrum, beta):
        weights = [mp.e ** (-beta * s) for s in spectrum]
        z = sum(weights)
        return [w / z for w in weights]

    def entropy(p):
        return -sum(x * mp.log(x) for x in p if x > 0)

    # Thermometer identification: on a nondegenerate spectrum the Gibbs
    # ratio between two levels pins beta; distinct betas give distinct
    # states.
    thermometer_checked = 0
    for _ in range(10):
        spectrum = sorted(
            {Fraction(rng.randint(0, 9), rng.randint(1, 3))
             for _ in range(4)}
        )
        if len(spectrum) < 2:
            continue
        spectrum = [mp.mpf(s.numerator) / s.denominator for s in spectrum]
        b1 = mp.mpf(rng.randint(1, 9)) / 4
        b2 = b1 + mp.mpf(rng.randint(1, 5)) / 4
        g1, g2 = gibbs(spectrum, b1), gibbs(spectrum, b2)
        gap = max(abs(a - c) for a, c in zip(g1, g2))
        require(
            gap > mp.mpf("1e-6"),
            "distinct betas give the same thermometer state",
        )
        recovered = -mp.log(g1[1] / g1[0]) / (spectrum[1] - spectrum[0])
        require(
            abs(recovered - b1) < mp.mpf("1e-40"),
            "two-level ratio does not recover beta",
        )
        thermometer_checked += 1

    # Additive contact: at fixed total energy, the entropy of the product
    # of two Gibbs branches is maximal at the equal-beta split.
    contact_checked = 0
    for _ in range(6):
        spec_a = [mp.mpf(rng.randint(0, 6)) / 2 for _ in range(4)]
        spec_b = [mp.mpf(rng.randint(0, 6)) / 2 for _ in range(4)]
        if len(set(spec_a)) < 2 or len(set(spec_b)) < 2:
            continue
        spec_a = [s - min(spec_a) for s in spec_a]
        spec_b = [s - min(spec_b) for s in spec_b]
        beta = mp.mpf(rng.randint(2, 8)) / 4
        pa, pb = gibbs(spec_a, beta), gibbs(spec_b, beta)
        u_total = (
            sum(p * s for p, s in zip(pa, spec_a))
            + sum(p * s for p, s in zip(pb, spec_b))
        )
        s_equal = entropy(pa) + entropy(pb)

        def split_entropy(beta_a, beta_b, target):
            qa, qb = gibbs(spec_a, beta_a), gibbs(spec_b, beta_b)
            u = (
                sum(p * s for p, s in zip(qa, spec_a))
                + sum(p * s for p, s in zip(qb, spec_b))
            )
            return entropy(qa) + entropy(qb), u

        # perturbed unequal-beta splits matched to the same total energy
        # by a one-dimensional search over the second branch
        for delta in (mp.mpf("0.3"), mp.mpf("-0.2")):
            beta_a = beta + delta
            lo, hi = mp.mpf("0.01"), mp.mpf(60)
            _, u_lo = split_entropy(beta_a, lo, u_total)
            _, u_hi = split_entropy(beta_a, hi, u_total)
            if not (u_hi <= u_total <= u_lo):
                continue
            for _ in range(220):
                mid = (lo + hi) / 2
                _, u = split_entropy(beta_a, mid, u_total)
                if u > u_total:
                    lo = mid
                else:
                    hi = mid
            s_uneq, u_uneq = split_entropy(beta_a, (lo + hi) / 2, u_total)
            require(
                abs(u_uneq - u_total) < mp.mpf("1e-20"),
                "contact energy match fails",
            )
            require(
                s_uneq <= s_equal + mp.mpf("1e-25"),
                "unequal-beta split beats the equal-beta split",
            )
            contact_checked += 1
    return {
        "statement": (
            "a nondegenerate finite thermometer identifies beta through "
            "its level ratios, and at fixed total energy the equal-beta "
            "split maximizes the additive-contact entropy; the Lean "
            "theorem gibbs_beta_injective carries the exact transitivity "
            "core"
        ),
        "thermometer_instances": thermometer_checked,
        "contact_instances": contact_checked,
    }


def clausius_certificate() -> dict[str, Any]:
    import mpmath

    mp = mpmath.mp
    mp.dps = 50
    rng = random.Random(20260808)

    def entropy(p):
        return -sum(x * mp.log(x) for x in p if x > 0)

    worst = None
    for _ in range(15):
        size = rng.randint(4, 8)
        energies = [mp.mpf(rng.randint(0, 8)) / 2 for _ in range(size)]
        beta = mp.mpf(rng.randint(1, 8)) / 4
        weights = [mp.e ** (-beta * e) for e in energies]
        z = sum(weights)
        tau = [w / z for w in weights]
        # a fibre structure unrelated to the energy, so heat flows
        labels = [rng.randrange(2) for _ in range(size)]
        tau_frac_like = list(tau)
        mass = {}
        for lab, w in zip(labels, tau_frac_like):
            mass[lab] = mass.get(lab, mp.mpf(0)) + w
        kernel = [
            [
                (tau[y] / mass[labels[x]]) if labels[y] == labels[x]
                else mp.mpf(0)
                for y in range(size)
            ]
            for x in range(size)
        ]
        raw = [mp.mpf(rng.randint(1, 30)) for _ in range(size)]
        p = [r / sum(raw) for r in raw]
        pushed = [
            sum(p[x] * kernel[x][y] for x in range(size))
            for y in range(size)
        ]
        delta_s = entropy(pushed) - entropy(p)
        heat = beta * (
            sum(q * e for q, e in zip(pushed, energies))
            - sum(q * e for q, e in zip(p, energies))
        )
        margin = delta_s - heat
        require(
            margin >= -mp.mpf("1e-40"),
            "Clausius inequality fails on a tau-preserving repair",
        )
        if worst is None or margin < worst:
            worst = margin
    # Landauer instance: a two-state record at inverse temperature beta
    # relaxed by the full resampling kernel from the uniform record; the
    # entropy drop times 1/beta lower-bounds the expelled energy.
    landauer_rows = []
    for beta_int in (1, 2, 5):
        beta = mp.mpf(beta_int) / 2
        energies = [mp.mpf(0), mp.mpf(1)]
        weights = [mp.e ** (-beta * e) for e in energies]
        z = sum(weights)
        tau = [w / z for w in weights]
        kernel = [[tau[0], tau[1]], [tau[0], tau[1]]]
        p = [mp.mpf(1) / 2, mp.mpf(1) / 2]
        pushed = [
            sum(p[x] * kernel[x][y] for x in range(2)) for y in range(2)
        ]
        s_drop = entropy(p) - entropy(pushed)
        expelled = (
            sum(a * e for a, e in zip(p, energies))
            - sum(a * e for a, e in zip(pushed, energies))
        )
        require(
            expelled >= s_drop / beta - mp.mpf("1e-40"),
            "Landauer bound fails",
        )
        landauer_rows.append(
            {
                "beta": mpf_str(beta),
                "entropy_drop": mpf_str(s_drop),
                "expelled_energy": mpf_str(expelled),
                "landauer_floor": mpf_str(s_drop / beta),
            }
        )
    return {
        "statement": (
            "on tau-preserving repair kernels with fibres transverse to "
            "the energy, the entropy change dominates beta times the "
            "heat, matching the Lean clausius theorem; the Landauer "
            "instances confirm that an entropy drop of c expels at "
            "least c/beta of energy"
        ),
        "instances": 15,
        "worst_margin": mpf_str(worst),
        "landauer": landauer_rows,
    }


def first_law_and_support_certificate() -> dict[str, Any]:
    # Exact first-law split on a rational matrix instance.
    n = 3
    rng = random.Random(20260809)

    def rmat():
        return [
            [Fraction(rng.randint(-4, 4), rng.randint(1, 3))
             for _ in range(n)]
            for _ in range(n)
        ]

    def mul(a, bm):
        return [
            [sum(a[i][k] * bm[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)
        ]

    def tr(a):
        return sum(a[i][i] for i in range(n))

    for _ in range(10):
        rho, drho, ham, dham = rmat(), rmat(), rmat(), rmat()
        left = tr(mul(
            [[rho[i][j] + drho[i][j] for j in range(n)] for i in range(n)],
            [[ham[i][j] + dham[i][j] for j in range(n)] for i in range(n)],
        )) - tr(mul(rho, ham))
        right = tr(mul(drho, ham)) + tr(mul(rho, dham)) + tr(
            mul(drho, dham)
        )
        require(left == right, "first-law split fails exactly")

    # Full support is preserved by the repair kernel, exactly.
    for _ in range(10):
        size = rng.randint(3, 7)
        labels = [rng.randrange(3) for _ in range(size)]
        pi = [Fraction(rng.randint(1, 9)) for _ in range(size)]
        kernel, _ = heat_bath(labels, pi)
        raw = [Fraction(rng.randint(1, 9)) for _ in range(size)]
        p = [r / sum(raw) for r in raw]
        pushed = [
            sum(p[x] * kernel[x][y] for x in range(size))
            for y in range(size)
        ]
        require(
            all(entry > 0 for entry in pushed),
            "repair extinguishes an atom",
        )
    return {
        "statement": (
            "the exact first-law split holds on rational matrix "
            "instances with its bilinear cross term, and the repair "
            "kernel preserves full support exactly, the finite-step "
            "unattainability mechanism"
        ),
        "instances": 20,
    }


def fluctuation_certificate() -> dict[str, Any]:
    """Exact fluctuation identities over the rationals.

    The exponential of the fluctuating entropy production
    sigma(x, y) = log(p x / pi x) - log(q y / pi y) is the rational
    ratio (p x / pi x) / (q y / pi y), so the integral fluctuation
    identity, the pointwise and level-set detailed fluctuation
    relations, and finite Onsager reciprocity all close exactly over
    the rationals with no numerical logarithm. Level sets of sigma
    coincide with level sets of the ratio because the logarithm is
    injective on the positives.
    """
    rng = random.Random(20260810)
    instances = 0
    level_classes = 0
    for _ in range(25):
        size = rng.randint(3, 8)
        blocks = rng.randint(1, 3)
        labels, pi_raw = random_instance(rng, size, blocks)
        total = sum(pi_raw)
        pi = [w / total for w in pi_raw]
        kernel, _ = heat_bath(labels, pi)
        raw = [Fraction(rng.randint(1, 12)) for _ in range(size)]
        tot = sum(raw)
        p = [r / tot for r in raw]
        q = [
            sum(p[x] * kernel[x][y] for x in range(size))
            for y in range(size)
        ]
        require(all(w > 0 for w in q), "pushforward lost support")
        ratio = [
            [(p[x] / pi[x]) / (q[y] / pi[y]) for y in range(size)]
            for x in range(size)
        ]
        ift = sum(
            p[x] * kernel[x][y] / ratio[x][y]
            for x in range(size)
            for y in range(size)
        )
        require(ift == 1, "integral fluctuation identity failed")
        for x in range(size):
            for y in range(size):
                require(
                    p[x] * kernel[x][y]
                    == ratio[x][y] * q[y] * kernel[y][x],
                    "pointwise detailed fluctuation relation failed",
                )
        groups: dict[Fraction, list[tuple[int, int]]] = {}
        for x in range(size):
            for y in range(size):
                groups.setdefault(ratio[x][y], []).append((x, y))
        for value, pairs in groups.items():
            forward = sum(p[x] * kernel[x][y] for x, y in pairs)
            reverse = sum(q[y] * kernel[y][x] for x, y in pairs)
            require(
                forward == value * reverse,
                "level-set fluctuation identity failed",
            )
            level_classes += 1
        f = [Fraction(rng.randint(-9, 9)) for _ in range(size)]
        g = [Fraction(rng.randint(-9, 9)) for _ in range(size)]
        lhs = sum(
            pi[x] * kernel[x][y] * f[x] * g[y]
            for x in range(size)
            for y in range(size)
        )
        rhs = sum(
            pi[x] * kernel[x][y] * g[x] * f[y]
            for x in range(size)
            for y in range(size)
        )
        require(lhs == rhs, "Onsager reciprocity failed")
        instances += 1
    # The mean entropy production equals the certified descent, checked
    # at high precision alongside the exact identities.
    import mpmath

    mp = mpmath.mp
    saved_dps = mp.dps
    try:
        mp.dps = 60

        def klp(qv, pv):
            out = mp.mpf(0)
            for qi, pi_ in zip(qv, pv):
                if qi > 0:
                    out += qi * mp.log(qi / pi_)
            return out

        mean_gap = None
        rng2 = random.Random(20260811)
        for _ in range(10):
            size = rng2.randint(3, 7)
            labels = [rng2.randrange(2) for _ in range(size)]
            pi_f = [Fraction(rng2.randint(1, 9)) for _ in range(size)]
            tot_pi = sum(pi_f)
            pi_n = [w / tot_pi for w in pi_f]
            kernel, _ = heat_bath(labels, pi_f)
            raw = [Fraction(rng2.randint(1, 9)) for _ in range(size)]
            tot_p = sum(raw)
            p = [r / tot_p for r in raw]
            q = [
                sum(p[x] * kernel[x][y] for x in range(size))
                for y in range(size)
            ]

            def mpf_of(fr):
                return mp.mpf(fr.numerator) / fr.denominator

            mean_sigma = mp.mpf(0)
            for x in range(size):
                for y in range(size):
                    if kernel[x][y] == 0:
                        continue
                    sigma = mp.log(mpf_of(p[x]) / mpf_of(pi_n[x])) \
                        - mp.log(mpf_of(q[y]) / mpf_of(pi_n[y]))
                    mean_sigma += mpf_of(p[x] * kernel[x][y]) * sigma
            descent = klp(
                [mpf_of(v) for v in p], [mpf_of(v) for v in pi_n]
            ) - klp(
                [mpf_of(v) for v in q], [mpf_of(v) for v in pi_n]
            )
            gap = abs(mean_sigma - descent)
            require(
                gap <= mp.mpf("1e-45"),
                "mean entropy production drifts from the descent",
            )
            if mean_gap is None or gap > mean_gap:
                mean_gap = gap
        mean_gap_str = mpmath.nstr(mean_gap, 8)
    finally:
        mp.dps = saved_dps
    return {
        "statement": (
            "the integral fluctuation identity, the pointwise and "
            "level-set detailed fluctuation relations, and finite "
            "Onsager reciprocity close exactly over the rationals for "
            "the repair kernel, and the mean entropy production "
            "reproduces the certified relative-entropy descent"
        ),
        "instances": instances,
        "level_set_classes": level_classes,
        "mean_descent_max_gap": mean_gap_str,
    }


def cap_first_law_certificate() -> dict[str, Any]:
    """Cap first law with exact remainder and the cap Clausius bound.

    With modular Hamiltonian K = -log pi split as K = 2 pi_circ B + Z
    for a fibre-measurable central charge Z, the identity
    S(p) - S(pi) = 2 pi_circ Delta<B> + Delta<Z> - D(p||pi) and the
    repair-step inequality 2 pi_circ Delta<B> <= Delta S are replayed
    at sixty digits; the central-charge mean is conserved by the
    kernel exactly.
    """
    import mpmath

    mp = mpmath.mp
    saved_dps = mp.dps
    try:
        mp.dps = 60
        rng = random.Random(20260812)
        two_pi = 2 * mp.pi

        def mpf_of(fr):
            return mp.mpf(fr.numerator) / fr.denominator

        def shannon_mp(vec):
            out = mp.mpf(0)
            for v in vec:
                if v > 0:
                    out -= v * mp.log(v)
            return out

        def klp(qv, pv):
            out = mp.mpf(0)
            for qi, pi_ in zip(qv, pv):
                if qi > 0:
                    out += qi * mp.log(qi / pi_)
            return out

        identity_max = mp.mpf(0)
        clausius_min = None
        conserved_checks = 0
        for _ in range(12):
            size = rng.randint(3, 7)
            labels = [rng.randrange(3) for _ in range(size)]
            pi_f = [Fraction(rng.randint(1, 9)) for _ in range(size)]
            tot = sum(pi_f)
            pi_n = [w / tot for w in pi_f]
            kernel, _ = heat_bath(labels, pi_f)
            # fibre-measurable central charge
            zvals = {lab: Fraction(rng.randint(-6, 6), rng.randint(1, 3))
                     for lab in set(labels)}
            z = [zvals[lab] for lab in labels]
            pi_mp = [mpf_of(v) for v in pi_n]
            bc = [
                (-mp.log(pi_mp[i]) - mpf_of(z[i])) / two_pi
                for i in range(size)
            ]
            raw = [Fraction(rng.randint(1, 9)) for _ in range(size)]
            tot_p = sum(raw)
            p = [r / tot_p for r in raw]
            q = [
                sum(p[x] * kernel[x][y] for x in range(size))
                for y in range(size)
            ]
            # exact conservation of the central-charge mean
            require(
                sum(q[i] * z[i] for i in range(size))
                == sum(p[i] * z[i] for i in range(size)),
                "central charge mean drifts under repair",
            )
            conserved_checks += 1
            p_mp = [mpf_of(v) for v in p]
            q_mp = [mpf_of(v) for v in q]
            # exact-remainder identity against the split
            for state in (p_mp, q_mp):
                lhs = shannon_mp(state) - shannon_mp(pi_mp)
                rhs = two_pi * (
                    sum(state[i] * bc[i] for i in range(size))
                    - sum(pi_mp[i] * bc[i] for i in range(size))
                ) + (
                    sum(state[i] * mpf_of(z[i]) for i in range(size))
                    - sum(pi_mp[i] * mpf_of(z[i]) for i in range(size))
                ) - klp(state, pi_mp)
                gap = abs(lhs - rhs)
                if gap > identity_max:
                    identity_max = gap
            require(
                identity_max <= mp.mpf("1e-45"),
                "cap first-law identity drifts",
            )
            # cap Clausius margin for the repair step
            margin = (shannon_mp(q_mp) - shannon_mp(p_mp)) - two_pi * (
                sum(q_mp[i] * bc[i] for i in range(size))
                - sum(p_mp[i] * bc[i] for i in range(size))
            )
            require(
                margin >= -mp.mpf("1e-45"),
                "cap Clausius inequality fails",
            )
            if clausius_min is None or margin < clausius_min:
                clausius_min = margin
        identity_str = mpmath.nstr(identity_max, 8)
        clausius_str = mpmath.nstr(clausius_min, 8)
    finally:
        mp.dps = saved_dps
    return {
        "statement": (
            "with the modular Hamiltonian split K = 2 pi B + Z over a "
            "fibre-measurable central charge, the exact-remainder cap "
            "first law S(p) - S(pi) = 2 pi Delta B + Delta Z - D and "
            "the repair-step cap Clausius bound 2 pi Delta B <= "
            "Delta S replay at sixty digits, with the central-charge "
            "mean conserved exactly by the kernel"
        ),
        "instances": 12,
        "central_conservation_checks": conserved_checks,
        "identity_max_gap": identity_str,
        "clausius_min_margin": clausius_str,
    }


def open_receipts() -> dict[str, Any]:
    return {
        "THERMO-GLOBAL": (
            "the weighted local third-axiom objective must carry its "
            "declared global representation for both the state and the "
            "transition object on the thermodynamic branch"
        ),
        "THERMO-COMMON-REFERENCE": (
            "the state and transition optimizers must descend from one "
            "common source-derived quotient weight; the transition "
            "matrix may not be manufactured from a desired equilibrium "
            "output"
        ),
        "THERMO-REALIZATION": (
            "the source-derived collar transition matrix must equal the "
            "conditional-resampling kernel or pass stochasticity, "
            "stationarity, and protected-charge preservation with "
            "detailed balance where microscopic reversibility is "
            "claimed; the committed collar-matrix probe measures these "
            "properties on the earned sixty-four-k quotient matrix and "
            "exhausts all fifteen nonempty field-subset projections of its "
            "four committed fields. The repair-load quotient gives an eight-state "
            "raw ergodic nonreversible H-theorem probe, but its declared "
            "record charge is constant, its stationary law is not "
            "identified with the state optimizer's source reference, "
            "and the fine chain's only recurrent restriction is a "
            "singleton freezeout state"
        ),
        "THERMO-ENERGY-CLOCK": (
            "one modular charge must be identified with physical energy "
            "and calibrated clock time; the cap branch supplies the "
            "finite modular split, the calibration stays open"
        ),
        "THERMO-LOW-T-REFINEMENT": (
            "a continuum third law needs refinement-uniform control of "
            "the ground multiplicity, gap, and state-count growth; the "
            "fixed-regulator theorem needs none of it"
        ),
    }


def build_receipt() -> dict[str, Any]:
    body: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "three_map_distinction": (
            "the strict-descent normalizer settles public facts, the "
            "conditional-resampling kernel relaxes unresolved degrees, "
            "and modular flow supplies the reversible ordering; only "
            "the kernel carries the entropy inequality"
        ),
        "lean_bindings": {
            "conditional_repair": (
                "Thermodynamics/FiniteConditionalRepair.lean: kl_nonneg, "
                "logSum, kl_push_le, heatBath_row_sum, "
                "heatBath_idempotent, heatBath_detailedBalance, "
                "heatBath_stationary, heatBath_fixes_fiberObservable, "
                "heatBath_row_optimal, heatBath_secondLaw, "
                "gibbs_pythagorean, gibbs_minimizer, excitedMass_le, "
                "excitedMass_lt_of_beta_large, gibbs_beta_injective, "
                "clausius, landauer, kl_eq_energy_sub_shannon, "
                "heatBath_preserves_pos, push_total, mixture_row_sum, "
                "mixture_stochastic, mixture_stationary, block_entropy_le"
            ),
            "stationary_realization": (
                "Thermodynamics/StationaryRealization.lean: "
                "stationary_secondLaw, constantObservable_fixed, "
                "directedLazy3_stationary, "
                "directedLazy3_not_detailedBalance, "
                "directedLazy3_secondLaw"
            ),
            "first_law": (
                "Thermodynamics/FirstLawIdentity.lean: firstLaw_split"
            ),
            "record_channel": (
                "EventAlgebra/PartitionPinchingCP.lean: kraus_complete, "
                "partitionPinching_kraus_form"
            ),
            "fluctuation": (
                "Thermodynamics/FluctuationTheorems.lean: sigmaEP, "
                "integral_fluctuation, crooks_pointwise, "
                "crooks_level_set, sigma_mean_eq_kl_descent, "
                "correlation_symm, heatBath_integral_fluctuation, "
                "heatBath_crooks, heatBath_correlation_symm"
            ),
            "cap_first_law": (
                "Thermodynamics/CapFirstLaw.lean: kl_self, "
                "shannon_ref_eq_modular_energy, cap_firstLaw_exact, "
                "modular_energy_split, cap_firstLaw_split, "
                "cap_clausius_of_central_conserved, heatBath_nonneg, "
                "push_heatBath_fixes_mean, heatBath_cap_clausius"
            ),
            "einstein_premise_link": (
                "Thermodynamics/EinsteinPremiseLink.lean: "
                "shannon_diff_eq_pairing_sub_kl, thermoFirstLawData, "
                "thermoFirstLawData_passes, "
                "thermo_first_law_on_simplex_tangent, "
                "repair_variation_mem_massZero, "
                "repair_variation_central_pairing_zero"
            ),
        },
        "strict_descent_control": strict_descent_control(),
        "kernel_algebra": kernel_algebra_certificate(),
        "entropy_checks": entropy_certificate(),
        "low_temperature": low_temperature_certificate(),
        "zeroth_law": zeroth_law_certificate(),
        "clausius": clausius_certificate(),
        "first_law_and_support": first_law_and_support_certificate(),
        "fluctuation": fluctuation_certificate(),
        "cap_first_law": cap_first_law_certificate(),
        "record_tower": pinching_tower_certificate(),
        "law_map": {
            "zeroth": (
                "equal intensive multipliers at the interior equilibrium "
                "of the additive-contact Gibbs branch; transitivity "
                "through any thermometer with two distinct levels"
            ),
            "first": (
                "exact bookkeeping split with explicit bilinear cross "
                "term, plus dynamical conservation of every "
                "fibre-constant charge under the repair kernel"
            ),
            "second": (
                "relative entropy to the common reference is "
                "nonincreasing under conditional repair; the modular "
                "form Delta S >= beta Q follows on the identified "
                "energy branch; the inequality is the mean of the "
                "fluctuating entropy production, whose exponential "
                "averages to one exactly and whose level sets obey "
                "the detailed fluctuation relation"
            ),
            "third": (
                "excited Gibbs mass below (d-g0)/g0 exp(-beta Delta) at "
                "every finite regulator, entropy limit log g0, and "
                "finite-step unattainability for faithful repair maps"
            ),
        },
        "open_receipts": open_receipts(),
        "no_comparison": (
            "the certificate reads no laboratory data and asserts no "
            "physical temperature; every physical attachment lives in "
            "the named open receipts"
        ),
    }
    body["receipt_sha256"] = tagged_sha256(canonical_json_bytes(body))
    return body


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    receipt = build_receipt()
    if args.write:
        RECEIPT_PATH.parent.mkdir(parents=True, exist_ok=True)
        RECEIPT_PATH.write_text(
            canonical_json_bytes(receipt).decode() + "\n"
        )
        print(f"wrote {RECEIPT_PATH}")
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
