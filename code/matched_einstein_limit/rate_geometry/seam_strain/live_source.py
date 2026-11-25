"""Keep both native marks with the written strain after pumping one mobile.

At each departing port p the three retained registers have Kraus matrix K_p.
Their 8-dimensional density is sum_p K_p (rail lottery ⊗ head state) K_p†.
Reference, ordinary and decorrelation use native formation and re-pairing. The
additional hidden direction solves a four-column linear system on the allowed
head algebra. Positivity gives a two-sided interval about the formed reference;
native production of every such perturbed density is a separate premise.
"""

from functools import lru_cache

import numpy as np

from . import hidden as qubit_hidden
from .source import measure as source


def input_witness():
    """Unit coefficient direction on XZ+ZX, center of the exact hidden tangent."""
    x, _, z = source.formation.PAULI
    return (np.kron(x, z) + np.kron(z, x)) / (4 * np.sqrt(2))


def apply_channel(kraus, heads, rail_weights):
    state = np.kron(np.diag(rail_weights), heads)
    return np.einsum("pai,ij,pbj->ab", kraus, state, kraus.conj(), optimize=True)


@lru_cache(maxsize=2)
def head_channel(coupled=True):
    """Kraus matrices with indices (discarded port, retained strain/head, incoming rail/head)."""
    board = source.formation.federation(1)
    route = source.formation.region(board)
    a, _, _ = source.system(board, route, -1, coupled=coupled)
    wave = next(
        w for t, w in source.fast_flow(a, source.launch(board, route) / 2) if t == 4
    )
    # Sum over every old-mobile port. Its position is never selected.
    kraus = wave.reshape(2, -1, 4, 8).transpose(1, 0, 2, 3).reshape(-1, 8, 8)
    kraus.setflags(write=False)
    return kraus


def propagate_heads(heads, coupled=True):
    """Apply the same four-tick channel using the fixed 2:1 rail lottery."""
    return apply_channel(head_channel(coupled), heads, source.RAILS["reference"])


@lru_cache(maxsize=2)
def postpump(coupled=True):
    """Four-tick 80-cell vacancy pump states in strain-major, head-minor order.

    The 'hidden' array is a traceless derivative at the reference, per unit
    correlation coefficient norm. Each of the other three arrays is a density.
    """
    board = source.formation.federation(1)
    rho, product = source.form(board, source.formation.region(board), -1)
    ch = head_channel(coupled)
    out = {
        name: apply_channel(ch, product if name == "decorrelated" else rho, weights)
        for name, weights in source.RAILS.items()
    }
    out["hidden"] = propagate_heads(input_witness(), coupled)
    for state in out.values():
        state.setflags(write=False)
    return out


def source_calibration():
    """Independent initial energy and modular contrasts for the incoming arms."""
    board = source.formation.federation(1)
    route = source.formation.region(board)
    rho, product = source.form(board, route, -1)
    k = -source.pairing.log_density(rho)
    return dict(
        ordinary_modular=float(np.log(2) / 3),
        decorrelated_modular=float(np.trace(k @ (product - rho)).real),
        hidden_modular=float(np.trace(k @ input_witness()).real),
        ordinary_bare_energy=1 / 3,
        ordinary_coupled_energy=41 / 120,
        decorrelated_bare_energy=0.0,
        decorrelated_coupled_energy=0.0,
        hidden_bare_energy=0.0,
        hidden_coupled_energy=0.0,
    )


def precise_wave(a, start):
    """Split each amplitude at bit 36, with total fixed-point scale 2**82."""
    from math import factorial, isqrt
    import sympy as sp
    from .source import certificate as cert

    radix = 1 << 36
    precision = 1 << 82
    degree = 26
    d, n = start.shape
    hi = np.column_stack((start * (precision // (2 * radix)), np.zeros_like(start)))
    lo = np.zeros_like(hi)
    row = int(abs(a).sum(axis=1).max())
    length = sp.Rational(row, 40 * cert.STEPS)
    assert length <= 1 and (a.T != a).nnz == 0
    rounding = sp.Rational(3 * degree * (1 + isqrt(2 * d)), 2 * precision)
    err = (
        4 * cert.STEPS * (rounding + 2 * length ** (degree + 1) / factorial(degree + 1))
    )
    assert err < sp.Rational(1, 10**15)
    for _ in range(4 * cert.STEPS):
        th, tl = hi.copy(), lo.copy()
        for k in range(1, degree + 1):
            den = 40 * cert.STEPS * k
            assert row * int(abs(th).max()) < 2**61 and (row + den) * radix < 2**61
            ph, pl = a @ th, a @ tl
            dh = np.column_stack((ph[:, n:], -ph[:, :n]))
            dl = np.column_stack((pl[:, n:], -pl[:, :n]))
            quotient, rem = np.divmod(dh, den)
            carry = (rem * radix + dl + den // 2) // den
            th = quotient + (carry >> 36)
            tl = carry & (radix - 1)
            lo += tl
            assert int(abs(hi).max()) + int(abs(th).max()) + 1 < 2**61
            hi += th + (lo >> 36)
            lo &= radix - 1
    return hi, lo, err


def rail_grams(hi, lo):
    """Exact Gram of the pump over the discarded position, preserving both heads."""
    mask = (1 << 18) - 1
    result = []
    for rail in (0, 1):

        def rearrange(w):
            # Groups are (real/imag, output strain, output head, input head).
            w = w.reshape(2, -1, 4, 2, 2, 4)
            return w[:, :, :, :, rail, :].transpose(1, 3, 0, 2, 4).reshape(-1, 64)

        u, v = rearrange(hi), rearrange(lo)
        digits = [v & mask, v >> 18, u & mask, (u >> 18) & mask, u >> 36]
        assert len(v) * max(int(abs(z).max()) for z in digits) ** 2 < 2**62
        gram = sum(
            (s.T @ t).astype(object) * (1 << (18 * (i + j)))
            for i, s in enumerate(digits)
            for j, t in enumerate(digits)
        )
        # G[b,j,a,i] = sum_p K[p,b,j]^* K[p,a,i]
        real = (gram[:32, :32] + gram[32:, 32:]).reshape(8, 4, 8, 4)
        imag = (gram[:32, 32:] - gram[32:, :32]).reshape(8, 4, 8, 4)
        result.append((real, imag))
    return result


def rational_read(grams, head, numerators):
    """Exact reduced matrix for a rational center (integer Gram, no eigensolver)."""
    import sympy as sp

    re, im = head.as_real_imag()
    # Put all entries on one shared denominator to keep integer contractions cheap.
    den = sp.ilcm(*(sp.denom(z) for part in (re, im) for z in part))
    r, i = (np.array(part * den, dtype=object) for part in (re, im))
    r = np.vectorize(int, otypes=[object])(r).T.reshape(16)
    i = np.vectorize(int, otypes=[object])(i).T.reshape(16)
    nr = np.zeros((8, 8), dtype=object)
    ni = nr.copy()
    for weight, (rg, ig) in zip(numerators, grams):
        # Q_ab = sum_ij G[b,j,a,i] rho[i,j].
        rg = rg.transpose(2, 0, 1, 3).reshape(8, 8, 16)
        ig = ig.transpose(2, 0, 1, 3).reshape(8, 8, 16)
        nr += weight * (rg @ r - ig @ i)
        ni += weight * (rg @ i + ig @ r)
    return (sp.Matrix(nr) + sp.I * sp.Matrix(ni)) / (3 * (1 << 164) * den)


def exact_witness(coupled_grams, col_error, formed, state_error):
    """Solve the three exact Bloch constraints with a positive modular slope."""
    import sympy as sp
    from .source import certificate as cert

    axes = qubit_hidden.correlation_axes()
    # The selected four columns are (ZZ, XZ, ZY, ZX). Fix the ZX coefficient
    # at one before Euclidean normalization. B contains the first three.
    chosen = [8, 2, 7, 6]
    dirs = [axes[j] for j in chosen]
    outputs = [rational_read(coupled_grams, n, (2, 1)) for n in dirs]

    def reduced(m):
        q = sp.Matrix(2, 2, lambda a, b: sum(m[4 * a + k, 4 * b + k] for k in range(4)))
        return sp.Matrix([2 * sp.re(q[0, 1]), -2 * sp.im(q[0, 1]), q[0, 0] - q[1, 1]])

    nominal = sp.Matrix.hstack(*(reduced(f) for f in outputs))
    b = nominal[:, :3]
    e = 4 * col_error + 4 * col_error**2
    inv = b.inv()
    inverse = inv.applyfunc(abs)
    inverse_norm = max(sum(inverse.row(j)) for j in range(3))
    eta = 3 * inverse_norm * e
    assert eta < sp.Rational(1, 100)
    u0 = sp.Matrix([0, 1, 0, 1])
    # Exact fixed ZX=1 solution differs from u0 by B^{-1} M u0 in other entries.
    errors = inverse * ((nominal * u0).applyfunc(abs) + sp.ones(3, 1) * 2 * e)
    # Componentwise bound accounts for unknown true B when solving for correction.
    maximum = max(errors) / (1 - eta)
    errors = errors + 3 * e * maximum * (inverse * sp.ones(3, 1))
    # ||normalize(u)-normalize(u0)||_2 <= 2||u-u0||_2/||u0||_2
    # and ||u-u0||_2 <= ||u-u0||_1. The rational factor 3/2 bounds sqrt(2).
    unit_error = sp.Rational(3, 2) * sum(errors)
    assert unit_error < sp.Rational(1, 10**7)
    log, elog = cert.prior.entropy_coordinates(
        formed, state_error, sp.Rational(1, 10**8), dirs
    )
    inv_sqrt = sp.Rational(float(1 / np.sqrt(2)))
    assert abs(2 * inv_sqrt**2 - 1) < sp.Rational(1, 2**51)
    mid = (log[1] + log[3]) * inv_sqrt
    # The four-head-state log has spectral norm < 21, and the Pauli/4 axes
    # make ||head(d)-head(d0)||_1 <= ||d-d0||_2.
    slope_err = 2 * elog + 21 * unit_error + sp.Rational(42, 2**51)
    assert mid - slope_err > sp.Rational(33, 1000)
    # A sharper dyadic center solves the actual fixed-point minor.
    u = (-inv * nominal[:, 3]).col_join(sp.Matrix([1]))
    center = np.array(u, dtype=float).ravel()
    center /= np.linalg.norm(center)
    proposal = sp.Matrix([sp.Rational(float(v)) for v in center])
    assert proposal[3] > 0
    norm_defect = abs(proposal.dot(proposal) - 1)
    rounding = sum(abs(z) for z in u - proposal / proposal[3])
    raw_error = 3 * inverse_norm * e * sum(abs(z) for z in u) / (1 - eta)
    center_error = norm_defect + 2 * (rounding + raw_error)
    assert center_error < sp.Rational(1, 10**13)
    refined = dict(
        coefficients=center.tolist(), l2_error=cert.enclosed(center_error, 0)[1]
    )
    return unit_error, cert.enclosed(mid, slope_err), outputs[1] + outputs[3], refined


@lru_cache(maxsize=1)
def formation_center():
    from .source import certificate as cert

    board = source.formation.federation(1)
    return cert.formed_state(
        board, source.formation.region(board), -1, cert.prior.weights()[7]
    )


@lru_cache(maxsize=2)
def pump_center(coupled=True):
    board = source.formation.federation(1)
    route = source.formation.region(board)
    a, _, _ = source.system(board, route, -1, coupled)
    hi, lo, colerr = precise_wave(a, source.launch(board, route))
    return hi, lo, colerr, rail_grams(hi, lo)


def certify():
    import sympy as sp
    from .source import certificate as cert

    board = source.formation.federation(1)
    (formed, prod), ferr, (mlo, mhi) = formation_center()
    dirs = qubit_hidden.correlation_axes()
    inv_sqrt = sp.Rational(float(1 / np.sqrt(2)))
    result = {}
    for coupled in (True, False):
        hi, lo, colerr, grams = pump_center(coupled)
        if coupled:
            direction_error, modular, raw_hidden, refined = exact_witness(
                grams, colerr, formed, ferr
            )
        else:
            raw_hidden = rational_read(grams, dirs[2] + dirs[6], (2, 1))
        exact = {
            "reference": rational_read(grams, formed, (2, 1)),
            "ordinary": rational_read(grams, formed, (1, 2)),
            "decorrelated": rational_read(grams, prod, (2, 1)),
            "hidden": raw_hidden * inv_sqrt,
        }
        error = {}
        for name, state in postpump(coupled).items():
            diff = exact[name] - cert.prior.sym_matrix(state)
            rr, ii = diff.as_real_imag()
            entry = sum(abs(z) for mat in (rr, ii) for z in mat)
            native = 4 * colerr + 4 * colerr**2
            if name == "hidden":
                native += direction_error + sp.Rational(2, 2**51)
            else:
                native += (2 if name == "decorrelated" else 1) * ferr
            error[name] = cert.enclosed(entry + native, 0)[1]
        # Direct per-column bound comparing the center from the exact replay
        # with the actual floating Kraus bank used by arbitrary subsequent reads.
        fast = (
            head_channel(coupled)
            .reshape(board.n, 2, 4, 8)
            .transpose(1, 0, 2, 3)
            .reshape(-1, 8)
        )
        # Both integer limbs convert to exact doubles. Reconstruction in one
        # addition costs at most 2^-52 componentwise (unit column plus tiny error).
        assert abs(hi).max() < 2**52
        rr, ii = np.split(hi / float(2**46) + lo / float(2**82), 2, axis=1)
        peak = max(np.max(abs(rr - fast.real)), np.max(abs(ii - fast.imag)))
        c_error = colerr + 125 * (sp.Rational(2 * float(peak)) + sp.Rational(1, 2**51))
        delta = 2 * c_error
        result["coupled" if coupled else "uncoupled"] = dict(
            nominal_state_trace_norm_errors=error,
            numerical_kraus_column_error=cert.enclosed(c_error, 0)[1],
            kraus_channel_trace_norm_error_per_unit_input=cert.enclosed(
                2 * delta + delta**2, 0
            )[1],
        )
    # Each |±| normalized head direction has op norm <= its HS norm = 1/2.
    # Thus the faithful reference lower bound 1e-8 permits e.g. steps 1e-8.
    return dict(
        kernel_basis=["ZZ", "XZ", "ZY", "ZX"],
        nominal_unit_coefficients=[0.0, float(inv_sqrt), 0.0, float(inv_sqrt)],
        refined_unit_center=refined,
        normalized_direction_l2_error_upper=cert.enclosed(direction_error, 0)[1],
        hidden_modular_nats=modular,
        ordinary_modular_nats="log(2)/3",
        decorrelation_modular_nats=[cert.enclosed(mlo, 0)[0], cert.enclosed(mhi, 0)[1]],
        certified_two_sided_step=1e-8,
        channels=result,
    )


def read_record():
    def matrices(coupled):
        return {
            name: dict(real=m.real.tolist(), imag=m.imag.tolist())
            for name, m in postpump(coupled).items()
        }

    return dict(
        cells=80,
        forming_charge=-1,
        pump_charge=-1,
        pump_ticks=4,
        retained_order=["strain", "left_head", "right_head"],
        discarded="old mobile position",
        ordinary_and_decorrelated="densities",
        hidden="derivative about reference",
        source=source_calibration(),
        coupled=matrices(True),
        uncoupled=matrices(False),
        bounds=certify(),
    )


if __name__ == "__main__":
    import json

    print(json.dumps(read_record(), indent=2, allow_nan=False))
