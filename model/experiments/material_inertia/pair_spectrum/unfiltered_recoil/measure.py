"""Enclose every final address, including interference with the unresolved spectrum."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import dispersion, measure as flight

HERE = Path(__file__).resolve().parent
SETTINGS = (("calibration", 0.125, 16384), ("held", 0.25, 32768))
PAD = 1e-12


def norm(matrix):
    return float(np.linalg.norm(matrix, 2))


def enclosure(
    basis,
    energies,
    rest_coeff,
    admission,
    residual,
    band_err,
    rest_err,
    coordinate,
    impulse,
    strength,
    dwell,
):
    rest = basis[:, :2]
    u = np.exp(-1j * dwell * (energies - energies[0]))
    # Each column stands for an independent initial rest vector. Taking the
    # operator norm below also covers every coherent vector in the doublet.
    signed = np.zeros((2, 2), complex)
    systematic = 0.0
    arms = []
    for sign in (-1, 1):
        launch = np.exp(sign * 1j * strength * impulse)[:, None] * rest
        coefficients = basis.T @ launch
        leak = norm(launch - basis @ coefficients) + band_err
        visible = basis @ (u[:, None] * coefficients)
        # Error between this vector and the true evolution of P*launch:
        # the initial projector error plus the Duhamel residual of the frame.
        evolution_err = band_err + dwell * residual + PAD
        read = coordinate[:, None] * visible
        cross = norm(read - basis @ (basis.T @ read)) + band_err + evolution_err
        arm_err = leak**2 + 2 * leak * cross + 2 * evolution_err
        systematic += arm_err
        arms.append(
            dict(
                sign=sign,
                outside_amplitude_bound=leak,
                cross_read_bound=cross,
                mean_error_bound=arm_err,
            )
        )
        signed += sign * (visible.conj().T @ read)
    signed = (signed + signed.conj().T) / 2
    # The true input rest vector differs by <= 2*rest_err from its normalized
    # projection into this numerical doublet. Each sign's X operator has norm 1.
    systematic += 8 * rest_err + PAD
    # On this doublet a coefficient error acts only on the nonconstant part
    # of the signed response. Center the 2x2 matrix to remove the identity.
    low, high = np.linalg.eigvalsh(signed)
    coeff_error = min(2.0, 2 * rest_err / np.sqrt(admission))
    systematic += (high - low) * coeff_error
    prediction = float(np.vdot(rest_coeff, signed @ rest_coeff).real)
    bounds = [prediction - systematic, prediction + systematic]
    return dict(
        strength=strength,
        dwell=dwell,
        arms=arms,
        displacement=bounds,
        mobility=[v / (2 * strength * dwell) for v in bounds],
        unconditional_read="all final addresses after initial admission",
    )


def run():
    board = flight.federation(0)
    excess, vacancy, _, _, h, cost = flight.defect.generators(board)
    source = json.loads((flight.HERE / "results.json").read_text())
    reduced = flight.PairReduction(h, cost)
    energy, basis, residual, _ = reduced.lowest(h, 24)
    if norm(basis.T @ basis - np.eye(24)) > PAD:
        raise ValueError("nonorthonormal Ritz frame")
    residual += PAD
    cut = (energy[1] + energy[2]) / 2
    if (
        reduced.inertia(cut)["upper_lines_below"] != 2
        or reduced.inertia(flight.CUT)["upper_lines_below"] != 24
    ):
        raise ValueError("spectral counts changed")
    band_err = (
        max(source["cluster_projector_error"], residual / (flight.CUT - energy[-1]))
        + PAD
    )
    # The parent Schur count also puts all lower bands below 30.
    delta = min(cut - energy[1], flight.CUT - energy[-1])
    rest_err = residual / (delta - residual) + PAD
    shell = dispersion.geometry(board)[0][1]
    outputs = []
    for launch in source["flights"]:
        if launch["admission"] < 0.001:
            continue
        owner = launch["owner"]
        axis = min(
            v // 12
            for edge in board.external
            for u, v in (edge, edge[::-1])
            if u // 12 == owner
        )
        f = shell[:, axis] / np.max(abs(shell[:, axis]))
        x = (f[excess // 12] + f[vacancy // 12]) / 2
        signed_load = f[excess // 12] - f[vacancy // 12]
        i = flight.defect.ordered_index(*launch["ports"], board.n)
        g = basis[i, :2]
        admission = float(g @ g)
        root = np.sqrt(admission)
        runs = {}
        for name, strength, dwell in SETTINGS:
            for key, impulse in (("common", x), ("signed_load", signed_load)):
                runs[f"{name}_{key}"] = enclosure(
                    basis,
                    energy,
                    g / root,
                    admission,
                    residual,
                    band_err,
                    rest_err,
                    x,
                    impulse,
                    strength,
                    dwell,
                )
        # Independently evolve both complete 57,360-component vectors at the
        # historical short dwell. The closure enclosure must include this
        # direct address reading even when its sign differs from the late one.
        short_time = 1 / 16
        short = enclosure(
            basis,
            energy,
            g / root,
            admission,
            residual,
            band_err,
            rest_err,
            x,
            x,
            0.25,
            short_time,
        )
        state = basis[:, :2] @ (g / root)
        both = state[:, None] * np.exp(0.25j * x[:, None] * np.array([-1, 1]))
        arrived = expm_multiply(-1j * short_time * h, both)
        np.testing.assert_allclose(
            (abs(arrived) ** 2).sum(axis=0), 1, atol=1e-12, rtol=0
        )
        direct = float(np.diff(x @ (abs(arrived) ** 2))[0])
        if (
            direct >= 0
            or not short["displacement"][0] <= direct <= short["displacement"][1]
        ):
            raise ValueError("complete short-time address reading failed")
        # Keep the clock apparatus and pulse fixed, halve all in-flight
        # generator rates. Its eigenspaces are identical; its phases halve.
        slow = enclosure(
            basis,
            energy / 2,
            g / root,
            admission,
            residual / 2,
            band_err,
            rest_err,
            x,
            x,
            0.25,
            32768,
        )
        if slow["displacement"][1] >= runs["held_common"]["displacement"][0]:
            raise ValueError("recoil cannot distinguish the slower dynamics")
        runs["half_rate_held_common"] = slow
        near = runs["calibration_common"]["displacement"]
        far = runs["held_common"]["displacement"]
        if near[0] <= 0:
            raise ValueError("calibration cannot predict a positive recoil")
        ratio = [far[0] / (4 * near[1]), far[1] / (4 * near[0])]
        outcome = (
            "pass"
            if ratio[0] >= 0.9 and ratio[1] <= 1.1
            else "reject"
            if ratio[0] > 1.1 or ratio[1] < 0.9
            else "inconclusive"
        )
        outputs.append(
            dict(
                owner=owner,
                axis=axis,
                initial_ports=launch["ports"],
                initial_admission_interval=[
                    (root - rest_err) ** 2,
                    (root + rest_err) ** 2,
                ],
                readings=runs,
                short_time_direct_displacement=direct,
                short_time_enclosure=short["displacement"],
                held_ratio=ratio,
                held_outcome=outcome,
            )
        )
    return dict(
        panels=outputs,
        residual_with_pad=residual,
        band_projector_bound=band_err,
        rest_projector_bound=rest_err,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = run()
    if args.output:
        with args.output.open("x") as out:
            json.dump(data, out, indent=2, allow_nan=False)
            out.write("\n")
        return
    receipt = (HERE / "results.json").read_bytes()
    if (
        hashlib.sha256(receipt).hexdigest()
        != "349a64bbf65a47129cd232d8f9b99957fa11edea960b9c0d40bfcbb482433eef"
    ):
        raise ValueError("the enrolled numerical receipt changed")
    saved = json.loads(receipt)
    for old, new in zip(saved["panels"], data["panels"], strict=True):
        if (
            old["owner"] != new["owner"]
            or old["axis"] != new["axis"]
            or old["initial_ports"] != new["initial_ports"]
        ):
            raise ValueError("source or coordinate changed")
        for name in old["readings"]:
            a = old["readings"][name]
            b = new["readings"][name]
            if (a["strength"], a["dwell"]) != (b["strength"], b["dwell"]):
                raise ValueError("changed impulse or exposure")
            # Different numerical eigensolvers must certify the same scale.
            # The bounds depend on the Ritz residual, so compare their
            # centres tightly and the enclosures by intersection and width.
            p, q = np.array(a["mobility"]), np.array(b["mobility"])
            if (
                abs(p.mean() - q.mean()) > 1e-10
                or max(p[0], q[0]) > min(p[1], q[1])
                or np.ptp(q) > 1.5 * np.ptp(p)
            ):
                raise ValueError(f"recoil enclosure changed for {name}")
        if (
            new["held_outcome"] != old["held_outcome"]
            or new["readings"]["held_common"]["mobility"][0] <= 0
        ):
            raise ValueError("held recoil gate changed")
    print(
        "Verified unfiltered long-flight recoil and held ratios on all three launches."
    )


if __name__ == "__main__":
    main()
