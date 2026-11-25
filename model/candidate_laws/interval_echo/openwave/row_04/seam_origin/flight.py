"""Propagate independently commissioned material with the filed source and reader."""

from functools import lru_cache
from math import log, sqrt

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.row_30_memory import experiment as gateway
from model.executable_model import federation
from model.executable_model.carrier import Carrier
from model.experiments.interval_echo_einstein_station.native_tidal_response import measure as tidal
from .. import open_pair
from ..commissioned_medium import experiment as transfer


UNITS = 10**11
GATE = 0.94
NUMERICAL_ALLOWANCE = 1e-7
SEEDS = range(178123, 178251)


def bins(home, total=1.0):
    """Round each probability independently, including the inclusive complement."""
    home = np.asarray(home)
    if (not np.isfinite(home).all() or home.min() < -1e-12 or
            home.sum() > total + 2e-11 or not 0 <= total <= 1 + 2e-11):
        raise ArithmeticError("Invalid projected record law or preparation mass")
    out = np.rint(np.r_[home, 1 - home.sum()] * UNITS).astype(np.int64).tolist()
    if min(out) < 0 or abs(sum(out)-UNITS) > len(out)/2:
        raise ArithmeticError("Invalid rounded record law")
    return out


def material(board, seed, fair):
    # formation returns the complete two-bit history as well as transport diagnostics.
    record = transfer.donor.formation(board, seed, 0, fair)[-1]
    hexwords = record["coin_pairs_packed_little_bit_first_hex"]
    words = np.frombuffer(bytes.fromhex(hexwords), dtype=np.uint8)
    bits = np.unpackbits(words, bitorder="little", count=2*len(board.seams))
    closed = {e for e, mark in zip(board.seams, bits.reshape(-1, 2).all(axis=1), strict=True) if mark}
    if len(closed) != record["written"]:
        raise ArithmeticError("Formation history and write count differ")
    return closed, hexwords


@lru_cache(maxsize=16)
def acquired(kappa, x_blocked, z_blocked):
    # The acquisition addresses X, Z and the original degree reference alone.
    # check_meter verifies both memory classes on the vacant-ended meter.
    closed = {e for e, cut in ((transfer.X, x_blocked), (transfer.Z, z_blocked)) if cut}
    return transfer.acquire(federation(0), closed, kappa)


def fly(board, home, closed, kappa, t, *, comparisons=False):
    exchange, contrast, launch, detector = transfer.generator(board, closed, home)
    rho, prep = acquired(kappa, transfer.X in closed, transfer.Z in closed)
    mass = rho.trace().real
    failed = prep["exhaustion"] + prep["failed_with_reset_unavailable"]
    if abs(mass + failed - 1) > 2e-12:
        raise ArithmeticError("Incomplete acquisition")
    amplitude = transfer.filed.flight_amplitudes(exchange, contrast, kappa, [t], launch, detector)[0]
    home_p = bins(transfer.home_bins(amplitude, rho), mass)
    # Poisson full swaps act on the acquired diagonal, also when a donor blocks a pulse.
    stochastic = expm_multiply((-t) * exchange, launch @ rho.diagonal().real)
    if abs(stochastic.sum() - mass) > 2e-11 or min(stochastic) < -1e-13:
        raise ArithmeticError("Swap flight lost population")
    result = {"acquired": home_p, "random_swaps": bins(stochastic[detector], mass),
              "source_exhaustion": round(prep["exhaustion"] * UNITS),
              "reset_unavailable": round(prep["failed_with_reset_unavailable"] * UNITS)}
    if comparisons:
        h = exchange + sparse.diags(kappa * contrast, dtype=float)
        direct = expm_multiply(-1j * t * h, launch)
        np.testing.assert_allclose(direct.conj().T @ direct, np.eye(2), atol=2e-11, rtol=0)
        # Compare effects, which also test the relative phase of the two source columns.
        effect = amplitude[:, :, None].conj() * amplitude[:, None, :]
        ref = direct[detector]
        np.testing.assert_allclose(effect, ref[:, :, None].conj()*ref[:, None, :], atol=2e-11, rtol=0)
        result["filed_input"] = bins(transfer.home_bins(amplitude, transfer.FILED_DENSITY),
                                     transfer.FILED_DENSITY.trace().real)
        result["optimal_XZ"] = round(float(np.linalg.eigvalsh(amplitude.conj().T @ amplitude)[-1]) * UNITS)
    return result


def historical(board, home):
    guarded = tidal.build_station(board, tidal.local_calibration()[0])[3]
    filed = set(open_pair.geometry()[0])
    masks = {"first_gateway": {gateway.geometry(board)[0]}, "tidal_sham": guarded,
             "blank_control": set(), "eleven_control": filed}
    return {name: {"written": [list(e) for e in sorted(closed)],
                   **fly(board, home, closed, 0.5, 37, comparisons=True)}
            for name, closed in masks.items()}


def constituents(board, home):
    a, b, c = open_pair.PATH
    closed = set(open_pair.geometry()[0])
    active = Carrier("filed neutral material", board.n, [e for e in board.seams if e not in closed])
    k = active.laplacian()/2
    start = np.zeros((board.n, 3))
    start[b, 0] = 1
    start[[a, c], 1] = [1/sqrt(2), -1/sqrt(2)]
    start[[a, c], 2] = 1/sqrt(2)
    u, v = np.array(home).T
    local = np.arange(12*(b//12), 12*(b//12+1))
    rows = []
    for t in (1, 7, 37, 150):
        probs = abs(expm_multiply(-1j*t*k, start))**2
        np.testing.assert_allclose(probs.sum(axis=0), 1, rtol=0, atol=3e-12)
        center, tip = probs[:, 0], probs[:, 1]
        rows.append({"exposure": t,
                     **{name: bins(probs[local, j]) for j, name in enumerate(("center", "minus", "plus"))},
                     "independent_center_and_tip": bins(center[u]*tip[v] + center[v]*tip[u])})
    return {"home_ports": local.tolist(), "appointments": rows}


def ensemble(board, home):
    fair = transfer.donor.native()
    rows = []
    for seed in SEEDS:
        closed, packed = material(board, seed, fair)
        rows.append({"seed": seed, "coins_hex": packed, "writes": len(closed),
                     **fly(board, home, closed, 0, 7, comparisons=seed == SEEDS.start)})
    mean = sum(1-r["acquired"][-1]/UNITS for r in rows)/len(rows)
    swap_mean = sum(1-r["random_swaps"][-1]/UNITS for r in rows)/len(rows)
    upper = min(1, mean + sqrt(log(100)/(2*len(rows))) + NUMERICAL_ALLOWANCE)
    return {"strength": 0, "exposure": 7, "banks": rows, "mean_home": round(mean, 10),
            "mean_swap_home": round(swap_mean, 10), "one_sided_upper_99": round(upper, 10),
            "historical_benchmark": GATE, "reject_iid_mean_at_least_benchmark": upper < GATE}
