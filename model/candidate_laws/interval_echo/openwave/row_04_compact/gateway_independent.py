"""Recompute all 56 held terminal probabilities from the scored seam array alone.

Run as a script with ANCHOR --openwave-checkout CHECKOUT --output independent.json.
Fetch the official OpenWave main into the trusted checkout before running. The
receiver reads twelve sharp home ports, accepted exterior and exhausted source.
No producer acquisition, detector, generator or compiled receipt is imported.
"""
import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import re
import subprocess

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix


BASE = "f8fd02b2e0d96cfa0d868259daec2c9822375d34"
TASK = "openwave/xperiments/m10_oph/research/tasks/m10_4_task_details.md"
TASK_HASH = "af312744a6842aa093c8fdfc3961053bb6a5dd29f407f4b1ba11bff742cf3640"
BOARD_HASH = "99870398e3d516f899c600e9bf25ba53c7f513142ba8facbe1cedd94b23f915b"
BOARD = (Path(__file__).parent/"../../../../executable_model/carriers/FEDERATION_LEVEL0.npz").resolve()


def first_registration(checkout, anchor):
    if not isinstance(anchor, str) or not re.fullmatch(r"[a-f0-9]{40}", anchor):
        raise ValueError("Expected a lowercase full commit SHA")

    def git(*arguments):
        return subprocess.run(["git", "-C", str(checkout), *arguments], capture_output=True,
                              check=True).stdout

    lineage = git("rev-list", "--first-parent", "refs/remotes/origin/main").decode().splitlines()
    if BASE not in lineage or anchor not in lineage[:lineage.index(BASE)]:
        raise ValueError("Anchor outside the eligible official first-parent history")
    # Walk snapshots in chronological order, including merge commits.
    for commit in reversed(lineage[:lineage.index(BASE)]):
        if git("ls-tree", commit, "--", TASK):
            data = git("show", commit+":"+TASK)
            if hashlib.sha256(data).hexdigest() == TASK_HASH:
                if commit != anchor:
                    raise ValueError("Another commit is the first exact registration")
                return
    raise ValueError("No exact registration in the fetched official history")


def settings_for(anchor):
    word = hashlib.sha256(b"OPH-M10.4-gateway-v1:" + anchor.encode("ascii")).digest()
    return dict(home=word[0] % 20, start=word[1] % 12,
                direction=(-1 if word[2] & 1 else 1), exposure=15+word[3] % 24)


def laplacian(size, edges):
    incidence = np.zeros((len(edges), size))
    rows = np.arange(len(edges))
    incidence[rows, edges[:, 0]] = 1.
    incidence[rows, edges[:, 1]] = -1.
    return incidence.T @ incidence


def source_weights(rejected, start, direction):
    """Launch reach aggregated over four complete sweeps, including failures."""
    cycle = (start + direction*np.arange(12)) % 12
    prefixes = np.cumprod(np.r_[1., rejected[cycle]])
    cycle_failure = prefixes[-1]
    weights = np.empty(12)
    weights[cycle] = prefixes[:-1] * sum(cycle_failure**k for k in range(4))
    return weights, cycle_failure**4


def probabilities(board, home, start, direction, exposure):
    if (type(home) is not int or not 0 <= home < 20 or type(start) is not int
            or not 0 <= start < 12 or type(direction) is not int or direction not in (-1, 1)
            or type(exposure) is not int or not 0 <= exposure <= 38):
        raise ValueError("Settings outside the registered board or replay time range")
    board_bytes = board.read_bytes()
    if hashlib.sha256(board_bytes).hexdigest() != BOARD_HASH:
        raise ValueError("Incorrect scored board")
    with np.load(io.BytesIO(board_bytes), allow_pickle=False) as bundle:
        edges = bundle["edges"]
    inside = np.arange(12*home, 12*home+12)
    outside = np.ones(240, dtype=bool)
    outside[inside] = False
    internal = edges[(edges[:, 0] // 12 == home) & (edges[:, 1] // 12 == home)]
    boundary = edges[edges[:, 0] // 12 != edges[:, 1] // 12].ravel()
    gateways = boundary[boundary // 12 == home] % 12
    if len(gateways) != 3 or len(set(gateways)) != 3 or len(internal) != 30:
        raise ValueError("Invalid source geometry")

    # On a sharp one-load launch all internal terms beyond L/2 contribute
    # a common phase. A clear signed gateway query zeros just its own row.
    energies, modes = eigh(laplacian(12, internal % 12)/2)
    step = (modes*np.exp(-1j*energies)) @ modes.T
    columns = np.eye(12, dtype=complex)
    rejected = np.zeros(12)
    for _ in range(1024):
        columns = step @ columns
        rejected += (abs(columns[gateways])**2).sum(axis=0)
        columns[gateways] = 0.
    if np.max(abs((abs(columns)**2).sum(axis=0)+rejected-1)) > 2e-10:
        raise ArithmeticError("Incomplete sharp attempt")
    weights, exhausted = source_weights(rejected, start, direction)
    release = columns*np.sqrt(weights)
    diagonal = (abs(release)**2).sum(axis=1)

    # The complete free board also differs from L/2 by a scalar. Each
    # column is the free evolution of a sharp home port on that full board.
    full = laplacian(240, edges)
    energies, modes = eigh(full/2)
    flight = (modes*np.exp(-1j*exposure*energies)) @ modes[inside].T

    # Each seam swaps at rate 1/2. Uniformization at rate 3 gives a
    # stochastic matrix, including self-loops at ports with degree five.
    # For t <= 38 the Poisson tail after term 512 is at most
    # exp(-114) * (e*114/513)**513 by the Chernoff bound.
    chain = csr_matrix(np.eye(240)-full/6)
    if chain.data.min() < 0 or np.max(abs(np.asarray(chain.sum(axis=0))-1)) > 2e-15:
        raise ArithmeticError("Invalid stochastic transition")
    state = np.zeros(240)
    state[inside] = diagonal
    weight = math.exp(-3*exposure)
    heat = weight*state
    total = weight
    for n in range(1, 513):
        state = chain @ state
        weight *= (3*exposure)/n
        heat += weight*state
        total += weight
    if abs(total-1) > 1e-14:
        raise ArithmeticError("Poisson summation error")

    bins = {}
    for name, spatial, failure in (
        ("coherent", (abs(flight @ release)**2).sum(axis=1), exhausted),
        ("dephased", (abs(flight)**2) @ diagonal, exhausted),
        ("unfiltered_uniform", (abs(flight)**2).mean(axis=1), 0.),
        ("stochastic_swap", heat, exhausted),
    ):
        vector = np.r_[spatial[inside], spatial[outside].sum(), failure]
        if not np.isfinite(vector).all() or vector.min() < 0 or abs(vector.sum()-1) > 2e-10:
            raise ArithmeticError("Invalid complete receiver distribution for "+name)
        bins[name] = vector.tolist()
    return bins


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("anchor")
    cli.add_argument("--openwave-checkout", type=Path, required=True,
                     help="trusted official OpenWave checkout with origin/main fetched")
    cli.add_argument("--board", type=Path, default=BOARD)
    cli.add_argument("--output", type=Path, required=True, help="new independent receipt")
    args = cli.parse_args()
    first_registration(args.openwave_checkout, args.anchor)
    settings = settings_for(args.anchor)
    receipt = dict(anchor=args.anchor, settings=settings)
    with args.output.open("x", encoding="utf-8", newline="\n") as destination:
        try:
            receipt["bins"] = probabilities(args.board, **settings)
        except (ArithmeticError, ValueError, OSError) as error:
            receipt["error"] = str(error)
            raise
        finally:
            json.dump(receipt, destination, indent=2, allow_nan=False)
            destination.write("\n")


if __name__ == "__main__":
    main()
