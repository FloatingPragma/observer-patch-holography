"""Complete terminal channel on the intact federation."""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate as echo
from model.experiments.interval_echo_space.instrument import local_tables

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SPEC = json.loads((HERE / "protocol.json").read_text())
TOL = SPEC["numeric_tolerance"]

# Whole-file changes reviewed at OPH 0ebb9ba9d6280fec258a7d65e2966c13352dfa8d.
REPLAY_UPDATES = {
    "model/candidate_laws/interval_echo/SPEC.md": "46781970157edaec4b30e55bb673c90102ffcf937aeb6974927e58b08e74dff4",
    "model/candidate_laws/interval_echo/candidate.py": "a318f231c29753da6d0bfaf54ceb83f4410ede7b054b1edc1360620692876944",
    "model/candidate_laws/interval_echo/EVALUATION.json": "4a31cbd583995fa9af119b4e4a4f5545c3656f5f26dfd25de49d6a9440509331",
}

# Additional metadata changes reviewed at main 70e58e24f3be9295926164d38626d58acb87b719.
REPLAY_70E58E24 = {
    "model/candidate_laws/interval_echo/SPEC.md": "c9954a39eaef1de8788e5c4a0e1de1cbd81e8d4dd162566c6aac59467d1af506",
    "model/candidate_laws/interval_echo/EVALUATION.json": "7d95fd5f067232b9f9f3aba664148865ee6773b90ffa52e92f882b06bdae80d9",
}

# The regional conservation, coherent-sweep and coupling evidence changes metadata only.
REPLAY_REGIONAL_CONSENSUS = {
    "model/candidate_laws/interval_echo/candidate.py": "c55b45990fedd5b79ebf5d8e206334b7e90c519e72811b5687149c7f58d9a4a4",
    "model/candidate_laws/interval_echo/EVALUATION.json": "768e20a009f611483bd51a027cc47593a13084e276fb35a561a9cde6953d94cc",
}

# Candidate and evaluation metadata now also record the separate comparison standard.
REPLAY_MODULAR_STANDARD = {
    "model/candidate_laws/interval_echo/candidate.py": "6afe9f5a24481beb0331a7f6930fce6d22c72cb6147d4491aa68b49edeb6107d",
    "model/candidate_laws/interval_echo/SPEC.md": "846dc2cbe46a714dd87d735d3a9d4d21e11bd022da8a6d106b3b351ce0c5c00a",
    "model/candidate_laws/interval_echo/EVALUATION.json": "661a1fbd55b9766acad1b657da24578b3a2ab35cad2016aad129a5d6d165cce6"
}


# Reviewed premise and evidence metadata for the regional gravitational and dark-sector instruments.
REPLAY_GRAVITY_CAP = {
    "model/candidate_laws/interval_echo/candidate.py": "c69c10a59f40c5e71f20854d0fec80faa5a3be16cf4f9232adae884dadb0b59e",
    "model/candidate_laws/interval_echo/EVALUATION.json": "c2817432c2218a97220d56d107f8c81748178f5c5eb8788de665ec291e46920f"
}


# Nested cap, boundary-volume and clock comparisons update the premise/evidence metadata.
REPLAY_NESTED_CAP = {
    "model/candidate_laws/interval_echo/candidate.py": "c666a31015ea2203ac919d55d98ddd5e47fe78076b3d87dd94bed2dd9990a503",
    "model/candidate_laws/interval_echo/EVALUATION.json": "6ad745ce68fa19b1c2fbdc1ea31948403e1d33bc13b8b761440141a031b2661a",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sources():
    frozen = SPEC["sources"]
    actual = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in frozen}
    reviewed = frozen | REPLAY_UPDATES
    radial = reviewed | REPLAY_70E58E24
    regional = radial | REPLAY_REGIONAL_CONSENSUS
    standard = regional | REPLAY_MODULAR_STANDARD
    gravity = standard | REPLAY_GRAVITY_CAP
    if actual in (frozen, reviewed, radial, regional, standard, gravity, gravity | REPLAY_NESTED_CAP):
        return
    mismatches = [name for name, expected in reviewed.items() if actual.get(name) != expected]
    raise ValueError(f"unverified replay source: {', '.join(mismatches)}")


def encode(board, station):
    contact = next(edge for edge in board.external if edge[0] < 12)
    _, a = contact
    b = min(p for p in board.neighbours[a] if p // 12 == a // 12 and board.degree[p] == 5)
    launch = min(set(board.neighbours[a]) & set(board.neighbours[b]))
    require(station["cells"] == board.cells and list(contact) == station["contact"]
            and [a, b] == station["ports"] and launch == station["launch"]
            and launch // 12 == a // 12 != 0, "station map")
    law = echo.native_law
    # Address one seam with common start |1,0> and the same contact area.
    pulse = expm(-1j * np.pi * law.seam_term(echo.KAPPA))[:, law.SRC.index((1, 0))]
    columns = np.zeros((board.n, 2), complex)
    columns[launch, :] = pulse[law.SRC.index((1, 0))]
    columns[[a, b], [0, 1]] = pulse[law.SRC.index((0, 1))]
    require(abs(abs(columns[a, 0]) - 1) < TOL and
            np.linalg.norm(np.delete(pulse, law.SRC.index((0, 1)))) < TOL, "sender exchange")
    return columns


def validate_law(p):
    require(p.shape == (2, 13) and p.dtype.kind == "f" and np.isfinite(p).all()
            and p.min() >= 0 and np.max(abs(p.sum(axis=1) - 1)) <= TOL,
            "incomplete local law")


def readout(h, launch, s):
    state = expm_multiply((-1j * s) * h, launch) if s else launch
    mass = abs(state)**2
    # Audit the raw complement and total before the inherited instrument's guard.
    residual = mass.sum(axis=0) - 1
    outside = mass[12:].sum(axis=0)
    require(np.max(abs(residual)) <= TOL and (1 - mass[:12].sum(axis=0)).min() >= 0,
            "raw normalization or complement")
    p = local_tables(state, np.arange(12)[None, :])[0].T
    validate_law(p)
    require(np.max(abs(outside - p[:, -1])) <= TOL, "outside probability")
    return p, residual.tolist()


def distance(p):
    validate_law(p)
    return float(abs(p[0] - p[1]).sum() / 2)
