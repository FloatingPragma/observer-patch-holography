"""Local detector cumulants from the programmed scalar and a freely released pair."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_04 import __main__ as release
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
PAIRS = tuple(combinations(range(4), 2))


def moments(weights, marks):
    weights = np.asarray(weights, dtype=float)
    marks = np.asarray(marks, dtype=float)
    if (
        weights.ndim != 1
        or marks.shape != (len(weights), 4)
        or not np.isfinite(weights).all()
        or np.any(weights < 0)
        or weights.sum() <= 0
        or not np.isin(marks, (0, 1)).all()
    ):
        raise ValueError(
            "Expected a complete nonnegative distribution on four binary marks"
        )
    p = weights / weights.sum()
    mean = p @ marks
    centered = marks - mean
    cov = (centered.T * p) @ centered
    wick = [cov[0, j] * cov[k, l] for j, k, l in ((1, 2, 3), (2, 1, 3), (3, 1, 2))]
    fourth = p @ np.prod(centered, axis=1)
    return dict(
        means=mean.tolist(),
        pairs=[cov[i, j] for i, j in PAIRS],
        centered_fourth=fourth,
        wick_pairings=wick,
        connected_fourth=fourth - sum(wick),
    )


def run():
    blind_path = HERE / "blind_chart.json"
    release_path = release.HERE / "results.json"
    counts = np.asarray(json.loads(blind_path.read_text())["calibration_counts"])
    if counts.shape != (8, 8) or not np.all(
        counts.sum(axis=1) == counts.sum(axis=1)[0]
    ):
        raise ValueError("The pooled scalar requires eight equally sampled sources")
    position = moments(counts.sum(axis=0), np.eye(8)[:, :4])
    data = json.loads(release_path.read_text())
    panels = []
    for board in data["boards"]:
        sector = PairSector(federation(board["level"]))
        addresses = board["home_internal_addresses"]
        indices = [sector.index(*edge) for edge in addresses]
        marks = [[port in edge for port in range(4)] for edge in addresses] + [
            [False] * 4
        ]
        source = np.zeros(sector.P)
        source[sector.index(*board["original_seam"])] = 1
        for arm, coupling in (("sharp", candidate.KAPPA), ("exchange_only", 0.0)):
            row = next(r for r in board["arms"][arm] if r["time"] == 4)
            terminal_counts = row["home_seam_or_quiet_counts"]
            if sum(terminal_counts) != data["shots_per_setting"]:
                raise ValueError("Incomplete pair detector records")
            (psi,) = release.evolve(sector.hamiltonian(coupling), coupling, source, [4])
            probabilities = np.abs(psi) ** 2
            if abs(probabilities.sum() - 1) > 1e-10:
                raise ArithmeticError("Full sector evolution lost probability")
            home = probabilities[indices]
            law = np.r_[home, 1 - home.sum()]
            panels.append(
                dict(
                    cells=board["cells"],
                    arm=arm,
                    exposure=4,
                    shots=sum(terminal_counts),
                    quiet=terminal_counts[-1],
                    recorded=moments(terminal_counts, marks),
                    population=moments(law, marks),
                )
            )
    return dict(
        source_sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (blind_path, release_path, HERE / "correlation_protocol.json")
        },
        ports=[0, 1, 2, 3],
        pair_order=PAIRS,
        scalar_shots=int(counts.sum()),
        scalar_position=position,
        pair_detector=panels,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
