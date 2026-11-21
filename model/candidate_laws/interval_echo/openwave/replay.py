"""Replay the independent Interval Echo experiments using their deposited inputs."""

import os
import subprocess
import sys


EXPERIMENTS = (
    "row_04",
    "row_04_compact.experiment",
    "row_04_identity.experiment",
    "row_04_conjugates.experiment",
    "row_04_large.experiment",
    "row_30.experiment",
    "row_30.patch_time",
    "row_30_free.experiment",
    "row_30_branch.experiment",
    "row_30_branch.gateways",
    "row_30_branch.classes",
    "row_30_branch.refinement",
    "row_30_branch.precision",
    "row_30_memory.experiment",
    "row_30_arrows.experiment",
    "row_30_arrows.phase",
    "row_30_arrows.spectrum",
    "row_30_arrows.clock",
    "row_30_arrows.phase_spectrum",
    "row_30_selective.experiment",
    "row_30_selective.staggered",
    "row_30_autonomous.experiment",
    "row_30_autonomous.midpoint",
    "row_30_autonomous.reference_control",
    "row_30_autonomous.uncertainty",
    "row_30_selective.channel_review",
    "row_30_route.experiment",
    "row_30_route.blank",
)


def main():
    env = dict(
        os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1"
    )
    for name in EXPERIMENTS:
        module = f"model.candidate_laws.interval_echo.openwave.{name}"
        print(f"Verifying {name}", flush=True)
        subprocess.run(
            [sys.executable, "-B", "-W", "error", "-m", module], check=True, env=env
        )
    print(f"All {len(EXPERIMENTS)} experiments reproduced.")


if __name__ == "__main__":
    main()
