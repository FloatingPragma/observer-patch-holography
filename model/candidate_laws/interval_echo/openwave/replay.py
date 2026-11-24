"""Replay the independent Interval Echo experiments using their deposited inputs."""

import os
import subprocess
import sys


EXPERIMENTS = (
    "row_04",
    "row_04_compact.experiment",
    "row_04_compact.anchored_vacancy",
    "row_04_compact.material_reference",
    "row_04_compact.localized_release",
    "row_01_charge.measure",
    "row_02_dressed.measure",
    "row_03_clock.measure",
    "row_04_identity.experiment",
    "row_04_conjugates.experiment",
    "row_04_large.experiment",
    "row_05_dipole.experiment",
    "row_05_dipole.exterior",
    "row_05_dipole.vector",
    "row_05_dipole.calibration",
    "row_05_dipole.transport",
    "row_05_dipole.orbital",
    "row_05_counter.experiment",
    "row_05_counter.dispersive",
    "row_05_counter.beat",
    "row_08_annihilation.experiment",
    "row_08_annihilation.converter",
    "row_08_annihilation.calorimeter",
    "row_08_annihilation.walls",
    "row_08_annihilation.wall_flight",
    "row_08_annihilation.bare_receivers",
    "row_08_annihilation.energy_budget",
    "row_11_spectrum.measure",
    "row_11_spectrum.bounds",
    "row_11_spectrum.transport",
    "row_12_charge.experiment",
    "row_12_charge.exterior",
    "row_12_flux.measure",
    "row_12_flux.dressing",
    "row_12_flux.cycle",
    "row_13_baryon.measure",
    "row_13_baryon.frames",
    "row_13_walls.energy",
    "row_20_electric.experiment",
    "row_20_electric.neutral",
    "row_20_electric.schedule",
    "row_20_electric.flight",
    "row_21_boost.measure",
    "row_22_confinement.measure",
    "row_23_coupling.experiment",
    "row_23_coupling.cutoff",
    "row_27_metric.experiment",
    "row_27_metric.exterior",
    "row_28_acceleration.experiment",
    "row_28_acceleration.signals",
    "row_28_acceleration.homogeneity",
    "row_28_acceleration.capacity",
    "row_28_acceleration.radar",
    "row_28_acceleration.cut_signal",
    "row_28_acceleration.cut_code",
    "row_28_acceleration.paired_probe",
    "row_28_acceleration.volume",
    "row_29_maxwell.experiment",
    "row_29_maxwell.survey",
    "row_29_circulation.experiment",
    "row_29_circulation.precision",
    "row_29_circulation.local_flight",
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
    "row_30_route.receiver",
    "row_30_formed.experiment",
    "row_30_repair.experiment",
    "row_30_local_routes.experiment",
    "row_30_local_routes.junction",
    "row_30_local_routes.filter_read",
    "row_30_optical.experiment",
    "row_30_optical.complement",
    "row_30_collective.formation",
    "row_30_collective.coupling_meter",
    "row_30_collective.experiment",
    "row_30_collective.holdout",
    "row_30_collective.uncertainty",
    "row_30_collective.packets",
    "row_30_collective.rest.experiment",
    "row_30_collective.decode_example",
    "row_30_collective.propagating_control",
    "row_31_orbitals.measure",
    "row_31_orbitals.mediator",
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
