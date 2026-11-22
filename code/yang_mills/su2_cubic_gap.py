"""Finite-volume min-max bound for pure SU(2) Hamiltonian lattice gauge theory."""

import argparse
from fractions import Fraction
import json
from pathlib import Path


def cubic_bound(side, magnetic):
    """Use electric Casimir units. The spatial graph is a periodic simple cubic lattice."""
    if type(side) is not int or side < 3:
        raise ValueError(
            "A simple periodic cubic graph requires at least three sites per axis"
        )
    if not isinstance(magnetic, Fraction) or magnetic < 0:
        raise ValueError("The plaquette coefficient must be a nonnegative rational")
    plaquettes = 3 * side**3
    # Every nonempty spin network contains a cycle. A fundamental loop attains this minimum.
    electric_gap = Fraction(3, 4) * min(4, side)
    vacuum_upper = magnetic * plaquettes
    gap = max(Fraction(0), electric_gap - vacuum_upper)
    return dict(
        side=side,
        links=3 * side**3,
        plaquettes=plaquettes,
        electric_gap=str(electric_gap),
        vacuum_energy_upper=str(vacuum_upper),
        spectral_gap_lower=str(gap),
        positive_bound=gap > 0,
    )


def run():
    coupling = Fraction(1, 512)
    sides = (4, 5, 6, 8, 9)
    return dict(
        group="SU(2)",
        spatial_dimension=3,
        time="continuous",
        energy_unit="coefficient of the single-link electric Casimir",
        magnetic_coefficient=str(coupling),
        representation_cutoff=None,
        electric_only=[cubic_bound(side, Fraction(0)) for side in sides],
        electric_and_magnetic=[cubic_bound(side, coupling) for side in sides],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(json.dumps(run(), indent=2) + "\n")
