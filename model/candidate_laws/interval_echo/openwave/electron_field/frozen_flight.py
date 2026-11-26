"""Matched no-ring control for the archived written-field source and receiver."""

import json
from pathlib import Path
import sys

import numpy as np


def run(experiments):
    # Use the archive's actual geometry and instruments after replay.py authenticates them.
    sys.path.insert(0, str(experiments / "interval_echo_vacuum_audit" / "repair"))
    from flight import Instrument, Sector, Sphere, federation
    from oracle import local_matrices
    from release import instruments

    instrument = Instrument()
    # Integer columns for the two coherent written heads and an equal-load (1,1) seam.
    heads = np.rint(np.sqrt(2) * instrument.basis).astype(int)
    loaded = np.eye(9, dtype=int)[:, [instrument.pairs.index((1, 1))]]
    columns = np.kron(heads, loaded)
    np.testing.assert_array_equal(instrument.h8 @ columns, 64 * columns)
    # A signed vacancy on a live seam *can* hop and reverse the head on the very
    # same native law. The physical matrix element is -4: h8=8H and ||head||²=2.
    load_basis = np.eye(9, dtype=int)
    before = np.kron(heads, load_basis[:, [instrument.pairs.index((0, 1))]])
    after = np.kron(heads, load_basis[:, [instrument.pairs.index((1, 0))]])
    np.testing.assert_array_equal(after.T @ instrument.h8 @ before, [[0, 0], [-64, 0]])
    source = local_matrices(instrument)
    detector = instrument.detector(1)
    np.testing.assert_allclose(detector, np.diag(np.diag(detector)), atol=1e-13, rtol=0)
    arrival = np.diag(detector).real
    deposit = json.loads((experiments / "interval_echo_electron_field" / "results.json").read_text())

    for level, record in zip((1, 2), deposit["coherent_ground_release"], strict=True):
        sphere = Sphere(federation(level), "depth-first within shell")
        sector = Sector(sphere)
        assert sphere.board.cells == record["cells"] and len(sector.words) == record["component_states"]
        source_heads = np.logical_xor(sphere.head[sphere.source], sector.flipped[:, 0]).astype(int)
        remote_heads = np.logical_xor([sphere.head[e] for e in sphere.detectors], sector.flipped).astype(int)
        # The original parking sequence leaves all non-unit loads in isolated reservoirs.
        # The signed source and receiver have disjoint seams for every enrolled station.
        for seam, arms, station in zip(sphere.detectors, sphere.arms, record["records"], strict=True):
            assert set(sphere.source).isdisjoint(seam)
            assert tuple(station["seam"]) == seam
            instruments(instrument, sphere.source, seam, arms)

        response = arrival[remote_heads]
        sham = response.mean(axis=0)
        np.testing.assert_allclose(sham, [s["sham_positive_probe_crosses"] for s in record["records"]],
                                   atol=2e-13, rtol=0)
        residuals = []
        for branches in source.values():
            # The click and miss prepare different source heads, but both outcomes are kept.
            # With scalar native flight, the remote diagonal effect can be summed on sharp
            # words even for the original coherent, equal-amplitude preparation.
            weights = sum(np.sum(np.abs(kraus[:, source_heads])**2, axis=0) for kraus in branches)
            signed = weights @ response / len(sector.words)
            residuals.append(float(np.max(np.abs(signed - sham))))
        assert max(residuals) < 2e-13
        far = record["farthest_representative"]
        assert abs(far["signed_odd_peak"]) > record["threshold"]
        print(f"{sphere.board.cells} cells, {len(sham)} detectors, both signs: "
              f"maximum native difference from sham {max(residuals):.3g}. "
              f"At {far['cell_steps']} steps the native probability is {sham[far['face']]:.9f}; "
              f"archived ring-on odd shift {far['signed_odd_peak']:.9g}.")


if __name__ == "__main__":
    run(Path(sys.argv[1]).resolve())
