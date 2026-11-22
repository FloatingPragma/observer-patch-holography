"""Read the deposited joint histograms without drawing a second acquisition."""
import argparse
from hashlib import sha256
import json

import numpy as np

from model.experiments.seam_program.gate import serialize
from . import dwell
from .dwell_checks import extended
from .dwell_source import Station, require

# Exact raw files in the primary and control acquisition commits.
ARCHIVES = {
    "dwell_counts.npz": "5f4289eb5f9cc3d54ab289bca7e529cede5c616e094c1955b79f5653ca4aa29a",
    "dwell_controls_counts.npz": "d31ec93aa81cdb350f66ea1909db20908d4055f8d3a9d4555ad56fb16f0ca7e3",
}


def records(path):
    require(sha256(path.read_bytes()).hexdigest() == ARCHIVES[path.name],
            f"Changed acquisition archive: {path.name}")
    return np.load(path, allow_pickle=False)


def control_plan(directory):
    """Check the registered control protocol and both propagation implementations."""
    manifest = (directory / "dwell_controls_sources.json").read_bytes().replace(b"\r\n", b"\n")
    require(sha256(manifest).hexdigest() ==
            "b590a65457a29a8e8bc41dae943e7728a1f4fc60b129b7c640cdb6a828d3bbf1",
            "Changed control source manifest")
    for name, digest in json.loads(manifest).items():
        source = (directory / name).read_bytes().replace(b"\r\n", b"\n")
        require(sha256(source).hexdigest() == digest, f"Changed control source: {name}")
    return json.loads((directory / "dwell_controls_protocol.json").read_text(encoding="utf-8"))


def verify():
    controls = control_plan(dwell.HERE)
    base = dwell.plan()
    station = Station(base["flight_level"], base["tag"])
    shots = base["shots"]
    appointments = base["last_tick"] + 1
    primary = {f"{arm}_{strength}": (arm, strength)
               for strength in base["strengths"] for arm in base["arms"]}
    interventions = {arm: (arm, controls["coupling"]) for arm in controls["controls"]}
    for stem, settings in (("dwell", primary), ("dwell_controls", interventions)):
        retained = json.loads((dwell.HERE / f"{stem}_result.json").read_text())
        panels = appointments * len(settings)
        error = np.sqrt(2 * np.log(8 * panels / base["alpha"]) / shots) + 1e-8
        with records(dwell.HERE / f"{stem}_counts.npz") as data:
            require(set(data.files) == set(settings), "Changed acquisition settings")
            for name, (arm, strength) in settings.items():
                counts = data[name]
                distortion = None
                if arm in ("always_ticking", "binding_only"):
                    distortion = 0.
                elif arm == "native":
                    bounds = [station.isolation(strength, parity) for parity in (0, 1)]
                    if all(bounds):
                        distortion = sum(np.arcsin((1 - b["weight"]) / b["weight"]) for b in bounds)
                observed = dwell.read(counts, shots, error, distortion)
                if stem == "dwell":
                    published = retained["arms"][name]["observed"]
                else:
                    terminal = counts[:, -1].sum(axis=-1, dtype=np.int64) / shots
                    real, imag = terminal[:, 0] - terminal[:, 1]
                    observed["last_principal_phase"] = float(np.angle(complex(real, imag)))
                    published = retained[name]
                require(serialize(observed) == published, f"Changed observer result: {name}")
    return base


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extended", action="store_true",
                        help="also check the recorded counts against bounded full-state propagation")
    args = parser.parse_args()
    protocol = verify()
    if args.extended:
        extended(protocol)
    print("Reproduced all eight observer results from the two deposited acquisitions.")
