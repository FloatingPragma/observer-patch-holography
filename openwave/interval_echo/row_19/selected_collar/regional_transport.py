"""Test the complete integrated net energy transport on the ordinary control flights."""

import argparse
from hashlib import sha256
import json

from .regional_source_calibration import HERE, calibrate
from .regional_rank import rank


def run():
    path = HERE / "distributed_source.json"
    record = json.loads(path.read_bytes())
    result = record["result"]
    ports = list(map(str, record["ports"]))
    sites = ports + ["complement"]

    def terminal(z):
        return [
            z[n][field] for field in ("energy", "connected", "marginal") for n in sites
        ]

    source = {}
    for h, z in result["ordinary_regional_sources"].items():
        _, sender, _ = h.split("/")
        source[h] = terminal(z) + [
            result["ordinary_total_energy"][h] if p == sender else [0.0, 0.0]
            for p in ports[:-1]
        ]
    shape = result["groups"]["shape"]
    material = result["groups"]["material"]
    held = terminal(result["held_sources"]) + [[0.0, 0.0] for _ in ports[:-1]]
    alternative = {"actual": (held, result["sources"]["connected"]["target"])}
    test = calibrate(source, result["columns"], alternative, shape + material)
    ranks = {}
    if test["status"] == "verified inverse":
        ranks = {
            name: rank(labels, test["ordinary"])
            for name, labels in (
                ("shape", shape),
                ("material", material),
                ("both", shape + material),
            )
        }
    return {
        "retained_launchers": ports[:-1],
        "dependent_launcher": ports[-1],
        "calibration": test,
        "residual_rank": ranks,
        "inputs": {
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (path, HERE / "REGIONAL_TRANSPORT.md", HERE / "regional_rank.json")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "regional_transport.json"
    if args.write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    print(
        result["calibration"]["status"],
        {g: v["rank_lower"] for g, v in result["residual_rank"].items()},
    )
