"""Ordinary calibration of independent link meters against held clock and material records."""

import argparse
from hashlib import sha256
import json

from . import distributed_source as dist
from .regional_source_calibration import calibrate
from .regional_rank import rank

HERE, lab, iv = dist.HERE, dist.lab, dist.iv


def sources(ports, local, total, origin, flow, spec):
    fields, transport, names = spec
    sites = ports + ["complement"]
    values = [local[site][key] for key in fields for site in sites]
    if transport:
        values.extend(total if site == origin else [0.0, 0.0] for site in ports[:-1])
    return values + [flow[n] for n in names]


def menus(currents):
    def names(meter, external):
        return [
            n
            for n in currents
            if n.startswith(meter + "/") and (external or not n.endswith("/outside"))
        ]

    energy = ["energy"]
    all_fields = ["energy", "connected", "marginal"]
    return {
        "energy_current": (energy, False, names("current", True)),
        "energy_current_work": (
            energy,
            False,
            names("current", True) + names("work", True),
        ),
        "net_source_internal_current": (all_fields, True, names("current", False)),
        "net_source_internal_current_work": (
            all_fields,
            True,
            names("current", False) + names("work", False),
        ),
    }


def run(extra=None):
    iv.dps = 35
    records = dist.read("distributed_source")
    flow = dist.read("directed_service")["columns"]
    prior = records["result"]
    ports = list(map(str, records["ports"]))
    shape = prior["groups"]["shape"]
    material = prior["groups"]["material"]
    target = prior["sources"]["connected"]["target"]
    programs = menus(next(iter(flow.values()))["ordinary"])
    mixed = {
        n: lab.six.bounds(
            sum(iv.mpf(z["overlap"][n]) for z in flow.values()) / len(flow)
        )
        for n in next(iter(flow.values()))["ordinary"]
    }
    rows = {}
    for key, spec in programs.items():
        source = {
            h: sources(
                ports,
                z,
                prior["ordinary_total_energy"][h],
                h.split("/")[1],
                flow[h]["ordinary"],
                spec,
            )
            for h, z in prior["ordinary_regional_sources"].items()
        }
        alternatives = {
            "programmed_overlap": (
                sources(ports, prior["held_sources"], None, None, mixed, spec),
                target,
            )
        }
        if extra is not None:
            for h, z in extra.items():
                alt = sources(
                    ports, z["regional"], z["total"], h.split("/")[1], z["flow"], spec
                )
                alternatives[h] = (alt, z["readings"])
        fit = calibrate(source, prior["columns"], alternatives, shape + material)
        ranks = {}
        if fit["status"] == "verified inverse" and fit["ordinary"]:
            ranks = {
                g: rank(ns, fit["ordinary"])
                for g, ns in (
                    ("shape", shape),
                    ("material", material),
                    ("both", shape + material),
                )
            }
        rows[key] = {
            "source_size": len(next(iter(source.values()))),
            "specification": list(spec),
            "calibration": fit,
            "residual_rank": ranks,
        }
    return {
        "fits": rows,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "distributed_source.json",
                "directed_service.json",
                "DIRECTED_SERVICE_PROTOCOL.md",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "current_geometry.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    for k, z in r["fits"].items():
        c = z["calibration"]
        print(
            k,
            c["status"],
            sum(bool(x["excluded"]) for x in c.get("ordinary", {}).values()),
            {h: len(x["excluded"]) for h, x in c.get("overlap", {}).items()},
            {g: v["rank_lower"] for g, v in z["residual_rank"].items()},
        )
