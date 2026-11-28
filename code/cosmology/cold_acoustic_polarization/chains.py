"""Summarize physical cold-to-baryon ratios in the Planck polarization-only chains."""

import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
URL = "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/COM_CosmoParams_base-plikHM_R3.01.zip"
CHANNELS = ("TE", "EE")


def extract(path):
    arrays = {}
    provenance = {"archive": URL, "channels": {}}
    with zipfile.ZipFile(path) as archive:
        for ch in CHANNELS:
            name = f"base_plikHM_{ch}_lowE"
            prefix = f"base/plikHM_{ch}_lowE/"
            hashes = {}

            def read(member):
                data = archive.read(prefix + member)
                hashes[prefix + member] = hashlib.sha256(data).hexdigest()
                return data

            properties = read(name + ".properties.ini").decode().splitlines()
            if "burn_removed=T" not in properties:
                raise ValueError("Primary burn-in convention changed")
            names = [
                s.split()[0] for s in read(name + ".paramnames").decode().splitlines()
            ]
            # CosmoMC column 0 is multiplicity; column 1 is minus log posterior.
            cols = (0, 2 + names.index("omegabh2"), 2 + names.index("omegach2"))
            ns = []
            for chain in range(1, 5):
                data = read(f"{name}_{chain}.txt")
                rows = np.loadtxt(io.BytesIO(data), usecols=cols)
                ns.append(len(rows))
                arrays[f"{ch}_{chain}"] = rows
            stats = read("dist/" + name + ".margestats").decode().splitlines()
            convergence = read("dist/" + name + ".converge").decode().splitlines()
            settings = dict(
                line.split(" = ", 1)
                for line in read(name + ".inputparams").decode().splitlines()
                if " = " in line
            )
            quoted = {}
            for par in ("omegabh2", "omegach2"):
                first = next(s.split() for s in stats if s.startswith(par + " "))
                autocorr = next(
                    s.split() for s in convergence if s.startswith(par + " ")
                )
                quoted[par] = {
                    "mean": float(first[1]),
                    "sd": float(first[2]),
                    "reported_effective_sample_size": int(autocorr[3]),
                }
            provenance["channels"][ch] = {
                "chain_lengths": ns,
                "sampled_density_prior_ranges": {
                    par: [float(x) for x in settings[f"param[{par}]"].split()[1:3]]
                    for par in ("omegabh2", "omegach2")
                },
                "H0_prior_range": [int(settings["H0_min"]), int(settings["H0_max"])],
                "sha256": hashes,
                "burn_removed": True,
                "source_marginals": quoted,
            }
    np.savez_compressed(HERE / "chain_density_columns.npz", **arrays)
    (HERE / "chain_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def weighted_samples(rows):
    weights, baryon, cold = rows.T
    if not np.all((weights > 0) & (weights == np.floor(weights))):
        raise ValueError("Expected integer positive chain multiplicities")
    return np.repeat(cold / baryon, weights.astype(int))


def measure():
    out = {}
    with np.load(HERE / "chain_density_columns.npz") as archive:
        for ch in CHANNELS:
            chains = [archive[f"{ch}_{i}"] for i in range(1, 5)]
            rows = np.vstack(chains)
            r = weighted_samples(rows)
            means = np.average(rows[:, 1:], weights=rows[:, 0], axis=0)
            variances = np.average(
                (rows[:, 1:] - means) ** 2, weights=rows[:, 0], axis=0
            )
            out[ch] = {
                "distinct_rows": len(rows),
                "total_weight": len(r),
                "density_means": means.tolist(),
                "density_sds": np.sqrt(variances).tolist(),
                "ratio_quantiles_025_16_50_84_975": np.quantile(
                    r, [0.025, 0.16, 0.5, 0.84, 0.975]
                ).tolist(),
                "min_max_sampled_ratio": [float(r.min()), float(r.max())],
                "fraction_sampled_at_or_below": {
                    str(q): float(np.mean(r <= q)) for q in (3 / 7, 1, 2)
                },
                "chain_medians": [
                    float(np.median(weighted_samples(v))) for v in chains
                ],
            }
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archive", type=Path, help="Rebuild packaged columns from the official ZIP"
    )
    args = parser.parse_args()
    if args.archive is not None:
        extract(args.archive)
    (HERE / "chain_results.json").write_text(
        json.dumps(measure(), indent=2, allow_nan=False) + "\n"
    )
