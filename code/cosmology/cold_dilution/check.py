"""Check the conserved-density shape against published binned dark-matter fits."""

import argparse
from decimal import Decimal
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import re
import tarfile
from urllib.request import urlopen

HERE = Path(__file__).parent
RECEIPT = HERE / "results.json"
SOURCES = {
    "1802.09541v2": (
        "varywletter.tex",
        "fe479a134df81a659b74720c2dce2aaa607981f068ae85fa22412ccc0f14cdce",
    ),
    "2004.09572v2": (
        "var-gdm-paper.tex",
        "4deab608381ef24d53dc90528471ea313894e89a08785497068c27ac0da7502e",
    ),
    "1601.05097v2": (
        "MCMC_I.tex",
        "e1d6ce72302d68f8c93fee6be4c31ae24dd5b194e7311dae9089c070097c0fc4",
    ),
}
DATASETS = ("PPS", "PPS+BAO", "PPS+HST")
MODELS = ("var-c", "var-wc", "var-w", "LCDM")


def read_sources(source_dir):
    texts = []
    for arxiv, (filename, expected) in SOURCES.items():
        if source_dir is None:
            with urlopen(f"https://arxiv.org/src/{arxiv}", timeout=40) as response:
                with tarfile.open(
                    fileobj=BytesIO(response.read()), mode="r:gz"
                ) as archive:
                    raw = archive.extractfile(filename).read()
        else:
            raw = (source_dir / filename).read_bytes()
        if sha256(raw).hexdigest() != expected:
            raise ValueError(f"Unexpected content in {filename}")
        texts.append(raw.decode("utf8"))
    return texts


def extract(kopp, ilic, thomas):
    number = r"[-+]?\d+(?:\.\d+)?"
    error_bar = rf"({number})\^\{{\+({number})\}}_\{{-({number})\}}"
    data = {}
    for i in range(8):
        (line,) = [
            s for s in kopp.splitlines() if re.match(r"\{\$w_" + str(i) + r"\s", s)
        ]
        cells = re.findall(error_bar, line)
        if len(cells) != len(DATASETS):
            raise ValueError("Missing marginalized interval")
        data[str(i)] = {
            key: dict(mean=mean, upper_error=upper, lower_error=lower)
            for key, (mean, upper, lower) in zip(DATASETS, cells, strict=True)
        }
    # Table IV's first row contains one group identified by font per model
    # in each data column. A match includes only each color's numeric entries.
    table = ilic.split(r"$-\ln(\Lcal_{\rm max}) $ &", 1)[1].split(r"\hline", 1)[0]
    likelihoods = {key: {} for key in ("PPS", "PPS+Lens+BAO")}
    for macro, model in zip(
        ("colorvarc", "colorvarwc", "colorvarw", "colorLCDM"), MODELS, strict=True
    ):
        values = re.findall(r"\\" + macro + r"\{\s*(\d+\.\d+)\s*\}", table)
        for dataset, value in zip(likelihoods, values, strict=True):
            likelihoods[dataset][model] = value
    # The results paragraph distinguishes this generalised-dark-matter
    # interval from the adjacent fit that fixes the perturbation parameters.
    (constant_w,) = re.findall(
        rf"\$({number})<w<({number})\$ at \$99\.7\\%\$ CL in the case of\s+\$\\Lambda\$-GDM",
        thomas,
    )
    (equality_fit,) = re.findall(
        r"equality \$\\omega\^\{\\rm eq\}_g\s*=\s*" + error_bar + r"\$ \(95\\% c\.l\.\)",
        kopp,
    )
    mean, upper, lower = equality_fit
    return {
        "kopp_table_I": data,
        "kopp_PPS_equality_density_95pct": dict(mean=mean, upper_error=upper, lower_error=lower),
        "ilic_table_IV_neg_log_max_likelihood": likelihoods,
        "thomas_PPS_lens_BAO_constant_w_99_7pct": list(constant_w),
    }


def evaluate(data):
    bins = {}
    for i, panels in data["kopp_table_I"].items():
        bins[i] = {}
        for dataset, cell in panels.items():
            center = Decimal(cell["mean"])
            lo = center - Decimal(cell["lower_error"])
            hi = center + Decimal(cell["upper_error"])
            bins[i][dataset] = {
                "marginal_95pct_w": [float(lo), float(hi)],
                "marginal_95pct_dilution_exponent": [
                    float(3 * (1 + lo)),
                    float(3 * (1 + hi)),
                ],
                "includes_cold_value": lo <= 0 <= hi,
            }
    comparisons = {}
    for dataset, values in data["ilic_table_IV_neg_log_max_likelihood"].items():
        nll = {k: Decimal(v) for k, v in values.items()}
        comparisons[dataset] = {}
        for label, restricted, free in (
            ("sound_and_viscosity_zero", "LCDM", "var-w"),
            ("sound_and_viscosity_vary_in_both", "var-c", "var-wc"),
        ):
            comparisons[dataset][label] = {
                "restricted": restricted,
                "free": free,
                "added_eos_parameters": 8,
                "twice_log_likelihood_gain": float(2 * (nll[restricted] - nll[free])),
            }
        # The larger models contain the restricted corners. These small
        # reversed improvements measure how the reported extrema respect nesting.
        comparisons[dataset]["adding_perturbation_parameters_at_fixed_w_zero"] = float(
            2 * (nll["LCDM"] - nll["var-c"])
        )
        comparisons[dataset]["adding_perturbation_parameters_at_free_w"] = float(
            2 * (nll["var-w"] - nll["var-wc"])
        )
    w_lo, w_hi = map(Decimal, data["thomas_PPS_lens_BAO_constant_w_99_7pct"])
    n_lo, n_hi = 3 * (1 + w_lo), 3 * (1 + w_hi)
    return {
        "marginal_intervals": bins,
        "profile_comparisons": comparisons,
        "constant_exponent_test": {
            "dilution_exponent_99_7pct": [float(n_lo), float(n_hi)],
            "density_ratio_for_doubling_99_7pct": [
                float(Decimal(2) ** (-n_hi)), float(Decimal(2) ** (-n_lo))
            ],
            "conserved_charge_prediction": {"exponent": 3, "density_ratio_for_doubling": 0.125},
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group()
    sources.add_argument(
        "--fetch",
        action="store_true",
        help="read the three versioned arXiv sources",
    )
    sources.add_argument(
        "--source-dir", type=Path, help="directory containing the three original TeX members"
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="create the first receipt from verified source tables",
    )
    args = parser.parse_args()
    if args.write and not (args.fetch or args.source_dir):
        parser.error("--write needs verified source tables")
    if args.fetch or args.source_dir:
        data = extract(*read_sources(args.source_dir))
    else:
        data = json.loads(RECEIPT.read_text())["published"]
    result = {
        "sources": {v: {"tex_member": f, "sha256": d} for v, (f, d) in SOURCES.items()},
        "published": data,
        **evaluate(data),
    }
    if args.write:
        with RECEIPT.open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
    elif result != json.loads(RECEIPT.read_text()):
        raise ValueError(
            "Published inputs or derived comparisons differ from the receipt"
        )
    print(json.dumps(result["profile_comparisons"], indent=2))
    print(json.dumps(result["constant_exponent_test"], indent=2))


if __name__ == "__main__":
    main()
