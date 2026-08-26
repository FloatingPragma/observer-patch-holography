#!/usr/bin/env python3
"""Download and hash-check the official DESI DR2 chains used by the audit."""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

from official_desi_dr2_chain_audit import DOWNLOAD_SPECS, SOURCE_COBAYA_ROOT, sha256


def source_url(spec: dict[str, object], index: int) -> str:
    return f"{SOURCE_COBAYA_ROOT}/{spec['model']}/{spec['directory']}/chain.{index}.txt"


def download_dataset(
    data_dir: Path, name: str, spec: dict[str, object], replace: bool
) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    for index, expected in enumerate(spec["sha256"], start=1):
        target = data_dir / f"{spec['slug']}_chain.{index}.txt"
        if target.exists() and sha256(target) == expected:
            print(f"verified existing {target}")
            continue
        if target.exists() and not replace:
            raise ValueError(
                f"existing file has the wrong hash: {target}; pass --replace to overwrite it"
            )
        partial = target.with_suffix(target.suffix + ".part")
        if partial.exists():
            partial.unlink()
        print(f"downloading {name} chain {index}: {source_url(spec, index)}")
        urllib.request.urlretrieve(source_url(spec, index), partial)
        actual = sha256(partial)
        if actual != expected:
            partial.unlink()
            raise ValueError(
                f"official SHA-256 mismatch for {name} chain {index}: {actual} != {expected}"
            )
        partial.replace(target)
        print(f"verified {target}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--dataset", action="append", choices=sorted(DOWNLOAD_SPECS))
    parser.add_argument("--replace", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    names = list(DOWNLOAD_SPECS) if args.dataset is None else args.dataset
    for name in names:
        download_dataset(args.data_dir, name, DOWNLOAD_SPECS[name], args.replace)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
