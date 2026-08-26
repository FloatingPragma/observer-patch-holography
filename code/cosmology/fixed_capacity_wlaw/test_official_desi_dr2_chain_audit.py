import importlib.util
import json
import sys
from pathlib import Path

import pytest


MODULE_PATH = Path(__file__).with_name("official_desi_dr2_chain_audit.py")
SPEC = importlib.util.spec_from_file_location("official_desi", MODULE_PATH)
assert SPEC and SPEC.loader
mod = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = mod
SPEC.loader.exec_module(mod)


def write_chain(path: Path, rows: list[tuple[float, float, float]]) -> None:
    lines = ["# weight w wa\n"]
    lines.extend(f"{weight} {w0} {wa}\n" for weight, w0, wa in rows)
    path.write_text("".join(lines), encoding="utf-8")


def write_base_chain(path: Path, rows: list[tuple[float, float, float]]) -> None:
    lines = ["# weight H0 omegal\n"]
    lines.extend(f"{weight} {h0} {omega_lambda}\n" for weight, h0, omega_lambda in rows)
    path.write_text("".join(lines), encoding="utf-8")


def test_cpl_full_interval_condition_uses_both_endpoints(tmp_path: Path) -> None:
    path = tmp_path / "chain.txt"
    write_chain(
        path,
        [
            (1, -0.9, -0.1),  # both endpoints pass
            (2, -0.9, -0.2),  # a=1 passes, a=1/3 fails
            (3, -1.1, 1.0),  # early endpoint passes, a=1 fails
            (4, -1.0, 0.0),  # boundary passes
        ],
    )
    summary = mod.read_chain(path).summary()
    assert summary["expanded_posterior_weight"] == 10
    assert summary["posterior_mass_w_ge_minus_one_for_0_le_z_le_2"] == 0.5
    assert summary["raw_rows_in_monotone_subset"] == 2


def test_weight_column_is_used(tmp_path: Path) -> None:
    path = tmp_path / "chain.txt"
    write_chain(path, [(9, -1.0, 0.0), (1, -1.2, 0.0)])
    summary = mod.read_chain(path).summary()
    assert summary["posterior_mass_w_ge_minus_one_for_0_le_z_le_2"] == 0.9


def test_missing_required_column_fails_closed(tmp_path: Path) -> None:
    path = tmp_path / "chain.txt"
    path.write_text("# weight w\n1 -1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="required column absent: wa"):
        mod.read_chain(path)


def test_nonpositive_weight_fails_closed(tmp_path: Path) -> None:
    path = tmp_path / "chain.txt"
    write_chain(path, [(0, -1.0, 0.0)])
    with pytest.raises(ValueError, match="weights must be positive"):
        mod.read_chain(path)


def test_hash_mutation_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "toy_chain.1.txt"
    write_chain(path, [(1, -1.0, 0.0)])
    monkeypatch.setitem(
        mod.DATASETS,
        "toy",
        {
            "slug": "toy",
            "directory": "toy",
            "sha256": ["0" * 64],
        },
    )
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        mod.audit_dataset(tmp_path, mod.DATASETS["toy"])


def test_base_lcdm_display_is_computed_sample_by_sample_with_weights(
    tmp_path: Path,
) -> None:
    path = tmp_path / "base.txt"
    rows = [(9.0, 60.0, 0.7), (1.0, 80.0, 0.6)]
    write_base_chain(path, rows)
    summary = mod.read_base_lcdm_chain(path).summary()
    direct = sum(
        weight * mod.lambda_lp2_from_base_lcdm_sample(h0, omega_lambda)
        for weight, h0, omega_lambda in rows
    ) / sum(weight for weight, _h0, _omega_lambda in rows)
    assert summary["Lambda_lP2"]["weighted_mean"] == pytest.approx(direct, rel=1e-15)
    assert summary["expanded_posterior_weight"] == 10.0
    assert summary["H0_OmegaLambda_weighted_correlation"] == pytest.approx(-1.0)


def test_base_lcdm_missing_or_unphysical_columns_fail_closed(tmp_path: Path) -> None:
    missing = tmp_path / "missing.txt"
    missing.write_text("# weight H0\n1 68\n", encoding="utf-8")
    with pytest.raises(ValueError, match="required column absent: omegal"):
        mod.read_base_lcdm_chain(missing)

    unphysical = tmp_path / "unphysical.txt"
    write_base_chain(unphysical, [(1.0, 68.0, -0.7)])
    with pytest.raises(ValueError, match="OmegaLambda must lie"):
        mod.read_base_lcdm_chain(unphysical)


def test_source_pins_and_committed_receipt_bind_current_producer() -> None:
    assert mod.SOURCE_MANIFEST_SHA256 == (
        "df78872aa8b2d3473a9e8de78f498180efd7cbcbeb18211ce4787fac52067ee5"
    )
    assert mod.BASE_LCDM_DATASET["sha256"] == [
        "00f3766f7a7b6370d21323886cd72869087b2b1346a04d729c8f3bc9e65ef698",
        "33b154eebdf4e9dca3b8f02ed2680120879d35c10b32fef42261a490104e1dc1",
        "d4717e7e5a13de851c86f24c87213faccef2b5f8747900274ab509d9dfa40aa2",
        "c827cd767a4864ca28aa15c902bda32004e803050d4be330e25aefddd78b5c36",
    ]
    receipt_path = (
        MODULE_PATH.with_name("runtime") / "official_desi_dr2_fz13_retrospective.json"
    )
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    assert receipt["schema"] == "oph.official_desi_dr2_fz13_retrospective.v2"
    assert receipt["producer"] == mod.producer_metadata()
    assert receipt["source"]["official_sha256_manifest_sha256"] == (
        mod.SOURCE_MANIFEST_SHA256
    )
    assert receipt["base_lcdm_capacity_display"]["combined"]["raw_rows"] > 0
