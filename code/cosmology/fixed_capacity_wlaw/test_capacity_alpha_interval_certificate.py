from __future__ import annotations

import copy
import importlib.util
import sys
from decimal import Decimal
from pathlib import Path

import pytest


HERE = Path(__file__).resolve().parent


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


producer = _load(
    "capacity_alpha_interval_certificate",
    "capacity_alpha_interval_certificate.py",
)
verifier = _load(
    "verify_capacity_alpha_interval_certificate",
    "verify_capacity_alpha_interval_certificate.py",
)


@pytest.fixture(scope="module")
def reduced_certificate() -> dict:
    return producer.build_certificate(
        mp_dps=35,
        iv_dps=35,
        su2_cutoff=24,
        su3_cutoff=16,
        log_alpha_half_width="0.00001",
    )


def test_interval_ift_certificate_excludes_zero(reduced_certificate: dict) -> None:
    branch = reduced_certificate["branch_certificate"]
    tangent = branch["local_derivatives"]["d_log_N_d_log_alpha"]
    reciprocal = branch["local_derivatives"]["d_log_alpha_d_log_N"]
    assert Decimal(tangent["lo"]) < Decimal(tangent["hi"]) < 0
    assert Decimal(reciprocal["lo"]) < Decimal(reciprocal["hi"]) < 0
    assert tangent["zero_excluded"] is True
    assert reciprocal["computed_only_after_zero_exclusion"] is True

    roots = branch["implicit_function_certificate"]["alpha_U_pixel_closure"]
    assert Decimal(roots["R_u_enclosure"]["hi"]) < 0
    assert roots["R_u_sign_definite"] is True
    for key in ("m_z_ift_for_R_u", "m_z_ift_for_input_seed"):
        assert Decimal(roots[key]["h_m_enclosure"]["hi"]) < 0
        assert roots[key]["h_m_sign_definite"] is True


def test_mean_value_bound_and_physical_boundary(reduced_certificate: dict) -> None:
    branch = reduced_certificate["branch_certificate"]
    mean_value = branch["mean_value_certificate"]
    assert mean_value[
        "valid_for_every_pair_of_alpha_values_in_certified_domain"
    ] is True
    assert Decimal(0) < Decimal(mean_value["abs_slope_lower"]) < Decimal(
        mean_value["abs_slope_upper"]
    )
    assert reduced_certificate["classification"][
        "mathematical_branch_differentiability"
    ] == "attained_on_certified_domain"
    assert reduced_certificate["classification"]["physical_epoch_evolution"] == (
        "undischarged"
    )
    assert all(
        premise["status"].startswith("undischarged")
        for premise in reduced_certificate["premises"].values()
    )
    domain = branch["certified_domain"]
    assert Decimal(
        domain["guaranteed_symmetric_log_alpha_inner_radius"]
    ) == Decimal("0.00001")
    assert Decimal(
        domain["outer_enclosure_max_abs_log_alpha_displacement"]
    ) >= Decimal("0.00001")


def test_independent_factorized_replay(reduced_certificate: dict) -> None:
    result = verifier.verify_certificate(reduced_certificate, replay=True)
    assert result["verified"] is True
    assert result["independent_interval_replay"] is True
    replay = result["replay"]
    assert Decimal(replay["tangent"]["hi"]) < 0
    assert Decimal(replay["R_u"]["hi"]) < 0
    assert Decimal(replay["h_m"]["hi"]) < 0


def test_uncertified_wide_box_fails_closed() -> None:
    with pytest.raises(
        producer.CapacityCertificateError,
        match="d log N / d log alpha interval contains zero",
    ):
        producer.certify_branch(
            mp_dps=35,
            iv_dps=35,
            su2_cutoff=24,
            su3_cutoff=16,
            log_alpha_half_width="0.001",
        )


def test_mutated_tangent_fails_closed(reduced_certificate: dict) -> None:
    mutated = copy.deepcopy(reduced_certificate)
    tangent = mutated["branch_certificate"]["local_derivatives"][
        "d_log_N_d_log_alpha"
    ]
    tangent["hi"] = "0"
    with pytest.raises(verifier.VerificationError, match="contains zero"):
        verifier.verify_certificate(mutated, replay=False)


def test_narrowed_enclosure_fails_independent_replay(
    reduced_certificate: dict,
) -> None:
    mutated = copy.deepcopy(reduced_certificate)
    tangent = mutated["branch_certificate"]["local_derivatives"][
        "d_log_N_d_log_alpha"
    ]
    tangent["lo"] = "-0.2102"
    tangent["hi"] = "-0.2101"
    tangent["encloses_direct_dual_evaluation"] = {
        "lo": "-0.2102",
        "hi": "-0.2101",
    }
    tangent["encloses_factorized_chain_rule_evaluation"] = {
        "lo": "-0.2102",
        "hi": "-0.2101",
    }
    mean_value = mutated["branch_certificate"]["mean_value_certificate"]
    mean_value["abs_slope_lower"] = "0.2101"
    mean_value["abs_slope_upper"] = "0.2102"
    with pytest.raises(verifier.VerificationError, match="does not contain replay"):
        verifier.verify_certificate(mutated, replay=True)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda artifact: artifact["premises"]["B1_same_capacity"].update(
                status="discharged"
            ),
            "physical premise statements or statuses changed",
        ),
        (
            lambda artifact: artifact["branch_certificate"][
                "implicit_function_certificate"
            ]["alpha_U_pixel_closure"].update(R_u_sign_definite=False),
            "R_u sign-definiteness missing",
        ),
        (
            lambda artifact: artifact["bindings"].update(
                interval_engine_sha256="0" * 64
            ),
            "stale or mutated bound file",
        ),
        (
            lambda artifact: artifact["bindings"].update(
                paper_math_sha256="0" * 64
            ),
            "stale or mutated bound file",
        ),
    ],
)
def test_status_denominator_and_binding_mutations_fail_closed(
    reduced_certificate: dict, mutation, message: str
) -> None:
    mutated = copy.deepcopy(reduced_certificate)
    mutation(mutated)
    with pytest.raises(verifier.VerificationError, match=message):
        verifier.verify_certificate(mutated, replay=False)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda artifact: artifact["branch_certificate"][
                "mean_value_certificate"
            ].update(abs_slope_lower="1", abs_slope_upper="2"),
            "mean-value slopes do not match tangent enclosure",
        ),
        (
            lambda artifact: artifact["branch_certificate"][
                "implicit_function_certificate"
            ].update(selected_branch_C1_on_domain=False),
            "selected branch C1 verdict missing",
        ),
        (
            lambda artifact: artifact["classification"].update(
                physical_attachment="attained"
            ),
            "classification or physical attachment changed",
        ),
        (
            lambda artifact: artifact.update(nonclaims=[]),
            "required physical nonclaims changed",
        ),
        (
            lambda artifact: artifact["premises"]["B2_epochwise_bridge"].update(
                statement="mathematics implies physical evolution"
            ),
            "physical premise statements or statuses changed",
        ),
        (
            lambda artifact: artifact["branch_certificate"][
                "implicit_function_certificate"
            ]["alpha_U_pixel_closure"].update(endpoint_signs_verified=False),
            "alpha_U endpoint-sign verdict missing",
        ),
        (
            lambda artifact: artifact["branch_certificate"][
                "implicit_function_certificate"
            ]["alpha_U_pixel_closure"].update(orientation="increasing"),
            "alpha_U orientation contradicts R_u",
        ),
        (
            lambda artifact: artifact["branch_certificate"]["local_derivatives"][
                "d_log_N_d_log_alpha"
            ].update(
                encloses_direct_dual_evaluation={"lo": "-1", "hi": "-0.9"}
            ),
            "direct-dual tangent subblock does not contain replay",
        ),
    ],
)
def test_rehashed_semantic_mutations_fail_without_replay(
    reduced_certificate: dict, mutation, message: str
) -> None:
    """Source hashes remain valid, so rejection comes from semantic checks."""

    mutated = copy.deepcopy(reduced_certificate)
    mutation(mutated)
    assert mutated["bindings"] == reduced_certificate["bindings"]
    with pytest.raises(verifier.VerificationError, match=message):
        verifier.verify_certificate(mutated, replay=False)


def test_coherently_inflated_inner_radius_fails_without_replay(
    reduced_certificate: dict,
) -> None:
    mutated = copy.deepcopy(reduced_certificate)
    branch = mutated["branch_certificate"]
    branch["configuration"]["requested_log_alpha_half_width"] = "0.1"
    branch["certified_domain"][
        "guaranteed_symmetric_log_alpha_inner_radius"
    ] = "0.1"
    branch["certified_domain"][
        "outer_enclosure_max_abs_log_alpha_displacement"
    ] = "0.1"
    branch["mean_value_certificate"][
        "guaranteed_symmetric_log_alpha_inner_radius"
    ] = "0.1"
    branch["mean_value_certificate"][
        "outer_enclosure_max_abs_log_alpha_displacement"
    ] = "0.1"
    assert mutated["bindings"] == reduced_certificate["bindings"]
    with pytest.raises(
        verifier.VerificationError,
        match="certified domain does not contain inner log ball",
    ):
        verifier.verify_certificate(mutated, replay=False)


@pytest.mark.parametrize("bad_value", ["NaN", "not-a-number", None])
def test_malformed_decimal_fields_fail_cleanly(
    reduced_certificate: dict, bad_value
) -> None:
    mutated = copy.deepcopy(reduced_certificate)
    mutated["branch_certificate"]["mean_value_certificate"][
        "abs_slope_lower"
    ] = bad_value
    with pytest.raises(verifier.VerificationError, match="decimal value"):
        verifier.verify_certificate(mutated, replay=False)
