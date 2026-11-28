"""Guard the source and the same-event original-start receiver comparison."""

import json

import pytest

from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.theorem_interfaces.interacting_chain.__main__ import (
    HERE,
    RECEIVER,
    joined_observer_records,
)


def records():
    return (
        json.loads((HERE / "results.json").read_bytes()),
        json.loads((HERE / "independent/records.json").read_bytes()),
    )


def test_local_spin_witness_on_preserved_banks():
    summary = joined_observer_records(*records())
    assert summary["cohorts"] == 16
    assert summary["simultaneous_raw_report_error_probability_upper"] < 0.01
    assert (
        min(
            summary[
                "minimum_raw_X_minus_Y_plus_Z_minus_admission_lower_by_wait"
            ].values()
        )
        > 0.0005
    )


@pytest.mark.parametrize(
    "alteration,reason",
    [
        ("length", "eight-site source"),
        ("address", "same terminal gate"),
        ("source", "source ancestry differs"),
        ("refusal", "omits or corrupts"),
    ],
)
def test_no_join_across_changed_domains(alteration, reason):
    baseline, independent = records()
    branch = independent["signed_branches"][0]["free_readouts"]["half_beat"][RECEIVER][
        0
    ]
    if alteration == "length":
        independent["length"] = 4
    elif alteration == "address":
        branch["spin_site"] = 4
    elif alteration == "source":
        independent["archived_source_sha256"] = {}
    else:
        branch["simulated_local_Pauli_comparison"]["banks"][0]["counts"][0] -= 1
    with pytest.raises(ValueError, match=reason):
        joined_observer_records(baseline, independent)
