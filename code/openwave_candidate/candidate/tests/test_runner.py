import copy
import pytest
from oph_openwave.runner import validate_result


def test_failure_flags_and_physical_promotions_are_rejected():
    base = {"scope": "finite control", "checks": {"identity": True},
            "mutations": {"broken_identity_rejected": True}, "full_criterion_met": False}
    validate_result("control", base)
    for field, value in (("checks", {}), ("checks", {"identity": 1}),
                         ("checks", {"identity": False}), ("mutations", {}),
                         ("mutations", {"broken": False}), ("scope", ""),
                         ("full_criterion_met", True)):
        bad = copy.deepcopy(base)
        bad[field] = value
        with pytest.raises(ValueError):
            validate_result("control", bad)
