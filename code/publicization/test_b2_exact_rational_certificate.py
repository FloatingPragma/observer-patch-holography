from b2_exact_rational_certificate import ONE, build_certificate


def test_b2_exact_rational_certificate() -> None:
    certificate = build_certificate()
    assert certificate["schema"] == "oph.b2-publicization.exact-rational.v1"
    assert certificate["basis_rows_checked"] == 9
    assert certificate["trace_rows_checked"] == 9
    assert certificate["partition_ranks"] == ["2", "1"]
    assert certificate["kraus_completeness"] == [
        str(ONE[i][j]) for i in range(3) for j in range(3)
    ]
