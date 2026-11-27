from tools.check_observables import completion


def test_established_weights_apply_to_every_observable():
    rows = [
        dict(id="conditional", kind="observable", status="Conditional", premises=["p", "q"]),
        dict(id="p", kind="premise", status="Demonstrated"),
        dict(id="q", kind="premise", status="Open"),
        dict(id="formal", kind="observable", status="Demonstrated", target_class="formalism"),
        dict(id="comparison", kind="observable", status="Partial", scope="comparison"),
        dict(id="open", kind="observable", status="Open"),
    ]
    assert completion(rows) == (0.75 + 1 + 0.25 + 0) / 4
    rows[0]["premises"] = []
    assert completion(rows) == (0.5 + 1 + 0.25 + 0) / 4
