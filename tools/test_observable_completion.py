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


def test_duplicate_premise_cannot_reweight_a_conditional_observable(tmp_path, monkeypatch, capsys):
    import json
    from tools import check_observables as ledger

    data = json.loads(ledger.LEDGER.read_bytes())
    row = next(r for r in data["rows"] if r["status"] == "Conditional" and len(r["premises"]) > 1)
    row["premises"].append(row["premises"][0])
    path = tmp_path / "ledger.json"
    path.write_text(json.dumps(data))
    monkeypatch.setattr(ledger, "LEDGER", path)
    assert ledger.main() == 1
    assert f"{row['id']}: premises contains repeated entries" in capsys.readouterr().err
