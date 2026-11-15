from spectral_apparatus import run

def test_same_H_complete_spectral_work_reference_ledger():
    r=run()
    assert all(r['checks'].values())
    assert r['dimension']==6
    assert len(r['complete_records'])==36*len({tuple(sorted(x['native_record'].items())) for x in r['complete_records']})
    assert any(x['source_energy_change']<-.01 and x['probability']>1e-9 for x in r['complete_records'])
    assert not r['full_criterion_met']
