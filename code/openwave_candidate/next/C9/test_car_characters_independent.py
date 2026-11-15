from review_car_characters import run


def test_entire_character_space_and_global_reference():
    r=run()
    assert all(r['checks'].values())
    assert r['global_ground_multiplicity']==2
    assert r['gap_above_ground_multiplet']>.28
    assert r['H_error']<1e-12 and r['all_K_error']<1e-12
    assert not r['full_criterion_met']
