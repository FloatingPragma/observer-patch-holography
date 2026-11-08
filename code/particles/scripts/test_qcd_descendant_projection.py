"""Current hadron projections reproduce their conditional parent receipts."""
import json

import build_postdiction_ledger as ledger
import build_source_only_mass_prediction_surface as surface


def test_current_postdiction_hadrons_are_the_official_parent_projection():
    actual=json.loads(ledger.DEFAULT_OUT.read_text())['sections']['hadrons']
    expected=ledger._hadron_rows(ledger._load('hadron_payload'),ledger._load('solver_standby'),
                                ledger._load('lambda_transmutation'),ledger._load('nucleon_external_ratio'))
    assert actual==expected
    nucleon=next(row for row in actual if row['id']=='nucleon_mass_external_ratio')
    assert nucleon['row_class']=='oph_plus_external_numeric_ratio'
    assert 'provenance is unverified' in nucleon['claim_boundary']


def test_current_mass_surface_hadrons_are_the_official_parent_projection():
    actual=json.loads(surface.DEFAULT_JSON_OUT.read_text())
    hadron=next(row for row in actual['families'] if row['family']=='hadrons')
    expected=surface._hadron_rows(ledger._load('lambda_transmutation'),ledger._load('nucleon_external_ratio'))
    assert hadron==expected
    assert hadron['rows'][0]['range_scope']=='sampled_threshold_corner_ranges'


def test_live_ward_source_grid_is_complete_and_unpromoted():
    path=ledger.PARTICLES/'hadron/ward_projected_payload/runtime/ward_projected_payload_source_bracket_current.json'
    value=json.loads(path.read_text())
    assert value['grid']['runs']==58
    assert value['promotion_allowed'] is False
    assert value['scoring_status']=='NOT_EVALUABLE_SOURCE_DIAGNOSTIC'
