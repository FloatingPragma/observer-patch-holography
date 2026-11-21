import json

from model.__main__ import main
from .capacity import HERE, profile


def test_whole_interval_domain_and_overflow_controls():
    rows = [profile(width) for width in (3, 4, 5)]
    assert rows == json.loads((HERE / 'capacity.json').read_text())['profiles']
    assert [r['outside_tolerance'] for r in rows] == [list(range(8, 32)), list(range(16, 32)), []]


def test_precision_trip_is_executable_through_model(capsys):
    assert main(['run', '--law', 'interval_echo', '--clock-trip', '8', '16']) == 0
    trip = json.loads(capsys.readouterr().out)
    assert (trip['age_bits'], trip['moving_cells'], trip['transfer_cells']) == (5, 14, 15)
    assert trip['reunion']['traveller'] == trip['home_retained_traveller']
