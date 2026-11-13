from fractions import Fraction as F
from .native import *
from .source import source_identity,calculation,support_norm


def test_complete_CAR_source_matches_independent_triquark_and_meson_frames():
    d=source_identity();assert d['triquark_floor']>1700 and d['meson_floor']>270
    assert sum(d['triquark_blocks'])==100 and sum(d['meson_blocks'])==48


def test_antiquark_join_CAR_hermiticity_and_global_singlet():
    column=pair_column(((0,0,1),(1,1,0)))
    assert not action_vector(colour_casimir,column,(0,1))
    for s in column:
        for target,value in base_action(s).items():
            assert base_action(target).get(s,0)==value
    # All actual outputs outside the original two cells remain represented.
    out=action_vector(base_action,column)
    assert any(any(mode[1]not in(0,1)for mode in s)for s in out)


def test_complete_autonomous_projection_and_positive_remainder():
    d=calculation()
    assert F('0.02761018')<d['a2']<F('0.02761019')
    assert d['full']!=d['components']['base']
    assert d['tau']*d['a']-d['tau']**2*d['B']/2==d['amplitude']>0
    assert 0<d['probability']<d['cab']/216
    assert support_norm(3)>support_norm(1)>0
