"""Complete records, covariance-independent exact probabilities, and finite-domain mutation controls."""
from fractions import Fraction as F
import pytest
from . import model as m
from .certificate import certificate


def test_every_native_event_preserves_load_and_decreases_combined_energy():
    for state in m.STATES:
        for seam in range(3):
            for coin in range(2):
                target,record=m.event(state,seam,coin)
                assert sum(target[:3])==sum(state[:3])
                assert m.energy(target)<=m.energy(state)
                assert record['no_operation']==(target==state)


def test_no_operation_records_keep_probability_and_operational_time():
    state=(2,2,2,3)
    assert m.kernel()[state]=={state:F(1)}
    assert m.evolve({state:F(1)},8)=={state:F(1)}
    records=[m.event(state,s,c)[1] for s in range(3) for c in range(2)]
    assert len(records)==6 and all(record['no_operation'] for record in records)


def test_class_only_projection_does_not_define_a_closed_abundance_law():
    # The same zero flux gives different birth probabilities at different load contrasts.
    assert all(s[3]==0 for s in m.kernel()[(0,0,0,0)])
    assert any(s[3]!=0 for s in m.kernel()[(0,0,5,0)])


def test_rejected_register_twist_does_not_delete_the_load_event():
    witnesses=[]
    for state in m.STATES:
        target,record=m.event(state,0,0)
        if record['current'] and not record['applied'] and target!=state:witnesses.append((state,target))
    assert witnesses
    assert all(a[3]==b[3] and a[:3]!=b[:3] for a,b in witnesses)


def test_domain_and_probability_mutations_rejected():
    for args in [((6,0,0,0),0,0),((0,0,0,6),0,0),((0,0,0,0),3,0),((0,0,0,0),0,2)]:
        with pytest.raises(ValueError):m.step(*args)
    with pytest.raises(ValueError):m.evolve({(0,0,0,0):F(1,2)})
    with pytest.raises(ValueError):m.evolve(m.initial_distribution(),-1)
    for steps in (-1,1.0,True):
        with pytest.raises(ValueError):m.evolve(m.initial_distribution(),steps)
    for distribution in ({(6,0,0,0):F(1)},{(0,0,0,0):1.0},{(0,0,0,0):True}):
        with pytest.raises(ValueError):m.evolve(distribution,0)
    for args in [((0.,0,0,0),0,0),((0,0,0,0),True,0),((0,0,0,0),0,1.0)]:
        with pytest.raises(ValueError):m.event(*args)



def test_population_curve_receipt_is_exact_and_recomputes():
    import json
    from pathlib import Path
    result=certificate()
    assert result['all_checks_pass']
    assert json.loads(Path(__file__).with_name('results.json').read_text())==result
    assert result['construction']['temperature'] is None


def test_persistent_abundance_solves_the_exact_backward_equation():
    values=m.persistent_probability()
    for state,row in m.kernel().items():
        assert values[state]==sum(weight*values[target] for target,weight in row.items())
        if m.extinct_or_frozen(state):assert values[state]==int(state[3]!=0)
    assert sum(weight*values[state] for state,weight in m.initial_distribution().items())>0


def test_uniform_epoch_bound_uses_every_initial_state():
    survival=m.survival_probabilities(32)
    assert len(survival)==1296 and max(survival.values())<=F(1,25)
    assert max(m.survival_probabilities(3).values())==F(107,108)
    assert all(value==0 for state,value in survival.items() if m.extinct_or_frozen(state))
