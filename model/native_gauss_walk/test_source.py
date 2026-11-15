from fractions import Fraction as F
from collections import defaultdict
import pytest
from . import source as s


@pytest.fixture(scope='module')
def primary():return s.exact_protocol()


def test_preparation_exact_native_strings_and_full_events():
    w=s.prepare();assert w.probability()==1 and len(w.amplitudes)==720
    for(c,e),a in w.amplitudes.items():
        assert e<60 and a==(1,0)
        assert not any(s.electric.gauss(s.state(c,1)))
        assert c==s.canonical_config(c[0])
        assert s.shortest_path(c[0])[0]==c[0][0]
    assert s.prepare(q=0).probability()==1


def test_all_reached_event_columns_and_padding_are_involutions(primary):
    configs={c for w in primary[1]for c,e in w.amplitudes}
    for c in configs:
      for e in range(64):
        target,p,j,failed=s.event_column(c,e,1,1)
        returned,p2,j2,failed2=s.event_column(target,e,1,1)
        assert returned==c and(p+p2)%6==0
        assert failed==failed2 and j==-j2
        assert s.electric.gauss(s.state(c,1))==s.electric.gauss(s.state(target,1))
    c=s.canonical_config((0,2));w=s.Wave({(c,63):(1,0)})
    assert s.step(w).amplitudes=={(c,63):(-1,0)}
    assert s.equivalent(s.step(s.step(w)),w)
    edge=s.electric.EDGE_INDEX[(0,5)];state=s.electric.pair_string();flux=list(state.flux);flux[edge]=-4
    c=(state.pair,state.tagged_port,tuple(flux));e=2*edge+1
    assert s.event_column(c,e,1,1)==(c,0,0,True)


def test_complete_joint_readouts_and_inverse_no_reset(primary):
    result,states=primary
    assert[result['records'][i]['nonzero_joint_arcs']for i in range(3)]==[720,720,7320]
    assert[r['electric_energy_mean']for r in result['records']]==['5/32','13/80','2793/16000']
    assert result['source_plus_electric_work_per_step']==['55/9216','1913719/165888000']
    for record in result['records']:
        rows=s.unpack_records(record)
        assert len(rows)==record['nonzero_joint_arcs']
        assert sum(r[-1]for r in rows)==record['joint_probability_denominator']
        assert all(s.norm(tuple(r[-2]))==r[-1]for r in rows)
        assert sum(map(F,record['all64_event_probabilities']))==1
        assert record['gauss_sector_probabilities']==[{'G':[0]*12,'probability':'1'}]
    assert result['records'][2]['canonical_dressing_readout']['complete_outside_dressing_probability']=='32821/81000'
    for inverse in result['inverse_controls']:
        assert inverse['full_state_exactly_restored'] and inverse['decoded']['canonical_dressing_probability']=='1'
        assert not inverse['source_link_tag_event_reset']


def test_electric_phase_first_observed_on_third_step(primary):
    result,states=primary;phase,phased=s.exact_protocol(electric_phase=True)
    assert phase['records'][2]['source_pair_probabilities']==result['records'][2]['source_pair_probabilities']
    assert phase['records'][2]['all64_event_probabilities']==result['records'][2]['all64_event_probabilities']
    a=s.complete_records(s.step(states[2]));b=s.complete_records(s.step(phased[2],electric_phase=True))
    tv=sum(abs(F(x['probability'])-F(y['probability']))for x,y in zip(a['source_pair_probabilities'],b['source_pair_probabilities']))/2
    assert tv==F(51876833,583200000)
    assert a['nonzero_joint_arcs']==46020 and b['nonzero_joint_arcs']==46020
    assert a['gauss_sector_probabilities']==b['gauss_sector_probabilities']==[{'G':[0]*12,'probability':'1'}]


def test_charge_off_and_undressed_controls_are_distinct():
    neutral,_=s.exact_protocol(q=0);bad,_=s.exact_protocol(dressed=False)
    assert all(r['electric_energy_mean']=='0'for r in neutral['records'])
    assert all(r['gauss_sector_probabilities']==[{'G':[0]*12,'probability':'1'}]for r in neutral['records'])
    assert any(r['G']!=[0]*12 and F(r['probability'])>0 for r in bad['records'][1]['gauss_sector_probabilities'])
    assert all(r['norm_squared']=='1'for r in bad['records'])


def test_arbitrary_reference_columns_keep_their_exact_gram():
    w=s.prepare();pairs=set(s.MATCHINGS[0]);a=s.Wave({k:v for k,v in w.amplitudes.items()if k[0][0]in pairs})
    b=s.Wave({k:s.rotate(v,1)for k,v in w.amplitudes.items()if k[0][0]not in pairs})
    def inner(x,y):
        out=(0,0)
        for k,v in x.amplitudes.items():
            z=y.amplitudes.get(k,(0,0));conj=(v[0]+v[1],-v[1])
            product=(conj[0]*z[0]-conj[1]*z[1],conj[0]*z[1]+conj[1]*z[0]+conj[1]*z[1])
            out=s.add(out,product)
        return tuple(F(v,720*x.denominator*y.denominator)for v in out)
    before=[[inner(x,y)for y in(a,b)]for x in(a,b)]
    x,y=a,b
    for _ in range(2):x=s.step(x,electric_phase=True);y=s.step(y,electric_phase=True)
    assert [[inner(v,u)for u in(x,y)]for v in(x,y)]==before
    for _ in range(2):x=s.step(x,inverse=True,electric_phase=True);y=s.step(y,inverse=True,electric_phase=True)
    assert s.equivalent(x,a) and s.equivalent(y,b)


def test_invalid_event_or_link_cannot_be_silently_dropped():
    c=s.canonical_config((0,2))
    with pytest.raises(ValueError):s.Wave({(c,64):(1,0)})
    with pytest.raises(ValueError):s.Wave({(c,0):(1.0,0)})
    with pytest.raises(ValueError):s.Wave({(c,0):(1,0)},0)
    bad=(c[0],c[1],(5,)+c[2][1:])
    with pytest.raises(ValueError):s.Wave({(bad,0):(1,0)})
