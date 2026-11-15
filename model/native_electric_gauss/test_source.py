from dataclasses import replace
from fractions import Fraction as F
import pytest
from . import source as s


def test_current_sign_tag_and_boundary_rejection():
    state=s.pair_string((0,1));before=s.gauss(state)
    hop=s.apply_native_current(state,(0,5),1,1)
    assert hop.state.pair==(1,5) and hop.state.tagged_port==5
    assert hop.realized_signed_current==1 and s.gauss(hop.state)==before
    undressed=replace(hop.state,flux=state.flux)
    assert s.gauss(undressed)!=before
    assert s.apply_native_current(hop.state,(0,5),1,1).state==state
    blocked=s.apply_native_current(s.pair_string((0,1),2),(0,5),1,2)
    assert blocked.state.tagged_port==5
    edge=s.EDGE_INDEX[(0,5)];f=list(state.flux);f[edge]=-1
    edge_state=replace(state,flux=tuple(f));failed=s.apply_native_current(edge_state,(0,5),1,1)
    assert failed.overflow and failed.realized_signed_current==0 and failed.state==edge_state
    bad=list(f);bad[edge]=1
    assert s.gauss(replace(hop.state,flux=tuple(bad)))!=s.gauss(edge_state)


def test_full_native_raw_slots_source_and_memory_agree():
    for pair in s.PAIRS:
      for seam in s.SEAMS:
        state=s.ElectricState(pair,pair[0],1,0,(0,)*30)
        rows=s.complete_raw_events(state,seam,limit=2)
        assert len(rows)==50 and sum(r['probability']for r in rows)==1
        for slot,branch in enumerate(s.native.raw_local_branches()):
            target=s.active.lifted_column(pair,0,seam,branch)
            row=rows[2*slot]
            if target is None:assert row['state']is None
            elif row['state']is not None:
                assert(row['state'].pair,row['state'].memory)==target
            assert rows[2*slot+1]['probability']==0


def test_complete_loop_records_work_and_quadratures():
    r=s.loop_experiment()
    assert[r['initial_electric_energy'],r['after_first_nonselective_energy'],r['after_second_nonselective_energy']]==['1/16','5/32','11/32']
    records=r['two_complete_reference_blocks']
    for block in(records[:4],records[4:]):
        assert sum(F(x['joint_probability'])for x in block)==1
        for row in block:
            assert sum(F(y['probability'])for y in row['complete_flux_source_memory_weights'])==F(row['joint_probability'])
    assert r['boundary_control_overflows_retained']>0
    assert r['boundary_word_inverse_exact']
    # Commit the actual cycle-flux record between blocks. Coherent memory
    # is retained within each flux block, so only path information is removed.
    first=s.loop_reference_measurement(s.character_preparation(),1,2)
    for state in first:
        groups={}
        for x,a in state.items():groups.setdefault(x.flux,{})[x]=a
        probabilities=[F(0),F(0)]
        for group in groups.values():
            for j,out in enumerate(s.loop_reference_measurement(group,1,2)):
                probabilities[j]+=s.probability(out)
        assert probabilities==[F(1,4),F(1,4)]


def test_charge_preparation_poisson_and_electric_phase():
    prep=s.preparation_certificate();assert prep['executed_preparations']==90
    rows=s.coulomb_comparison()['rows'];assert len(rows)==11
    assert rows[0]['integer_ground_electric_energy']=='1/16'
    assert rows[0]['graph_poisson_relaxed_energy']=='11/480'
    assert all(F(r['strict_integer_excess'])>0 for r in rows)
    phase=s.electric_phase_readout()
    assert phase['records'][1]['complete_pointer_probabilities']==['0','1']
    assert all(F(r['field_work_before_phase'])+F(r['field_work_after_inverse'])==0 for r in phase['records'])
    assert phase['zero_kappa_negative_control_probabilities']==['1','0']
    assert s.maxwell_identities()['cycle_divergence_exact']==[0]*12


def test_complete_instrument_with_arbitrary_reference_superposition():
    a=s.pair_string();b=s.apply_native_current(a,(0,5)).state
    # Distinct Gauss states can carry arbitrary reference amplitudes. This
    # exact superposition includes a nontrivial root phase and a boundary.
    state={a:s.ONE,b:s.ROOTS[1]}
    for limit in(1,2):
        image=s.loop_operator(state,limit=limit)
        assert s.loop_operator(image,inverse=True,limit=limit)==state
        for phase in range(6):
            outcomes=s.loop_reference_measurement(state,phase,limit)
            assert sum(s.probability(o)for o in outcomes)==s.probability(state)
    with pytest.raises(ValueError):s.apply_native_current(a,(0,2))
    with pytest.raises(ValueError):s.pair_string((0,2))


def test_all_cut_charge_records_and_coherent_half_control():
    r=s.cut_charge_certificate()
    assert r['regions_per_state']==4096 and r['native_source_states_checked']==8
    assert len(r['initial_full_region_charge_table'])==4096
    assert set(r['initial_full_region_charge_table'])=={-1,0,1}
    assert r['loop_deformation_cut_tables_exact']
    control=r['coherent_half_expectation_control']
    assert control['regional_expectations']==['1/2','1/2']
    assert control['joint_records']==[{'charges':[0,1],'probability':'1/2'},{'charges':[1,0],'probability':'1/2'}]
    assert control['sharp_half_half_record_probability']=='0'
    assert sum(F(x['probability'])for x in control['all12_reference_memory_records'])==1
    with pytest.raises(ValueError):s.cut_flux(s.pair_string(),1<<12)
