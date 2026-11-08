from fractions import Fraction as F
from math import comb
import pytest
from . import bell


@pytest.mark.parametrize('level',(0,1))
@pytest.mark.parametrize('choi',(False,True))
def test_actual_source_complete_parity_instrument(level,choi):
    row,out=bell.branch_run(level,choi)
    assert row['complete_raw_slots']==25
    assert row['nonzero_raw_slots']==6
    assert row['Bell_success_mass']==F(1,2)
    assert len(out)==25
    assert bell.charge_block_check(level)


@pytest.mark.parametrize('attempts',(0,1,2,5))
def test_bounded_restoration_retains_all_nonzero_paths(attempts):
    _,rows=bell.branch_run();_,edge,_=bell.layout()
    for row in rows:
        if row.density and row.record.current:
            mask=next(iter(row.density))[0]
            paths=bell.restore_basis(mask,edge,1,attempts)
            assert sum(w for _,_,w,_ in paths)==1
            assert sum(w for _,_,w,s in paths if s=='RESTORATION_TIMEOUT')==F(1,2**attempts)
            assert all(len(h)<=attempts for h,_,_,_ in paths)


def test_retry_complete_law_and_resource_cost():
    for rounds in (1,2,8,16):
        for restoration in (0,1,16):
            row=bell.retry_budget(rounds,restoration)
            assert row['Bell_success_mass_given_two_eligible_cells']+row['restoration_timeout_mass']+row['round_exhaustion_mass']==1
            assert row['maximum_native_attempts']==22+rounds+(rounds-1)*restoration
            assert row['maximum_local_source_pulses']==44*rounds
    assert bell.retry_budget(8,16)['Bell_success_mass_given_two_eligible_cells']>F(996,1000)


def test_fixed_number_source_handoff_rejects_leakage_and_no_direct_sum_schema():
    with pytest.raises(ValueError):
        bell.native.carrier({'kind':'federation','level':0,'total':'all'})
    with pytest.raises(ValueError):
        bell.native.carrier({'kind':'federation_direct_sum','level':0})
    carrier=bell.native.carrier({'kind':'federation','level':0,'total':2})
    with pytest.raises(ValueError):
        bell.native._validate_density({(0,0,0,0):F(1)},carrier)


def test_declared_local_beams_do_not_use_cross_cell_targets():
    words=bell.service_words()
    assert len(words)==4
    assert all(len(set(row['local_ports']))==2 and max(row['local_ports'])<12 for row in words)
    assert all(row['resources']['native_pulses']==22 for row in words)


def test_full_column_beam_and_inverse_with_complete_failure_bound():
    rows=bell.beam_execution_certificate()
    assert len(rows)==2
    assert all(row['complete_input_columns']==4 and row['pulses']==22 for row in rows)
    assert all(F(row['error']['complete_halfdiamond_upper'])<F(1,10**9) for row in rows)


def test_restoration_rejects_nonrestorable_or_malformed_inputs():
    _,edge,_=bell.layout()
    with pytest.raises(ValueError):bell.restore_basis(0,edge,1,2)
    with pytest.raises(ValueError):bell.restore_basis(True,edge,1,2)
    with pytest.raises(ValueError):bell.restore_basis(0,(0,1),0,2)


def test_source_dicke_separation_is_exact():
    sizes=[comb(12,2),12*12]
    assert max(sizes)==144
    assert F(144,comb(240,2))==F(6,1195)
    assert 1-F(144,comb(240,2))==F(1189,1195)


def test_C_no_current_is_not_E_even_projection():
    graph,edge,spare=bell.layout();masks=bell.code_masks(edge,spare)
    command={'op':'repair','selector':{'kind':'fixed','seam':list(edge)},'law':'C','repeat':1}
    program=bell.native.program(federation={'level':0,'total':2},biased=False,commands=[command])
    odd=masks[1]
    rows=bell.native._outcomes(program,command,{(odd,0,odd,0):F(1)})
    survived=sum(bell.native.trace(r.density) for r in rows if r.record.law=='C' and r.record.current==0)
    assert survived==F(1,2)


def test_literal_angle_evaluator_cannot_silently_replace_pi():
    from . import literal_angles
    row=literal_angles.certificate()
    assert row['pulses']==22
    assert len(row['nonzero_number_leakage_entries'])==8
    for witness in row['nonzero_number_leakage_entries']:
        lo,hi=map(F,witness['interval'])
        assert lo<=hi and (lo>0 or hi<0)


def test_actual_full_M6_literal_control_leakage_intervals_against_independent_precision():
    import mpmath as mp
    from . import literal_controls as lc
    row=lc.certificate()
    assert row['complete_complement_rows']==80
    assert len(row['successful_work_number_leakage_entries'])==8
    with mp.workdps(100):
        number=lambda q:mp.mpf(q.numerator)/q.denominator
        output=mp.eye(4)
        for event in bell.coupling.words()['beam'][1]:
            controls=lc.source.literal_controls(event['kind'],event['angle'])
            exact=lc.source.specification()
            generator=mp.matrix(6)
            for i in range(6):
                for j in range(6):
                    generator[i,j]=sum(number(c)*(number(g[i][j][0][0])+number(g[i][j][0][1])*mp.sqrt(5)
                                      +1j*(number(g[i][j][1][0])+number(g[i][j][1][1])*mp.sqrt(5)))
                                      for c,g in zip(controls,exact,strict=True))
            full=mp.expm(generator)
            if event['kind']=='J':
                order=(0,2,1,3) if event['targets']==[0,1] else (0,1,2,3)
                gate=mp.matrix([[full[i,j]for j in order]for i in order])
            else:
                pair=mp.matrix([[full[i,j]for j in range(2)]for i in range(2)])
                gate=mp.matrix(4)
                bit=event['targets'][0]
                for i in range(4):
                    for j in range(4):
                        if (i>>(1-bit))%2==(j>>(1-bit))%2:
                            gate[i,j]=pair[(i>>bit)%2,(j>>bit)%2]
            output=gate*output
        for witness in row['successful_work_number_leakage_entries']:
            value=output[witness['output'],witness['input']]
            value=mp.re(value)if witness['component']=='real' else mp.im(value)
            lo,hi=map(F,witness['interval'])
            assert number(lo)<=value<=number(hi)
