from copy import deepcopy
from fractions import Fraction as F
import json
import math
import numpy as np
import pytest

from . import intervals as iv, synthesis as s, packet


@pytest.fixture(scope='module')
def short():
    word=s.pair_word(0,65);s.certify_word(word)
    return word


@pytest.fixture(scope='module')
def general():
    name,vector,strength=packet.examples()[1]
    word=s.rank_one_word(vector,strength);s.certify_word(word)
    return word,vector,strength


def test_rational_pi_and_rounding_domain():
    lower,upper=iv.pi_interval()
    assert lower<upper and upper-lower<F(1,10**100)
    assert max(abs(F(math.pi)-lower),abs(F(math.pi)-upper))<F(1,10**15)
    assert max(abs(F(iv.pi_grid(),iv.AQ)-lower),abs(F(iv.pi_grid(),iv.AQ)-upper))<F(1,iv.AQ)
    assert iv.binary_exact(17,2**50)==float(F(17,2**50))
    with pytest.raises(ArithmeticError):iv.binary_exact(2**54+1,2**55)


@pytest.mark.parametrize('value',[F(),F(1,3),F(-7,8),F(4),F(-4)])
def test_directed_trigonometric_enclosures(value):
    sine,cosine=iv.sincos(value)
    # The interval proof is integer arithmetic; libm supplies only an
    # independent low-precision sanity check here.
    for interval,expected in ((sine,math.sin(float(value))),(cosine,math.cos(float(value)))):
        assert interval[0]<=interval[1]
        assert abs(float(F(sum(interval),2*iv.Q))-expected)<2e-15
        assert interval[1]-interval[0]<10**7
    with pytest.raises(ValueError):iv.sincos(F(5))


@pytest.mark.parametrize('left,right',[(F(9),F(16)),(F(1),F(1)),(F(),F(7)),
                                      (F(7),F()),(F(1,2**60),F(1))])
def test_split_grid_and_endpoint_controls(left,right):
    angle=F(iv.split_angle(left,right),iv.AQ)
    expected=math.atan2(math.sqrt(float(right)),math.sqrt(float(left)))
    assert abs(float(angle)-expected)<=1.01/iv.AQ
    if left and right:
        _,c=iv.sincos(2*angle);wanted=(left-right)/(left+right)
        assert max(abs(F(v,iv.Q)-wanted) for v in c)<F(3,iv.AQ)
    with pytest.raises(ValueError):iv.split_angle(F(-1),F(1))


@pytest.mark.parametrize('x,y',[(F(3),F(4)),(F(-3),F(4)),(F(-3),F(-4)),
                              (F(3),F(-4)),(F(),F(-1)),(F(-1),F())])
def test_complex_quadrants_preserve_phase(x,y):
    angle=F(iv.phase_angle(x,y),iv.AQ)
    assert abs(float(angle)-math.atan2(float(y),float(x)))<=2.01/iv.AQ


def test_unequal_weight_effect_is_constructed_on_full_carrier(short):
    assert len(short)==80
    array=s.compiler.execute(short,9)
    wanted=np.zeros(66,complex);wanted[0]=.6;wanted[65]=.8
    row=array[s.tomo.codes()[0],list(s.tomo.codes())]
    assert np.linalg.norm(row-wanted)<1e-8
    effect=np.outer(row.conj(),row)
    assert np.linalg.norm(effect-np.outer(wanted,wanted))<1e-8
    # Unequal amplitudes cannot be a positive sum of balanced one-pair atoms.
    assert F(9,25)!=F(16,25) and F(12,25)**2==F(9,25)*F(16,25)
    assert array.shape==(512,512)


def test_general_all_coordinate_rank_one_and_source_budget(general):
    word,vector,strength=general
    norm=sum((x*x+y*y for x,y in vector),F())
    wanted=np.array([complex(float(x),-float(y)) for x,y in vector])*np.sqrt(float(strength/norm))
    actual=s.selected_row(word)[list(s.tomo.codes())]
    assert np.linalg.norm(actual-wanted)<1e-8
    assert len(word)<=23579<65536
    assert s.interpreter.source_word(word,9)['declared_resources']['native_pulses']==len(word)
    assert F(s.bounds(len(word),general=True)['complete_halfdiamond_upper'])<F(1,10**6)


def test_zero_strength_keeps_the_zero_effect_without_renormalization():
    vector=tuple((F(i==7),F()) for i in range(66))
    word=s.rank_one_word(vector,F())
    actual=s.selected_row(word)[list(s.tomo.codes())]
    assert np.linalg.norm(actual)<1e-8
    assert len(word)>23000
    assert sum(e['angle']==0 for e in word)>0


def test_unit_strength_uses_full_complex_target_without_extra_qubits():
    _,vector,_=packet.examples()[0]
    word=s.rank_one_word(vector,F(1),output=65)
    row=s.selected_row(word,65)[list(s.tomo.codes())]
    wanted=np.zeros(66,complex);wanted[0]=.6;wanted[65]=-.8j
    assert np.linalg.norm(row-wanted)<1e-8
    assert len(word)<=23579 and max(max(e['targets']) for e in word)==8


def test_final_address_routing_is_required(general):
    word,vector,strength=general
    omitted=word[:12291+11270]
    assert len(omitted)<len(word)
    correct=s.selected_row(word)[list(s.tomo.codes())]
    wrong=s.selected_row(omitted)[list(s.tomo.codes())]
    assert np.linalg.norm(correct-wrong)>.1


def test_generic_host_stop_and_full_outcome_catalogue(general):
    word,_,_=general
    density=np.array([[F(2,3)]],dtype=object)
    out=s.interpreter.execute_word(s.interpreter.source_word(word,9),density)
    assert out['status']=='RESOURCE_STOP' and out['pulses_executed']==0
    assert np.array_equal(out['blocks'][('RESOURCE_STOP',)],density)
    catalogue=s.outcome_catalogue(word)
    assert catalogue['terminal_padding_kraus_rows']==446
    assert catalogue['observed_terminal_padding_abort_flags']==1
    assert catalogue['native_tree_raw_words']==25**11
    assert catalogue['nonzero_native_tree_words']==66*2**11
    assert catalogue['processor_first_failure_labels']==sum(6-2**len(e['targets']) for e in word)
    assert catalogue['native_attempts']==11


def test_selected_record_keeps_entangled_reference_and_padding(short):
    rng=np.random.default_rng(221)
    factor=rng.normal(size=(66,2))+1j*rng.normal(size=(66,2));factor/=np.linalg.norm(factor)
    blocks=s.coarse_output(short,factor)
    assert len(blocks)==67 and all(b.shape==(2,2) for b in blocks)
    expected=.6*factor[0]+.8*factor[65]
    ideal=np.outer(expected,expected.conj())
    assert np.linalg.norm(blocks[0]-ideal,'nuc')<2*float(F(s.bounds(len(short),general=False)['complete_halfdiamond_upper']))
    assert np.trace(blocks[-1]).real>1e-6
    assert abs(sum(np.trace(b).real for b in blocks)-1)<1e-8


def test_fine_native_and_padding_rows_remain_addressable(short):
    run=s.detector.run(0,(0,)*11)
    weight,terminal,row,records=s.fine_success_row(short,run['slots'])
    assert weight==run['weight'] and terminal==run['final'] and records==run['records']
    assert np.allclose(row,s.selected_row(short,0)[list(s.tomo.codes())])
    slots=list(run['slots']);slots[3]=next(i for i,x in enumerate(s.detector.kernel(s.detector.tree()[3])) if not x.pairs)
    weight,terminal,row,records=s.fine_success_row(short,slots)
    assert terminal is None and not np.any(row) and len(records)==11
    p=next(i for i in range(512) if i not in s.tomo.codes())
    assert s.padding_kraus_row(short,p).shape==(66,)
    with pytest.raises(ValueError):s.padding_kraus_row(short,s.tomo.codes()[0])


def test_first_failure_export_is_unnormalized_and_retains_reference(short,monkeypatch):
    factor=np.zeros((66,2),complex);factor[0]=(1,2j)
    actual=s.processor_failure(short,0,2,factor)
    assert actual.shape==(512,2) and np.linalg.norm(actual)<1e-8
    original=s.compiler.native_pulse
    def leaky(kind,angle):
        result=original(kind,angle);result[2,:2]+=1/1000
        return result
    monkeypatch.setattr(s.compiler,'native_pulse',leaky)
    changed=s.processor_failure(short,0,2,factor)
    assert np.linalg.norm(changed)>1e-4 and np.linalg.norm(changed)<.01
    assert np.allclose(changed[:,1],2j*changed[:,0])
    s.certify_primitive.cache_clear()
    with pytest.raises(ArithmeticError,match='six-level'):
        s.certify_primitive(short[0]['kind'],short[0]['angle'])
    s.certify_primitive.cache_clear()


def test_walsh_coefficient_mutation_is_rejected(monkeypatch):
    original=s._walsh
    def broken(values):
        actual=list(original(values));actual[0]+=1/16
        return tuple(actual)
    monkeypatch.setattr(s,'_walsh',broken)
    with pytest.raises(ArithmeticError,match='angle changed'):
        s._state_word([F(1)]+[F()]*511)


def test_input_and_budget_guards():
    good=tuple((F(i==0),F()) for i in range(66))
    for bad in ([],[(F(),F())]*66,[(0.0,F())]*66):
        with pytest.raises(ValueError):s.rank_one_word(bad)
    with pytest.raises(ValueError):s.rank_one_word(good,F(2))
    with pytest.raises(ValueError):s.rank_one_word(good,True)
    with pytest.raises(ValueError):s.rank_one_word(good,F(1,2**65))
    with pytest.raises(ValueError):s.rank_one_word(good,output=True)
    with pytest.raises(ValueError):s.bounds(65536,general=True)
    with pytest.raises(ValueError):s.pair_word(2,1)


def test_complete_bound_includes_failure_and_nonnormalized_gain():
    value=s.bounds(23579,general=True)
    e=F(value['full_word_operator_error_upper']);failure=F(value['all_first_processor_failure_mass_upper'])
    assert failure>0 and F(value['complete_halfdiamond_upper'])==(e*(2+e)+failure)/2
    assert F(value['complete_map_gain_upper'])>1
    assert F(value['fixed_binary_angle_target_error_upper'])==F(11,2**41)


def test_receipt_replay():
    assert packet.verify()


def test_typed_receipt_rejects_boolean_integer_substitution(tmp_path,monkeypatch):
    original=json.loads(packet.OUTPUT.read_text());bad=deepcopy(original)
    bad['scope']['fixed_precision_rank_one_effects_constructed']=1
    path=tmp_path/'forged.json';path.write_text(packet.canonical(bad))
    monkeypatch.setattr(packet,'produce',lambda:original)
    with pytest.raises(ValueError,match='differs'):packet.verify(path)
