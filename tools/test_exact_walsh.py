"""Arithmetic equivalence and source-guard tests for host-only compilation speed."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import random
from types import SimpleNamespace

import pytest

from tools import exact_walsh as w


def direct(values):
    return tuple(sum((x*(-1)**((i&j).bit_count()) for j,x in enumerate(values)),0)
                 for i in range(len(values)))


@pytest.mark.parametrize('n', [1,2,4,8,16,32,64,128,256,512])
def test_exact_signed_integer_and_fraction_transforms(n):
    rng=random.Random(511+n)
    integers=tuple(rng.randrange(-2**70,2**70) for _ in range(n))
    rational=tuple(F(x,(i%7)+1) for i,x in enumerate(integers))
    for values in (integers,rational):
        assert w.transform(values)==direct(values)
        assert w.transform(w.transform(values))==tuple(n*x for x in values)


def test_every_512_coordinate_matches_the_character_matrix():
    for j in range(512):
        basis=tuple(int(i==j) for i in range(512))
        assert w.transform(basis)==tuple((-1)**((i&j).bit_count()) for i in range(512))


@pytest.mark.parametrize('value', [[],[1,2,3],[True,0],[1.0,0],['1',0],[complex(1),0]])
def test_transform_rejects_outside_exact_domain(value):
    with pytest.raises(ValueError): w.transform(value)


def source_module():
    from model.native_effect_synthesis import synthesis
    return synthesis


def adapter(source):
    return w.SynthesisAdapter(source,Path(source.__file__).resolve().parents[2])


def test_forward_preserves_hex_zero_and_error_interfaces():
    source=source_module(); new=w.replacements(source)['_walsh']; rng=random.Random(613)
    assert source._walsh([])==new([])==()
    for n in (1,2,4,8,16,32,64,128,256,512):
        values=tuple(rng.randrange(-source.iv.AQ,source.iv.AQ) for _ in range(n))
        assert tuple(x.hex() for x in source._walsh(values))==tuple(x.hex() for x in new(values))
    for values in ([True],[1.0],[1,2,3],[F(1)]):
        with pytest.raises(ValueError) as old: source._walsh(values)
        with pytest.raises(ValueError) as changed: new(values)
        assert str(old.value)==str(changed.value)
    with pytest.raises(ArithmeticError): new([2**100+1])


def test_full_emitted_phase_and_amplitude_words_and_native_work_are_identical():
    source=source_module(); rng=random.Random(716)
    weights=[F((i%13)+1) for i in range(512)]; total=sum(weights)
    weights=[x/total for x in weights]
    phases=[rng.randrange(-source.iv.AQ,source.iv.AQ) for _ in range(512)]
    old_state=source._state_word(weights); old_phase=source._phase_word(phases)
    old_zero=source._phase_word([0]*512)
    slot=adapter(source);slot.install()
    try:
        state=source._state_word(weights);phase=source._phase_word(phases)
        zero=source._phase_word([0]*512)
        slot.verify()
    finally: slot.restore()
    encode=lambda value:json.dumps(value,sort_keys=True,allow_nan=False,separators=(',',':'))
    assert encode(state)==encode(old_state)
    assert encode(phase)==encode(old_phase)
    assert len(state)==12291 and len(phase)==11270
    # Zero-angle/no-op requests remain actual words and are not pruned.
    assert encode(zero)==encode(old_zero) and len(zero)==11270
    assert sum(x['angle']==0 for x in zero)>0
    assert hashlib.sha256(encode(state+phase).encode()).digest()==hashlib.sha256(encode(old_state+old_phase).encode()).digest()


def test_full_rank_one_word_is_identical():
    source=source_module()
    vector=tuple((F(i%5-2),F((3*i)%7-3)) for i in range(66))
    old=source.rank_one_word(vector,F(2,3),output=17)
    slot=adapter(source);slot.install()
    try:
        actual=source.rank_one_word(vector,F(2,3),output=17)
        slot.verify()
    finally: slot.restore()
    assert json.dumps(actual,sort_keys=True)==json.dumps(old,sort_keys=True)


def test_forward_mutation_still_fails_the_inverse_checks(monkeypatch):
    source=source_module();slot=adapter(source);slot.install()
    try:
        honest=source._walsh
        def corrupt(values):
            result=list(honest(values));result[0]+=1.0;return tuple(result)
        with monkeypatch.context() as m:
            m.setattr(source,'_walsh',corrupt)
            with pytest.raises(ArithmeticError,match='angle changed'):
                source._state_word([F(1)]+[F()]*511)
            with pytest.raises(ArithmeticError,match='lost a phase'):
                source._phase_word([0]*512)
        slot.verify()
    finally: slot.restore()


def test_preexisting_function_substitution_is_rejected(monkeypatch):
    source=source_module()
    monkeypatch.setattr(source,'_walsh',lambda values:())
    with pytest.raises(ValueError,match='already substituted'): adapter(source)


def test_foreign_source_origin_is_rejected(tmp_path):
    fake=SimpleNamespace(__file__=str(tmp_path/'synthesis.py'))
    with pytest.raises(ValueError,match='requires renewed review'):
        w.SynthesisAdapter(fake,tmp_path)


def test_late_function_change_is_rejected_and_restored(monkeypatch):
    source=source_module();old=source._phase_word;slot=adapter(source);slot.install()
    try:
        with monkeypatch.context() as m:
            m.setattr(source,'_phase_word',old)
            with pytest.raises(ValueError,match='changed during replay'): slot.verify()
    finally: slot.restore()
    assert source._phase_word is old
