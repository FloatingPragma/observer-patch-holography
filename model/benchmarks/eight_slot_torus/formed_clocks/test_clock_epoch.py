import math

import numpy as np
import pytest

from . import clock_epoch as epoch


def test_independent_native_root_word_and_injectivity():
    bank = epoch.likelihoods()
    angle = math.fsum(math.pi*(1-math.cos(math.pi*(j+1)/128))/4 for j in range(35,128))
    totals = []
    for lap in range(4):
        totals.append((1+np.real((-1j)**lap*np.exp(-2j*angle)))/8)
    assert min(totals)>.0157
    assert min(np.diff(sorted(totals)))>.048
    for i,t in enumerate((*epoch.CALIBRATION,*epoch.VISITS)):
        s = totals[(t//8192)%4]
        angle = (36-t)*math.pi/4096
        word = np.array([s*(1+v)/4 for v in (math.cos(angle),math.sin(angle),
                                               -math.cos(angle),-math.sin(angle))])
        np.testing.assert_allclose(bank[0,i,:4],word,atol=2e-14,rtol=0)
        np.testing.assert_allclose(bank[0,i,4:],[.5-s,.5],atol=2e-14,rtol=0)
    assert len({(lap,(36-age)%8192) for lap in range(4)
                for age in range(128,8192,2)})==16128


def test_epoch_decoder_rejects_incomplete_denominator():
    base = np.zeros((4,6),np.int64)
    with pytest.raises(ValueError):
        epoch.decode(base,base[0])


def test_published_epoch_words_and_competing_controls(tmp_path):
    result = epoch.experiment()
    assert result["failed_gates"]==[]
    assert result["launches"]==36*2**38
    assert result["uniform_budget"]["mature_even_epochs"]==16128
    for menu,number in (("native",1),("saturated",4),("dephased",4032)):
        assert [r["candidate_count"] for r in result["words"][menu]]==[number]*8
    damaged = tmp_path/"epoch.npz"
    record = bytearray(epoch.ARCHIVE.read_bytes())
    record[-12] ^= 1
    damaged.write_bytes(record)
    with pytest.raises(ValueError,match="changed"):
        epoch.experiment(damaged)
