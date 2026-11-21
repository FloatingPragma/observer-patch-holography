"""Position-space check of every unfiltered port, birthplace and age allowance."""
import io
import math

import numpy as np

from . import clock_factory, clock_records, clock_stream, clock_unfiltered
from .test_clock_loader_literal import _routing


def test_two_complete_intervals_on_all_original_starts():
    shape = (3,3,3,3)
    sites = np.array(list(np.ndindex(shape)))
    difference = ((sites[:,None,:]-sites[None,:,:])%3).transpose(2,0,1)
    offsets = np.ravel_multi_index(difference,shape)
    with np.load(io.BytesIO(clock_records.snapshot(clock_stream.ARCHIVE,clock_stream.PUBLISHED))) as data:
        hist = data["writes"].reshape(193,-1)
        unmarked = data["unmarked"]
    # Independent direct evolution of the matched diagonal. The native
    # first writes supply distinct input histories at every source site.
    wave = np.zeros((2,len(sites),8),complex)
    wave[:,0] = .25
    perm = _routing(shape).indices
    gates = clock_factory.uniform_schedule()
    phase = np.exp(1j*np.array([128,129])*math.pi/4096)[:,None,None]
    impulse = []
    for age in range(194):
        if age in gates:
            for r in (0,1):
                wave[r] = wave[r] @ gates[age][r].T
        impulse.append(np.stack([np.sum(abs(wave[0]+p*wave[1])**2,axis=-1)/4
                                 for p in (1,1j,-1,-1j)],axis=-1))
        bright = wave.mean(axis=-1,keepdims=True)
        coin = phase.conj()*bright - phase*(wave-bright)
        routed = coin.reshape(2,-1)[:,perm].reshape(wave.shape)
        theta = math.pi*(1-math.cos(math.pi*(age+1)/128))/4 if age<128 else math.pi/2
        e = np.exp(-2j*theta)
        wave = ((1+e)*coin+(1-e)*routed)/2
    for summary,words,residual,bound in clock_unfiltered.experiment():
        t,u = summary["station_ticks"]
        old_hand, full_hand, young_power = [], [], []
        for i,tick in enumerate((t,u)):
            cohorts = np.stack([np.einsum("o,sof->sf",hist[b],impulse[tick-b][offsets])
                                for b in range(tick)])
            all_ports = cohorts.sum(axis=0)
            h = all_ports @ np.array([1,1j,-1,-1j])
            full_hand.append(h)
            old_hand.append(cohorts[:t-127].sum(axis=0) @ np.array([1,1j,-1,-1j]))
            young_power.append(cohorts[t-127:].sum(axis=(0,2)))
            local = all_ports.sum(axis=-1)
            same_public = (sites[:,:3,None]==sites.T[None,:3,:]).all(axis=1)
            same_elsewhere = same_public @ local-local
            other = (~same_public) @ local
            expected = np.column_stack([all_ports,same_elsewhere,other,np.full(len(sites),unmarked[tick])])
            np.testing.assert_allclose(words[i].reshape(-1,7),expected,atol=2e-12,rtol=0)
        turn = np.exp(-2j*math.pi/4096)
        np.testing.assert_allclose(old_hand[1],turn*old_hand[0],atol=2e-12,rtol=0)
        np.testing.assert_allclose(residual.ravel(),full_hand[1]-turn*full_hand[0],atol=2e-12,rtol=0)
        np.testing.assert_allclose(bound.ravel(),np.mean(young_power,axis=0),atol=2e-12,rtol=0)
        assert np.all(abs(residual) <= bound+2e-12)
