"""Literal one-time loading, repeated finite ages and an explicitly occupied capture head."""
import math

import numpy as np
import pytest

from .. import formation
from . import clock_cycle, clock_factory
from .test_clock_loader_literal import _launch, _routing


SHAPE = (5,5,3,3)
ROOT = (2,2,1,1)
LAP = 8192
MU = np.array([128,129])*math.pi/4096
EQUAL = np.full(8,1/math.sqrt(8))


def _literal_step(wave, perm, age):
    bright = wave.mean(axis=-1,keepdims=True)
    coin = (np.exp(-1j*MU[:,None,None])*bright
            -np.exp(1j*MU[:,None,None])*(wave-bright))
    flipped = coin.reshape(2,-1)[:,perm].reshape(wave.shape)
    if age >= 127:
        return flipped
    theta = math.pi*(1-math.cos(math.pi*(age+1)/128))/4
    e = np.exp(-2j*theta)
    return ((1+e)*coin+(1-e)*flipped)/2


def _birth(route):
    refused,pair = _launch(SHAPE,ROOT,route)
    assert abs(refused["unmarked"]-.5) < 2e-14
    # These two half-weight rail columns may carry distinct untouched reference states.
    return np.stack([pair.diagonal().reshape(-1,8)/math.sqrt(2)]*2)


def _pulse(wave, gates, age):
    if age in gates:
        return np.einsum("rij,rnj->rni",gates[age],wave)
    return wave


def _head(wave):
    site = np.ravel_multi_index(ROOT,SHAPE)
    bright = wave[:,site] @ EQUAL
    keep = wave.copy()
    keep[:,site] += (math.cos(math.pi/4)-1)*bright[:,None]*EQUAL
    return keep, abs(bright)**2/2


def test_real_space_dark_storage_from_a_single_ordinary_birth():
    route = _routing(SHAPE)
    perm = route.indices
    wave = _birth(route)
    gates = clock_factory.uniform_schedule()
    reference = None
    checkpoints = {128+j*LAP:j for j in range(5)}
    for tick in range(max(checkpoints)+1):
        # The commissioning apparatus is disarmed after service 35.
        if 1 <= tick <= 35:
            wave = _pulse(wave,gates,tick)
        if tick in checkpoints:
            n = checkpoints[tick]
            if not n:
                reference = wave.copy()
            reflected = reference.reshape(2,-1)[:,perm].reshape(reference.shape)
            expected = ((1+(-1j)**n)*reference+(1-(-1j)**n)*reflected)/2
            np.testing.assert_allclose(wave,expected,atol=2e-11,rtol=0)
            chart_reference = reference.reshape(2,*SHAPE,8)
            native = clock_cycle.transfer(chart_reference,128,n*LAP)
            np.testing.assert_allclose(native.reshape(wave.shape),wave,atol=2e-11,rtol=0)
        if tick and tick%LAP == 0:
            for width in (64,128,4096):
                assert tick%width == 0
            # The original site wire has an occupied head. A routed wire
            # can arrive elsewhere, so also check zero bright weight everywhere.
            assert np.max(abs(wave @ EQUAL)) < 1e-11
            wave,release = _head(wave)
            assert max(release) < 1e-20
        if tick < max(checkpoints):
            wave = _literal_step(wave,perm,tick%LAP)
    # The full state returns, including arbitrary rail/reference phases.
    np.testing.assert_allclose(wave,reference,atol=2e-11,rtol=0)
    np.testing.assert_array_equal(formation.birth_ages(2*LAP+1,np.array([0,1,LAP-1,LAP+1])),
                                  [1,0,2,0])


def test_head_on_a_bright_or_unaligned_material_input():
    field = np.zeros((2,*SHAPE,8),complex)
    field[(slice(None),*ROOT)] = np.array([1/math.sqrt(3),1j*math.sqrt(2/3)])[:,None]*EQUAL
    flat = field.reshape(2,-1,8)
    expected, release = _head(flat)
    kept, measured = clock_cycle.return_probe(field,ROOT)
    np.testing.assert_allclose(kept.reshape(flat.shape),expected,atol=1e-14,rtol=0)
    np.testing.assert_allclose(measured,release,atol=1e-14,rtol=0)
    np.testing.assert_allclose(release,[1/6,1/3],atol=1e-14,rtol=0)
    other = (0,0,0,0)
    kept,measured = clock_cycle.return_probe(field,other)
    np.testing.assert_array_equal(kept,field)
    np.testing.assert_array_equal(measured,[0,0])


def test_an_age_table_rearmed_at_the_wrap_changes_the_carried_pair():
    route = _routing(SHAPE)
    perm = route.indices
    gates = clock_factory.uniform_schedule()
    wave = _birth(route)
    for age in range(35):
        wave = _pulse(_literal_step(wave,perm,age),gates,age+1)
    for age in range(35,LAP):
        wave = _literal_step(wave,perm,age)
    reference,rearmed = wave.copy(),wave.copy()
    for age in range(35):
        reference = _literal_step(reference,perm,age)
        rearmed = _pulse(_literal_step(rearmed,perm,age),gates,age+1)
    # Compare local native capture-active power, independently of global phase.
    assert np.linalg.norm(reference @ EQUAL) < 1e-11
    assert np.linalg.norm(rearmed @ EQUAL) > 1e-5


def test_local_digit_is_reversible_on_coherent_lap_blocks():
    rng = np.random.default_rng(481)
    wave = rng.normal(size=(5,2,*SHAPE,8))+1j*rng.normal(size=(5,2,*SHAPE,8))
    route = _routing(SHAPE).indices
    plates = clock_factory.uniform_schedule()
    for age in (0,26,34,LAP-1):
        moved = np.empty_like(wave)
        reconstructed = np.empty_like(wave)
        for lap in range(5):
            field,new_age,new_lap = clock_cycle.program_step(wave[lap],age,lap,5,plates)
            assert new_age==(age+1)%LAP
            assert new_lap==(lap+int(age==LAP-1))%5
            moved[new_lap] = field
        for new_lap,field in enumerate(moved):
            lap = (new_lap-int(age==LAP-1))%5
            literal = _literal_step(wave[lap].reshape(2,-1,8),route,age)
            if new_lap==0:
                literal = _pulse(literal,plates,new_age)
            np.testing.assert_allclose(field.reshape(literal.shape),literal,rtol=0,atol=4e-14)
            back = field.reshape(2,-1,8)
            if new_lap==0 and new_age in plates:
                back = np.einsum("rji,rnj->rni",plates[new_age].conj(),back)
            flip = back.reshape(2,-1)[:,route].reshape(back.shape)
            theta = math.pi*(1-math.cos(math.pi*min(age+1,128)/128))/4
            e = np.exp(2j*theta)
            back = ((1+e)*back+(1-e)*flip)/2
            bright = back.mean(axis=-1,keepdims=True)
            back = (np.exp(1j*MU[:,None,None])*bright
                    -np.exp(-1j*MU[:,None,None])*(back-bright))
            reconstructed[lap] = back.reshape(2,*SHAPE,8)
        np.testing.assert_allclose(reconstructed,wave,atol=5e-14,rtol=0)


def test_register_arithmetic_large_and_small_timestamp_types():
    for dtype in (np.int8,np.uint8,np.int32,np.uint64):
        np.testing.assert_array_equal(formation.birth_ages(10**60+1,np.array([0,1,2],dtype=dtype)),
                                      [(10**60+1-j)%LAP for j in range(3)])


def test_formation_run_uses_separate_cyclic_bound_ages(monkeypatch):
    # Use a short age register to cross two wraps on the ordinary entry point.
    seen, batches = [], []
    native_angle, native_coin = formation.hop_angle, formation.contact_coin
    def angle(age, maturation):
        seen.append(age)
        return native_angle(age,maturation)
    def coin(fields, phase):
        batches.append(tuple(seen))
        seen.clear()
        return native_coin(fields,phase)
    monkeypatch.setattr(formation,"AGE_PERIOD",4)
    monkeypatch.setattr(formation,"hop_angle",angle)
    monkeypatch.setattr(formation,"contact_coin",coin)
    formation.formation_run(size=8,updates=10,orders=(1,-1),support=1,law_key="torus")
    assert len(batches)==10 and len(batches[0])>0
    previous = ()
    for batch in batches:
        expected = tuple((j+1)%4 for j in previous)
        assert batch[:len(previous)]==expected
        assert batch[len(previous):]==(0,)*(len(batch)-len(previous))
        previous = batch


def test_complete_return_and_controller_census():
    result = clock_cycle.experiment()
    assert result["failures"]==[]
    for row in result["eligible_returns"]:
        assert row["wire_at_occupied_head"]==[64,4096]
    rows = result["records_conditional_on_first_write"]
    assert [r["elapsed"] for r in rows]==[128+j*8192 for j in range(5)]
    for row in rows:
        ports = row["root_ports"]
        assert len(ports)==4 and min(ports)>=0 and row["born_root_miss"]>=0
        np.testing.assert_allclose(sum(ports)+row["born_root_miss"],1,atol=3e-13,rtol=0)
    np.testing.assert_allclose(rows[0]["root_ports"],rows[4]["root_ports"],atol=2e-13,rtol=0)


def test_capture_preserves_digit_and_reference_through_discharge_and_rebirth():
    # Five program digits, two rails, three untouched reference components.
    rng = np.random.default_rng(35)
    initial = rng.normal(size=(5,2,3,9))+1j*rng.normal(size=(5,2,3,9))
    a,b = clock_cycle.capture_vertex(initial[...,0],initial[...,1:])
    m = np.eye(9)
    c = 1/math.sqrt(2)
    m[0,0],m[0,1:],m[1:,0] = c,-c*EQUAL,c*EQUAL
    m[1:,1:] += (c-1)*np.outer(EQUAL,EQUAL)
    np.testing.assert_allclose(np.concatenate([a[...,None],b],axis=-1),initial@m.T,atol=3e-15,rtol=0)
    back, matched = clock_cycle.capture_vertex(a,b,inverse=True)
    np.testing.assert_allclose(back,initial[...,0],atol=3e-15,rtol=0)
    np.testing.assert_allclose(matched,initial[...,1:],atol=3e-15,rtol=0)


def test_passive_darkness_check_preserves_rare_branch_directions():
    field = np.zeros((2,*SHAPE,8),complex)
    field[(0,*ROOT)] = EQUAL
    with pytest.raises(ValueError, match="brightness"):
        clock_cycle.transfer(1e-200*field,128,LAP)
    # This input is dark at the root, but bright at its two routed endpoints.
    field[(0,*ROOT)] = [1j,-1j,0,0,0,0,0,0]
    with pytest.raises(ValueError, match="brightness"):
        clock_cycle.transfer(1e-200*field,128,LAP)
    field[0] = clock_factory.packet(SHAPE,ROOT)
    actual = clock_cycle.transfer(1e-200*field,128,LAP)/1e-200
    np.testing.assert_allclose(actual,clock_cycle.transfer(field,128,LAP),rtol=0,atol=1e-14)


def test_material_controllers_never_drop_an_extra_rail():
    field = np.zeros((3,*SHAPE,8),complex)
    field[(2,*ROOT)] = EQUAL
    with pytest.raises(ValueError,match="two finite rail fields"):
        clock_factory.advance(field,0,{})
    with pytest.raises(ValueError,match="two finite rail fields"):
        clock_cycle.program_step(field,LAP-1,0,5,{})
