"""Source, instruments, spectral bounds and failure checks for a write-formed pair."""
import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1 import reader
from .dwell import counts_for, plan, read, returns
from .dwell_checks import PhysicalFlight
from .dwell_source import Station, local_blocks
from .dwell_scaling import coefficient, memory_charges, secular


@pytest.mark.parametrize("kappa",[8,0,-8])
def test_native_source_and_inverse_echo_receiver(kappa):
    U=law.echo()
    native,control=local_blocks(kappa)
    seed=np.eye(99)[:,law.SRC.index((2,0))]
    raw_signs=np.column_stack([U[:,law.SRC.index(word)] for word in [(2,0),(0,2)]])
    expected=raw_signs@np.array([1,-1j])/np.sqrt(2)
    second=U@expm(-1j*1.5*np.pi*native)@seed
    assert abs(np.vdot(second,expected))==pytest.approx(1,abs=2e-13)
    sign=np.array([1.,-1.])/np.sqrt(2)
    assert np.linalg.norm(native@raw_signs)<1e-12
    np.testing.assert_allclose(control@(raw_signs@sign),(4*kappa+1)*raw_signs@sign,atol=1e-12)
    # Read the restored double's endpoint. In all other binary words both sign
    # branches have a quiet outcome, including any relative phase between them.
    labels=tuple(((c,s),word) for c,s in law.LABELS for word in law.SRC)
    for word in ((0,0),(0,1),(1,0),(1,1)):
        before=np.zeros((99,2))
        for j,s in enumerate((2,-2)):
            for c in (0,1):
                before[9+labels.index(((c,s),word)),j]=1/np.sqrt(2)
        after=U@before
        blocks=after.reshape(11,9,2)
        effects=[blocks[:,law.SRC.index(d),:].T@blocks[:,law.SRC.index(d),:] for d in ((2,0),(0,2))]
        for j,effect in enumerate(effects):
            wanted=np.diag(np.eye(2)[j]) if word==(1,1) else np.zeros((2,2))
            np.testing.assert_allclose(effect,wanted,atol=1e-14)
    # The public agreement reader can commit either of these sharp flags.
    assert reader.gate(1,0)>>6==29 and reader.gate(0,0)>>6==1
    # After the first hit, the two exterior units are on the tag. Echoing
    # (port 0, empty port 5) and recording (1,1) distinguishes their old sign.
    basis=U[:,[law.SRC.index((n,0)) for n in (2,0)]]
    pair_effect=np.arange(law.DIM)%9==law.SRC.index((1,1))
    for pointer,mask in ((1,pair_effect),(0,~pair_effect)):
        public=(reader.gate(pointer,0)>>8)&7
        expected=np.diag([1,0] if public==7 else [0,1])
        np.testing.assert_allclose(basis[mask].T@basis[mask],expected,atol=1e-14)


def test_local_reference_and_stationary_memory_charge():
    ranks=memory_charges()
    assert ranks=={"0_1":2,"1_1":2,"1_-1":2,"2_1":2,"2_-1":2}
    cfg=plan()
    for level in cfg["geometry_levels"]:
        b=Station(level,cfg["tag"])
        assert len(b.cuts)==18 and b.clock_distance>=5
        assert b.isolation(8,0)["survival"]>.9766
    is_written=np.diag(np.r_[np.zeros(9),np.ones(90)])
    native,_=local_blocks(8)
    np.testing.assert_allclose(native@is_written,is_written@native,atol=1e-13)
    # The supplied echo can create the mark, unlike the free generator.
    assert np.linalg.norm(law.echo()@is_written-is_written@law.echo())>1


def test_projected_operators_and_parameter_free_coefficient():
    b=Station(0,(0,1))
    b.flight()
    literal=PhysicalFlight(b)
    assert str(coefficient(b))=="29/180"
    for k in (-8,0,8):
        for j in (0,1):
            for arm in ("native","always_ticking"):
                assert (b.generator(k,j,arm)-literal.generator(k,j,arm)).count_nonzero()==0
    # Preserve the native binding surface while allowing the occupied hand to tick.
    delta=literal.generator(8,1,"binding_only")-literal.generator(8,0,"native")
    assert (delta-sparse.eye(b.sector.P)).count_nonzero()==0


def test_polynomial_phase_and_secular_crossing():
    H=np.array([[-12.,-.5,-.5],[-.5,-3.,-.5],[-.5,-.5,8.]])
    times=np.array([0.,.3,2.,96.])
    value,x=secular(sparse.csr_matrix(H),0,8.)
    eig,basis=np.linalg.eigh(H)
    assert value==pytest.approx(eig[0],abs=1e-12)
    assert 1/(1+x@x)==pytest.approx(basis[0,0]**2,abs=1e-12)
    sharp=np.array([1.,0.,0.])
    direct=np.array([sharp@expm(-1j*t*H)@sharp for t in times])
    np.testing.assert_allclose(returns(sparse.csr_matrix(H),0,times),direct,atol=5e-12)
    np.testing.assert_equal(returns(sparse.csr_matrix(H),0,np.zeros(2)),np.ones(2))


def test_quiet_and_impossible_histograms():
    c=np.zeros((2,3,3,4),np.int32)
    c[:,:,2,0]=64
    assert read(c,64,.01,0.)["rate_interval"] is None
    for bad in [c[:,:,:,:2],c*.5,c*0]:
        with pytest.raises(ValueError,match="invalid joint"):
            read(bad,64,.01,0.)
    c[0,0,1,:2]=[-1,1]
    with pytest.raises(ValueError,match="invalid joint"):
        read(c,64,.01,0.)
    with pytest.raises(ValueError,match="invalid joint"):
        read(c*0,0,.01,0.)
    # No signal return in either branch gives all the weight to quiet.
    hands=np.ones((2,4))/4
    out=counts_for(np.zeros(2),np.zeros(2),hands,64,np.random.default_rng(5))
    assert out[:,:,:2].sum()==0 and np.all(out.sum(axis=(2,3))==64)
    for clock in (hands*2,hands[:,0],hands*np.nan,np.tile([-1,1,1,0],(2,1))):
        with pytest.raises(ValueError,match="clock probabilities"):
            counts_for(np.zeros(2),np.ones(2),clock,64,np.random.default_rng(5))
    for bad in (False,0,2**31,1.5):
        with pytest.raises(ValueError,match="trial count"):
            counts_for(np.zeros(2),np.ones(2),hands,bad,np.random.default_rng(5))
    for amp in (np.ones(3),np.full(2,1.001),np.full(2,np.inf)):
        with pytest.raises(ValueError,match="return amplitudes"):
            counts_for(np.zeros(2),amp,hands,64,np.random.default_rng(5))
