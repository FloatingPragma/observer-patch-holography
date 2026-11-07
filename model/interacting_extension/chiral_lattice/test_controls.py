"""Independent finite stencils, measurement probabilities and cutoff controls."""
import itertools

import numpy as np
import pytest
from scipy.linalg import expm

from . import certificate, light_branch as light, rotor_response as flux
from .source_data import load


def test_independent_ring_and_complete_gauss_basis():
    for charge in (-4,-3,1,2,6):
        model=flux.build_flux_ring(charge,3)
        labels=model["labels"]
        expected=np.zeros_like(model["H"])
        for j,(x,ell) in enumerate(labels):
            expected[j,j]=sum((ell+charge*(e<x))**2 for e in range(4))/2
            for i,(y,m) in enumerate(labels):
                hop=(y==(x+1)%4 and m==ell+charge*(x==3))
                loop=(y==x and m==ell+1)
                if hop:expected[i,j]-=1;expected[j,i]-=1
                if loop:expected[i,j]-=.5;expected[j,i]-=.5
        np.testing.assert_array_equal(model["H"],expected)
        assert not np.any(model["Gauss"])
        receipt=flux.operator_receipt(model)
        assert receipt["T_unitarity_defect_norm"]==1
        assert receipt["W_unitarity_defect_norm"]==1


def test_record_law_uses_complete_density_instrument():
    model=flux.build_flux_ring(-3,3)
    p=model["position_projectors"][0]
    times=(.07,.11,.03)
    actual=flux.complete_binary_records(model["H"],model["initial_state"],p,times)
    expected=[]
    for history in itertools.product((0,1),repeat=3):
        rho=np.outer(model["initial_state"],model["initial_state"])
        for outcome,t in zip(history,times,strict=True):
            u=expm(-1j*model["H"]*t);effect=p if outcome else np.eye(len(p))-p
            rho=effect@u@rho@u.conj().T@effect
        expected.append(np.trace(rho).real)
    np.testing.assert_allclose([r["probability"] for r in actual["records"]],expected,atol=2e-14)
    assert sum(expected)==pytest.approx(1,abs=2e-14)


def test_cutoff_bound_for_superposition_and_reference():
    """A two-column purification tests a stronger input than the localized witness."""
    q=6;small=flux.build_flux_ring(q,5);large=flux.build_flux_ring(q,17)
    rng=np.random.default_rng(3107)
    psi=np.zeros((len(small["labels"]),2),complex)
    for i,(_,ell) in enumerate(small["labels"]):
        if abs(ell)<=1:psi[i]=rng.normal(size=2)+1j*rng.normal(size=2)
    psi/=np.linalg.norm(psi)
    inclusion=np.zeros((len(large["labels"]),len(small["labels"])))
    lookup={v:i for i,v in enumerate(large["labels"])}
    for j,v in enumerate(small["labels"]):inclusion[lookup[v],j]=1
    t=.04
    error=np.linalg.norm(expm(-1j*large["H"]*t)@inclusion@psi-inclusion@expm(-1j*small["H"]*t)@psi)
    a=flux.dyson_bound(5,t,q,initial_flux_support=1)["state_vector_error_bound"]
    b=flux.dyson_bound(17,t,q,initial_flux_support=1)["state_vector_error_bound"]
    assert error<=a+b


def test_infinite_ground_bracket_contains_larger_cutoff_values():
    for q in (-4,1,6):
        exact=np.linalg.eigvalsh(flux.build_flux_ring(q,14)["H"])[0]
        for cutoff in (1,3,5):
            model=flux.build_flux_ring(q,cutoff)
            finite=np.linalg.eigvalsh(model["H"])[0]
            bracket=flux.infinite_ground_bracket(model,finite)
            assert bracket["infinite_ground_lower_bound"]-1e-12<=exact<=finite+1e-12


def test_free_cone_all_corners_and_operator_error_band():
    gamma=load()["gamma"];gs,_=gamma;beta=gs[3]
    for corner in itertools.product((0,np.pi),repeat=3):
        h,_,data=light.spatial_block(corner,1,gamma)
        assert data["energy"]==pytest.approx(0 if not any(corner) else 2,abs=1e-14)
        np.testing.assert_allclose(h@h,data["energy"]**2*np.eye(4),atol=3e-14)
    rng=np.random.default_rng(6701)
    for _ in range(50):
        p=rng.normal(size=3);a=.25/(1+np.linalg.norm(p))
        h,_,d=light.spatial_block(p,a,gamma)
        expected=sum(1j*beta@gs[i]*p[i] for i in range(3))
        assert np.linalg.norm(h-expected,2)<=a*(p@p)+2e-14
        u=1-np.cos(a*p)
        assert d["R"]**2==pytest.approx(1+2*sum(u[i]*u[j] for i in range(3) for j in range(i)),abs=2e-14)


def test_isolated_one_bond_sign_change_is_not_spatial_propagation():
    source=load();gamma,g5=source["gamma"]
    plus=(np.eye(4)+gamma[0])/2;minus=(np.eye(4)-gamma[0])/2
    chirality=np.kron(np.eye(2),g5)
    for dimension in (1,2,3,6):
        s=np.sqrt(dimension/6)
        hop=np.block([[np.zeros((4,4)),s*plus],[s*minus,np.zeros((4,4))]])
        np.testing.assert_allclose(hop@hop,0,atol=1e-15)
        for b in (-1,1):
            kernel=chirality@(b*np.eye(8)-hop)
            eigen,v=np.linalg.eigh(kernel)
            overlap=np.eye(8)+chirality@(v*np.sign(eigen))@v.conj().T
            expected=2+4*b/np.sqrt(4+dimension/6)
            np.testing.assert_allclose(np.linalg.svd(overlap,compute_uv=False)**2,expected,atol=1e-14)
            # On the 7,080-vector zero-transport complement, the sign alone
            # changes all overlap eigenvalues from two to zero.
            assert 1+np.sign(b)==(2 if b==1 else 0)


@pytest.mark.parametrize("kwargs",[{"q":1.5},{"q":True},{"K":2.5},{"K":True},{"L":3.5},
                                    {"electric":0},{"electric":float("nan")},{"loop":float("inf")}])
def test_fractional_truncation_and_nonfinite_couplings_rejected(kwargs):
    args=dict(q=1,K=2);args.update(kwargs)
    with pytest.raises(ValueError):flux.build_flux_ring(**args)


@pytest.mark.parametrize("kwargs",[{"K":1,"initial_flux_support":2},{"K":-1},{"K":1,"total_time":-1},
                                    {"K":1,"total_time":float("nan")}])
def test_cutoff_estimate_requires_its_hypotheses(kwargs):
    args=dict(K=1,total_time=.1);args.update(kwargs)
    with pytest.raises(ValueError):flux.dyson_bound(**args)


def test_receipt_replays_without_platform_residuals():
    import json
    actual=certificate.report()
    assert actual==json.loads(certificate.OUTPUT.read_text())
    assert "maximum_projection_error" not in json.dumps(actual)
