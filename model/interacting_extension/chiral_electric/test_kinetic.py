"""Independent projector and positive-parent controls for the electric join."""
import json

import numpy as np
import pytest
from scipy.linalg import block_diag, eigh

from . import certificate, covariant_kinetic as kinetic


def rotating_frame(theta):
    v=np.array([[np.cos(theta)],[np.sin(theta)]])
    dv=np.array([[-np.sin(theta)],[np.cos(theta)]])
    return {"V":v,"dV":dv,"P":v@v.T,"dP":dv@v.T+v@dv.T}


def test_connection_and_geometric_term_on_exact_rotating_line():
    d=rotating_frame(.37);v,dv,p,dp=(d[k] for k in ("V","dV","P","dP"))
    c=p@dp-dp@p
    np.testing.assert_allclose(dp+c@p-p@c,0,atol=2e-16)
    np.testing.assert_allclose(dv+c@v,0,atol=2e-16)
    np.testing.assert_allclose(dv.T@(np.eye(2)-p)@dv,[[1]],atol=2e-16)
    assert np.linalg.norm(dp-c@p+p@c)>1


def test_finite_difference_of_geometric_kinetics_is_indefinite():
    grid=kinetic.grid_model(rotating_frame,5,halfwidth=.7,inertia=1)
    eigen=eigh(grid["compressed"]-grid["parallel"],eigvals_only=True)
    assert eigen[0]<-.1 and eigen[-1]>.1
    # Compression is positive on every retained vector.
    assert eigh(grid["compressed"],eigvals_only=True)[0]>0


def test_compression_frame_covariance_and_complete_parent():
    grid=kinetic.grid_model(rotating_frame,5,halfwidth=.7,inertia=1)
    qs=block_diag(*(np.array([[np.exp(1j*j/3)]]) for j in range(5)))
    changed=grid["V"]@qs
    np.testing.assert_allclose(changed.conj().T@grid["parent"]@changed,
                               qs.conj().T@grid["compressed"]@qs,atol=1e-13)
    # The two-dimensional fibers retain a five-dimensional complement.
    assert grid["parent"].shape==(10,10)
    result,_=kinetic.complete_penalty(grid,count=4)
    for row in result["runs"]:
        assert row["maximum_difference_from_compression"]<=row["analytic_Schur_bound_float64"]+1e-12


def test_random_positive_parent_schur_bracket():
    rng=np.random.default_rng(620620)
    for _ in range(20):
        x=rng.normal(size=(13,13))+1j*rng.normal(size=(13,13))
        k=x.conj().T@x
        a=eigh(k[:5,:5],eigvals_only=True)
        penalty=2*np.linalg.norm(k,2)
        full=k+penalty*np.diag([0]*5+[1]*8)
        ev=eigh(full,eigvals_only=True)[:5]
        error=a-ev
        bound=np.linalg.norm(k[:5,5:],2)**2/(penalty-a[-1])
        assert error.min()>-1e-11 and error.max()<=bound+1e-11


@pytest.mark.parametrize("kwargs",[{"points":True},{"points":0},{"points":2.5},
                                    {"points":3,"halfwidth":0},{"points":3,"inertia":float("nan")}])
def test_grid_domain_rejected(kwargs):
    with pytest.raises(ValueError):kinetic.grid_model(rotating_frame,**kwargs)


def test_source_projector_and_complete_receipt():
    actual=certificate.report()
    assert actual==json.loads(certificate.OUTPUT.read_text())
    assert actual["complete_parent"]["dimension"]==576
    assert actual["scope"]["interval_arithmetic_certificate"] is False
