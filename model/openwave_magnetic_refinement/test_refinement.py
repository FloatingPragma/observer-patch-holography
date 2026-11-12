from fractions import Fraction as F
import numpy as np
import pytest
from scipy.linalg import eigh
from .certified_probe import rational_residual,scalar_mass_gap,uniform_field_bound
from .experiment import HERE,check_witness
from model.openwave_intrinsic_magnetic.solver import dirac


def test_exact_rational_source_residual_and_corrupted_state():
    n=3;a=F(1);N=9;B=F(1,100)
    h,_=dirac(n,float(a),N,np.zeros((N,)*3+(3,)),external_B=float(B))
    w,v=eigh(h.toarray());index=np.flatnonzero(w>0)[0]
    cert=rational_residual(v[:,index],w[index],n,a,B)
    assert F(cert['analytic_matrix_residual_upper'])<F(1,10**11)
    bad=v[:,index].copy();bad[0]+=.01;bad/=np.linalg.norm(bad)
    changed=rational_residual(bad,w[index],n,a,B)
    assert F(changed['analytic_matrix_residual_upper'])>F(1,1000)


def test_all_mode_mass_gap_and_uniform_nonlinear_bound():
    bound=scalar_mass_gap(17,F(1,3),F(1,100))
    assert F(bound['fifth_absolute_D_eigenvalue_lower'])>F(8162,1000)
    dynamic=uniform_field_bound(17,F(1,3),63,F(1,10**14))
    assert F(dynamic['contraction_upper'])<F(1,10**8)
    assert F(dynamic['uniform_kernel_change_upper'])*6000<F(1,10**5)


def test_witness_residual_detects_changed_field(tmp_path):
    source=HERE/'witnesses/stationary_n5.npz';check_witness(source)
    data=np.load(source,allow_pickle=False);copy={k:data[k].copy() for k in data.files}
    copy['field'][10,10,10,0]+=.01
    changed=tmp_path/'changed.npz';np.savez_compressed(changed,**copy)
    with pytest.raises(AssertionError):check_witness(changed)
