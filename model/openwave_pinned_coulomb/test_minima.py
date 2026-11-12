import numpy as np
import sympy as s
import pytest
from .source import native_compression,exact_rows,exact_minimum_difference,profile
from model.openwave_current_coulomb.source import cochains
from model.coherent_defect_wave.chiral_joint.periodic_family.model import yukawa_words
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory


def test_full_native_hopping_compression_and_local_gauge_invariance():
    result=native_compression();assert result['all_actual_hopping_outputs_leave_pinned_density']
    y=representation()[2][11].diagonal();sy=higgs_memory()[1][11].diagonal()
    for channel,row,col,ops,coefficient in yukawa_words():
        assert sy[row]-sy[col]+sum((1 if create else -1)*y[i%15] for i,create in ops)==0


def test_exact_minimum_differences_independent_green_inverse():
    for n in (2,4):
        d,g=cochains(n);independent=np.linalg.pinv(np.array(d.T*d,dtype=float))
        for cell in ((0,0,0),(1,0,0),(1,1,1)):
            index=(cell[0]*n+cell[1])*n+cell[2]
            assert abs(float(exact_minimum_difference(n,cell))-(independent[0,0]-independent[0,index]))<1e-14
    assert exact_minimum_difference(2,(1,0,0))==s.Rational(7,48)
    assert exact_minimum_difference(4,(1,0,0))==s.Rational(21,128)
    assert exact_minimum_difference(4,(1,0,0),2)==s.Rational(21,64)
    with pytest.raises(ValueError):exact_minimum_difference(2,(1,0,0),-1)


def test_all_native_Gauss_string_witnesses():
    assert all(row['all_twelve_Gauss_constraints'] for row in exact_rows())


def test_unknown_local_minimum_and_loop_multiplicity_cancel():
    # Independent tensor spectrum with a hidden charged level below the
    # visible trial. It changes the absolute minimum, not its separation.
    h0=np.array([[2.,.7],[.7,-3.]])
    hchi=np.array([[-11.,.3,0],[.3,-8.,.2],[0,.2,-12.]])
    local=np.kron(h0,np.eye(3))+np.kron(np.eye(2),hchi)
    differences=[]
    for loops,cost in ((1,0.),(2,7/48),(5,21/128)):
        full=np.kron(local,np.eye(loops))+cost*np.eye(len(local)*loops)
        differences.append(np.linalg.eigvalsh(full)[0])
    assert abs((differences[1]-differences[0])-7/48)<1e-13
    assert abs((differences[2]-differences[0])-21/128)<1e-13
    # Restoring a path-length penalty changes the static difference.
    assert abs((differences[1]+1-differences[0])-7/48)>.9


def test_static_force_sign_and_normalization():
    row=profile(64)['rows'][0]
    assert row['finite_difference_force']<0
    assert .9<row['potential_difference_amplitude_ratio']<1.1
