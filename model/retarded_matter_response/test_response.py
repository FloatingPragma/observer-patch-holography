from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import null_space, solve

from model.cube_eta_protection import cube
from model.interacting_vacuum import ground
from model.retarded_matter_response import response as r


def test_exact_all_channel_source_moments_and_memory():
    value=r.certificate()
    channels={(c['reference'],c['channel']):c for c in value['channels']}
    assert channels['ground','down']['memory_Gram_diagonal']==['1/4']*8
    assert channels['charged','down']['memory_Gram_diagonal']==['0','2/9','2/9','1/4','2/9','1/4','1/4','1/4']
    assert all(channels[ref,'up']['free_resolvent_by_exact_full_cyclic_closure'] for ref in ('ground','charged'))
    assert value['free_down_control']['memory_Gram_diagonal']==['0']*8
    assert value['common_metric_or_Einstein_premise_completed'] is False


def test_restricted_matrix_is_actual_full_source_with_no_discarded_couplings():
    source=cube.operators().hubbard_H.real
    for up,down in ((4,1),(5,1),(3,0),(4,0),(5,0),(6,0)):
        masks=r.sector_masks(up,down)
        block=r.source_block(up,down)
        assert (source[masks][:,masks]-block).nnz==0
        assert source[:,masks].nnz==block.nnz


def test_charged_hartree_change_is_the_actual_added_orbital_density():
    adjacency=ground.band_data()[0]
    # sqrt2 P+ delta0 = (delta0+A delta0/sqrt3)/sqrt2.
    density=[F(int(x==0),2)+F(int(adjacency[x,0])**2,6) for x in range(8)]
    difference=r.first_numerator('charged','down')-r.first_numerator('ground','down')
    assert all(F(int(difference[x,x]),48)==density[x] for x in range(8))
    assert sum(density)==1
    assert np.count_nonzero(difference-np.diag(np.diag(difference)))==0


def test_complete_addition_removal_sign_and_free_comparator():
    for reference in ('ground','charged'):
        data=r.dilation(reference,'up')
        a,b=data.vectors
        v=np.sqrt(data.norm_square)*(a+np.sqrt(3)*b)/36
        generator=data.h16.toarray()/16+np.diag(data.sqrt3_shift*np.sqrt(3)/16)
        z=complex(1,2)
        actual=v.T@solve(z*np.eye(len(v))-generator,v)
        expected=np.linalg.inv(z*np.eye(8)-ground.band_data()[0]/16)
        assert np.linalg.norm(actual-expected,2)<2e-14
        # Forgetting the removal sign is a genuine failed control.
        wrong=generator.copy();split=data.addition_dimension
        wrong[split:,split:]*=-1
        assert np.linalg.norm(v.T@solve(z*np.eye(len(v))-wrong,v)-expected,2)>1e-3


def test_schur_memory_matches_direct_complete_sector_resolvent():
    # Independent numerical diagnostic, separate from exact admission.
    data=r.dilation();a,b=data.vectors
    v=(a+np.sqrt(3)*b)/36
    k=data.h16.toarray()/16+np.diag(data.sqrt3_shift*np.sqrt(3)/16)
    q=null_space(v.T)
    first=v.T@k@v;w=q.T@k@v;compressed=q.T@k@q
    z=complex(.1,.4)
    direct=v.T@solve(z*np.eye(560)-k,v)
    memory=w.T@solve(z*np.eye(552)-compressed,w)
    schur=np.linalg.inv(z*np.eye(8)-first-memory)
    assert np.linalg.norm(direct-schur,2)<3e-13
    assert np.linalg.eigvalsh(-direct.imag).min()>0
    assert np.linalg.eigvalsh(-memory.imag).min()>0
    assert np.linalg.norm(direct-np.linalg.inv(z*np.eye(8)-first),2)>.1


@pytest.mark.parametrize('reference,spin',[('vacuum','down'),('ground','physical'),('', 'up')])
def test_reference_domain(reference,spin):
    with pytest.raises(ValueError):r.dilation(reference,spin)


def test_overflow_rejected_before_integer_multiplication():
    with pytest.raises(ArithmeticError):
        r._matmul(np.array([[-2**63]],dtype=np.int64),np.array([[2]],dtype=np.int64))
    with pytest.raises(ValueError):r.dilation(interaction=1)
    with pytest.raises(ValueError):r.dilation().moment(5)
