#!/usr/bin/env python3
"""Finite entangled-pair spectral candidate and its normalization boundary.

Reference: A. D. Wissner-Gross, A Parameter-Free Candidate for the
Fine-Structure Constant from Entangled Qudits, preprint (2026).
The photon/Thomson identification is a separate physical conjecture.
"""
from __future__ import annotations
from collections import Counter
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from evidence_contract import public_receipt
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def run():
    dimensions=(3,2,3,2)
    labels=[]
    for i,j in combinations(range(4),2):
        for a,b in product(range(dimensions[i]**2-1),range(dimensions[j]**2-1)):
            label=[-1]*4;label[i],label[j]=a,b;labels.append(tuple(label))
    assert len(labels)==169
    index={x:i for i,x in enumerate(labels)}
    exchange=np.zeros((169,169),dtype=np.int64)
    X=np.eye(172,dtype=np.int64);Y=X.copy()
    for col,label in enumerate(labels):
        row=index[(label[2],label[3],label[0],label[1])]
        exchange[row,col]=1
        positions=[i for i,x in enumerate(label) if x>=0]
        target=X if dimensions[positions[0]]==dimensions[positions[1]] else Y
        target[:169,col]=-exchange[:,col]
    X[169:,169:]=-np.eye(3,dtype=np.int64)
    Y[169:,169:]=-np.eye(3,dtype=np.int64)
    I=np.eye(172,dtype=np.int64)
    assert np.array_equal(X@X,I) and np.array_equal(Y@Y,I)
    assert np.array_equal(X@Y,Y@X)
    Ktwice=(2*I-X)@(3*I-Y)
    ranks={};projectors=[]
    for x,y in product((1,-1),repeat=2):
        Q=(I+x*X)@(I+y*Y)
        assert np.array_equal(Q@Q,4*Q)
        gain=(2-x)*(3-y)//2
        assert np.array_equal(Ktwice@Q,2*gain*Q)
        assert int(np.trace(Q))%4==0
        ranks[str(gain)]=int(np.trace(Q))//4
        projectors.append(Q)
    assert ranks=={'1':79,'2':48,'3':42,'6':3}
    assert np.array_equal(sum(projectors),4*I)
    normalization=(3**2-1)*(2**2-1)*(3*2)**4
    assert normalization==31104
    with mp.workdps(80):
        alpha=sum(mp.mpf(n)*mp.sqrt(int(gain)) for gain,n in ranks.items())/normalization
        inverse=1/alpha
        pixel=(1+mp.sqrt(5))/2+mp.sqrt(mp.pi)*alpha
        numbers={'alpha':mp.nstr(alpha,70),'inverse_alpha':mp.nstr(inverse,70),
                 'P_after_declared_same_cell_attachment':mp.nstr(pixel,70)}
    # Bell projection on two flat 3x2 pairs keeps its complementary outcome.
    pair=np.zeros((3,2));pair[0,0]=pair[1,1]=1/np.sqrt(2)
    initial=np.einsum('ab,cd->abcd',pair,pair)
    omega=np.eye(3)/np.sqrt(3)
    far=np.einsum('ac,abcd->bd',omega,initial)
    success=float(np.vdot(far,far).real)
    assert abs(success-1/6)<1e-14
    A=np.array([[1.,0.,0.],[0.,1.,0.]])/np.sqrt(2)
    B=np.array([[1.,0.],[0.,.5],[0.,np.sqrt(3)/2]])/np.sqrt(2)
    assert np.allclose(A@A.T,np.eye(2)/2) and np.allclose(B.T@B,np.eye(2)/2)
    joined=A@B;joined/=np.linalg.norm(joined)
    spectrum=np.linalg.eigvalsh(joined@joined.T)
    assert np.allclose(spectrum,[.2,.8])
    artifact=HERE/'results/entangled_alpha_operators.npz'
    np.savez_compressed(artifact,exchange=exchange,X=X,Y=Y,twice_K=Ktwice,
                        four_times_joint_projectors=np.array(projectors))
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    result={'status':'PASS_EXACT_ENTANGLED_PAIR_SPECTRAL_CANDIDATE_AND_SCOPED_NORMALIZATION_CONTROLS',
        'producer_sha256':sha(__file__),'source_sha256':bound_pins({}),'artifacts':{artifact.name:sha(artifact)},
        'primary_reference':'https://alexwg.org/publications/ChainTheoryAlpha.pdf',
        'attribution':{'author':'A. D. Wissner-Gross',
            'title':'A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits',
            'year':2026,'borrowed_construction':'Entangled-pair exchange operators, spectral weights and normalization giving the exact scalar.',
            'OPH_contribution':'Independent reconstruction and complete-outcome, support-alignment and normalized-readout controls; physical electromagnetic identification remains conditional.'},
        'empirical_inputs':[],'operator_dimension':36,'mixed_traceless_response_dimension':24,
        'product_records':169,'extra_exchange_odd_records':3,'nonzero_feature_dimension':172,
        'full_feature_normalization':normalization,'gain_multiplicities':ranks,
        'exact_integer_involution_and_projector_checks':True,
        'spectral_expression':'(79+48 sqrt(2)+42 sqrt(3)+3 sqrt(6))/31104',**numbers,
        'complete_Bell_join':{'success_probability':'1/6','failure_probability':'5/6'},
        'misaligned_flat_support_control':{'conditional_reduced_spectrum':['1/5','4/5'],
            'conclusion':'Flat factors alone do not force a flat contracted pair. Aligned supports or an adjusted contraction are extra hypotheses.'},
        'invariant_readout_control':{'example':'I4/4 and the normalized Bell projector are invariant under U tensor conjugate(U), both read I4 as 1, and read that projector as 1/4 and 1.',
            'matching_invariance_and_identity_normalization_select_trace':False},
        'physical_scope':{'finite_scalar_evaluated':True,'Thomson_photon_residue_identified':False,
            'Maxwell_kinetic_normalization_selected':False,'same_cell_P_alpha_bridge_proved':False,
            'candidate_overrides_existing_fixed_point_maps':False}}
    return public_receipt(result)


if __name__=='__main__':
    result=run()
    (HERE/'results/entangled_alpha.json').write_text(json.dumps(result,indent=2)+'\n', encoding='utf-8', newline='\n')
    print(result['status'])
