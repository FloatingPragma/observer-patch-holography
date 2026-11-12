"""Independent spectral, Clebsch--Gordan and projector-tail controls."""
from pathlib import Path
import hashlib
import json

import numpy as np
import scipy.linalg as la
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TARGET=ROOT/'model/openwave_energy_resolved_transport'


def spin_control():
    # A three-dimensional triplet Gram, obtained by deleting the invariant
    # Bell vector from C^2 tensor C^2, supplies independent CG maps.
    cg=s.Matrix([[0,1,0,0],[0,0,1,0],[1/s.sqrt(2),0,0,-1/s.sqrt(2)]])*2/s.sqrt(3)
    maps=[cg[:,0:2],cg[:,2:4]]
    pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    assert cg.rank()==3
    for i in range(2):
        for j in range(2):
            gram=maps[i].H*maps[j]
            expected=int(i==j)*s.eye(2)-sum((p[j,i]*p for p in pauli),s.zeros(2))/3
            assert s.simplify(gram-expected)==s.zeros(2)
    for axis,p in enumerate(pauli):
        effects=s.zeros(4)
        for spin in range(2):
            amplitude=s.zeros(6,4)
            for a in range(2):
                for b in range(2):
                    amplitude+=p[a,b]*s.kronecker_product(maps[a].H*maps[spin],maps[b])/(4*s.sqrt(2))
            effects+=amplitude.H*amplitude
        bond=sum(((1 if j==axis else -1)*s.kronecker_product(q,q)/4
                  for j,q in enumerate(pauli)),s.zeros(4))
        assert s.simplify(effects-s.eye(4)/12-2*bond/27)==s.zeros(4)
        assert effects!=s.eye(4)/12+bond/27  # wrong-spin-coefficient control
    return {'triplet_Gram_rank':3,'direct_CG_first_hop_axes':3,
            'wrong_exchange_coefficient_rejected':True}


def native_residue_control(receipt):
    from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source
    from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.model import transfer
    def array(matrix):return np.array(matrix.tolist(),complex)
    sb=source(1);g=array(sb['gram']);h=array(sb['H'])
    e,u=la.eigh(g@h,g);ground=u[:,0]
    x=s.Symbol('x')
    roots=sorted(float(s.re(z)) for z in s.nroots(x**4-2*x**3-212*x**2+213*x+11114))
    pole=roots[0];weights=[];dimensions=[]
    for spin in range(2):
        tb=transfer('ec',spin);tg=array(tb['gram']);th=array(tb['H'])
        values,vectors=la.eigh(tg@th,tg)
        hits=np.where(abs(values-pole)<1e-9)[0]
        assert len(hits)==1
        amplitude=vectors[:,hits[0]].conj()@tg@array(tb['map'])@ground
        weights.append(float(abs(amplitude)**2));dimensions.append(th.shape[0])
    assert abs(weights[1]-2*weights[0])<2e-12
    w=sum(weights)/2
    interval=[float(s.Rational(v)) for v in receipt['atomic_pole']['pump_pole_weight_interval']]
    assert abs(w-sum(interval)/2)<2e-12
    return {'generalized_native_eigenproblem_dimensions':dimensions,
            'direct_up_down_residues':weights,'pump_weight':w,
            'residue_ratio_error':abs(weights[1]-2*weights[0]),
            'floating_calculation_is_independent_diagnostic_not_interval_certificate':True}


def projector_controls():
    # Independent finite examples preserve only the hypotheses used by the
    # theorem; no global atomic-band isolation or detector commutation.
    rng=np.random.default_rng(20843);cases=0;largest_identity=0.;positive_tails=0
    for size in (3,4,6):
        h0=np.diag(np.arange(size,dtype=float));f=np.zeros((size,size));f[1,1]=1
        for coupling in (1e-4,0.003,0.02):
            for repeat in range(5):
                raw=rng.normal(size=(size,size))+1j*rng.normal(size=(size,size))
                t=(raw+raw.conj().T)/2;t[1,1]=0
                t/=la.norm(t,2);h=h0+coupling*t
                energies,vectors=la.eigh(h);mask=abs(energies)<0.2
                e=vectors[:,mask]@vectors[:,mask].conj().T
                x0=np.eye(size)[:,0];beta=coupling/0.8
                assert la.norm(f@e,2)<=beta+1e-13
                exact=f@e@x0
                rhs=-coupling*f@t@x0+coupling*f@e@t@x0+coupling*f@t@(np.eye(size)-e)@x0
                error=la.norm(exact-rhs);largest_identity=max(largest_identity,error)
                assert error<2e-14
                assert la.norm(exact)>1e-8;positive_tails+=1
                assert la.norm(f@x0)==0  # no-filter control
                assert la.norm(exact-coupling*f@t@x0)>la.norm(exact-rhs)*1000
                weights=rng.random(size);weights[1]=0;weights*=0.5/sum(weights)
                y=np.sqrt(weights);selected=np.sqrt(weights[0])*x0
                distances=[0.2]+[j-0.2 for j in range(1,size)]
                remainder=coupling*(beta+coupling/distances[0])*np.sqrt(weights[0])
                remainder+=beta*coupling*sum(np.sqrt(weights[j])/distances[j] for j in range(1,size))
                assert la.norm(f@e@y+coupling*f@t@selected)<=remainder+1e-13
                for j in range(1,size):
                    xj=np.eye(size)[:,j]
                    assert la.norm(f@e@xj)<=beta*coupling/distances[j]+1e-13
                c=coupling*(1/np.sqrt(2)+sum(j*np.sqrt(weights[j])/distances[j] for j in range(1,size)))
                for time in (0.03,1.,30.):
                    u=(vectors*np.exp(-1j*time*energies))@vectors.conj().T
                    change=la.norm(f@(u-np.eye(size))@e@y)
                    assert change<=min(beta*time*c,np.sqrt(2)*beta)+1e-13
                cases+=1
    # Idempotence FE=FEE supplies the second v factor, including the test
    # atomic vector that lies in the detector range.
    return {'independent_Hermitian_examples':cases,'time_comparisons':3*cases,
            'largest_commutator_identity_error':largest_identity,
            'positive_time_zero_dressing_examples':positive_tails,
            'wrong_leading_sign_rejected':True,'no_filter_control_is_zero':True}


def record_control():
    q,a,b=s.symbols('q a b');r=s.symbols('r:6');u=s.symbols('u:6')
    records=[1-q,q*(1-a),q*a/2,q*(a/2-b),*(q*x for x in r),*(q*x for x in u),q*(b-sum(r)-sum(u))]
    assert len(records)==17 and s.expand(sum(records))==1
    values={q:s.Rational(3,4),a:s.Rational(4,5),b:s.Rational(1,3)}
    values.update({x:s.Rational(1,100) for x in r+u})
    assert all(x.subs(values)>=0 for x in records)
    return {'complete_outcomes':17,'exact_probability_sum':'1','positive_example_checked':True}


def report():
    paths=[TARGET/name for name in ('response.py','PROOF.md','results.json')]+[Path(__file__)]
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    receipt=json.loads((TARGET/'results.json').read_text())
    for name,digest in receipt['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    result={'schema':'oph.independent_energy_filtered_arrival_review.v1','input_sha256':pins,
            'spin':spin_control(),'native_residues':native_residue_control(receipt),
            'projector_controls':projector_controls(),'records':record_control(),
            'scope':'Positive t=0 energy-filter dressing is supported. These checks do not establish propagation, a free electron or a relativistic mass shell.',
            'verdict':'PASS_INDEPENDENT_ENERGY_FILTERED_ARRIVAL_REVIEW'}
    assert all(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest for name,digest in pins.items())
    return result


if __name__=='__main__':
    result=report();path=HERE.parent/'energy_review.json'
    path.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(result['verdict'])
