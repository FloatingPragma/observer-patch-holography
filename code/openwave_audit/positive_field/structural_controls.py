"""Finite noncommuting-field controls of the block argument, not a field cutoff."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import hashlib
import json
import numpy as np
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'inputs/color_source.py'
spec=importlib.util.spec_from_file_location('colour_exact_source',SOURCE)
source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)

def run():
    sx=np.array([[0,1],[1,0]],complex);sz=np.diag([1,-1])
    links={0:expm(.37j*sx),1:expm(-.29j*sz)}
    assert np.linalg.norm(links[0]@links[1]-links[1]@links[0])>.1
    cases=[]
    for L in [2,3]:
      for kappa in [F(1,2),F(2)]:
        d=source.data(L,kappa,1);D=d['dimension'];labels=d['labels'];ix={n:i for i,n in enumerate(d['counts'])}
        plain=np.zeros((D,D),complex);cov=np.zeros((2*D,2*D),complex)
        for col,occ in enumerate(d['counts']):
          for x in range(L-1):
            for sw in range(4):
              left=4*x+sw;right=4*(x+1)+sw
              diagonal=(occ[left]+occ[right])/8
              plain[col,col]+=diagonal
              cov[2*col:2*col+2,2*col:2*col+2]+=diagonal*np.eye(2)
              for dest,origin,W in [(right,left,links[sw%2]),(left,right,links[sw%2].conj().T)]:
                target,z=source.act(occ,dest,origin)
                if target is None:continue
                row=ix[target];v=-source.value(z)/8
                plain[row,col]+=v
                cov[2*row:2*row+2,2*col:2*col+2]+=v*W
        assert np.linalg.norm(cov-cov.conj().T)<1e-13
        assert np.linalg.eigvalsh(cov)[0]>-1e-12
        H0=source.dense(d)[:D,:D]
        matter=np.kron(H0-plain,np.eye(2))+cov
        P=[i for i,s in enumerate(d['span']) if s==0]
        Q=[i for i in range(D) if i not in P]
        p=[2*i+j for i in P for j in range(2)];q=[2*i+j for i in Q for j in range(2)]
        cross=matter[np.ix_(q,p)]
        expected=np.zeros((len(p),len(p)))
        for z,i in enumerate(P):
          x=labels[i][0]//4;degree=int(x>0)+int(x<L-1)
          expected[2*z:2*z+2,2*z:2*z+2]=3*degree/64*np.eye(2)
        gram_error=np.linalg.norm(cross.conj().T@cross-expected,2)
        assert gram_error<1e-13
        # A positive charge energy and a field energy that does not commute with links.
        B=np.zeros((L,L-1))
        for x in range(L-1):B[x,x]=-1;B[x+1,x]=1
        green=np.linalg.pinv(B@B.T)
        coulomb=[]
        for occ in d['counts']:
          charge=np.array([(2*(occ[4*x]+occ[4*x+2])-occ[4*x+1]-occ[4*x+3])/3 for x in range(L)])
          coulomb.append(float(charge@green@charge)/8)
        assert max(abs(coulomb[i]) for i in P)<1e-14
        field=np.diag([0,1.]);H=matter+np.kron(np.diag(coulomb),np.eye(2))+np.kron(np.eye(D),field)
        endpoint=np.zeros(D,complex)
        endpoint[labels.index((0,1,3))]=1/np.sqrt(3)
        endpoint[labels.index((1,1,2))]=-np.sqrt(2/3)
        delta=3/2000;psi=np.kron(endpoint,[np.sqrt(1-delta),np.sqrt(delta)])
        assert abs(np.vdot(psi,H@psi).real-(57/16+delta))<1e-12
        outputs=[]
        for weak in [0,.25,-.25]:
          Hw=H+weak*np.kron(np.diag([float(x) for x in d['dipole'][:D]]),np.eye(2))
          floor=float(np.linalg.eigvalsh(Hw[np.ix_(q,q)])[0])
          assert floor>= (281/48 if weak==0 else 91/16)-1e-12
          values,vectors=np.linalg.eigh(Hw);coeff=vectors.conj().T@psi;maximum=0
          for time in np.linspace(0,20,81):
            state=vectors@(np.exp(-1j*time*values)*coeff)
            maximum=max(maximum,float(np.vdot(state[q],state[q]).real))
          assert maximum<(.07 if weak==0 else .08)
          outputs.append({'probe':weak,'actual_separated_floor':floor,'sampled_max_breakup':maximum})
        cases.append({'L':L,'kappa':str(kappa),'physical_matter_dimension':D,
                      'cross_gram_error':gram_error,'positive_hopping_minimum':float(np.linalg.eigvalsh(cov)[0]),'outputs':outputs})
    # Removing the separated-energy floor permits complete resonant breakup.
    b=np.sqrt(3/32);Hbad=np.array([[57/16,b],[b,57/16]])
    resonant=expm(-1j*np.pi/(2*b)*Hbad)@np.array([1.,0.])
    badprob=float(abs(resonant[1])**2);assert badprob>.999999999999
    return {'schema':'oph.color_positive_field.structural_controls.v1',
      'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
      'cases':cases,'negative_energy_floor_control_breakup':badprob,
      'field_domain':'Two-dimensional noncommuting bounded-field control only; no approximation to the actual Maxwell field is inferred.',
      'actual_Maxwell_field_executed':False,'full_row_award':False}

if __name__=='__main__':
    result=run();HERE.joinpath('STRUCTURAL_CONTROLS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'cases':len(result['cases']),'negative_control_breakup':result['negative_energy_floor_control_breakup'],'pass':True}))
