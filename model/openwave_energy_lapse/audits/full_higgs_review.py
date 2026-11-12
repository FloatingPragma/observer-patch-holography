"""Independent native tensors, full-CAR floor and actual colour-singlet checks."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,itertools
import numpy as np
import sympy as s
from model.common_source import native_matter_attachment as matter
ROOT=Path(__file__).resolve().parents[3]

def car(n,j):
    a=s.zeros(2**n)
    for mask in range(2**n):
        if mask>>j&1:a[mask^(1<<j),mask]=(-1)**((mask&((1<<j)-1)).bit_count())
    return a

def run():
    rep=matter.representation();indices=list(rep['selected'])
    weak=rep['operators']['weak_C2'].extract(indices,indices);Q=(rep['operators']['T3']+rep['operators']['q']/6).extract(indices,indices)
    vertices={tuple(map(str,x['origin'])):x for x in matter.invariant_higgs_lines()}
    channels={r['existing_action_channel']:s.sympify(r['same_current_tensor_scale_exact'])*vertices[tuple(r['charge_origin'])]['components'][1]for r in matter.action_intertwiner()[1]}
    C=sum(channels.values(),s.zeros(15));P=C**2;charged=s.diag(*(int(Q[i,i]!=0)for i in range(15)))
    assert P==charged and P.rank()==14 and Q*(s.eye(15)-P)==s.zeros(15)
    assert Q.T*C+C*Q==s.zeros(15)
    for a,b in itertools.combinations(channels.values(),2):assert a*b==s.zeros(15)
    bad=[]
    for name in('up','down'):
        reduced=C-channels[name];indices_zero=[i for i in range(15)if reduced[:,i]==s.zeros(15,1)and Q[i,i]!=0];assert len(indices_zero)==6;bad.append(indices_zero)
    # Independent full four-mode exterior matrices, with both spin pairs.
    c=[np.asarray(car(4,j),complex)for j in range(4)];N=[a.conj().T@a for a in c]
    V=c[0].conj().T@c[3].conj().T-c[1].conj().T@c[2].conj().T;V+=V.conj().T
    minmargin=0.;count=0
    for mass in(.1,.7,1.3):
        for momentum in(0.,-.3,.8):
            H=momentum*(N[0]+N[3]-N[1]-N[2])+mass*V;values,vectors=np.linalg.eigh(H);ground=values[0]
            for charge in(-1.,2/3,-1/3):
                q=charge*(N[0]+N[1]-N[2]-N[3]);assert np.linalg.norm(H@q-q@H)<1e-13
                absq=np.diag(abs(np.diag(q)));margin=np.linalg.eigvalsh(H-ground*np.eye(16)-mass*absq)[0]
                assert margin>-2e-14;minmargin=min(minmargin,float(margin));count+=1
                assert abs(np.vdot(vectors[:,0],q@vectors[:,0]))<1e-13
    # The native weak-doublet colour action, retained for two up and one down.
    selected=[i for charge in(s.Rational(2,3),s.Rational(-1,3))for i in range(15)if Q[i,i]==charge and weak[i,i]!=0]
    assert len(selected)==6
    generators=[]
    for i in range(3):
        for j in range(i+1,3):
            generators.extend([matter.unit(5,i,j)+matter.unit(5,j,i),s.I*(matter.unit(5,i,j)-matter.unit(5,j,i))])
    generators.extend([matter.unit(5,0,0)-matter.unit(5,1,1),matter.unit(5,1,1)-matter.unit(5,2,2)])
    small=[matter.exterior_derivative(g).extract(indices,indices).extract(selected,selected)for g in generators]
    a=[car(6,j)for j in range(6)];number=sum((v.T*v for v in a),s.zeros(64));charge=s.zeros(64)
    for j,index in enumerate(selected):charge+=Q[index,index]*a[j].T*a[j]
    sector=[mask for mask in range(64)if number[mask,mask]==3 and charge[mask,mask]==1];assert len(sector)==9
    second=[]
    for g in small:
        big=sum((g[i,j]*a[i].T*a[j]for i in range(6)for j in range(6)),s.zeros(64));second.append(big)
    casimir=sum((g*g for g in second),s.zeros(64));kernel=casimir.extract(sector,sector).nullspace();assert len(kernel)==1
    singlet=s.zeros(64,1)
    for i,mask in enumerate(sector):singlet[mask]=kernel[0][i]
    assert all(g*singlet==s.zeros(64,1)for g in second)
    assert charge*singlet==singlet and number*singlet==3*singlet
    # Direct 180-row spectral comparison is independent of producer kinetic.py.
    Y=np.diag([.13,.47,.83]);kernel_mass=np.kron(np.asarray(C,float),Y)
    sig=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1])]
    eps=np.array([[0,1],[-1,0]]);delta=np.kron(eps,kernel_mass)
    errors=[]
    for p in[(0.,0.,0.),(.23,-.37,.41),(1.,1.,0.)]:
        h=-np.kron(sum((v*t for v,t in zip(p,sig)),np.zeros((2,2),complex)),np.eye(45))
        H=np.block([[h,delta],[delta.conj().T,h.T]])
        expected=np.kron(np.eye(4),np.kron(np.asarray(P,float),Y@Y))+sum(v*v for v in p)*np.eye(180)
        error=np.linalg.norm(H@H-expected);assert error<1e-13;errors.append(float(error))
    path=ROOT/'model/openwave_full_higgs_floor/results.json';receipt=json.loads(path.read_text())
    assert receipt['exact_native_symbol']['native_internal_rank']==14
    refs=set(receipt['source_sha256'])|{'model/openwave_full_higgs_floor/results.json',str(Path(__file__).relative_to(ROOT))}
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v for p,v in receipt['source_sha256'].items())
    return {'verdict':'PASS','independent_full_CAR_floor_members':count,'minimum_numerical_operator_margin':minmargin,
      'native_internal_charge_support_rank':14,'neutral_kernel_dimension':1,
      'removed_quark_channel_charged_zero_controls':bad,'native_colour_singlet_sector_dimension':len(sector),
      'unique_native_colour_singlet_charge':'1','native_colour_singlet_quasiparticle_count':3,
      'all_eight_colour_generators_annihilate_singlet':True,'independent180row_squared_symbol_errors':errors,
      'scope':'Shared family coefficients and the analyzer are declared. The exact charged floor retains neutral zeros and spatial valleys; equal quark/lepton family spectra do not fit observed quark masses.',
      'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(refs)}}
if __name__=='__main__':
    r=run();Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(r['verdict'])
