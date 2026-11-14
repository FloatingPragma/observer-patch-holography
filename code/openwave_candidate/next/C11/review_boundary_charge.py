"""Independent complete cut/boundary charge projector audit."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).parent;BASE=HERE.parents[1];sys.path.insert(0,str(BASE/'next/C10'))
from charged_reservoir import BoundarySource


def run():
    paths=[Path(__file__),HERE/'boundary_charge.py',HERE/'BOUNDARY_CHARGE.json',BASE/'next/C10/charged_reservoir.py']
    pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    data=json.loads((HERE/'BOUNDARY_CHARGE.json').read_text());m=BoundarySource(6);numbers=[mask.bit_count() for mask,_ in m.states];checked=0;outcomes=0;zero=0
    for packet in data['complete_regions']:
        region=set(packet['region']);r=len(region);gauss=sum(z<3 for z in region);groups={};cuts=[]
        for i,(mask,f) in enumerate(m.states):
            q=sum((mask>>z)&1 for z in region);n=numbers[i];boundary=3-n;cut=0
            for e,(a,b) in enumerate(m.construction.seams):
                cut+=f[2*e]*(int(a in region)-int(3+e in region));cut+=f[2*e+1]*(int(b in region)-int(3+e in region))
            cut%=6;cuts.append(cut);groups.setdefault((n,cut),[]).append(i)
            if (q+cut+boundary*int(0 in region)-gauss)%6:raise ValueError('dressedGaussidentityfailed')
            checked+=1
        seen=[0]*m.n
        for record in packet['outcomes']:
            n=record['total_number'];cut=record['cut_residue'];expected=groups.get((n,cut),[])
            if record['basis_indices']!=expected:raise ValueError('projectorwrong')
            admitted=[q for q in range(7) if q<=r and 0<=n-q<=6-r and (q+cut+(3-n)*int(0 in region)-gauss)%6==0]
            if len(admitted)>1 or record['regional_charge']!=(admitted[0] if admitted else None):raise ValueError('decoderwrong')
            if record['boundary']!=3-n:raise ValueError('boundarywrong')
            for i in expected:seen[i]+=1
            outcomes+=1;zero+=not expected
        if seen!=[1]*m.n:raise ValueError('notacompleteinstrument')
    # Isolate the added boundary measurement itself: same cut AND local charge,
    # even total fermion parity, but distinct B labels.
    region={3};a=next(i for i,(mask,_) in enumerate(m.states) if mask==0);b=next(i for i,(mask,_) in enumerate(m.states) if mask==6)
    controls=[]
    for i in (a,b):
        mask,f=m.states[i];q=int(bool(mask&(1<<3)));cut=(-f[0]-f[1])%6
        controls.append({'mask':mask,'number':mask.bit_count(),'parity':mask.bit_count()%2,'regional_charge':q,'cut_residue':cut,'boundary':3-mask.bit_count()})
    if controls[0]['cut_residue']!=controls[1]['cut_residue']:raise ValueError('controlcutdiffers')
    coherent=np.array([[.5,.5j],[-.5j,.5]]);dephased=np.diag(np.diag(coherent));distance=float(np.linalg.norm(coherent-dephased,ord='nuc')/2)
    # General Pauli interval width by explicit admissible occupation counts.
    for d in range(1,13):
        for n in range(d+1):
            for r in range(d+1):
                possible=[q for q in range(d+1) if q<=r and 0<=n-q<=d-r]
                low=max(0,n-d+r);high=min(n,r)
                if possible!=list(range(low,high+1)) or high-low> d//2:raise ValueError('Pauliintervalproofcontrolfailed')
    checks={'all24576basis_cut_relations':checked==24576,'all2688projectorbranches_equal_and_complete':outcomes==2688,'same_even_parity_same_cut_Q_boundary_disturbance':abs(distance-.5)<1e-12 and controls[0]['parity']==controls[1]['parity'],'general_Pauli_interval_bound':True}
    if not all(checks.values()):raise ValueError(checks)
    if any(hashlib.sha256((BASE/p).read_bytes()).hexdigest()!=pin for p,pin in pins.items()):raise ValueError('sourceschanged')
    return {'schema':'oph.c11.boundary_charge_independent_review.v1','checks':checks,'zero_outcomes_retained':zero,'boundary_only_disturbance_control':{'region':[3],'states':controls,'trace_distance':distance},
      'proof_scope':'Complete joint(B,cut) versus(N,Q) diagonalprojectors areidentical withknownPauliinterval; tensorswith everyexternalreference. Q-only instrumentdiffers fromaddedBmeasurement, evenwithoutcrossparity orchangedcut. No physicalelectricdictionary/topologicalelectronwinding/fractionalcharge exclusion beyonddeclaredoccupancy.',
      'source_sha256':pins,'full_criterion_met':False,'official_rows_changed':[]}

if __name__=='__main__':
 r=run();HERE.joinpath('BOUNDARY_CHARGE_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
