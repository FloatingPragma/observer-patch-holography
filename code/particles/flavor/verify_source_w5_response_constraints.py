#!/usr/bin/env python3
"""Independent rational/group-character verifier of the source W5 receipt.

No producer import. Geometry comes from an independent distance-isometry
enumeration; color cumulants use moment recursion rather than derivatives.
The source-boundary audit checks declared imports/file literals, not a sandbox.
"""
from __future__ import annotations
import argparse,ast,hashlib,json,math
from collections import Counter,deque
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
REPO=Path(__file__).resolve().parents[3]
PRODUCER_REL='code/particles/flavor/source_w5_response_constraints.py'
SPEC_REL='code/particles/flavor/source_w5_response_constraints_spec.json'
SPEC_SHA256='71a112358ff9e24efa843cfebec96678b0c9d90e428467ce080f5da8642bcb00'
SOURCE_FILES=['code/a5_closure/manifests/echosahedral_federation_reference.json','code/a5_closure/manifests/a3_scheduler_kernel_reference.json','code/a5_closure/manifests/record_counting_mechanism_reference.json','code/a5_closure/source_repair_generator_certificate.py']
DEFAULT_INPUT=REPO/'code/particles/runs/flavor/source_w5_response_constraints.json'

def need(ok,why):
    if not ok:raise ValueError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def encode(A):return [[str(x) for x in row] for row in A]
def outer(x):return [[a*b for b in x] for a in x]
def add(A,B):return [[x+y for x,y in zip(a,b)] for a,b in zip(A,B)]
def mulscale(A,c):return [[x*c for x in row] for row in A]
def mm(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def zeros(n):return [[F(0) for _ in range(n)] for _ in range(n)]
def field_add(x,y):return (x[0]+y[0],x[1]+y[1])
def field_mul(x,y):return (x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def field_div(x,y):
    d=y[0]**2-5*y[1]**2
    return field_mul(x,(y[0]/d,-y[1]/d))

def audit_source(text):
    tree=ast.parse(text)
    allowed={'argparse','hashlib','itertools','json','sys','collections','pathlib','sympy','source_repair_generator_certificate','__future__'}
    allowed_files={*SOURCE_FILES,PRODUCER_REL,SPEC_REL,'code/particles/runs/flavor/source_w5_response_constraints.json',*[Path(p).name for p in SOURCE_FILES]}
    for node in ast.walk(tree):
        if isinstance(node,ast.Import):need(all(n.name in allowed for n in node.names),'undeclared producer import')
        if isinstance(node,ast.ImportFrom):need(node.module in allowed,'undeclared producer import')
        if isinstance(node,ast.Constant) and isinstance(node.value,str) and node.value.endswith(('.json','.csv','.tsv','.npz')):
            need(node.value in allowed_files,'undeclared producer data path')

@lru_cache(maxsize=4)
def reconstruct_geometry(manifest_bytes):
    d=json.loads(manifest_bytes);ports=d['carrier']['ports'];index={p:i for i,p in enumerate(ports)}
    faces=[tuple(index[p] for p in f) for f in d['carrier']['oriented_faces']]
    edges={tuple(sorted((a,b))) for face in faces for a,b in zip(face,face[1:]+face[:1])}
    need(edges=={tuple(sorted((index[a],index[b]))) for a,b in d['carrier']['edges']},'source incidence')
    adjacency=[set() for _ in ports]
    for i,j in edges:adjacency[i].add(j);adjacency[j].add(i)
    distances=[]
    for i in range(12):
        dist=[None]*12;dist[i]=0;q=deque([i])
        while q:
            j=q.popleft()
            for k in adjacency[j]:
                if dist[k] is None:dist[k]=dist[j]+1;q.append(k)
        distances.append(dist)
    def cyc(face):return min(tuple(face[i:]+face[:i]) for i in range(3))
    oriented={cyc(list(f)) for f in faces}
    group=[]
    def visit(partial):
        i=len(partial)
        if i==12:
            if all(cyc([partial[j] for j in f]) in oriented for f in faces):group.append(tuple(partial))
            return
        for j in set(range(12))-set(partial):
            if all(distances[i][k]==distances[j][partial[k]] for k in range(i)):visit(partial+[j])
    visit([])
    need(len(group)==60,'independent proper-isometry count')
    anti=[row.index(3) for row in distances]
    pairs=sorted({tuple(sorted((i,anti[i]))) for i in range(12)});axis={j:i for i,p in enumerate(pairs) for j in p}
    pair_group=[tuple(axis[g[p[0]]] for p in pairs) for g in group]
    face=set(faces[0]);opp={anti[i] for i in face}
    H=[g for g in group if {g[i] for i in face} in (face,opp)]
    hp=[tuple(axis[g[p[0]]] for p in pairs) for g in H]
    def order(g):
        p=tuple(range(6))
        for n in range(1,7):
            p=tuple(g[p[i]] for i in range(6))
            if p==tuple(range(6)):return n
        raise ValueError('element order')
    orders=[order(g) for g in hp]
    chars=[sum(i==j for i,j in enumerate(g))-1 for g in pair_group]
    hc=[sum(i==j for i,j in enumerate(g))-1 for g in hp]
    irreps={'trivial':[1]*6,'sign':[-1 if n==2 else 1 for n in orders],'standard':[2 if n==1 else 0 if n==2 else -1 for n in orders]}
    multiplicity=lambda v:{k:int(sum(F(x*y,6) for x,y in zip(v,ch))) for k,ch in irreps.items()}
    P=[[F(int(i==j))-F(1,6) for j in range(6)] for i in range(6)]
    PH=[[sum(F(int(g[j]==i),6) for g in hp)-F(1,6) for j in range(6)] for i in range(6)]
    source_gram=[(F(1),F(0)),(F(0),F(1,5)),(F(0),F(-1,5)),(F(-1),F(0))]
    norm=(F(0),F(0))
    for i in face:
        for j in face:norm=field_add(norm,source_gram[distances[i][j]])
    line=[]
    for pair in pairs:
        z=(F(0),F(0))
        for j in face:z=field_add(z,source_gram[distances[pair[0]][j]])
        z=field_add(field_div(field_mul(z,z),norm),(F(-1,3),F(0)))
        need(z[0]==0 and abs(z[1])==F(2,15),'face-axis line')
        line.append(('' if z[1]>0 else '-')+'2*sqrt(5)/15')
    gram=[]
    for g in hp:
        inv=[g.index(i) for i in range(6)];row=[]
        for h in hp:row.append(F((sum(inv[h[i]]==i for i in range(6))-1)**2))
        gram.append(row)
    invgram=[[F(int(i==j),24)-F(1,720) for j in range(6)] for i in range(6)]
    need(mm(gram,invgram)==[[F(int(i==j)) for j in range(6)] for i in range(6)],'channel dual')
    quotient=Counter(tuple(sorted((axis[i],axis[j]))) for i,j in edges)
    # Native Laplacian descends on pair totals, with no choice of subgroup.
    L6=[[F(6*int(i==j)-1) for j in range(6)] for i in range(6)]
    for i,pair in enumerate(pairs):
        for vertex in range(12):
            image=5*int(vertex in pair)-sum(vertex in adjacency[p] for p in pair)
            need(image==L6[i][axis[vertex]],'native pair-total descent')
    increments=[]
    for i,j in sorted(edges):
        v=[F(int(k==axis[j])-int(k==axis[i])) for k in range(6)];increments.extend([v,[-x for x in v]])
    covariance=zeros(6)
    for v in increments:covariance=add(covariance,mulscale(outer(v),F(1,60)))
    need(covariance==mulscale(P,F(2,5)),'source alphabet covariance')
    profiles={'balanced_three_and_three':[1,1,1,-1,-1,-1],'dipole':[1,-1,0,0,0,0],'one_high_five_low':[5,-1,-1,-1,-1,-1]}
    moments={}
    for name,a in profiles.items():
        norm2=sum(x*x for x in a)
        m={k:sum(sum(x*y for x,y in zip(a,v))**k for v in increments)/60/F(norm2)**(k//2) for k in (2,4,6)}
        moments[name]={'second':str(m[2]),'fourth_cumulant':str(m[4]-3*m[2]**2),'sixth_cumulant':str(m[6]-15*m[4]*m[2]+30*m[2]**3)}
    dipole=[1]*12;dipole[0]=2;dipole[min(adjacency[0])]=0
    witnesses=[]
    for N in ([1]*12,dipole,list(range(12))):
        steps=[];missing=zeros(6)
        for i,j in sorted(edges):
            v=[F(int(k==axis[j])-int(k==axis[i])) for k in range(6)]
            if N[i]-N[j]>=2:steps.append(v)
            elif N[j]-N[i]>=2:steps.append([-x for x in v])
            else:steps.append([F(0)]*6);missing=add(missing,mulscale(outer(v),F(1,30)))
        mean=[sum(v[i] for v in steps)/30 for i in range(6)];M=zeros(6)
        for v in steps:M=add(M,mulscale(outer(v),F(1,30)))
        Cov=add(M,mulscale(outer(mean),F(-1)))
        need(add(M,missing)==mulscale(P,F(2,5)),'accepted bound factorization')
        witnesses.append({'counts':N,'mean':list(map(str,mean)),'raw_second':encode(M),'covariance':encode(Cov),'bound_remainder':encode(missing),'accepted_seams':sum(any(v) for v in steps)})
    return {
     'geometry':{'ports':12,'seams':len(edges),'A5_order':len(group),'antipodal_pairs':[list(p) for p in pairs],'axis_pair_seam_multiplicity':sorted(quotient.values()),'D3_face_ports':sorted(face),'D3_order':len(H),'D3_element_orders':{str(k):v for k,v in Counter(orders).items()}},
     'carrier':{'dimension':5,'D3_line_complement_dimension':4,'A5_character_counts':{str(k):v for k,v in Counter(chars).items()},'character_norm':str(sum(F(c*c,60) for c in chars)),'D3_multiplicities':multiplicity(hc),'pair_centering_projector':encode(P),'quadrupole_Gram':encode(mulscale(P,F(4,5))),'D3_invariant_line':line,'D3_line_projector':encode(PH),'invariant_Hermitian_functional_dimensions':{'D3':sum(x*x for x in multiplicity(hc).values()),'A5':int(sum(F(c*c,60) for c in chars)),'A5_after_normalization':0},'normalized_native_projector_weights':['1/5','4/5']},
     'operators':{'EndW5_D3_multiplicities':multiplicity([x*x for x in hc]),'native_W5_laplacian_eigenvalues':{'6':5},'native_pair_laplacian':encode(L6),'D3_transposition_W5_spectrum':{'0':1,'3':4},'D3_transposition_EndW5_spectrum':{'0':5,'3':16,'6':4},'D3_transposition_regular_spectrum':{'0':1,'3':4,'6':1},'channel_Gram':encode(gram),'channel_Gram_inverse':encode(invgram),'central_rate_eigenvalue_coefficients':[[6,0],[3,3]],'faithful_group_heat_channel_identity':True},
     'proposal':{'undirected_seam_count':30,'symmetrized_oriented_count':len(increments),'distinct_pair_total_increments':len({tuple(v) for v in increments}),'symmetrized_covariance':encode(covariance),'raw_moment_coefficients':{'second_S2':'2/5','fourth_S4':'2/5','fourth_S2_squared':'1/5','fourth_cumulant_S4':'2/5','fourth_cumulant_S2_squared':'-7/25'},'symmetrized_unit_direction_examples':moments,'native_accepted_step_witnesses':witnesses}}

def reconstruct_color():
    # SU(n) completeness contracted over the middle indices gives C_F.
    n=3;CF=F(n*n-1,2*n);p=F(1,5)
    moment={k:p*CF**k for k in range(1,5)};cum={}
    for k in range(1,5):cum[k]=moment[k]-sum(F(math.comb(k-1,j-1))*cum[j]*moment[k-j] for j in range(1,k))
    return {'Casimir':encode([[CF*int(i==j) for j in range(n)] for i in range(n)]),'A_rank':3,'A_nonzero_eigenvalue':str(CF),'normalized_trace':str(p*CF),'normalized_HS_square':str(p*CF*CF),'spectral_values':['0',str(CF)],'spectral_probabilities':[str(1-p),str(p)],'scalar_cumulants':[str(cum[k]) for k in range(1,5)],'twirled_rank_r_first_trace_coefficient':'1','twirled_rank_r_HS_square_coefficient':str(F(1,3)),'conditional_Gaussian_coefficients':{'d432_rank1_untwirled':str(F(1,432)),'d432_rank1_twirl_before_response':str(F(1,3*432)),'d15_rank4_untwirled':str(F(4,15)),'d15_rank4_twirl_before_response':str(F(4,3*15)),'d15_Casimir_line_response':str(CF*CF/5)}}

def verify(payload,repo=REPO):
    repo=Path(repo);specpath=repo/SPEC_REL
    need(sha(specpath)==SPEC_SHA256,'specification pin')
    spec=json.loads(specpath.read_text())
    need(json.dumps(payload.get('specification'),sort_keys=True)==json.dumps(spec,sort_keys=True),'receipt scientific scope or specification')
    a3=json.loads((repo/SOURCE_FILES[1]).read_text());repair=json.loads((repo/SOURCE_FILES[2]).read_text())
    need(a3['a3_selection']['reference']=='uniform_on_moves' and a3['a3_selection']['selected_kernel_probability_per_seam']=='1/30','source proposal law')
    need(repair['mechanism']['repair_move']=='conservative unit transfer across a seam with difference magnitude at least two','source admission law')
    audit_source((repo/PRODUCER_REL).read_text())
    expected={'schema':'oph.source_w5_response_constraints.v1','status':'EXACT_CONSTRAINTS_WITH_EXPLICIT_RESPONSE_HYPOTHESES','specification':spec,'source_pins':{p:sha(repo/p) for p in SOURCE_FILES+[PRODUCER_REL,SPEC_REL]},**reconstruct_geometry((repo/SOURCE_FILES[0]).read_bytes()),'color':reconstruct_color(),'exact_checks':dict.fromkeys(['source_geometry','A5_irreducibility','canonical_quadrupole_intertwiner','D3_restriction','operator_carrier','positive_central_jump_family','native_pair_laplacian','proposal_moment_identities','accepted_covariance_factorizations','color_Casimir','color_twirl','scalar_MGF_cumulants'],True)}
    need(set(payload)==set(expected),'receipt fields')
    for k in expected:need(json.dumps(payload[k],sort_keys=True)==json.dumps(expected[k],sort_keys=True),f'independent {k} reconstruction mismatch')
    return True

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path,nargs='?',default=DEFAULT_INPUT);a=p.parse_args()
    verify(json.loads(a.input.read_text()));print('source W5 response constraints: independent rational verification passed');return 0
if __name__=='__main__':raise SystemExit(main())
