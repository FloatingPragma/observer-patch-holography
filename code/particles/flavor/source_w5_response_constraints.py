#!/usr/bin/env python3
"""Exact native W5 response constraints and explicitly conditional readouts.

This constructs the source carrier from oriented incidence, not quark data.
The D3 background/jump law and color-response interpretation remain supplied.
Native accepted-step covariance bounds follow from sums of positive outer
products; symmetrized alphabet moments assert neither IID histories nor a clock.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,sys
from collections import Counter,deque
from pathlib import Path
import sympy as s
REPO=Path(__file__).resolve().parents[3]
PRODUCER_REL='code/particles/flavor/source_w5_response_constraints.py'
SPEC_REL='code/particles/flavor/source_w5_response_constraints_spec.json'
DEFAULT_OUTPUT=REPO/'code/particles/runs/flavor/source_w5_response_constraints.json'
SOURCE_FILES=['code/a5_closure/manifests/echosahedral_federation_reference.json', 'code/a5_closure/manifests/a3_scheduler_kernel_reference.json', 'code/a5_closure/manifests/record_counting_mechanism_reference.json', 'code/a5_closure/source_repair_generator_certificate.py']
SPEC_SHA256='71a112358ff9e24efa843cfebec96678b0c9d90e428467ce080f5da8642bcb00'
sys.path.insert(0,str(REPO/'code/a5_closure'))
import source_repair_generator_certificate as src

def require(ok,why):
    if not ok:raise ValueError(why)
def sha256(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def matrix(M):return [[str(s.simplify(x)) for x in row] for row in M.tolist()]
def read_spec(path=None):
    path=REPO/SPEC_REL if path is None else Path(path)
    require(sha256(path)==SPEC_SHA256,'source-response specification differs from reviewed law')
    return json.loads(path.read_text())
def build():
    spec=read_spec()

    a3=json.loads((REPO/SOURCE_FILES[1]).read_text())
    require(a3['a3_selection']['reference']=='uniform_on_moves','uniform source reference')
    require(a3['a3_selection']['selected_kernel_probability_per_seam']=='1/30','native seam probability')
    repair=json.loads((REPO/SOURCE_FILES[2]).read_text())
    require(repair['mechanism']['repair_move']=='conservative unit transfer across a seam with difference magnitude at least two','native admission')
    manifest=REPO/'code/a5_closure/manifests/echosahedral_federation_reference.json'
    doc=json.loads(manifest.read_text());ports,faces,edges=src.data_from_faces(doc)
    G=src.proper_autos(12,edges,faces);assert len(G)==60
    adj=[set() for _ in ports]
    for i,j in edges:adj[i].add(j);adj[j].add(i)
    dist=[]
    for i in range(12):
     d={i:0};q=deque([i])
     while q:
      j=q.popleft()
      for k in adj[j]:
       if k not in d:d[k]=d[j]+1;q.append(k)
     dist.append(d)
    anti=[next(j for j in range(12) if dist[i][j]==3) for i in range(12)]
    pairs=sorted({tuple(sorted((i,anti[i]))) for i in range(12)})
    axis={i:k for k,pair in enumerate(pairs) for i in pair};assert len(pairs)==6
    pair_G=[tuple(axis[g[pair[0]]] for pair in pairs) for g in G]
    assert len(set(pair_G))==60
    P=s.eye(6)-s.ones(6)/6

    def perm(g):
     return s.Matrix(6,6,lambda i,j:int(i==g[j]))
    B=[perm(g)*P for g in pair_G]
    chars=[s.trace(b) for b in B]
    assert sum(c*c for c in chars)==60 # irreducible A5 carrier, not just dimension coincidence
    assert Counter(chars)=={s.Integer(5):1,s.Integer(1):15,s.Integer(-1):20,s.Integer(0):24}
    face=frozenset(faces[0]);opposite=frozenset(anti[i] for i in face)
    Hix=[i for i,g in enumerate(G) if frozenset(g[j] for j in face) in (face,opposite)]
    assert len(Hix)==6
    H=[B[i] for i in Hix]

    def order(g):
     z=tuple(range(6));n=0
     while True:
      z=tuple(g[z[i]] for i in range(6));n+=1
      if z==tuple(range(6)):return n
    orders=[order(pair_G[i]) for i in Hix]
    assert Counter(orders)=={1:1,2:3,3:2}
    sign=[-1 if n==2 else 1 for n in orders]
    std=[2 if n==1 else 0 if n==2 else -1 for n in orders]
    chi=[s.trace(b) for b in H]
    mult={name:sum(a*b for a,b in zip(chi,c))/6 for name,c in [('trivial',[1]*6),('sign',sign),('standard',std)]}
    assert mult=={'trivial':1,'sign':0,'standard':2}
    PH=sum(H,s.zeros(6))/6
    assert PH*PH==PH and PH.rank()==1
    # Canonical quadrupole intertwiner T(b)=sum b_i(n_i n_i^T-I/3),
    # reconstructed from the incidence Gram, without choosing a spatial frame.
    port_gram=s.Matrix(12,12,lambda i,j: [1,s.sqrt(5)/5,-s.sqrt(5)/5,-1][dist[i][j]])
    frame_gram=s.Matrix(6,6,lambda i,j:s.simplify(port_gram[pairs[i][0],pairs[j][0]]**2-s.Rational(1,3)))
    assert frame_gram==s.Rational(4,5)*P
    norm_face=s.simplify(sum(port_gram[i,j] for i in face for j in face))
    line=s.Matrix([s.simplify(sum(port_gram[pair[0],j] for j in face)**2/norm_face-s.Rational(1,3)) for pair in pairs])
    assert sum(line)==0 and PH==s.simplify(line*line.T/line.dot(line))
    assert set(line)=={-2*s.sqrt(5)/15,2*s.sqrt(5)/15}
    # Operator conjugation carrier: its characters square those of W5.
    chi_end=[c*c for c in chi]
    end_mult={name:sum(a*b for a,b in zip(chi_end,c))/6 for name,c in [('trivial',[1]*6),('sign',sign),('standard',std)]}
    assert end_mult=={'trivial':5,'sign':4,'standard':8}
    EE=s.kronecker_product(P,P)
    Ad=[s.kronecker_product(b,b) for b in H]
    # Actual operator superoperator Gram: every nonidentity restricted W5 character is +/-1.
    channel_gram=s.Matrix(6,6,lambda i,j:s.trace(H[i].T*H[j])**2)
    assert channel_gram==24*s.eye(6)+s.ones(6)
    assert channel_gram.inv()==s.eye(6)/24-s.ones(6)/720
    LF=3*P-sum((h for h,n in zip(H,orders) if n==2),s.zeros(6))
    LE=3*EE-sum((h for h,n in zip(Ad,orders) if n==2),s.zeros(36))
    assert LF.eigenvals()=={0:2,3:4}
    # Exactly form the three central conjugation projectors; ambient complement discarded.
    E0=sum(Ad,s.zeros(36))/6
    Es=sum((x*h for x,h in zip(sign,Ad)),s.zeros(36))/6
    E2=EE-E0-Es
    assert all(e*e==e for e in (E0,Es,E2))
    assert [s.trace(e) for e in (E0,Es,E2)]==[5,4,16]
    assert LE==6*Es+3*E2
    c2,c3=s.symbols('c2 c3',nonnegative=True)
    Lcentral=c2*LE+c3*(2*EE-sum((h for h,n in zip(Ad,orders) if n==3),s.zeros(36)))
    assert Lcentral==6*c2*Es+3*(c2+c3)*E2
    # Faithful random-unitary image of the group heat law does not require
    # a regular submodule inside F. Its operator heat trace is a different readout.
    r,u=s.symbols('r u',real=True)
    probabilities=[(1+4*r+u)/6 if n==1 else (1-u)/6 if n==2 else (1-2*r+u)/6 for n in orders]
    channel=sum((p*h for p,h in zip(probabilities,Ad)),s.zeros(36))
    assert (channel-(E0+r*E2+u*Es)).applyfunc(s.simplify)==s.zeros(36)
    Hperms=[pair_G[i] for i in Hix];lookup={g:i for i,g in enumerate(Hperms)}
    regular=s.zeros(6)
    for j,g in enumerate(Hperms):
        regular[j,j]=3
        for h,n in zip(Hperms,orders):
            if n==2:regular[lookup[tuple(h[g[k]] for k in range(6))],j]-=1
    assert regular.eigenvals()=={0:1,3:4,6:1}
    # Full native 30-seam graph does not select H or LE.
    L12=s.diag(*[len(z) for z in adj])-s.Matrix(12,12,lambda i,j:int(j in adj[i]))
    Totals=s.Matrix(6,12,lambda i,j:int(j in pairs[i]))
    L6=6*s.eye(6)-s.ones(6)
    assert Totals*L12==L6*Totals and L6*P==6*P
    quotient=Counter(tuple(sorted((axis[i],axis[j]))) for i,j in edges)
    assert len(quotient)==15 and set(quotient.values())=={2}
    # Each oriented proposal is one native transfer. Pair TOTALS need no paired simultaneous move.
    increments=[]
    for i,j in sorted(edges):
     for source,target in ((i,j),(j,i)):
      v=Totals[:,target]-Totals[:,source];increments.append(v)
    assert len(increments)==60 and len({tuple(v) for v in increments})==30
    cov=sum((v*v.T for v in increments),s.zeros(6))/60
    assert cov==s.Rational(2,5)*P
    # Exact first four moments and fourth cumulant for every centered direction.
    a=s.symbols('a0:5');av=s.Matrix([*a,-sum(a)])
    second=s.expand(sum((av.dot(v))**2 for v in increments)/60)
    fourth=s.expand(sum((av.dot(v))**4 for v in increments)/60)
    S2=sum(z*z for z in av);S4=sum(z**4 for z in av)
    assert s.expand(second-s.Rational(2,5)*S2)==0
    assert s.expand(fourth-s.Rational(2,5)*S4-s.Rational(1,5)*S2**2)==0
    k4=s.expand(fourth-3*second**2)
    assert s.expand(k4-s.Rational(2,5)*S4+s.Rational(7,25)*S2**2)==0
    profiles={
     'balanced_three_and_three':s.Matrix([1,1,1,-1,-1,-1])/s.sqrt(6),
     'dipole':s.Matrix([1,-1,0,0,0,0])/s.sqrt(2),
     'one_high_five_low':s.Matrix([5,-1,-1,-1,-1,-1])/s.sqrt(30)}
    moments={}
    for key,v in profiles.items():
     m={k:s.simplify(sum(v.dot(z)**k for z in increments)/60) for k in (2,4,6)}
     moments[key]={'second':str(m[2]),'fourth_cumulant':str(s.simplify(m[4]-3*m[2]**2)),'sixth_cumulant':str(s.simplify(m[6]-15*m[4]*m[2]+30*m[2]**3))}
    assert moments['one_high_five_low']['fourth_cumulant']=='0' and moments['one_high_five_low']['sixth_cumulant']=='-48/125'
    # A physically typed color quadratic response reproduces one old rational
    # without pretending a rank-four sharp color-invariant projector exists.
    E=lambda i,j:s.Matrix(3,3,lambda a,b:int((a,b)==(i,j)))
    Ta=[(E(0,1)+E(1,0))/2,(-s.I*E(0,1)+s.I*E(1,0))/2,
        (E(0,0)-E(1,1))/2,(E(0,2)+E(2,0))/2,
        (-s.I*E(0,2)+s.I*E(2,0))/2,(E(1,2)+E(2,1))/2,
        (-s.I*E(1,2)+s.I*E(2,1))/2,
        (E(0,0)+E(1,1)-2*E(2,2))/(2*s.sqrt(3))]
    assert all(s.trace(a*b)==s.Rational(int(i==j),2) for i,a in enumerate(Ta) for j,b in enumerate(Ta))
    Casimir=sum((a*a for a in Ta),s.zeros(3)).applyfunc(s.simplify)
    assert Casimir==s.Rational(4,3)*s.eye(3)
    A=s.kronecker_product(Casimir,PH)
    assert s.trace(A)/15==s.Rational(4,15)
    assert s.trace(A*A)/15==s.Rational(16,45)
    assert A*A==s.Rational(4,3)*A and A.rank()==3
    assert s.trace(PH)/5==s.Rational(1,5) and s.trace(P-PH)/5==s.Rational(4,5)
    # Finite exact color twirl: diagonal phase and cyclic shift are both SU(3).
    omega=(-1+s.sqrt(3)*s.I)/2
    Z=s.diag(1,omega,omega**2);C=s.Matrix([[0,0,1],[1,0,0],[0,1,0]])
    X=s.Matrix(3,3,s.symbols('x0:9'))
    twirled=sum(((C**i)*(Z**j)*X*(Z**j).conjugate().T*(C**i).T for i in range(3) for j in range(3)),s.zeros(3))/9
    twirled=twirled.applyfunc(s.simplify)
    assert twirled==s.trace(X)*s.eye(3)/3

    w=s.symbols('w', real=True)
    scalar_mgf=(4+s.exp(4*w/3))/5
    cumulants=[s.simplify(s.diff(s.log(scalar_mgf),w,k).subs(w,0)) for k in range(1,5)]
    assert cumulants==[s.Rational(4,15),s.Rational(64,225),s.Rational(256,1125),s.Rational(1024,50625)]
    xs=s.symbols('z0:6', nonzero=True)
    assert s.expand(sum(xs[i]/xs[j] for i in range(6) for j in range(6) if i!=j)-sum(xs)*sum(1/x for x in xs)+6)==0
    witnesses=[]
    neighbor=min(adj[0]);dipole=[1]*12;dipole[0]=2;dipole[neighbor]=0
    for N in ([1]*12,dipole,list(range(12))):
        steps=[];missing=s.zeros(6)
        for i,j in sorted(edges):
            v=Totals[:,j]-Totals[:,i]
            if N[i]>=N[j]+2:step=v
            elif N[j]>=N[i]+2:step=-v
            else:step=s.zeros(6,1);missing+=v*v.T/30
            steps.append(step)
        mean=sum(steps,s.zeros(6,1))/30
        rawsecond=sum((z*z.T for z in steps),s.zeros(6))/30
        covariance=rawsecond-mean*mean.T
        assert s.Rational(2,5)*P-rawsecond==missing
        assert covariance==sum(((z-mean)*(z-mean).T for z in steps),s.zeros(6))/30
        witnesses.append({'counts':N,'mean':list(map(str,mean)),'raw_second':matrix(rawsecond),'covariance':matrix(covariance),'bound_remainder':matrix(missing),'accepted_seams':sum(z!=s.zeros(6,1) for z in steps)})
    assert witnesses[0]['accepted_seams']==0
    pins={name:sha256(REPO/name) for name in SOURCE_FILES+[PRODUCER_REL,SPEC_REL]}
    return {
      'schema':'oph.source_w5_response_constraints.v1','status':'EXACT_CONSTRAINTS_WITH_EXPLICIT_RESPONSE_HYPOTHESES',
      'specification':spec,'source_pins':pins,
      'geometry':{'ports':12,'seams':len(edges),'A5_order':len(G),'antipodal_pairs':[list(x) for x in pairs],'axis_pair_seam_multiplicity':sorted(quotient.values()),'D3_face_ports':sorted(face),'D3_order':len(H),'D3_element_orders':{str(k):v for k,v in Counter(orders).items()}},
      'carrier':{'dimension':5,'D3_line_complement_dimension':4,'A5_character_counts':{str(k):v for k,v in Counter(chars).items()},'character_norm':str(sum(c*c for c in chars)/60),'D3_multiplicities':{k:int(v) for k,v in mult.items()},'pair_centering_projector':matrix(P),'quadrupole_Gram':matrix(frame_gram),'D3_invariant_line':list(map(str,line)),'D3_line_projector':matrix(PH),'invariant_Hermitian_functional_dimensions':{'D3':5,'A5':1,'A5_after_normalization':0},'normalized_native_projector_weights':['1/5','4/5']},
      'operators':{'EndW5_D3_multiplicities':{k:int(v) for k,v in end_mult.items()},'native_W5_laplacian_eigenvalues':{'6':5},'native_pair_laplacian':matrix(L6),'D3_transposition_W5_spectrum':{'0':1,'3':4},'D3_transposition_EndW5_spectrum':{'0':5,'3':16,'6':4},'D3_transposition_regular_spectrum':{'0':1,'3':4,'6':1},'channel_Gram':matrix(channel_gram),'channel_Gram_inverse':matrix(channel_gram.inv()),'central_rate_eigenvalue_coefficients':[[6,0],[3,3]],'faithful_group_heat_channel_identity':True},
      'proposal':{'undirected_seam_count':30,'symmetrized_oriented_count':len(increments),'distinct_pair_total_increments':len({tuple(v) for v in increments}),'symmetrized_covariance':matrix(cov),'raw_moment_coefficients':{'second_S2':'2/5','fourth_S4':'2/5','fourth_S2_squared':'1/5','fourth_cumulant_S4':'2/5','fourth_cumulant_S2_squared':'-7/25'},'symmetrized_unit_direction_examples':moments,'native_accepted_step_witnesses':witnesses},
      'color':{'Casimir':matrix(Casimir),'A_rank':int(A.rank()),'A_nonzero_eigenvalue':'4/3','normalized_trace':'4/15','normalized_HS_square':'16/45','spectral_values':['0','4/3'],'spectral_probabilities':['4/5','1/5'],'scalar_cumulants':list(map(str,cumulants)),'twirled_rank_r_first_trace_coefficient':'1','twirled_rank_r_HS_square_coefficient':'1/3','conditional_Gaussian_coefficients':{'d432_rank1_untwirled':'1/432','d432_rank1_twirl_before_response':'1/1296','d15_rank4_untwirled':'4/15','d15_rank4_twirl_before_response':'4/45','d15_Casimir_line_response':'16/45'}},
      'exact_checks':{'source_geometry':True,'A5_irreducibility':True,'canonical_quadrupole_intertwiner':True,'D3_restriction':True,'operator_carrier':True,'positive_central_jump_family':True,'native_pair_laplacian':True,'proposal_moment_identities':True,'accepted_covariance_factorizations':True,'color_Casimir':True,'color_twirl':True,'scalar_MGF_cumulants':True}}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT);a=p.parse_args()
    value=build();raw=json.dumps(value,sort_keys=True,indent=2)+'\n'
    if a.check:require(a.output.read_text()==raw,'source W5 response receipt is stale')
    else:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(raw)
    print('source W5 response constraints: '+('exact rebuild matches' if a.check else 'receipt written'))
    return 0
if __name__=='__main__':raise SystemExit(main())
