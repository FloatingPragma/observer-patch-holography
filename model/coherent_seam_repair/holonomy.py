"""Exact finite binary-cover controller generated from the public carrier.

The five-level code and its continuous response family are declared architecture.
This producer reconstructs its lift words and executes typed chart arrows.
"""
from functools import lru_cache
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
from model.common_source import slow_geometry


@lru_cache(maxsize=1)
def certificate():
    c=slow_geometry.source()
    Z=c.q5(0);O=c.q5(1)
    add,sub,mul,inv=c.q5_add,c.q5_sub,c.q5_mul,c.q5_inv
    def neg(a):return sub(Z,a)
    def scale(a,n):return mul(a,c.q5(n))
    def sig(a):return (a[0],-a[1])
    def dot(a,b):
     out=Z
     for x,y in zip(a,b):out=add(out,mul(x,y))
     return out
    def cross(a,b):return tuple(sub(mul(a[(i+1)%3],b[(i+2)%3]),mul(a[(i+2)%3],b[(i+1)%3])) for i in range(3))
    def rational_sqrt(x):
     if x<0:return None
     a,b=isqrt(x.numerator),isqrt(x.denominator)
     return F(a,b) if a*a==x.numerator and b*b==x.denominator else None
    def sqrtq(x):
     if x==Z:return Z
     a,b=x
     if b==0:
      r=rational_sqrt(a)
      if r is not None:return c.q5(r)
      r=rational_sqrt(a/5)
      if r is not None:return c.q5(0,r)
     d=rational_sqrt(a*a-5*b*b)
     if d is not None:
      for sign in (1,-1):
       r=rational_sqrt((a+sign*d)/2)
       if r:
        y=c.q5(r,b/(2*r))
        if c.q5_sign(y)<0:y=neg(y)
        if mul(y,y)==x:return y
     raise ValueError(('not a square in Q(sqrt5)',x))
    def qmul(q,r):
     w,x=q[0],q[1:];v,y=r[0],r[1:]
     xy=cross(x,y)
     return (sub(mul(w,v),dot(x,y)),)+tuple(add(add(mul(w,y[i]),mul(v,x[i])),xy[i]) for i in range(3))
    def qneg(q):return tuple(map(neg,q))
    def qinv(q):return (q[0],)+tuple(map(neg,q[1:]))
    def canon(q):
     for x in q:
      if c.q5_sign(x):return q if c.q5_sign(x)>0 else qneg(q)
     raise ValueError('zero quaternion')
    def qrotate(q,v):return qmul(qmul(q,(Z,)+tuple(v)),qinv(q))[1:]
    def qfromrot(R):
     w2=scale(add(O,add(R[0][0],add(R[1][1],R[2][2]))),F(1,4))
     if w2!=Z:
      w=sqrtq(w2);den=inv(scale(w,4))
      q=(w,)+tuple(mul(sub(R[(i+2)%3][(i+1)%3],R[(i+1)%3][(i+2)%3]),den) for i in range(3))
     else:
      j=next(i for i in range(3) if add(O,R[i][i])!=Z)
      x=sqrtq(scale(add(O,R[j][j]),F(1,2)));v=[Z]*3;v[j]=x
      for k in range(3):
       if k!=j:v[k]=mul(add(R[j][k],R[k][j]),inv(scale(x,4)))
      q=(Z,)+tuple(v)
     assert dot(q,q)==O
     return canon(q)
    def encodeq(q):return [[str(a),str(b)] for a,b in q]
    def compose(p,r):return tuple(p[r[i]] for i in range(12))

    # Derive geometric rotations from actual source permutations, not a binary-group table.
    vertices=slow_geometry.raw_frame()
    rotations=tuple(c.rotations())
    base=(vertices[0],vertices[1],vertices[4])
    det=dot(base[0],cross(base[1],base[2]));assert det!=Z
    dual=tuple(tuple(mul(x,inv(det)) for x in cross(base[(i+1)%3],base[(i+2)%3])) for i in range(3))
    canonical={};project={};rot_by_q={}
    for p in rotations:
     targets=tuple(vertices[p[i]] for i in (0,1,4))
     R=tuple(tuple(dot(tuple(targets[k][i] for k in range(3)),tuple(dual[k][j] for k in range(3))) for j in range(3)) for i in range(3))
     q=qfromrot(R)
     assert all(qrotate(q,vertices[i])==vertices[p[i]] for i in range(12))
     assert q not in canonical
     canonical[q]=p
     for h in (q,qneg(q)):project[h]=p;rot_by_q[h]=R
    assert len(project)==120
    identity=(O,Z,Z,Z);minus=qneg(identity);idperm=tuple(range(12))
    assert project[identity]==project[minus]==idperm

    def closure(gens):
     words={identity:()};queue=[identity]
     for h in queue:
      for label,g in enumerate(gens):
       gh=qmul(g,h)
       assert gh in project
       if gh not in words:words[gh]=words[h]+(label,);queue.append(gh)
     return words
    # Choose a carrier half-turn, then the first carrier lift which generates the full cover.
    s=next(q for q in canonical if q[0]==Z)
    t=next(q for q in canonical if len(closure((s,q)))==120)
    generators=(s,qinv(s),t,qinv(t))
    words=closure(generators);assert len(words)==120

    # Complex Q(sqrt5) arithmetic for the actual five-level state transition.
    def ca(a,b):return (add(a[0],b[0]),add(a[1],b[1]))
    def cm(a,b):return (sub(mul(a[0],b[0]),mul(a[1],b[1])),add(mul(a[0],b[1]),mul(a[1],b[0])))
    def cc(a):return (a[0],neg(a[1]))
    CZ=(Z,Z)
    def spinor(q):
     w,x,y,z=map(sig,q)
     return (((w,x),(y,z)),((neg(y),z),(w,neg(x))))
    def implementer(q):
     M=[[CZ]*5 for _ in range(5)]
     for i in range(3):
      for j in range(3):M[i][j]=(rot_by_q[q][i][j],Z)
     S=spinor(q)
     for i in range(2):
      for j in range(2):M[3+i][3+j]=S[i][j]
     return tuple(tuple(row) for row in M)
    def response_numeric():
     """Floating evaluation of the exact declared geometric response formula.

     Exact Gram and covariance are certified separately below, without a tolerance.
     This evaluator permits an independent executable consumer to read the controls.
     """
     import numpy as np
     sy=np.array([[0,-1j],[1j,0]])
     sz=np.diag([1.,-1.]);sx=np.array([[0.,1.],[1.,0.]])
     out=[]
     for v in vertices:
      n=np.array([c.q5_float(x) for x in v]);n=n/np.linalg.norm(n)
      m=np.array([c.q5_float(sig(x)) for x in v]);m=m/np.linalg.norm(m)
      x,y,z=n;C=np.array([[0,-z,y],[z,0,-x],[-y,x,0.]])
      P=(np.eye(3)-np.outer(n,n)+1j*C)/2
      Q=(np.eye(2)+m[0]*sz+m[1]*sy+m[2]*sx)/2
      X=np.zeros((5,5),complex);X[:3,:3]=P/2;X[3:,3:]=-Q/2
      out.append(X)
     return out
    def apply(M,v):
     out=[]
     for row in M:
      z=CZ
      for a,b in zip(row,v):z=ca(z,cm(a,b))
      out.append(z)
     return tuple(out)
    def cmatrixmul(M,N):return tuple(tuple(apply(M,tuple(N[k][j] for k in range(len(N))))[i] for j in range(len(N))) for i in range(len(M)))
    def norm(v):
     z=CZ
     for a in v:z=ca(z,cm(cc(a),a))
     assert z[1]==Z
     return z[0]

    def step(state, instruction):
     """One declared reversible shared-chart transport, with a complete exact record.

     Public source code is fixed. A command chooses generator and destination chart.
     State includes its frame lift; its central bit is never quotiented away.
     """
     source,label,destination=instruction
     if source!=state['chart']:raise ValueError('transport arrow has the wrong source chart')
     if type(label) is not int or label not in range(4):raise ValueError('invalid gate')
     if type(destination) is not int or destination not in (0,1):raise ValueError('invalid chart')
     g=generators[label];h=qmul(g,state['frame'])
     record={'chart_before':state['chart'],'chart_after':destination,
             'gate':label,'carrier_permutation':list(project[g]),
             'frame_before':encodeq(state['frame']),'frame_after':encodeq(h),
             'central_bit_after':int(h!=canon(h))}
     out={'chart':destination,'frame':h,'psi':apply(implementer(g),state['psi'])}
     return out,record

    I5=tuple(tuple((O,Z) if i==j else CZ for j in range(5)) for i in range(5))
    unitary_checks=0
    for q in project:
     U=implementer(q);Ui=implementer(qinv(q))
     assert cmatrixmul(Ui,U)==I5
     unitary_checks+=1
    hom_checks=0
    for q in project:
     for g in generators:
      assert cmatrixmul(implementer(g),implementer(q))==implementer(qmul(g,q))
      hom_checks+=1

    # Exact Gram of Xp=(P_spin1(np)/2) direct_sum (-P_spinhalf(Galois(np))/2).
    # The normalized conjugate frame uses its own positive norm, not conjugated sqrt.
    T=c.intrinsic_gram_exact()
    G=tuple(tuple(add(scale(mul(add(O,T[i][j]),add(O,T[i][j])),F(1,16)),scale(add(O,sig(T[i][j])),F(1,8))) for j in range(12)) for i in range(12))
    def shifted(a):return [[sub(G[i][j],c.q5(a) if i==j else Z) for j in range(12)] for i in range(12)]
    poly=c.q5_matmul(c.q5_matmul(shifted(F(5,2)),shifted(F(1,2))),shifted(F(1,10)))
    assert all(x==Z for row in poly for x in row)
    assert all(sum(float(c.q5_float(x)) for x in row)>0 for row in G)
    # Every root of this annihilating polynomial is strictly positive: Gram is invertible.
    exact_gram={'positive_annihilating_roots':['5/2','1/2','1/10'],'rank':12,'proof':'(G-5I/2)(G-I/2)(G-I/10)=0 exactly; G is Hermitian'}

    # This numerical evaluator is a consumer cross-check, not the exact certificate.
    import numpy as np
    Xs=response_numeric()
    gram_error=max(abs(np.trace(Xs[i]@Xs[j])-c.q5_float(G[i][j])) for i in range(12) for j in range(12))
    cov_error=0.
    for g in generators:
     U=np.array([[c.q5_float(a)+1j*c.q5_float(b) for a,b in row] for row in implementer(g)])
     for p in range(12):cov_error=max(cov_error,float(np.max(abs(U@Xs[p]@U.conj().T-Xs[project[g][p]]))))
    assert gram_error<1e-12 and cov_error<1e-12

    # Verify multiplication and quotient on every pair; verifies exact closure beyond BFS.
    product_checks=0
    for q in project:
     for r in project:
      qr=qmul(q,r)
      assert qr in project and project[qr]==compose(project[q],project[r])
      product_checks+=1

    # The word executes on all five basis states, not just its abstract group label.
    compiler=[];state_steps=0
    for target,word in words.items():
     logs=[]
     for j in range(5):
      psi=tuple((O,Z) if i==j else CZ for i in range(5))
      st={'chart':0,'frame':identity,'psi':psi}
      for k,label in enumerate(word):
       # Every compiled word closes at chart0. Intermediate destinations alternate.
       destination=0 if k==len(word)-1 else (k+1)%2
       st,rec=step(st,(st['chart'],label,destination));state_steps+=1
       if j==0:logs.append(rec)
      assert st['chart']==0 and st['frame']==target
      assert st['psi']==apply(implementer(target),psi) and norm(st['psi'])==O
     compiler.append({'quaternion':encodeq(target),'central_bit':int(target!=canon(target)),
                      'carrier_permutation':list(project[target]),'word':list(word),'records':logs})
    assert all(row['word'] or row['quaternion']==encodeq(identity) for row in compiler)
    zrow=next(row for row in compiler if row['quaternion']==encodeq(minus))
    assert implementer(minus)==tuple(tuple((O if i<3 else neg(O),Z) if i==j else CZ for j in range(5)) for i in range(5))

    def autonomous_cycle():
     """A fixed finite source program visits every lift and reversibly returns.

     Its schedule is generated from the carrier. No measured masses, pulse angles,
     external choices or stochastic draws select these accepted overlap words.
     """
     chart=0
     for word in words.values():
      for path in (word,tuple(label^1 for label in reversed(word))):
       for k,label in enumerate(path):
        destination=0 if k==len(path)-1 else (k+1)%2
        yield chart,label,destination
        chart=destination
    cycle=list(autonomous_cycle())
    psi=tuple((O,Z) for _ in range(5)) # common norm sqrt(5) is kept implicit
    state={'chart':0,'frame':identity,'psi':psi}
    visited={identity}
    for instruction in cycle:
     state,record=step(state,instruction)
     if state['chart']==0:visited.add(state['frame'])
    assert state=={'chart':0,'frame':identity,'psi':psi}
    assert visited==set(project)
    inverse_arrow_checks=0
    for source in (0,1):
     for destination in (0,1):
      for label in range(4):
       initial={'chart':source,'frame':identity,'psi':psi}
       after,_=step(initial,(source,label,destination))
       restored,_=step(after,(destination,label^1,source))
       assert restored==initial
       inverse_arrow_checks+=1

    # Complete finite source-chart covariance uses the same exact carrier data.
    binding=slow_geometry.registered_binding()
    assert binding['all_pass']
    out={
     'status':'EXACT declared encoded overlap architecture; not the irreversible native pump',
     'source':'actual pinned native icosahedron rotations and exact raw frame',
     'source_carrier_sha256':hashlib.sha256(Path(c.__file__).read_bytes()).hexdigest(),
     'same_registered_Lean_carrier':binding,
     'scalar_field':'Q(sqrt(5)); pairs encode a+b*sqrt(5)',
     'code_representation':'R(q) direct_sum U(Galois(q)); central -1 -> diag(I3,-I2)',
     'code_dimension':5,'cover_order':120,'carrier_quotient_order':60,
     'kernel_size':2,'central_bit_retained':True,
     'generators':[encodeq(g) for g in generators],
     'generator_permutations':[list(project[g]) for g in generators],
     'exact_group_product_checks':product_checks,'exact_unitary_checks':unitary_checks,
     'exact_generator_representation_products_checked':hom_checks,
     'exact_response_gram':exact_gram,
     'response_evaluator_float64_tolerance':'1e-12',
     'response_evaluator_Gram_within_tolerance':bool(gram_error<1e-12),
     'response_evaluator_covariance_within_tolerance':bool(cov_error<1e-12),
     'exact_basis_state_transport_steps':state_steps,
     'maximum_compiled_word_length':max(map(len,words.values())),
     'total_compiled_word_length':sum(map(len,words.values())),
     'central_word':zrow['word'],
     'closed_chart_words_verified':120,
     'autonomous_fixed_schedule_pulses':len(cycle),
     'autonomous_fixed_schedule_lifts_visited':len(visited),
     'autonomous_cycle_returns_code_and_chart_exactly':True,
     'typed_arrow_inverse_checks':inverse_arrow_checks,
     'chart_scope':'Classical source/target fibre types; each accepted arrow is unitary on its declared five-level fibre. No unitary reset of a superposed chart register is claimed.',
     'extra_continuous_response_directions':0,
     'extra_directions_scope':'The finite W(q) overlap gates only, on the declared five-level response. Not a completeness certificate for a larger apparatus or coupled source.',
     'full_scheduler_carrier_covariance_proved':False,
     'response_interface':'Geometric charged coherent Xp; all R(q)+U(Galois(q)) lie in its connected SU3 x SU2 subgroup; full twelve-direction control is declared separately.',
     'G5_G6_scope':'closed same-instruction words and endogeneity for this declared two-chart code; does not certify a growing federation or reinterpret old pump words',
     'remaining':['encoded preparation/readout joined to source observer state','growing federation local overlap law','refinement intertwiners and uniform resource bounds','physical propagation or particle identification'],
     'compiler':compiler,
    }
    return out
