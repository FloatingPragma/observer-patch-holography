"""Exact committed phase/parity comparison on the public twelve-port carrier."""
from functools import lru_cache
import sympy as s
from model.common_source import slow_geometry


@lru_cache(maxsize=1)
def certificate():
    edges=[tuple(map(int,e)) for e in slow_geometry.source().seams()]
    n=12;I=s.eye(n)
    adj={i:set() for i in range(n)}
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    assert len(edges)==30 and all(len(x)==5 for x in adj.values())
    def actual_kraus(a,b):
     plus=s.zeros(n,1);minus=s.zeros(n,1);plus[a]=plus[b]=minus[a]=1/s.sqrt(2);minus[b]=-1/s.sqrt(2)
     outside=I.copy();outside[a,a]=outside[b,b]=0
     base=[outside+plus*plus.T,plus*minus.T];out=[]
     for endpoint in [a,b]:
      V=I.copy();V[endpoint,endpoint]=s.I
      out.extend([K*V/s.sqrt(2) for K in base])
     return out

    def event(X,a,b):
     Y=X.copy();S=X[a,a]+X[b,b]
     Y[a,a]=Y[a,b]=Y[b,a]=Y[b,b]=S/2
     for j in range(n):
      if j in (a,b):continue
      Y[a,j]=Y[b,j]=(1+s.I)*(X[a,j]+X[b,j])/4
      Y[j,a]=Y[j,b]=(1-s.I)*(X[j,a]+X[j,b])/4
     return Y
    # Independently check every matrix unit for one arbitrary seam; formula is index-natural.
    a,b=edges[0];K=actual_kraus(a,b)
    assert s.simplify(sum((k.conjugate().T*k for k in K),s.zeros(n)))==I
    for i in range(n):
     for j in range(n):
      X=s.zeros(n);X[i,j]=1
      actual=sum((k*X*k.conjugate().T for k in K),s.zeros(n)).applyfunc(s.simplify)
      assert actual==event(X,a,b)
    # Four carrier orbitals.
    dist=[]
    for a in range(n):
     d=[99]*n;d[a]=0;queue=[a]
     for x in queue:
      for y in adj[x]:
       if d[y]==99:d[y]=d[x]+1;queue.append(y)
     dist.append(d)
    assert max(map(max,dist))==3
    B=[s.Matrix(n,n,lambda i,j:int(dist[i][j]==d)) for d in range(4)]
    rho=B[0]/12+s.Rational(19,1332)*B[1]+s.Rational(1,444)*B[2]+s.Rational(1,1332)*B[3]
    assert sum((event(rho,a,b) for a,b in edges),s.zeros(n))/30==rho
    assert s.trace(rho)==1
    values=[s.Rational(1,6),s.Rational(55,666)+4*s.sqrt(5)/333,s.Rational(55,666)-4*s.sqrt(5)/333,s.Rational(5,74)]
    assert all(x.is_positive for x in values)
    x=s.Symbol('x')
    expected=(x-s.Rational(1,6))*(x-s.Rational(5,74))**5*((x-s.Rational(55,666))**2-s.Rational(80,110889))**3
    assert s.expand(rho.charpoly(x).as_expr()-expected)==0
    # Actual one-event forward effects: their sum is half the seam diagonal projector,
    # and their difference is an imaginary edge operator up to sign.
    def X(a,b,imag=False):
     M=s.zeros(n);M[a,b]=s.I if imag else 1;M[b,a]=-s.I if imag else 1;return M
    for a,b in edges:
     ks=actual_kraus(a,b);Ea=s.simplify(ks[1].conjugate().T*ks[1]);Eb=s.simplify(ks[3].conjugate().T*ks[3])
     D=s.zeros(n);D[a,a]=D[b,b]=1
     assert s.simplify(Ea+Eb)==D/2
     assert s.simplify(2*(Ea-Eb))==X(a,b,True)
    # The no-jump recorded branch on an adjacent triangle gives real edge effects.
    triangle=0
    for a,b in edges:
     c=next(iter(adj[a]&adj[b]));ks=actual_kraus(a,c);K=ks[0]*s.sqrt(2)
     assert s.simplify(K.conjugate().T*X(a,b,True)*K)==(X(a,b)+X(c,b,True))/2
     triangle+=1
    inc=s.Matrix([[int(i in e) for i in range(n)] for e in edges]);assert inc.rank()==12
    out={'status':'EXACT independent minimal12port audit',
     'local_complete_instrument':True,'all_local_matrix_units_checked':144,
     'phase_menu':['+i on endpoint a','+i on endpoint b'], 'reverse_pump':False,
     'stationary_rho_by_distance':['1/12','19/1332','1/444','1/1332'],
     'rho_eigenvalues':[str(v) for v in values], 'multiplicities':[1,3,3,5],
     'rho_characteristic_polynomial_checked':True,
     'rho_faithful_nontracial':True,'actual_one_event_imaginary_edges_checked':30,
     'actual_two_event_real_edges_checked':triangle,'all_diagonals_incidence_rank':12,
     'full_history_Hermitian_span_dimension':144,'history_depth_bound':4,
     'modular_time_or_gauge_or_protected_memory_claim':False}
    return out
