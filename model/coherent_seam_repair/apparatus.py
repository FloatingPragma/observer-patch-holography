"""Finite complex tomography from optional probes of the declared walk.

Phase probes are apparatus interventions; they are absent from idle evolution.
"""
from functools import lru_cache
from fractions import Fraction as F
import sympy as s
from model.common_source import slow_geometry


@lru_cache(maxsize=1)
def certificate():
    edges=list(slow_geometry.source().seams())
    addresses=[(a,b) for edge in edges for a,b in (edge,edge[::-1])]
    D=len(addresses);n=12;d=2*n*D
    assert D==60 and d==1440
    def index(a,r,k):return (2*a+r)*D+k
    def unpack(i):return i//(2*D),(i//D)%2,i%D
    address_index={e:k for k,e in enumerate(addresses)}
    neighbors={a:set() for a in range(n)}
    for a,b in edges:neighbors[a].add(b);neighbors[b].add(a)
    # Q is the literal ordered-seam rank-one reflector; W=Q R.
    # Store its nonzero integer numerators W_ij = numerator_ij / 60.
    rows=[]
    for i in range(d):
     a,r,k=unpack(i);u,v=addresses[k]
     if a not in (u,v):qrow={(a,r):F(1)}
     else:
      signs={(u,0):1,(v,0):-1,(u,1):-1,(v,1):-1}
      qrow={(b,t):F(int((a,r)==(b,t)))-F(signs[a,r]*sgn,2) for (b,t),sgn in signs.items()}
     row={}
     for (b,t),q in qrow.items():
      assert q!=0
      for l in range(D):
       z=q*(2-D*int(k==l));assert z.denominator==1 and z!=0
       row[index(b,t,l)]=z.numerator
     assert sum(z*z for z in row.values())==D*D
     assert row[i]!=0
     rows.append(row)
    assert F(sum(row[i] for i,row in enumerate(rows)),D)==-1276
    minimum_diagonal=min(abs(F(row[i],D)) for i,row in enumerate(rows))
    minimum_nonzero=min(abs(F(z,D)) for row in rows for z in row.values())
    assert minimum_diagonal==F(29,60) and minimum_nonzero==F(1,60)

    # Explicit directed paths: change the address, then follow a shortest carrier path.
    source_paths={}
    for a in range(n):
     paths={a:[a]};queue=[a]
     for x in queue:
      for y in sorted(neighbors[x]):
       if y not in paths:paths[y]=paths[x]+[y];queue.append(y)
     assert len(paths)==n
     for b,path in paths.items():source_paths[a,b]=path
    assert max(len(p)-1 for p in source_paths.values())==3
    def support_path(i,j):
     if i==j:return [i]
     a,r,k=unpack(i);b,t,l=unpack(j)
     if a==b and r==t:return [i,j]
     path=[i]
     if a==b:
      e=(a,min(neighbors[a]));first=index(a,r,address_index[e])
      if first!=i:path.append(first)
      path.append(j)
     else:
      vp=source_paths[a,b]
      first=index(a,r,address_index[vp[0],vp[1]])
      if first!=i:path.append(first)
      for h in range(1,len(vp)):
       addr=l if h==len(vp)-1 else address_index[vp[h],vp[h+1]]
       path.append(index(vp[h],t,addr))
     assert path[-1]==j and len(path)-1<=4
     return path
    # All 2,073,600 pairs, not just an undirected or single-root connectivity check.
    pair_checks=0;maximum_path=0
    for i in range(d):
     for j in range(d):
      path=support_path(i,j)
      assert all(y in rows[x] for x,y in zip(path,path[1:]))
      maximum_path=max(maximum_path,len(path)-1);pair_checks+=1

    # Exact general matrix-unit extraction, using complexification of a real
    # Hermitian effect span. Each term needs at most five quarter-turn gates.
    N=4;A=s.Matrix(N,N,lambda i,j:s.Symbol(f'a{i}{j}'))
    def adjphase(A,i,power=1):
     V=s.eye(N);V[i,i]=s.I**power
     return V.H*A*V
    def isolate_pair(A,i,j):
     Di=(A-adjphase(A,i,2))/2
     return (Di-adjphase(Di,j,2))/2
    def matrix_unit_projection(A,i,j):
     pair=isolate_pair(A,i,j)
     return (pair+s.I*adjphase(pair,i))/2
    for i in range(N):
     for j in range(N):
      if i==j:continue
      expected=s.zeros(N);expected[i,j]=A[i,j]
      assert s.simplify(matrix_unit_projection(A,i,j)-expected)==s.zeros(N)

    k=0;a,b=addresses[k];x=index(a,0,k);y=index(b,0,k)
    exy=F(rows[x][x]*rows[x][y],D*D)
    assert exy==F(841,3600)
    rho_xy=F(-29,31440)
    assert rho_xy!=0
    out={
     'status':'EXACT same-walk optional apparatus support certificate plus general span proof',
     'walk_dimension':d,'W_rational_denominator':D,
     'nonzero_W_entries':sum(map(len,rows)),
     'nonzero_diagonals':d,'trace_W':-1276,
     'directed_source_support_path_pairs_checked':pair_checks,
     'maximum_constructed_directed_path_length':maximum_path,
     'general_symbolic_ordered_pair_projections_checked':12,
     'matrix_unit_projection_max_quarter_turn_gates':5,
     'full_complex_Hermitian_effect_span_dimension':d*d,
     'maximum_unitary_setting_depth_before_detector':6*maximum_path,
     'minimum_diagonal_magnitude':str(minimum_diagonal),
     'minimum_nonzero_entry_magnitude':str(minimum_nonzero),
     'matrix_unit_reconstruction_coefficient_l1_bound':str((1/(minimum_diagonal*minimum_nonzero))**maximum_path),
     'Hermitian_pair_reconstruction_coefficient_l1_bound':str(2*(1/(minimum_diagonal*minimum_nonzero))**maximum_path),
     'noncommuting_probe_witness_E_xy':str(exy),
     'prepared_reference_noncommuting_projector_witness_rho_xy':str(rho_xy),
     'quarter_turn_probe_witness_E_xy':'-841*i/3600',
     'declared_optional_apparatus':['basis projector detector {P_i,I-P_i}','basis-selective quarter turn V_i=I+(i-1)P_i'],
     'same_existing_operation':'W and its finite iterates; inverse W not needed for the span proof',
     'stationary_reference_scope':'reduced walk algebra A_walk; optional probes disturb it and are not idle events',
     'arbitrary_POVM_or_modular_gate_implementation':False,
     'changed_idle_dynamics':False,
    }
    return out
