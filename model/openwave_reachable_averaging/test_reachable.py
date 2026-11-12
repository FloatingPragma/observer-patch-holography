from collections import defaultdict
from fractions import Fraction as F
from .reachable import normal_order,wick,column,bounds,certificate
from model.openwave_current_averaging.averaging import act,parent,receiver,reduce_polynomial

def test_exact_normal_order_on_all_small_fock_words():
 from itertools import product
 for raw in product(range(4),repeat=4):
  ops=tuple((i,c)for i,c in zip(raw,(True,False,True,False)))
  for mask in range(16):
   target,sign=act(mask,ops);expected={}if not sign else {target:sign};got=defaultdict(int)
   for(I,J),a in normal_order(ops).items():
    k,z=act(mask,tuple((i,True)for i in I)+tuple((j,False)for j in J))
    if z:got[k]+=a*z
   assert {k:v for k,v in got.items()if v}==expected

def test_full_Wick_columns_match_original_source_exactly():
 for state in(parent(),receiver()):
  result=defaultdict(lambda:defaultdict(F))
  for(I,J),poly in wick().items():
   k,z=act(state,tuple((i,True)for i in I)+tuple((j,False)for j in J))
   if z:
    for term,c in poly.items():result[k][term]+=z*c
  result={k:q for k,v in result.items()if(q:=reduce_polynomial(v))}
  assert result==column(state)

def test_number_resolved_full_source_bound_and_failed_strong_control():
 b=bounds();assert set(b['primitive_by_number'])=={3,5,7};assert b['minimum_gap']>F(48,10000)
 r=certificate();assert r['coupling_scenarios'][2]['nontrivial_comparison'];assert not r['coupling_scenarios'][3]['nontrivial_comparison']
 assert F(r['coupling_scenarios'][2]['prepared_domain_uniform_state_error_upper'])<F(1,318)
