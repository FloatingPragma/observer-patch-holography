"""Exact native CAR columns and particle-number-resolved averaging bounds."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb,isqrt
import numpy as np
from model.openwave_current_averaging.averaging import words,source,act,reduce_polynomial,parent,receiver
from model.native_flavor_orbit.intervals import I,sqrt_rational

GRID=2**36

def upper(q):return F(-(-F(q).numerator*GRID//F(q).denominator),GRID)
def root(q):return sqrt_rational(q,36).hi

def multiply_poly(poly,scale):return {k:v*scale for k,v in poly.items()}

@lru_cache(None)
def coefficient_variables():
 d=source();q=(1-d['energy_boxes'][0])/2;out=[]
 for E in d['energy_boxes'][3:]:
  for sign in(1,-1):
   z=(E+sign*q)/(2*E);out.append(I(sqrt_rational(z.lo).lo,sqrt_rational(z.hi).hi))
 return out

@lru_cache(None)
def coefficient_interval(key):
 variables=coefficient_variables();sqrt2=sqrt_rational(2);result=I(0)
 for powers,r,c in key:
  term=I(c)*(sqrt2 if r else I(1))
  for x,p in zip(variables,powers):term*=x**p
  result+=term
 return result

def interval(poly):return coefficient_interval(tuple(sorted((p,r,c)for(p,r),c in poly.items())))
def bound(poly):
 a=interval(poly);return upper(max(abs(a.lo),abs(a.hi)))

@lru_cache(None)
def column(mask):
 out=defaultdict(lambda:defaultdict(F))
 for ops,(powers,r,c),_ in words()[0]:
  k,sign=act(mask,ops)
  if sign:out[k][powers,r]+=sign*c
 return {k:q for k,v in out.items()if(q:=reduce_polynomial(v))}

@lru_cache(None)
def normal_order(ops):
 # Original left-to-right CAR word. Swap every adjacent a_i a_j^*;
 # the contraction delta_ij retains all reference-induced lower degrees.
 for j in range(len(ops)-1):
  if not ops[j][1] and ops[j+1][1]:
   out=defaultdict(int)
   swapped=ops[:j]+(ops[j+1],ops[j])+ops[j+2:]
   for k,v in normal_order(swapped).items():out[k]-=v
   if ops[j][0]==ops[j+1][0]:
    for k,v in normal_order(ops[:j]+ops[j+2:]).items():out[k]+=v
   return {k:v for k,v in out.items()if v}
 creators=[i for i,c in ops if c];ann=[i for i,c in ops if not c]
 if len(set(creators))!=len(creators)or len(set(ann))!=len(ann):return {}
 inversions=sum(a>b for xs in(creators,ann)for j,a in enumerate(xs)for b in xs[j+1:])
 return {(tuple(sorted(creators)),tuple(sorted(ann))):(-1)**inversions}

@lru_cache(None)
def wick():
 out=defaultdict(lambda:defaultdict(F))
 for ops,(powers,r,c),_ in words()[0]:
  for pair,sign in normal_order(ops).items():out[pair][powers,r]+=c*sign
 return {k:q for k,v in out.items()if(q:=reduce_polynomial(v))}

def matrix_norms():
 groups=defaultdict(list)
 for(creates,ann),poly in wick().items():groups[len(creates),len(ann)].append((creates,ann,bound(poly)))
 result={};rows=[]
 for(a,b),entries in sorted(groups.items()):
  if a==b==0:continue
  r=defaultdict(F);c=defaultdict(F);squares=F(0)
  for i,j,z in entries:r[i]+=z;c[j]+=z;squares+=z*z
  schur=max(r.values())*max(c.values());norm=upper(root(min(schur,squares)));result[a,b]=norm
  rows.append({'creators':a,'annihilators':b,'nonzero_coefficients':len(entries),'coefficient_operator_norm_upper':str(norm),'coefficient_Frobenius_squared_upper':str(squares),'coefficient_Schur_squared_upper':str(schur)})
 scalar=wick().get(((),()),{});scalarbox=interval(scalar);return result,rows,scalarbox

def energy_vector(mask):
 d=source();return tuple(sum(v[k]for i,v in enumerate(d['energy_vectors'])if mask>>i&1)for k in range(6))
def gap(mask):
 v=list(energy_vector(mask));v[1]-=2
 z=sum((a*b/2 for a,b in zip(v,source()['energy_boxes'])),I(0));assert z.lo>0 or z.hi<0
 return min(abs(z.lo),abs(z.hi))

def bounds():
 P=(parent(),receiver());squares=defaultdict(F);details=[];min_gap=None
 for mask in P:
  col=column(mask);nons=0
  for target,poly in col.items():
   if target in P:continue
   q=gap(target);min_gap=q if min_gap is None else min(min_gap,q);v=bound(poly);squares[target.bit_count()]+=4*v*v/(q*q);nons+=1
  details.append({'prepared_mask':str(mask),'exact_nonzero_current_column_entries':len(col),'nonresonant_entries':nons})
 A={n:upper(root(s))for n,s in squares.items()};total=upper(root(sum((v*v for v in A.values()),F(0))))
 norms,normrows,scalar=matrix_norms();centre=F(343,4);constant=upper(max(abs(scalar.lo-centre),abs(scalar.hi-centre)))
 outgoing=defaultdict(F)
 for n,z in A.items():
  outgoing[n]+=constant*z
  for(a,b),norm in norms.items():
   if n<b:continue
   k=n-b+a;factor=upper(root(F(comb(n,b)*comb(k,a))))
   outgoing[k]+=norm*factor*z
 product_norm=upper(root(sum((v*v for v in outgoing.values()),F(0))))
 tau=F(44,35);coefficient=total+tau*(product_norm+F(5,4)*total)
 print('COLUMN_COUNTS',details,'GAP',float(min_gap),'A', {n:float(v)for n,v in A.items()},'NORMS',[(r['creators'],r['annihilators'],float(F(r['coefficient_operator_norm_upper'])))for r in normrows],'SCALAR',scalar.strings(),'K',float(product_norm),'ERROR/G',float(coefficient),flush=True)
 return dict(columns=details,minimum_gap=min_gap,primitive_by_number=A,primitive=total,norm_rows=normrows,scalar=scalar,centred_constant=constant,product_by_number=outgoing,product=product_norm,error_coefficient=coefficient)

if __name__=='__main__':bounds()

def certificate():
 b=bounds();scenarios=[]
 for exponent in(72,44,30,6):
  G=F(1,2**exponent);error=G*b['error_coefficient'];scenarios.append({'G':str(G),'binary_exponent':exponent,'coupling_multiple_of_2_minus72':str(F(2**72,2**exponent)),'prepared_domain_uniform_state_error_upper':str(error),'receiver_first_peak_probability_interval':[str(max(F(0),F(16,25)-2*error)),str(min(F(1),F(16,25)+2*error))],'receiver_first_peak_probability_decimal_interval':[float(max(F(0),F(16,25)-2*error)),float(min(F(1),F(16,25)+2*error))],'entire_orthogonal_source_remainder_probability_upper':str(min(F(1),error*error)),'nontrivial_comparison':error<F(1,10),'first_peak_source_time_interval':['2*pi/(5G)','2*pi/(5G)'],'source_energy_measurement_work_absolute_upper':str(2*G*41952)})
 assert F(scenarios[1]['prepared_domain_uniform_state_error_upper'])<F(1,5000000)
 assert F(scenarios[2]['prepared_domain_uniform_state_error_upper'])<F(1,318)
 assert scenarios[-1]['nontrivial_comparison']is False
 from model.openwave_current_averaging.averaging import exact_orbit,small_block_records
 return {'schema':'oph.prepared_native_current_averaging.v1','verdict':'PASS_NUMBER_RESOLVED_FULL_SOURCE_AVERAGING','full_native_CAR_modes':270,'full_native_current_words':len(__import__('model.openwave_weak_decay.decay',fromlist=['interaction']).interaction()[0]),
 'source_unchanged_except_supplied_G':True,'source_truncation_or_added_penalty':False,
 'prepared_quantum_domain':'span(parent,receiver) tensor any external reference; all full-source output complements retained',
 'exact_current_columns':b['columns'],'actual_column_nonresonant_gap_lower':str(b['minimum_gap']),'actual_column_gap_lower_diagnostic':float(b['minimum_gap']),
 'complete_Wick_scalar':b['scalar'].strings(),'centering_constant':'343/4','exact_resonant_block_centred_norm':'5/4','complete_Wick_coefficient_matrix_norms':b['norm_rows'],
 'oscillatory_primitive_Hilbert_Schmidt_bound_by_particle_number':{str(k):str(v)for k,v in sorted(b['primitive_by_number'].items())},
 'complete_primitive_operator_bound':str(b['primitive']),
 'current_times_primitive_bound_by_output_particle_number':{str(k):str(v)for k,v in sorted(b['product_by_number'].items())},
 'complete_current_times_primitive_operator_bound':str(b['product']),
 'prepared_domain_error_rule':'G*A+G^2*t*(K+5*A/4); 0<=Gt<=2pi/5',
 'uniform_first_peak_error_coefficient':str(b['error_coefficient']),
 'coupling_scenarios':scenarios,'exact_resonant_orbit':exact_orbit(),'executed_resonant_complete_records':small_block_records(),
 'all_reference_states_all_particle_numbers_and_failures_retained':True,
 'global_all_input_unitary_bound_claimed':False,
 'physical_weak_coupling_lifetime_or_irreversibility_claimed':False,
 'scope':'Exact current normal form and rational prepared-domain full-source bounds. Finite number sectors bound two applications of the current; they do not replace the source dynamics or truncate its evolving state. Architecture and calibration are permitted inputs. No physical weak-rate fit is asserted.'}
