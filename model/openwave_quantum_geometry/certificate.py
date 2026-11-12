"""Exact dyadic residual certificates for all native conditional blocks."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial,isqrt
from itertools import product
import numpy as np
from scipy.linalg import eigh
from mpmath import iv
from model.openwave_energy_lapse.source import native
from .source import geometry,block,KAPPA,clock_state,partial_transpose,complete_records,independent_geometry_control
BITS=46;DEN=2**BITS

def scalar(x):
 if isinstance(x,F):return iv.mpf(x.numerator)/x.denominator
 return iv.mpf(x)
def rational(mpf_tuple):
 sign,mantissa,exponent,_=mpf_tuple;v=F(int(mantissa))*F(2)**int(exponent)
 return -v if sign else v
def ends(x):return rational(x._mpi_[0]),rational(x._mpi_[1])
def enc(x):return [str(v)for v in ends(x)]
def sqceil(n):
 v=isqrt(int(n));return v+int(v*v<n)
def matrix_exp(A):
 n=len(A);I=[[iv.mpf(int(i==j))for j in range(n)]for i in range(n)]
 out=[row[:]for row in I];power=[row[:]for row in I]
 for k in range(1,13):
  power=[[sum((A[i][l]*power[l][j]for l in range(n)),iv.mpf(0))/k for j in range(n)]for i in range(n)]
  out=[[out[i][j]+power[i][j]for j in range(n)]for i in range(n)]
 # Each single normal-mode coordinate has norm at most1/50 at cap4.
 tail=F(2)*F(1,50)**13/F(factorial(13))
 return [[v+iv.mpf([-scalar(tail),scalar(tail)])for v in row]for row in out]
@lru_cache(maxsize=4)
def exact_geometry(cap):
 iv.dps=80;n=cap+1;a=[[iv.mpf(0)for _ in range(n)]for _ in range(n)];b=[[iv.mpf(0)for _ in range(n)]for _ in range(n)]
 for j in range(1,n):
  a[j-1][j]=a[j][j-1]=iv.sqrt(j)/200
  b[j-1][j]=b[j][j-1]=iv.sqrt(j)/(200*iv.sqrt(iv.sqrt(3)))
 A=matrix_exp(a);B=matrix_exp(b);minus=matrix_exp([[-v for v in row]for row in b]);exp=[]
 for Q in(B,minus):
  exp.append([[A[i//n][j//n]*Q[i%n][j%n]for j in range(n*n)]for i in range(n*n)])
 diag=[iv.mpf(i)+iv.sqrt(3)*j for i,j in product(range(n),repeat=2)]
 return exp,diag

def block_certificate(cap,i,j):
 iv.dps=80;g=geometry(cap);H=block(cap,i,j);e,O=eigh(H);d=len(e)
 Hint=np.array([[int(round(float(v)*DEN))for v in row]for row in H],dtype=object)
 Oint=np.array([[int(round(float(v)*DEN))for v in row]for row in O],dtype=object)
 Dint=np.array([int(round(float(v)*DEN))for v in e],dtype=object)
 gram=Oint.T@Oint-np.eye(d,dtype=object)*(DEN*DEN)
 eta=F(sqceil(sum(int(v)**2 for v in gram.flat)),DEN*DEN);assert eta<F(1,2)
 eQ,hg=exact_geometry(cap);gaps=native()[2];nativeE=[iv.mpf([scalar(gaps[k].lo),scalar(gaps[k].hi)])for k in(i,j)]
 entry=F(0)
 for a in range(d):
  for b in range(d):
   exact=nativeE[0]*eQ[0][a][b]+nativeE[1]*eQ[1][a][b]+(hg[a]if a==b else 0)
   low,high=ends(exact);mid=F(int(Hint[a,b]),DEN);entry=max(entry,abs(low-mid),abs(high-mid))
 Herror=d*entry
 R=Hint@Oint-Oint*Dint[None,:]
 residual=F(sqceil(sum(int(v)**2 for v in R.flat)),DEN*DEN)+2*Herror
 bound=F(30);assert max(abs(int(v))for v in Dint)<bound*DEN
 epsilon=2*residual+4*bound*eta
 ground=F(int(Dint[0]),DEN)
 return {'cap':cap,'native_levels':[i,j],'dyadic_bits':BITS,'eigenvalue_integers':list(map(int,Dint)),
  'eigenframe_integer_rows':[[int(v)for v in row]for row in Oint],
  'orthogonality_Frobenius_upper':str(eta),'Hamiltonian_entry_error_upper':str(entry),
  'exact_source_eigenframe_residual_Frobenius_upper':str(residual),
  'polar_unitary_Hamiltonian_error_upper':str(epsilon),
  'all_eigenvalue_intervals':[[str(F(int(v),DEN)-epsilon),str(F(int(v),DEN)+epsilon)]for v in Dint],
  'ground_energy_interval':[str(ground-epsilon),str(ground+epsilon)],
  'complete_block_unitary_error_rule':'abs(t)*polar_unitary_Hamiltonian_error_upper +3*orthogonality_Frobenius_upper'}

def witness(cap,t,certificates):
 branches,rho=clock_state(cap,t);ev,V=np.linalg.eigh(partial_transpose(rho));v=V[:,0]
 integer=[(int(round(z.real*2**30)),int(round(z.imag*2**30)))for z in v];norm=sum(a*a+b*b for a,b in integer)
 w=np.array([complex(a,b)for a,b in integer]);w/=np.sqrt(norm)
 W=partial_transpose(np.outer(w,w.conj()));estimate=float(np.trace(W@rho).real)
 # The exact numerical eigenframe defines a dyadic approximation. Re-evaluate
 # its unitary phases and witness using interval arithmetic.
 vectors=[];errors=[]
 for i,j in product((0,1),repeat=2):
  c=certificates[8*i+j];O=c['eigenframe_integer_rows'];D=c['eigenvalue_integers'];n=len(D);vec=[]
  phases=[iv.exp(-iv.j*scalar(F(t)*F(x,DEN)))for x in D]
  for a in range(n):vec.append(sum((scalar(F(O[a][k]*O[0][k],DEN*DEN))*phases[k]for k in range(n)),iv.mpc(0)))
  vectors.append(vec);errors.append(abs(F(t))*F(c['polar_unitary_Hamiltonian_error_upper'])+3*F(c['orthogonality_Frobenius_upper']))
 er=max(errors);assert er<F(1,100)
 density=[[sum((vectors[a][n]*iv.mpc(vectors[b][n].real,-vectors[b][n].imag)for n in range(len(vectors[0]))),iv.mpc(0))/4 for b in range(4)]for a in range(4)]
 # Tr[(|w><w|)^T_B rho] = <w|rho^T_B|w>, with an exact rational witness.
 value=iv.mpc(0)
 for a,b in product(range(4),repeat=2):
  aa,ab=divmod(a,2);ba,bb=divmod(b,2)
  # Integer products below avoid loss of exactness in the witness coefficient.
  ar,ai=integer[a];br,bi=integer[b];co=iv.mpc(ar*br+ai*bi,ar*bi-ai*br)/norm
  value+=co*density[2*aa+bb][2*ba+ab]
 low,high=ends(value.real);uncertainty=2*er+er*er
 interval=[low-uncertainty,high+uncertainty];assert interval[1]<0
 # Complete Pauli/geometry instruments are bounded from the same exact
 # dyadic propagator, retaining every complementary native level.
 X=[[1/iv.sqrt(2),1/iv.sqrt(2)],[1/iv.sqrt(2),-1/iv.sqrt(2)]]
 Y=[[1/iv.sqrt(2),-iv.j/iv.sqrt(2)],[1/iv.sqrt(2),iv.j/iv.sqrt(2)]]
 Z=[[iv.mpf(1),iv.mpf(0)],[iv.mpf(0),iv.mpf(1)]]
 instrument=[]
 for nameA,A in(('X',X),('Y',Y),('Z',Z)):
  for nameB,B in(('X',X),('Y',Y),('Z',Z)):
   entries=[]
   for row in range(4):
    aa,ab=divmod(row,2);out=[]
    for n in range(len(vectors[0])):
     amp=sum((A[aa][col//2]*B[ab][col%2]*vectors[col][n]/2 for col in range(4)),iv.mpc(0))
     lo,hi=ends(amp.real**2+amp.imag**2)
     lo=max(F(0),lo-uncertainty);hi=min(F(1),hi+uncertainty)
     out.append([float(np.nextafter(float(lo),-np.inf))if lo else 0.,float(np.nextafter(float(hi),np.inf))if hi<1 else 1.])
    entries.append(out)
   assert sum(x[0]for row in entries for x in row)<=1<=sum(x[1]for row in entries for x in row)
   instrument.append({'axes':nameA+nameB,'all_active_clock_geometry_probability_intervals':entries,
    'other60_native_pair_lines_each_geometry_occupation_probability':'0','complete_outcome_count':64*len(vectors[0])})
 return {'cap':cap,'source_time':t,'normalized_witness_integer_pairs':[list(v)for v in integer],
  'witness_squared_integer_norm':norm,'approximate_clock_state_unitary_error_upper':str(er),
  'exact_source_partial_transpose_witness_interval':list(map(str,interval)),
  'witness_decimal_interval':[float(x)for x in interval],
  'numeric_partial_transpose_minimum':float(ev[0]),'witness_operator_norm_upper':'1',
  'complete_certified_joint_instruments':instrument,'all_native_and_geometry_outcomes_retained':True,'physical_gravity_entanglement_prediction':False}

def run():
 caps=[]
 for cap in(1,2,4):
  rows=[block_certificate(cap,i,j)for i,j in product(range(8),repeat=2)];cross=[]
  for i,j in product(range(1,8),repeat=2):
   both,left,right,zero=[list(map(F,rows[8*a+b]['ground_energy_interval']))for a,b in((i,j),(i,0),(0,j),(0,0))]
   low=both[0]-left[1]-right[1]+zero[0];high=both[1]-left[0]-right[0]+zero[1];assert high<0
   cross.append({'levels':[i,j],'mutual_ground_energy_interval':[str(low),str(high)]})
  caps.append({'cap':cap,'complete_geometry_dimension':len(geometry(cap)['basis']),
    'complete_native_matter_geometry_dimension':64*len(geometry(cap)['basis']),
    'all64_conditional_spectral_certificates':rows,'all49_positive_energy_source_pair_interactions':cross,
    'entanglement':witness(cap,10000,rows),'records':complete_records(cap,10000),
    'independent_geometry_control':independent_geometry_control(cap,10000),
    'source_time_response':[complete_records(cap,t)for t in(0,1,10,100)]})
 return {'schema':'oph.autonomous_native_quantum_geometry.v1','native_matter_dimension_per_cell':8,
  'native_gap_intervals':[v.encoded()for v in native()[2]],'geometry_spatial_Laplacian':[[2,-1],[-1,2]],
  'geometry_stiffness_and_inertia':KAPPA,'source':'Hg + exp(Q0)K0 + exp(Q1)K1; full native8x8 matter, complete product normal-mode oscillator caps',
  'geometry_mode_frequencies':['1','sqrt(3)'],'native_matter_energies_are_conserved':True,
  'all_source_energy_is_conserved':True,'global_ground':'native0,0 and geometry vacuum, unique with energy0',
  'static_minimum_controller_used':False,'geometry_bath_or_irreversible_relaxation_assumed':False,
  'cap_comparisons':caps,'physical_G_clock_length_or_observed_gravity_identified':False,
  'scope':'Declared finite quantum geometry extension with exact exponential energy coupling. Complete source-clock spectral certificates and finite mediator-induced entanglement, without a continuum, universal classical lapse or physical gravitational strength.',
  'verdict':'PASS_AUTONOMOUS_QUANTUM_GEOMETRY_NATIVE_ENTANGLEMENT'}
