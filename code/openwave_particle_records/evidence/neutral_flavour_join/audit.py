#!/usr/bin/env python3
"""Post-run exact arithmetic audit of stored matrices; no source execution."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import gzip,hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'reverse-engineering-reality').is_dir())

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def q(x):return Q(float(x))
def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul_scalar(a,c):return (a[0]*c,a[1]*c) if c>=0 else(a[1]*c,a[0]*c)
def bound(a):return max(abs(a[0]),abs(a[1]))
def z5(text):
 v=s.expand(s.sympify(text));b=v.coeff(s.sqrt(5));a=s.expand(v-b*s.sqrt(5))
 assert a.is_Rational and b.is_Rational
 return Q(int(a.p),int(a.q)),Q(int(b.p),int(b.q))
def interval(pair,root):return add((pair[0],pair[0]),mul_scalar(root,pair[1]))
def real_array(x):return [[q(v) for v in row] for row in x]
def raw_diff(x,target):
 if isinstance(x,dict):real=x['real'];imag=x['imag']
 else:real=x;imag=[[0]*len(real[0]) for _ in real]
 return sum(abs(q(real[i][j])-target[i][j])+abs(q(imag[i][j])) for i in range(len(real)) for j in range(len(real[0])))
def main():
 for line in (HERE/'AUDIT_PREREG.txt').read_text().splitlines():
  h,p=line.split('  ',1);assert sha(HERE/p)==h
 d=json.loads(gzip.decompress((HERE/'RECORDS.json.gz').read_bytes()));fr=d['frame']
 old=json.loads((ROOT/'openwave/row_10_neutrino_oscillations/work/2026-11-14_codex_neutral_even_records/SOURCE_RECEIPT.json').read_text())
 previous=next(x for x in old['policies'] if x['policy']=='charged_oriented')
 assert fr['Y_exact']==previous['record_kernel']
 D=10**60;r=isqrt(5*D*D);root=Q(r,D),Q(r+1,D)
 assert root[0]**2<5<root[1]**2
 brackets=[tuple(map(Q,row)) for row in previous['spectrum']['eigenvalue_intervals']]
 X=s.Symbol('x');poly=s.Poly(s.sympify(previous['spectrum']['characteristic_polynomial']),X)
 coeff=[z5(str(v)) for v in poly.all_coeffs()]
 signs=[]
 for i,(lo,hi) in enumerate(brackets):
  sig=[]
  for value in(lo,hi):
   total=(Q(0),Q(0))
   for c in coeff:total=add(mul_scalar(total,value),interval(c,root))
   sig.append(1 if total[0]>0 else -1 if total[1]<0 else 0)
  assert sig==([-1,1] if i in(0,2) else[1,-1]);signs.append(sig)
 assert all(brackets[i][1]<brackets[i+1][0] for i in range(2)) and brackets[0][0]>0
 Y=[[interval(z5(v),root) for v in row] for row in fr['Y_exact']]
 F=real_array(fr['frame']);eigen=fr['charged_eigenvalues'];columns=[]
 for j in range(3):
  v=[F[i][j] for i in range(3)];lam=q(eigen[j]);n2=sum(x*x for x in v)
  assert abs(n2-1)<Q(1,2)
  residual=[]
  for i in range(3):
   total=(-lam*v[i],-lam*v[i])
   for k in range(3):total=add(total,mul_scalar(Y[i][k],v[k]))
   residual.append(bound(total))
  residual_bound=sum(residual)
  gap=min(min(abs(lam-a),abs(lam-b)) for k,(a,b) in enumerate(brackets) if k!=j)
  assert gap>Q(1,100)
  # ||v||>=1/2. Projection onto other exact eigenlines is <=r/g.
  # After normalizing, phase-aligned distance is <=|norm^2-1|+4r/g.
  delta=abs(n2-1)+4*residual_bound/gap
  assert delta<Q(1,10**10)
  columns.append({'residual_l1_upper':str(residual_bound),'separation_lower':str(gap),'eigenvector_distance_upper':str(delta)})
 eigenframe_error=2*sum(Q(x['eigenvector_distance_upper']) for x in columns)
 forward=raw_diff(fr['forward'],F)+eigenframe_error
 transpose=[list(row) for row in zip(*F,strict=True)]
 inverse=raw_diff(fr['backward'],transpose)+eigenframe_error
 ef=Q(2,10**10);ex=Q(1,10**10);ep=Q(1,10**11)
 assert forward<ef and inverse<ef
 basis=[]
 for a,m in enumerate(fr['basis_full_matrices']):
  phase=(1j)**a.bit_count();error=Q(0)
  for i in range(8):
   for j in range(8):
    v=phase if i==j^a else 0
    error+=abs(q(m['real'][i][j])-q(v.real))+abs(q(m['imag'][i][j])-q(v.imag))
  assert error<ex;basis.append(str(error))
 P=[[Q(0) for _ in range(8)] for _ in range(8)]
 for source in range(8):
  a,b,c=source&1,(source>>1)&1,source>>2;a^=c;b^=a;P[a+2*b+4*c][source]=1
 pointer=raw_diff(fr['pointer_matrix'],P);assert pointer<ep
 timeout=Q(d['reset']['timeout_bound']);dp=ep*(2+ep)/2;df=ef*(2+ef)/2;dx=ex*(2+ex)/2
 budgets=[]
 for k in range(7):
  reads=6+3*k;frames=2+2*k;xs=6
  gain=(1+ef)**(2*frames)*(1+ex)**(2*xs)*(1+ep)**(2*reads)
  total=gain*(reads*(timeout+dp)+xs*dx+frames*df)
  assert total<Q(1,10**6)
  budgets.append({'intermediate_full_family_reads':k,'native_bit_read_services':reads,'frame_words':frames,'X_words_max':xs,'complete_halfdiamond_upper':str(total),'decimal_diagnostic':float(total)})
 reset_p=Q(d['complete_apparatus']['exact_reset_success']);cycle_p=Q(255,256)**12
 out={'schema':'oph.neutral_three_tag_artifact_audit.v1','audit_prereg_sha256':sha(HERE/'AUDIT_PREREG.txt'),
   'records_sha256':sha(HERE/'RECORDS.json.gz'),'exact_field_root_brackets_checked':signs,'eigenvector_certificates':columns,
   'full8_forward_to_exact_source_frame_upper':str(forward),'full8_inverse_to_exact_source_frame_upper':str(inverse),
   'basis_source_word_errors_upper':basis,'pointer_source_word_error_upper':str(pointer),'complete_apparatus_budgets':budgets,
   'all_six_auxiliary_resets_success_diagnostic':float(reset_p**6),'ideal_unconditional_cycle3_correct_return_diagnostic':float(reset_p**6*cycle_p),
   'matrix_entries_treated_as_exact_binary_rationals':True,'new_native_execution':False,'hardware_accuracy_authenticated':False,'all_checks_pass':True}
 (HERE/'APPARATUS_AUDIT.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'all_checks_pass':True,'full8_source_frame_bound':float(max(forward,inverse)),'pointer_bound':float(pointer),'terminal_apparatus_bound':budgets[0]['decimal_diagnostic'],'six_intermediate_reads_bound':budgets[-1]['decimal_diagnostic'],'unconditional_cycle3_return':float(reset_p**6*cycle_p)},indent=2))
if __name__=='__main__':main()
