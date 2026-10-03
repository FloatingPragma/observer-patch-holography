"""Independently verify the conditional positive-capacity top coordinate band.

No producer, PaperMathContext or repository RGE engine is imported. Finite
heat sums reconstruct the source packet. Radau integrates the linear z=1/y²
equation independently of the producer's quadrature. Coupled controls use
squared non-GUT-hypercharge variables and a separately transcribed beta.
Symbolic checks are algebraic checks within the accompanying mathematical
argument, not formal verification of the complete analytic theorem.
"""
from __future__ import annotations
import argparse
import ast
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
SPEC=HERE/'top_positive_capacity_band_spec.json'
OUT=REPO/'code/particles/runs/calibration/top_positive_capacity_band.json'
REVIEWED_SPEC_SHA256='fda35fc4b193454fbaf7738f747b02ba29972d13feefd3a69becd5701192977d'
PINS={
 'code/particles/calibration/top_positive_capacity_band_spec.json',
 'code/P_derivation/codata_2022_alpha_fixture.json',
 'code/particles/calibration/sm_two_loop_rge_engine.py',
 'code/P_derivation/paper_math.py',
 'code/particles/calibration/top_positive_capacity_band.py',
}
pi=math.pi;k=1/(16*pi*pi)


def require(ok,message):
 if not ok:raise ValueError(message)


def strict_load(path):
 def pairs(items):
  out={}
  for key,value in items:
   require(key not in out,'duplicate JSON key');out[key]=value
  return out
 return json.loads(Path(path).read_text(),object_pairs_hook=pairs,
  parse_constant=lambda value:(_ for _ in ()).throw(ValueError('nonfinite JSON')))


def number(value):
 require(isinstance(value,(int,float)) and not isinstance(value,bool),'invalid numeric type')
 require(math.isfinite(value),'nonfinite value')
 return float(value)


def near(value,expected,atol=2e-6,rtol=2e-11):
 require(abs(number(value)-expected)<=atol+rtol*abs(expected),'numerical mismatch')


def keys(value,expected,label):
 require(isinstance(value,dict) and set(value)==set(expected),label+' schema changed')


def reviewed_spec():
 require(hashlib.sha256(SPEC.read_bytes()).hexdigest()==REVIEWED_SPEC_SHA256,'unreviewed specification')
 return strict_load(SPEC)


def source_ancestry_check():
 """Bound the actual producer's parsed-input and import surface."""
 tree=ast.parse((HERE/'top_positive_capacity_band.py').read_text())
 imports=set()
 for node in ast.walk(tree):
  if isinstance(node,ast.Import):imports.update(a.name for a in node.names)
  elif isinstance(node,ast.ImportFrom):imports.add(node.module)
 require(imports=={'argparse','decimal','pathlib','hashlib','json','math','sys','numpy',
                  'scipy.integrate','scipy.optimize','sm_two_loop_rge_engine','paper_math'},
         'unreviewed producer import ancestry')
 readers=[]
 for node in ast.walk(tree):
  if isinstance(node,ast.Call) and isinstance(node.func,ast.Attribute) and node.func.attr=='read_text':
   readers.append(ast.unparse(node.func.value))
 require(sorted(readers)==['ALPHA','SPEC'],'unreviewed numeric input path')
 classes=[n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='GaugeOnlyContext']
 require(len(classes)==1,'gauge-only source class missing')
 cls=classes[0]
 aliases={t.id:ast.unparse(n.value) for n in cls.body if isinstance(n,ast.Assign)
          for t in n.targets if isinstance(t,ast.Name)}
 require(all(aliases.get(n)=='forbidden' for n in ['diagonal_quark_masses','charged_lepton_masses',
  'alpha_external_from_d10','structured_thomson_running','structured_thomson_running_asymptotic','solve_closure']),
  'flavor-dependent source path no longer forbidden')
 helper=next((n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='_derive_stage5_integer_vectors'),None)
 require(helper is not None and len(helper.body)==1 and isinstance(helper.body[0],ast.Return)
  and isinstance(helper.body[0].value,ast.Dict) and not helper.body[0].value.keys,
  'unused flavor initialization no longer suppressed')
 calls={ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n,ast.Call)}
 require('ctx.p_from_inverse_alpha' in calls and 'ctx.build_d10_from_p' in calls,
         'independent-P gauge path missing')


def exact_identity_checks():
 """Exact polynomial/calculus identities; no sampled proof replacement."""
 import sympy as s
 t,tU,tE,u,wU,wE=s.symbols('t tU tE u wU wE',real=True)
 tb=(wU*tU+wE*tE)/(wU+wE)
 cost=lambda x:wU*(x-tU)**2+wE*(x-tE)**2
 require(s.cancel(cost(t)-cost(tb)-(wU+wE)*(t-tb)**2)==0,'quadratic identity')
 theta=s.symbols('theta',real=True);psi=s.Function('psi')
 thetaU,thetaE,b,a=s.symbols('thetaU thetaE b a',real=True)
 forward=(theta-thetaU)*s.diff(psi(theta),theta)-psi(theta)+psi(thetaU)
 backward=(thetaU-theta)*s.diff(psi(thetaU),thetaU)-psi(thetaU)+psi(theta)
 require(s.simplify(s.diff(forward,theta)-(theta-thetaU)*s.diff(psi(theta),theta,2))==0,'forward KL identity')
 require(s.simplify(s.diff(backward,theta)-s.diff(psi(theta),theta)+s.diff(psi(thetaU),thetaU))==0,'reverse KL identity')
 require(s.expand(b*(a+b*t-(a+b*tU))-b*b*(t-tU))==0,'affine natural-parameter identity')
 a,b,c,kk=s.symbols('a b c k',positive=True)
 Cp=s.Rational(17,12)*kk*s.Rational(41,3)*a*a+s.Rational(9,4)*(-kk*s.Rational(19,3)*b*b)+8*(-14*kk*c*c)
 require(s.expand(Cp-kk*(s.Rational(697,36)*a*a-s.Rational(57,4)*b*b-112*c*c))==0,'C derivative identity')
 require(s.simplify((a+s.sqrt(3)*b)**2-(a*a+2*a*b+3*b*b)-2*(s.sqrt(3)-1)*a*b)==0,'square-root bound identity')
 require(4-s.Rational(11,12)-9*s.sqrt(3)/8>0,'sensitivity coefficient sign')
 require(s.Rational(697,36)-28<0,'barrier coefficient sign')
 # z=1/y² converts y'=k*y*(9*y²/2-C) to a linear ODE.
 y,C=s.symbols('y C',positive=True)
 require(s.expand(-2*y**-3*kk*y*(s.Rational(9,2)*y*y-C)-(2*kk*C/y**2-9*kk))==0,'linearization identity')
 return {'quadratic_completion':True,'forward_KL_derivative':True,'reverse_KL_derivative':True,
         'affine_parameter':True,'sensitivity_and_barrier_algebra':True,'linearized_yukawa':True}


def independent_packet(spec):
 fixture=strict_load(REPO/spec['alpha_fixture'])['inverse_fine_structure_constant']
 require(fixture=={'value':spec['inverse_alpha'],'standard_uncertainty':spec['inverse_alpha_standard_uncertainty']},'independent alpha changed')
 P=(1+math.sqrt(5))/2+math.sqrt(pi)/float(fixture['value']);E=float(spec['E_star_display_GeV'])
 U=math.exp(-2*pi)*P**(1/6);cell=P**(-.5)
 def heat(t,kind):
  terms=[]
  if kind==2:
   for n in range(41):
    d=n+1;cas=n*(n+2)/4;terms.append((d*math.exp(-t*cas),math.log(d)))
  else:
   for i in range(41):
    for j in range(41):
     d=(i+1)*(j+1)*(i+j+2)/2;cas=(i*i+j*j+i*j+3*i+3*j)/3
     terms.append((d*math.exp(-t*cas),math.log(d)))
  return math.fsum(w*l for w,l in terms)/math.fsum(w for w,l in terms)
 def low(au):
  v=cell*math.exp(-pi/(2*au))
  def aa(z):return [1/(1/au+b/(2*pi)*(math.log(U)-z)) for b in [33/5,1,-3]]
  def f(z):
   a1,a2,a3=aa(z);return math.log(v/2*math.sqrt(4*pi*(a2+3*a1/5)))-z
  z=brentq(f,math.log(v)-2,math.log(v)+1,xtol=1e-13)
  return v,math.exp(z),aa(z)
 def f(au):
  _,_,(_,a2,a3)=low(au);return heat(4*pi*pi*a2,2)+heat(4*pi*pi*a3,3)-P/4
 au=brentq(f,.04,.042,xtol=5e-16)
 v,z,(a1,a2,a3)=low(au)
 return dict(P=P,E_star=E,alpha_U=au,mu_U=U*E,E_cell=cell*E,v=v*E,mz=z*E,
  low_g_squared=[4*pi*3*a1/5,4*pi*a2,4*pi*a3],alpha_inverse=fixture['value'])


def beta(t,x):
 """Buttazzo App. B in x=(gY²,g2²,g3²,yt²,lambda), t=log(mu)."""
 a,b,c,h,l=x
 da=2*a*a*(k*41/6+k*k*(199*a/18+9*b/2+44*c/3-17*h/6))
 db=2*b*b*(-k*19/6+k*k*(3*a/2+35*b/6+12*c-3*h/2))
 dc=2*c*c*(-7*k+k*k*(11*a/6+9*b/2-26*c-2*h))
 dh=2*h*(k*(9*h/2-17*a/12-9*b/4-8*c)+k*k*(
  -12*h*h-12*l*h+6*l*l+h*(131*a/16+225*b/16+36*c)
  +1187*a*a/216-3*a*b/4+19*a*c/9-23*b*b/4+9*b*c-108*c*c))
 dl=k*(24*l*l-6*h*h+3*(2*b*b+(b+a)**2)/8+l*(12*h-9*b-3*a))
 dl+=k*k*(-312*l**3-144*l*l*h+l*l*(108*b+36*a)
  +l*h*(-3*h+80*c+45*b/2+85*a/6)+l*(-73*b*b/8+629*a*a/24+39*a*b/4)
  +h*h*(30*h-32*c-8*a/3)+h*(-9*b*b/4-19*a*a/4+21*a*b/2)
  +305*b**3/16-379*a**3/48-289*b*b*a/48-559*b*a*a/48)
 return np.array([da,db,dc,dh,dl])


@lru_cache(maxsize=1)
def independent_expectation():
 spec=reviewed_spec();p=independent_packet(spec)
 low=np.array(p['low_g_squared']);tz=math.log(p['mz']);tU=math.log(p['mu_U']);tE=math.log(p['E_cell']);tL=math.log(50.)
 bvec=np.array([41/6,-19/6,-7]);cvec=np.array([17/12,9/4,8])
 def den(t):return 1/low-2*k*bvec*(t-tz)
 def gauge(t):return 1/den(t)
 def F(t):
  a,b,c=gauge(t);return ((a*a+2*a*b+3*b*b)/16)**.25
 def proxy(m,c):
  a=c/(4*pi*pi);return m*(1+4*a/3+(13.4434-1.0414*5)*a*a+(190.595-26.655*5+.653*25)*a**3)
 one=[];two=[]
 for u in spec['anchors']['control_fractions']:
  tb=(1-u)*tU+u*tE
  def dz(t,z):return 2*k*float(cvec@gauge(t))*z-9*k
  sol=solve_ivp(dz,(tb,tL),[1/F(tb)**2],method='Radau',rtol=2e-12,atol=2e-14,dense_output=True)
  require(sol.success,'linearized Yukawa integration failed')
  mass=brentq(lambda m:p['v']/math.sqrt(2*sol.sol(math.log(m))[0])-m,50.,p['v'],xtol=1e-11)
  a,b,c=gauge(tb);ap,bp,cp=2*k*bvec*gauge(tb)**2
  X=a*a+2*a*b+3*b*b
  seed=F(tb)*((2*a*ap+2*ap*b+2*a*bp+6*b*bp)/(4*X)-k*(4.5*F(tb)**2-float(cvec@gauge(tb))))
  one.append(dict(weight_E_over_sum=u,mu_boundary=math.exp(tb),y_boundary=F(tb),mass_coordinate=mass,
                  QCD_pole_proxy=proxy(mass,gauge(math.log(mass))[2]),boundary_sensitivity_seed=seed))
  def weak(g):
   hh=brentq(lambda h:beta(0,[*g,h,0.])[4],.05**2,2**2,xtol=1e-14)
   a,b,c=g
   rr=np.roots([30*k*k,-6*k+k*k*(-32*c-8*a/3),k*k*(-9*b*b/4-19*a*a/4+21*a*b/2),
    k*3*(2*b*b+(b+a)**2)/8+k*k*(305*b**3/16-379*a**3/48-289*b*b*a/48-559*b*a*a/48)])
   require(all(abs(x.imag)<1e-9 for x in rr),'critical root inventory')
   require(abs(hh-min(x.real for x in rr if x.real>0))<1e-10,'nonperturbative critical root')
   return hh
  def evolve(g,tol):
   ss=solve_ivp(beta,(tb,tL),[*g,weak(g),0.],method='Radau',rtol=tol,atol=tol/100,dense_output=True)
   require(ss.success,'independent coupled integration failed');return ss
  guess=gauge(tb)
  for tol in [2e-10,2e-12]:
   def residual(g):return evolve(g,tol).sol(tz)[:3]-low
   shot=root(residual,guess,tol=1e-10)
   require(shot.success and np.max(abs(residual(shot.x)))<1e-9,'independent shooting failed')
   guess=shot.x
  sol=evolve(guess,2e-12)
  mass=brentq(lambda m:math.sqrt(sol.sol(math.log(m))[3])*p['v']/math.sqrt(2)-m,50.,p['v'],xtol=1e-11)
  a,b,c,h,l=sol.sol(math.log(mass))
  two.append(dict(weight_E_over_sum=u,mu_boundary=math.exp(tb),mass_coordinate=mass,
                  QCD_pole_proxy=proxy(mass,c),Higgs_tree_proxy=p['v']*math.sqrt(2*l)))
 gL,gU,gE=gauge(tL),gauge(tU),gauge(tE)
 checks=dict(positive_inverse_gauge_minimum=float(min(np.min(den(tL)),np.min(den(tE)))),
  c_over_b_on_anchor_endpoints=[gU[2]/gU[1],gE[2]/gE[1]],
  c_over_a_on_domain_endpoints=[gL[2]/gL[0],gE[2]/gE[0]],
  positive_sensitivity_coefficient=4-11/12-9*math.sqrt(3)/8,
  negative_Cprime_coefficient=697/36-28,
  uniform_lower_y_at_50=(3*gE[1]**2/16)**.25,required_y_at_50=50*math.sqrt(2)/p['v'],
  uniform_upper_y_at_mu_U=math.sqrt(2*float(cvec@gU)/9),required_y_at_mu_U=p['mu_U']*math.sqrt(2)/p['v'])
 return dict(packet=p,one=one,two=two,checks=checks)


def verify_payload(receipt):
 spec=reviewed_spec()
 keys(receipt,['schema','specification','scope','packet','proof_checks','one_loop_monotone_band_endpoints',
  'one_loop_controls','coupled_two_loop_controls','numerical_interpretation','source_pins'],'receipt')
 require(receipt['schema']=='oph.top_positive_capacity_band.v1','wrong schema')
 require(receipt['specification']==spec,'specification changed')
 require(receipt['scope']==spec['scope'],'scope promotion')
 require(receipt['numerical_interpretation']==spec['uncertainty'],'uncertainty promotion')
 keys(receipt['source_pins'],PINS,'source pins')
 for path,digest in receipt['source_pins'].items():
  require(digest==hashlib.sha256((REPO/path).read_bytes()).hexdigest(),'source pin drift')
 source_ancestry_check()
 ref=independent_expectation()
 keys(receipt['packet'],ref['packet'],'source packet')
 for key,value in ref['packet'].items():
  actual=receipt['packet'][key]
  if isinstance(value,str):require(actual==value,'alpha calibration changed')
  elif isinstance(value,list):
   require(isinstance(actual,list) and len(actual)==len(value),'gauge packet shape')
   for a,b in zip(actual,value):near(a,b,atol=0,rtol=2e-11)
  else:near(actual,value,atol=0,rtol=2e-11)
 keys(receipt['proof_checks'],ref['checks'],'analytic hypothesis margins')
 for key,value in ref['checks'].items():
  actual=receipt['proof_checks'][key]
  if isinstance(value,list):
   require(isinstance(actual,list) and len(actual)==len(value),'margin shape')
   for a,b in zip(actual,value):near(a,b,atol=0,rtol=2e-11)
  else:near(actual,value,atol=0,rtol=2e-11)
 pc=receipt['proof_checks']
 require(pc['positive_inverse_gauge_minimum']>0 and min(pc['c_over_b_on_anchor_endpoints'])>.5
  and min(pc['c_over_a_on_domain_endpoints'])>.5,'analytic gauge-ratio hypotheses fail')
 require(pc['positive_sensitivity_coefficient']>0 and pc['negative_Cprime_coefficient']<0,'analytic coefficient sign fails')
 require(pc['uniform_lower_y_at_50']>pc['required_y_at_50'] and pc['uniform_upper_y_at_mu_U']<pc['required_y_at_mu_U'],'uniform self-scale bracket fails')
 rows=receipt['one_loop_controls']
 require(isinstance(rows,list) and len(rows)==5,'complete one-loop control inventory required')
 for actual,expected in zip(rows,ref['one']):
  keys(actual,[*expected,'self_scale_residual','quadrature_error_estimate_in_inverse_y_squared','independent_ODE_mass','ODE_difference'],'one-loop control')
  require(actual['weight_E_over_sum']==expected['weight_E_over_sum'],'capacity fraction changed')
  for key,value in expected.items():near(actual[key],value,atol=0 if key=='mu_boundary' else 2e-7,rtol=2e-11)
  near(actual['independent_ODE_mass'],expected['mass_coordinate'],atol=2e-7)
  near(actual['ODE_difference'],actual['independent_ODE_mass']-actual['mass_coordinate'],atol=1e-12,rtol=0)
  require(abs(number(actual['ODE_difference']))<1e-6,'quadrature/ODE mismatch')
  require(abs(number(actual['self_scale_residual']))<1e-8,'self-scale residual')
  require(0<=number(actual['quadrature_error_estimate_in_inverse_y_squared'])<1e-8,'quadrature error diagnostic')
 band=receipt['one_loop_monotone_band_endpoints']
 require(isinstance(band,list) and len(band)==2,'band endpoint inventory')
 for value,row in zip(band,[ref['one'][0],ref['one'][-1]]):near(value,row['mass_coordinate'],atol=2e-7)
 require(band==[rows[0]['mass_coordinate'],rows[-1]['mass_coordinate']],'band detached from limiting controls')
 rows=receipt['coupled_two_loop_controls']
 require(isinstance(rows,list) and len(rows)==5,'complete two-loop controls required')
 for actual,expected in zip(rows,ref['two']):
  keys(actual,['weight_E_over_sum','mu_boundary','rows','tolerance_mass_difference'],'coupled control')
  require(actual['weight_E_over_sum']==expected['weight_E_over_sum'],'capacity control altered')
  near(actual['mu_boundary'],expected['mu_boundary'],atol=0)
  require(isinstance(actual['rows'],list) and len(actual['rows'])==2,'both tolerance rows required')
  for row,tol in zip(actual['rows'],spec['transport']['two_loop_rtol']):
   keys(row,['rtol','mass_coordinate','QCD_pole_proxy','Higgs_tree_proxy','gauge_anchor_max_residual','criticality_residual','self_scale_residual'],'coupled tolerance row')
   require(row['rtol']==tol,'transport tolerance changed')
   for key in ['mass_coordinate','QCD_pole_proxy','Higgs_tree_proxy']:near(row[key],expected[key],atol=2e-6)
   require(0<=number(row['gauge_anchor_max_residual'])<1e-8,'gauge anchor residual')
   require(abs(number(row['criticality_residual']))<1e-12,'criticality residual')
   require(abs(number(row['self_scale_residual']))<1e-8,'self-scale residual')
  difference=abs(actual['rows'][0]['mass_coordinate']-actual['rows'][1]['mass_coordinate'])
  near(actual['tolerance_mass_difference'],difference,atol=1e-12,rtol=0)
  require(difference<1e-5,'coupled numerical convergence')
 exact_identity_checks()
 return {'verified':True,'one_loop_controls':5,'two_loop_controls':5,
         'independent_scalar_linear_ODE_and_coupled_replay':True,
         'formal_theorem_proof_claimed':False,'physical_matching_or_two_loop_interval_certified':False}


if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('receipt',type=Path,nargs='?',default=OUT)
 args=parser.parse_args();print(json.dumps(verify_payload(strict_load(args.receipt)),sort_keys=True))
