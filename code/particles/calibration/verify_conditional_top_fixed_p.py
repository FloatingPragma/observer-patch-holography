"""Independent numerical verification of the conditional top coordinate.

No producer, PaperMathContext, legacy transport or repository RGE import.
Heat sums are reconstructed, beta functions use squared non-GUT hypercharge,
and Radau replaces the producer's DOP853/RK4 integrators. This verifies the
specified approximation, not a physical matching or an uncertainty enclosure.
"""
from __future__ import annotations
import argparse
from decimal import Decimal, localcontext
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SPEC_PATH=HERE/'conditional_top_fixed_p_spec.json'
OUT_PATH=HERE.parent/'runs/calibration/conditional_top_fixed_p.json'
REVIEWED_SPEC_SHA256='842c2ae23168141de92f2348e0480baf1500963850683da29358a883a7281c6b'
PIN_PATHS=(
 'code/P_derivation/codata_2022_alpha_fixture.json',
 'code/P_derivation/paper_math.py',
 'code/particles/calibration/derive_d11_criticality_boundary_scan.py',
 'code/particles/calibration/sm_two_loop_rge_engine.py',
 'code/particles/calibration/conditional_top_fixed_p_spec.json',
 'code/particles/calibration/conditional_top_fixed_p.py',
)
pi=math.pi;k=1/(16*pi*pi)


def require(ok,message):
 if not ok:raise ValueError(message)


def number(value):
 require(isinstance(value,(int,float,str)) and not isinstance(value,bool),'not a numerical scalar')
 try:result=float(value)
 except (ValueError,TypeError,OverflowError) as exc:raise ValueError('invalid numerical scalar') from exc
 require(math.isfinite(result),'nonfinite numerical scalar')
 return result


def near(value,expected,absolute=2e-9,relative=2e-10):
 actual=number(value);expected=float(expected)
 require(abs(actual-expected)<=absolute+relative*abs(expected),f'numerical mismatch: {actual} versus {expected}')


def keys(value,expected,name):
 require(isinstance(value,dict) and set(value)==set(expected),name+' schema mismatch')


def strict_load(path):
 def pairs(items):
  result={}
  for key,value in items:
   require(key not in result,'duplicate JSON key');result[key]=value
  return result
 return json.loads(Path(path).read_text(),object_pairs_hook=pairs,
  parse_constant=lambda text:(_ for _ in ()).throw(ValueError('nonfinite JSON token')))


def reviewed_spec():
 require(hashlib.sha256(SPEC_PATH.read_bytes()).hexdigest()==REVIEWED_SPEC_SHA256,'unreviewed specification')
 return strict_load(SPEC_PATH)


@lru_cache(maxsize=1)
def independent_expectation():
 spec=reviewed_spec()
 P=(1+math.sqrt(5))/2+math.sqrt(pi)/float(spec['inverse_alpha'])
 E=float(spec['E_star_display_GeV'])
 U=math.exp(-2*pi)*P**(1/6);cell=P**(-.5)
 def heat(t,kind,cut=40):
  terms=[]
  if kind==2:
   for n in range(cut+1):
    dim=n+1;c=n*(n+2)/4;w=dim*math.exp(-t*c);terms.append((w,math.log(dim)))
  else:
   for i in range(cut+1):
    for j in range(cut+1):
     dim=(i+1)*(j+1)*(i+j+2)/2;c=(i*i+j*j+i*j+3*i+3*j)/3
     w=dim*math.exp(-t*c);terms.append((w,math.log(dim)))
  return math.fsum(w*d for w,d in terms)/math.fsum(w for w,d in terms)

 def low_packet(au):
  v=cell*math.exp(-pi/(2*au))
  def aa(logmz):return [1/(1/au+b/(2*pi)*(math.log(U)-logmz)) for b in [33/5,1,-3]]
  def eq(logmz):
   a1,a2,a3=aa(logmz)
   return math.log(v/2*math.sqrt(4*pi*(a2+3*a1/5)))-logmz
  z=brentq(eq,math.log(v)-2,math.log(v)+1,xtol=1e-13)
  a1,a2,a3=aa(z)
  return v,math.exp(z),a1,a2,a3

 def hclose(au):
  _,_,_,a2,a3=low_packet(au)
  return heat(4*pi*pi*a2,2)+heat(4*pi*pi*a3,3)-P/4
 au=brentq(hclose,.04,.042,xtol=5e-16)
 vbar,zbar,a1,a2,a3=low_packet(au)
 V=vbar*E;Z=zbar*E;B=math.sqrt(U*cell)*E
 low=np.array([4*pi*3*a1/5,4*pi*a2,4*pi*a3])

 # x=(gY^2,g2^2,g3^2,yt^2,lambda); d/d(log mu).
 # These are transformed Buttazzo app.B coefficients in non-GUT hypercharge.
 def beta(t,x):
  a,b,c,h,l=x
  da=2*a*a*(k*41/6+k*k*(199*a/18+9*b/2+44*c/3-17*h/6))
  db=2*b*b*(-k*19/6+k*k*(3*a/2+35*b/6+12*c-3*h/2))
  dc=2*c*c*(-7*k+k*k*(11*a/6+9*b/2-26*c-2*h))
  dh=2*h*(k*(9*h/2-17*a/12-9*b/4-8*c)+k*k*(
   -12*h*h-12*l*h+6*l*l+h*(131*a/16+225*b/16+36*c)
   +1187*a*a/216-3*a*b/4+19*a*c/9-23*b*b/4+9*b*c-108*c*c))
  dl=k*(24*l*l-6*h*h+3*(2*b*b+(b+a)**2)/8+l*(12*h-9*b-3*a))
  dl+=k*k*(-312*l**3-144*l*l*h+l*l*(108*b+36*a)
   +l*h*(-3*h+80*c+45*b/2+85*a/6)
   +l*(-73*b*b/8+629*a*a/24+39*a*b/4)
   +h*h*(30*h-32*c-8*a/3)+h*(-9*b*b/4-19*a*a/4+21*a*b/2)
   +305*b**3/16-379*a**3/48-289*b*b*a/48-559*b*a*a/48)
  return np.array([da,db,dc,dh,dl])

 def weak_h(g):return brentq(lambda h:beta(0,[*g,h,0])[4],.05**2,2**2,xtol=1e-14)

 def evolve(g,tol=2e-10):
  sol=solve_ivp(beta,(math.log(B),math.log(40)),[*g,weak_h(g),0.],method='Radau',rtol=tol,atol=tol/100,dense_output=True)
  assert sol.success
  return sol

 def residual(g,tol=2e-10):return evolve(g,tol).sol(math.log(Z))[:3]-low
 # One-loop boundary estimate uses this audit's independently solved low packet.
 guess=1/(1/low-np.array([41/6,-19/6,-7])/(8*pi*pi)*math.log(B/Z))
 shot=root(residual,guess,tol=2e-10)
 require(shot.success,'independent shooting failed')
 sr=root(lambda g:residual(g,2e-12),shot.x,tol=1e-10)
 require(sr.success,'independent refined shooting failed')
 sol=evolve(sr.x,2e-12)
 mass=brentq(lambda m:math.sqrt(sol.sol(math.log(m))[3])*V/math.sqrt(2)-m,40.,V,xtol=1e-11)
 x=sol.sol(math.log(mass))
 def qcd(m,alpha):
  a=alpha/pi
  return m*(1+4*a/3+(13.4434-1.0414*5)*a*a+(190.595-26.655*5+.653*25)*a**3)
 ay,aw,ass=sr.x
 coefficients=[30*k*k,-6*k+k*k*(-32*ass-8*ay/3),
  k*k*(-9*aw*aw/4-19*ay*ay/4+21*ay*aw/2),
  k*3*(2*aw*aw+(aw+ay)**2)/8+k*k*(305*aw**3/16-379*ay**3/48-289*aw*aw*ay/48-559*aw*ay*ay/48)]
 polynomial_roots=sorted(float(r) for r in np.roots(coefficients))
 coupled={'boundary_gauges':[math.sqrt(sr.x[0]*5/3),math.sqrt(sr.x[1]),math.sqrt(sr.x[2])],
  'boundary_top_yukawa':math.sqrt(weak_h(sr.x)),
  'critical_polynomial_roots_yt_squared':polynomial_roots,
  'top_running_mass_coordinate_GeV':mass,'top_QCD_pole_proxy_GeV':qcd(mass,x[2]/(4*pi)),
  'Higgs_tree_proxy_GeV':V*math.sqrt(2*x[4]),'alpha_s_self_scale':x[2]/(4*pi),
  'lambda_self_scale':x[4],'yukawa_self_scale':math.sqrt(x[3])}
 # Independent two-state hybrid integrations keep analytic one-loop gauge paths.
 hybrids={}
 def gauge(t):return 1/(1/low-np.array([41/6,-19/6,-7])/(8*pi*pi)*(t-math.log(Z)))
 for loops in [1,2]:
  gb=gauge(math.log(B))
  h0=math.sqrt((2*gb[1]**2+(gb[1]+gb[0])**2)/16) if loops==1 else weak_h(gb)
  def hybrid_beta(t,z):
   a,b,c=gauge(t);h,l=z
   if loops==2:return beta(t,[a,b,c,h,l])[3:]
   return [2*h*k*(9*h/2-17*a/12-9*b/4-8*c),
     k*(24*l*l-6*h*h+3*(2*b*b+(b+a)**2)/8+l*(12*h-9*b-3*a))]
  hs=solve_ivp(hybrid_beta,(math.log(B),math.log(40)),[h0,0],method='Radau',rtol=2e-12,atol=2e-14,dense_output=True)
  require(hs.success,'independent hybrid integration failed')
  mm=brentq(lambda m:math.sqrt(hs.sol(math.log(m))[0])*V/math.sqrt(2)-m,40.,V,xtol=1e-11)
  hh,ll=hs.sol(math.log(mm));al=gauge(math.log(mm))[2]/(4*pi)
  hybrids[str(loops)]={'boundary_scale_gev':B,'y_t_u':math.sqrt(h0),'y_t_mt':math.sqrt(hh),
   'lambda_mt':ll,'alpha_s_mt':al,'top_running_mass_coordinate_GeV':mm,
   'top_QCD_pole_proxy_GeV':qcd(mm,al),'Higgs_tree_proxy_GeV':V*math.sqrt(2*ll)}
 alphaY=3*a1/5;alphaEM=1/(1/a2+1/alphaY)
 return {'E':E,'P':P,'D10':{'n_c':3,'mu_u':U,'alpha_u':au,'mz_run':zbar,'v':vbar,
   'alpha1_mz':a1,'alpha2_mz':a2,'alpha3_mz':a3,'alpha_y_mz':alphaY,
   'alpha_em_mz':alphaEM,'alpha_em_inv_mz':1/alphaEM,'sin2w_mz':alphaEM/a2},
  'scales':{'mu_U_gauge_unification':U*E,'E_cell_pixel_energy':cell*E,'E_star':E,
   'log_midpoint_half_turn':B,'v_transmutation_gev':V,'mz_run_gev':Z},
  'hybrid':hybrids,'coupled':coupled,'heat_residual':hclose(au)}


def verify_payload(receipt):
 spec=reviewed_spec()
 keys(receipt,['schema','specification','scope','source_pins','P_from_independent_alpha','D10',
   'source_scales_GeV','hybrid_rows','fully_coupled_two_loop_same_low_gauge_anchor','numerical_diagnostics'],'receipt')
 require(receipt['schema']=='oph.conditional_top_fixed_p.v1','wrong receipt schema')
 require(receipt['specification']==spec,'conditional specification changed')
 require(receipt['scope']==spec['scope'],'scope promotion or missing qualification')
 keys(receipt['source_pins'],PIN_PATHS,'source pins')
 for path,digest in receipt['source_pins'].items():
  require(digest==hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),'stale source pin: '+path)
 fixture=strict_load(ROOT/spec['alpha_fixture'])['inverse_fine_structure_constant']
 require(fixture=={'value':spec['inverse_alpha'],'standard_uncertainty':spec['inverse_alpha_standard_uncertainty']},'alpha fixture mismatch')
 # Independent high-precision elementary P identity; source root is not fitted.
 with localcontext() as context:
  context.prec=65
  pi_decimal=Decimal('3.1415926535897932384626433832795028841971693993751058209749445923078')
  p=(1+Decimal(5).sqrt())/2+pi_decimal.sqrt()/Decimal(spec['inverse_alpha'])
  for value in [receipt['P_from_independent_alpha'],receipt['D10']['p']]:
   try:delta=abs(Decimal(value)-p)
   except Exception as exc:raise ValueError('invalid P') from exc
   require(delta.is_finite() and delta<Decimal('1e-55'),'independent P identity failed')
 ref=independent_expectation()
 keys(receipt['D10'],list(ref['D10'])+['p'],'D10')
 for name,value in ref['D10'].items():near(receipt['D10'][name],value,absolute=0,relative=1e-10)
 keys(receipt['source_scales_GeV'],ref['scales'],'source scales')
 for name,value in ref['scales'].items():near(receipt['source_scales_GeV'][name],value,absolute=0,relative=1e-10)
 keys(receipt['hybrid_rows'],['1','2'],'hybrid branches')
 for branch,rows in receipt['hybrid_rows'].items():
  require(isinstance(rows,list) and len(rows)==2,'missing hybrid convergence row')
  for row,steps in zip(rows,spec['transport']['hybrid_RK4_steps']):
   keys(row,list(ref['hybrid'][branch])+['top_running_coordinate_over_Estar','top_QCD_pole_proxy_over_Estar','RK4_steps'],'hybrid row')
   require(row['RK4_steps']==steps,'hybrid step count changed')
   for name,value in ref['hybrid'][branch].items():near(row[name],value,absolute=2e-6,relative=2e-10)
   near(row['top_running_coordinate_over_Estar'],number(row['top_running_mass_coordinate_GeV'])/ref['E'],absolute=0,relative=1e-12)
   near(row['top_QCD_pole_proxy_over_Estar'],number(row['top_QCD_pole_proxy_GeV'])/ref['E'],absolute=0,relative=1e-12)
 rows=receipt['fully_coupled_two_loop_same_low_gauge_anchor']
 require(isinstance(rows,list) and len(rows)==2,'missing coupled convergence row')
 extras=['relative_tolerance','low_anchor_gauge_residual','top_running_coordinate_over_Estar',
  'top_QCD_pole_proxy_over_Estar','top_fixed_point_residual_GeV','beta_lambda_boundary']
 for row,tolerance in zip(rows,spec['transport']['coupled_rtol']):
  keys(row,list(ref['coupled'])+extras,'coupled row')
  require(row['relative_tolerance']==tolerance,'coupled tolerance changed')
  for name,value in ref['coupled'].items():
   if isinstance(value,list):
    require(isinstance(row[name],list) and len(row[name])==len(value),'root/gauge inventory changed')
    for actual,expected in zip(row[name],value):near(actual,expected,absolute=2e-9,relative=2e-10)
   else:near(row[name],value,absolute=2e-6,relative=2e-10)
  root_values=list(map(number,row['critical_polynomial_roots_yt_squared']))
  require(root_values[0]<0<root_values[1]<4<root_values[2],'weak-root branch changed')
  near(number(row['boundary_top_yukawa'])**2,root_values[1],absolute=1e-12,relative=1e-12)
  require(len(row['low_anchor_gauge_residual'])==3,'wrong residual dimension')
  require(max(abs(number(v)) for v in row['low_anchor_gauge_residual'])<1e-8,'gauge shooting residual')
  require(abs(number(row['beta_lambda_boundary']))<1e-12,'criticality residual')
  require(abs(number(row['top_fixed_point_residual_GeV']))<1e-8,'top fixed-point residual')
  near(row['top_running_mass_coordinate_GeV'],number(row['yukawa_self_scale'])*ref['scales']['v_transmutation_gev']/math.sqrt(2),absolute=2e-7,relative=1e-10)
  near(row['Higgs_tree_proxy_GeV'],math.sqrt(2*number(row['lambda_self_scale']))*ref['scales']['v_transmutation_gev'],absolute=2e-7,relative=1e-10)
  for label,num in [('top_running_coordinate_over_Estar','top_running_mass_coordinate_GeV'),('top_QCD_pole_proxy_over_Estar','top_QCD_pole_proxy_GeV')]:
   near(row[label],number(row[num])/ref['E'],absolute=0,relative=1e-12)
 diagnostics=receipt['numerical_diagnostics']
 keys(diagnostics,['hybrid_top_step_differences_GeV','coupled_top_tolerance_difference_GeV','hybrid_to_coupled_top_shift_GeV','interpretation'],'diagnostics')
 require(diagnostics['interpretation']==spec['uncertainty'],'numerical error promoted to theory bound')
 keys(diagnostics['hybrid_top_step_differences_GeV'],['1','2'],'hybrid differences')
 for branch,rr in receipt['hybrid_rows'].items():
  near(diagnostics['hybrid_top_step_differences_GeV'][branch],abs(number(rr[1]['top_running_mass_coordinate_GeV'])-number(rr[0]['top_running_mass_coordinate_GeV'])),absolute=1e-12,relative=0)
 near(diagnostics['coupled_top_tolerance_difference_GeV'],abs(number(rows[1]['top_running_mass_coordinate_GeV'])-number(rows[0]['top_running_mass_coordinate_GeV'])),absolute=1e-12,relative=0)
 near(diagnostics['hybrid_to_coupled_top_shift_GeV'],number(receipt['hybrid_rows']['2'][-1]['top_running_mass_coordinate_GeV'])-number(rows[-1]['top_running_mass_coordinate_GeV']),absolute=1e-12,relative=0)
 return {'verified':True,'independent_integrator':'Radau','independent_variables':'gY_squared,g2_squared,g3_squared,yt_squared,lambda',
  'branches_checked':3,'top_running_mass_coordinate_GeV':ref['coupled']['top_running_mass_coordinate_GeV'],
  'physical_matching_or_theory_uncertainty_certified':False}


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('receipt',nargs='?',type=Path,default=OUT_PATH)
 args=parser.parse_args()
 print(json.dumps(verify_payload(strict_load(args.receipt)),indent=2))


if __name__=='__main__':main()
