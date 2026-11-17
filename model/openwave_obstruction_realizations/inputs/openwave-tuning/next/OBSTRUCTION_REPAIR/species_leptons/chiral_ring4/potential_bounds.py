"""Interval complete-band lower bounds for every infinite-rotor sector."""
import json,math
from pathlib import Path
import numpy as np
import mpmath as mp
mp.iv.dps=40
ROOT=Path(__file__).parent
GRID=1024;L=4;S=2;M=16
lower=np.full(M+1,np.inf);upper=np.full(M+1,np.inf);arg=np.zeros(M+1,dtype=int)
AMPS=[0,.25,-.25,.5,-.5,1.,-1.,2.,-2.]
barta=np.full((M+1,len(AMPS)),np.inf)
for j in range(GRID):
 theta=2*mp.iv.pi*j/GRID;hi=[];lo=[]
 trial_terms=[]
 for amplitude in AMPS:
  correction=(mp.iv.mpf(amplitude)*mp.iv.cos(theta)-mp.iv.mpf(amplitude)**2*mp.iv.sin(theta)**2)*L/4
  trial_terms.append(np.nextafter(float(correction.a),-np.inf))
 for q in range(L):
  k=(theta+2*mp.iv.pi*q)/L
  W=2*mp.iv.sin(k/2)**2
  ep2=2*W+mp.iv.mpf(1)/2+mp.iv.sqrt(1+4*W**2)/2
  em2=(3*W**2+2*W)/ep2
  for e in [8*mp.iv.sqrt(ep2),8*mp.iv.sqrt(em2)]:
   hi.append(np.nextafter(float(e.b),np.inf));lo.append(np.nextafter(float(e.a),-np.inf))
 hi.sort(reverse=True);lo.sort(reverse=True)
 hc=0.;lc=0.
 for n in range(9):
  if n:
   hc=np.nextafter(hc+hi[n-1],np.inf);lc=np.nextafter(lc+lo[n-1],-np.inf)
  a=-hc;b=-lc
  for N in {n,M-n}:
   if a<lower[N]:lower[N]=a;arg[N]=j
   upper[N]=min(upper[N],b)
   for ia,amplitude in enumerate(AMPS):
    term=trial_terms[ia]
    barta[N,ia]=min(barta[N,ia],np.nextafter(a+term,-np.inf))
rows=[]
for N in range(M+1):
 n=min(S,N,M-N)
 c=8*mp.iv.sqrt(((n+1)**2)//4)
 # Nearest grid point is at distance pi/GRID; derivative norm<=2c.
 correction=float((2*c*mp.iv.pi/GRID).b)
 bound=np.nextafter(lower[N]-correction,-np.inf)
 barta_rows=[]
 for ia,amplitude in enumerate(AMPS):
  err=((2*c+L/4*(abs(amplitude)+amplitude**2))*mp.iv.pi/GRID)
  berr=np.nextafter(float(err.b),np.inf)
  barta_rows.append({'a':amplitude,'grid_lower':float(barta[N,ia]),
   'Lipschitz_correction_upper':float(berr),'infinite_rotor_lower':float(np.nextafter(barta[N,ia]-berr,-np.inf))})
 rows.append({'N':N,'matter_grid_min_lower':float(lower[N]),'grid_min_upper':float(upper[N]),
  'global_Lipschitz_correction_upper':float(np.nextafter(correction,np.inf)),
  'infinite_rotor_sector_lower_m0':float(bound),'grid_argmin_index':int(arg[N]),
  'all_momenta':L,'all_one_particle_bands':M,'all_Barta_tests':barta_rows,
  'best_Barta_infinite_rotor_lower':max(x['infinite_rotor_lower']for x in barta_rows)})
result={'L':L,'S':S,'CAR_modes':M,'holonomy_grid':GRID,'interval_decimal_digits':40,
 'interval_backend':'mpmath.iv','electric_energy_positive':True,
 'reference_charge_energy':0,'all_number_sectors':rows,
 'nonzero_mass_uniform_lower_shift_bound':'L*abs(m), except N0 andN16 where zero'}
(ROOT/'POTENTIAL_BOUNDS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(rows,indent=2))
