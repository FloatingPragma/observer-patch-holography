"""Independent native reference, flux graph and autonomous current audit."""
from pathlib import Path
from itertools import product
from collections import defaultdict
import hashlib,json
import sympy as s
import numpy as np
from scipy.linalg import expm
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as pw
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.clebsch import coordinate
from model.openwave_autonomous_coulomb.source import basis,matrices,REFERENCE_MASK
ROOT=Path(__file__).resolve().parents[3];TARGET=ROOT/'model/openwave_autonomous_coulomb'

def ref_current(matrix,mask):
 out=defaultdict(lambda:s.Integer(0))
 for spin in range(2):
  for i,j in matrix.todok():
   a,b=15*spin+i,15*spin+j
   if not(mask>>b&1):continue
   n=mask^(1<<b);sign=(-1)**((mask&((1<<b)-1)).bit_count())
   if n>>a&1:continue
   sign*=(-1)**((n&((1<<a)-1)).bit_count());n|=1<<a
   out[n]+=sign*matrix[i,j]
 return {k:s.simplify(v)for k,v in out.items()if s.simplify(v)!=0}

def main():
 gauges=representation()[2]
 for j,T in enumerate(gauges):
  actual=ref_current(T,REFERENCE_MASK)
  assert actual==({REFERENCE_MASK:-1}if j==11 else{})
 U=coordinate(6,0,0);labels=(0,46,105);Y=pw.endpoint_current(11,0);right=pw.endpoint_current(11,1)
 assert [Y[j,j]for j in labels]==[0,-1,1]
 assert (Y*U-U*Y+U)==s.zeros(119)
 assert right[:,labels]==-Y[:,labels]
 P=s.diag(*[int(i in labels)for i in range(119)])
 assert U*P==P*U and U.H*P==P*U.H
 # Independently enumerate all flux circulation sectors, without package basis.
 states=[]
 for spin,y,f in product(range(2),range(4),product((-1,0,1),repeat=4)):
  if all(int(z==y)-int(z==0)+f[z]-f[(z-1)%4]==0 for z in range(4)):states.append((spin,y,f))
 assert len(states)==18
 where={v:i for i,v in enumerate(states)};H=s.zeros(18);pot=[]
 # Exact three-dimensional N4 Green sum from finite Fourier characters.
 for spin,y,f in states:
  value=s.Rational(0)
  for k in product(range(4),repeat=3):
   if k==(0,0,0):continue
   lap=4*sum((s.sin(s.pi*n/4)**2 for n in k),s.Rational(0))
   value+=(1-s.cos(s.pi*k[2]*y/2))/(64*lap)
  pot.append(s.simplify(value))
 V=s.diag(*pot)
 for col,(spin,y,f)in enumerate(states):
  H[col,col]=118+pot[col]
  for edge in range(4):
   coeff=(1 if edge<3 else-1)*s.I*s.sqrt(2)/8*(1 if spin==0 else-1)
   for direct in (True,False):
    start=(edge+1)%4 if direct else edge;end=edge if direct else(edge+1)%4
    if y!=start:continue
    target=list(f);target[edge]+=(-1 if direct else 1)
    if target[edge]not in (-1,0,1):continue
    H[where[spin,end,tuple(target)],col]+=coeff if direct else s.conjugate(coeff)
 assert H==H.H
 rows=basis();order=[where[spin,y,tuple(map(int,f))]for _,spin,y,f in rows]
 referenceH,referenceV=matrices()
 assert H.extract(order,order)==referenceH and V.extract(order,order)==referenceV
 # All positions, spins and loops enter; direct current derivative includes boundary sign.
 currents=[];theta=s.Symbol('theta',real=True)
 for edge in range(4):
  J=s.zeros(18)
  for a,(sa,ya,fa)in enumerate(states):
   for b,(sb,yb,fb)in enumerate(states):
    if a!=b and sum(x!=y for x,y in zip(fa,fb))==1 and fa[edge]!=fb[edge]:
     J[a,b]=s.I*(int(ya==edge)-int(yb==edge))*H[a,b]
  assert J==J.H;currents.append(J)
 for y in range(4):
  Q=s.diag(*[int(state[1]==y)for state in states])
  assert s.I*(H*Q-Q*H)==currents[(y-1)%4]-currents[y]
 left=next(i for i,(spin,y,f)in enumerate(states)if spin==0 and y==1 and f==(0,0,0,1))if (0,1,(0,0,0,1))in where else None
 # Select any connected outward bond1->2 and its lower-kinetic coherent state.
 left=next(i for i,(spin,y,f)in enumerate(states)if spin==0 and y==1)
 right_index=next(i for i,(spin,y,f)in enumerate(states)if spin==0 and y==2 and H[i,left]!=0)
 tau=s.sqrt(2)/8;p=s.zeros(18,1);p[left]=1/s.sqrt(2);p[right_index]=-H[right_index,left]/(s.sqrt(2)*tau)
 force=s.simplify((p.H*s.I*(V*currents[1]-currents[1]*V)*p)[0]);assert force==-19*s.sqrt(2)/5120
 A=np.array(H-118*s.eye(18),complex);initial=np.array(p,complex).ravel();numerical=[]
 for t in (.1,1.,3.,8.):
  v=expm(-1j*t*A)@initial;assert abs(np.vdot(v,v)-1)<1e-14
  assert abs(np.vdot(v,A@v)-np.vdot(initial,A@initial))<1e-14
  J=np.array(currents[1],complex);dt=1e-5
  j=lambda time:np.vdot(expm(-1j*time*A)@initial,J@(expm(-1j*time*A)@initial)).real
  residual=abs((j(t+dt)-j(t-dt))/(2*dt)-np.vdot(v,1j*(A@J-J@A)@v).real)
  assert residual<2e-10;numerical.append(residual)
 receipt=json.loads((TARGET/'results.json').read_text())
 for path,h in receipt['source_sha256'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,path
 paths=[p for p in TARGET.iterdir()if p.suffix in ('.py','.md','.tex','.json')]
 result={'verdict':'PASS','independent_reference_full_CAR_current_checks':12,'independent_PW119_reducing_projection':True,
   'independent_complete_flux_spin_position_rows':18,'native_matrix_entries_compared':324,
   'independent_exact_three_dimensional_Green_sum':True,'all_four_native_continuity_identities':True,
   'exact_potential_current_acceleration':str(force),'independent_full_unitary_current_derivative_residuals':numerical,
   'scientific_scope':'Autonomous motion of the declared zero-Yukawa, selected-hop native source with positive number energy and longitudinal Green interaction. All circulation sectors and original source/reference factors are retained; no transfer of the unit-Yukawa interacting ground, physical alpha or transverse radiation is inferred.',
   'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}}
 (Path(__file__).resolve().parents[1]/'autonomous_coulomb_review.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print('PASS_INDEPENDENT_AUTONOMOUS_NATIVE_COULOMB')
if __name__=='__main__':main()
