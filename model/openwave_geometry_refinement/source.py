"""Fixed cap-four source and its nested occupation compressions."""
from functools import lru_cache
from itertools import product
from fractions import Fraction as F
import math
import numpy as np
from scipy.linalg import eigh
from mpmath import iv
from model.openwave_dynamic_wave.source import shape,packet,packet_integers,energies,KAPPA
from model.openwave_dynamic_wave.certificate import scalar,ends,interval,rat,number
from model.openwave_energy_lapse.source import native
from model.openwave_quantum_geometry.certificate import sqceil
CAP=4
DEN=2**42
ROOTDEN=2**46
UNIT=F(1,2**53)
def upper(x,bits=80):
 x=F(x);return F((x.numerator*2**bits+x.denominator-1)//x.denominator,2**bits)
def gamma(n):return n*UNIT/(1-n*UNIT)
def fupper(x):return F(float(np.nextafter(float(x),np.inf)))
def frob_integer(A,den):return F(sqceil(sum(int(x)**2 for x in A.flat)),den)
def midpoint(x):return float(sum(ends(x))/2)
def error_from_interval(x,v):return max(abs(z-F(v))for z in ends(x))

@lru_cache(None)
def exact_geometry():
 iv.dps=80;sh=shape();M,S=sh['M_exact'],sh['S_exact']
 mu=[number(M[0,0]+M[0,1]),number(M[0,0]-M[0,1])];omega=[iv.sqrt(number(S[0,0]+S[0,1])/mu[0]),iv.sqrt(number(S[0,0]-S[0,1])/mu[1])]
 c=[1/iv.sqrt(4*KAPPA*omega[i]*mu[i])for i in range(2)]
 ev=[-iv.sqrt(5+iv.sqrt(10)),-iv.sqrt(5-iv.sqrt(10)),iv.mpf(0),iv.sqrt(5-iv.sqrt(10)),iv.sqrt(5+iv.sqrt(10))]
 columns=[]
 for x in ev:
  v=[iv.mpf(1),x]
  for k in range(1,4):v.append((x*v[-1]-iv.sqrt(k)*v[-2])/iv.sqrt(k+1))
  norm=iv.sqrt(sum(z*z for z in v));columns.append([z/norm for z in v])
 W1=[[columns[j][i]for j in range(5)]for i in range(5)];W1f=np.array([[midpoint(x)for x in row]for row in W1]);err1=5*max(error_from_interval(W1[i][j],W1f[i,j])for i,j in product(range(5),repeat=2))
 W=np.kron(W1f,W1f);werr=2*err1+err1*err1+25*UNIT*(1+err1)**2
 values=[[c[0]*ev[i]+c[1]*ev[j],c[0]*ev[i]-c[1]*ev[j]]for i,j in product(range(5),repeat=2)]
 profiles=[[number(sh['F_exact'][x,0])*a+number(sh['F_exact'][x,1])*b for x in range(49)]for a,b in values]
 qmax=max(max(abs(a),abs(b))for v in values for a,b in map(ends,v));assert qmax<F(19,100)
 assert werr<F(1,10**12)
 assert ends(iv.sqrt(8)*iv.exp(scalar(F(19,50))))[1]<F(9,2)
 assert ends(2*scalar(native()[2][7].hi)*iv.exp(scalar(F(19,100))))[1]<8
 assert ends(4*(omega[0]+omega[1]))[1]<9
 assert all(v>=0 for v in sh['F_exact']) and all(sum(sh['F_exact'][i,j]for j in range(2))<=1 for i in range(49))
 states=list(product(range(5),repeat=2));hg=[i*omega[0]+j*omega[1]for i,j in states];Hgf=np.array([midpoint(x)for x in hg]);hgerr=max(error_from_interval(x,y)for x,y in zip(hg,Hgf))
 return {'omega':omega,'c':c,'values':values,'profiles':profiles,'W':W,'W_error':werr,'Hg':Hgf,'Hg_error':hgerr,'states':states,'qmax':qmax}

@lru_cache(None)
def field_certificate():
 g=exact_geometry();L=shape()['L'];roots=[];cert=[]
 for profile in g['profiles']:
  Kexact=[[int(L[i,j])*iv.exp(2*(profile[i]+profile[j]))if L[i,j]else iv.mpf(0)for j in range(49)]for i in range(49)]
  K=np.array([[midpoint(x)for x in row]for row in Kexact]);e,V=eigh(K);A=(V*np.sqrt(e))@V.T;AI=np.rint(A*ROOTDEN).astype(np.int64);KI=np.rint(K*ROOTDEN).astype(np.int64)
  delta=max(error_from_interval(Kexact[i][j],F(int(KI[i,j]),ROOTDEN))for i,j in product(range(49),repeat=2))
  X=AI.astype(object);res=frob_integer(X@X-KI.astype(object)*ROOTDEN,ROOTDEN**2)+49*delta
  positive=min(F(int(AI[i,i])-sum(abs(int(AI[i,j]))for j in range(49)if i!=j),ROOTDEN)for i in range(49));assert positive>0
  # lambda_min(K)>=8/(7+1)^2 exp(-4*.19)>1/25.
  assert ends(iv.exp(-scalar(F(19,25)))/8)[0]>F(1,25)
  err=res/(positive+F(1,5));roots.append(np.array(AI,float)/ROOTDEN)
  cert.append({'dyadic_positive_root':AI.tolist(),'positive_Gershgorin_lower':str(positive),'squared_root_residual_upper':str(res),'source_root_error_upper':str(err)})
 # The einsum uses25 products of three factors per entry, with a conservative200-operation bound.
 W=g['W'];A=np.stack(roots);full=np.einsum('ak,bk,kij->aibj',W,W,A,optimize=False).reshape(1225,1225)
 error=max(F(x['source_root_error_upper'])for x in cert)+10*g['W_error']+1225*gamma(200)*5*(1+g['W_error'])**2
 return {'H':full,'error':error,'root_certificates':cert}

@lru_cache(None)
def vacuum(left=1,right=0):
 g=exact_geometry();W=g['W'];ee=[interval(native()[2][i])for i in(left,right)];diag=[ee[0]*iv.exp(a)+ee[1]*iv.exp(b)for a,b in g['values']];v=np.array([midpoint(x)for x in diag]);derr=max(error_from_interval(x,y)for x,y in zip(diag,v))
 H=np.einsum('ak,bk,k->ab',W,W,v,optimize=False)+np.diag(g['Hg']);err=g['Hg_error']+derr+20*g['W_error']+25*gamma(200)*10*(1+g['W_error'])**2
 J=np.rint((H+H.T)/2*DEN).astype(np.int64);err+=F(25,2*DEN)+25*gamma(3)*20
 return {'H':np.array(J,float)/DEN,'integer_H':J,'error':err}

@lru_cache(None)
def block(left=1,right=0):
 v=vacuum(left,right);f=field_certificate();H=np.kron(v['H'],np.eye(49))+f['H'];J=np.rint((H+H.T)/2*DEN).astype(np.int64);err=v['error']+f['error']+F(1225,2*DEN)+1225*gamma(4)*32
 assert max(sum(abs(int(x))for x in row)for row in J)<32*DEN
 return {'H':np.array(J,float)/DEN,'integer_H':J,'error':err}

def indices(cap,wave=True):
 if cap not in range(5):raise ValueError('nested cap0 through4 only')
 cells=[i for i,(a,b)in enumerate(exact_geometry()['states'])if max(a,b)<=cap]
 return np.array([49*i+x for i in cells for x in range(49)]if wave else cells,dtype=int)
def compression(left=1,right=0,cap=4,wave=True):
 source=block(left,right)if wave else vacuum(left,right);ix=indices(cap,wave);return {'H':source['H'][np.ix_(ix,ix)],'integer_H':source['integer_H'][np.ix_(ix,ix)],'error':source['error'],'indices':ix}
def initial(cap=4,wave=True):
 v=np.zeros((cap+1)**2*(49 if wave else 1),complex)
 if wave:v[:49]=packet()
 else:v[0]=1
 return v
