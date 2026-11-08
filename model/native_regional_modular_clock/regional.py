"""Exact regional restrictions of the declared fair transfer reference.

Only the 66-state source chart is
operationally attached here; growing references are mathematical states.
"""
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from math import comb
import numpy as np


def require(condition,message):
 if not condition:raise ArithmeticError(message)

def choose(n,k):
 if type(n)is not int or type(k)is not int or n<0:raise ValueError('exact binomial arguments required')
 return comb(n,k) if 0<=k<=n else 0

def local_tracial_eigenvalue(V,q,v,k):
 if any(type(x)is not int for x in(V,q,v,k))or not(0<=q<=V and 0<=k<=v<=32 and v<=V):
  raise ValueError('exact capacity and bounded local restriction required')
 if k>q or v-k>V-q:return F(0)
 result=F(1)
 for i in range(k):result*=q-i
 for i in range(v-k):result*=V-q-i
 for i in range(v):result/=V-i
 return result

def matrix_sector(v,q):
 if type(v)is not int or type(q)is not int or not 0<=q<=v<=32 or choose(v,q)>2048:
  raise ValueError('bounded exact sector: at most 32 ports and 2048 basis states')
 return choose(v,q)

def masks(v,q):
 matrix_sector(v,q)
 return tuple(sum(1<<i for i in a) for a in combinations(range(v),q))

def graph(v,edges):
 if type(v) is not int or v<2 or type(edges) is not tuple or not edges:
  raise ValueError('finite graph with at least one literal edge required')
 if any(type(e) is not tuple or len(e)!=2 or any(type(i) is not int for i in e) or not 0<=e[0]<e[1]<v for e in edges) or len(set(edges))!=len(edges):
  raise ValueError('distinct ordered integer edges required')
 return edges

def laplacian(v,edges,q):
 graph(v,edges)
 if type(q) is not int or not 0<=q<=v:raise ValueError('exact sector required')
 basis=masks(v,q);where={m:i for i,m in enumerate(basis)};out=np.zeros((len(basis),len(basis)),dtype=object)
 for j,m in enumerate(basis):
  for a,b in edges:
   if ((m>>a)^(m>>b))&1:
    out[j,j]+=1;out[where[m^(1<<a)^(1<<b)],j]-=1
 return basis,out

def reference(v,edges,q):
 basis,G=laplacian(v,edges,q);M=len(edges);d=len(basis)
 P=np.array([[F(int(i==j))-F(G[i,j],2*M) for j in range(d)]for i in range(d)],dtype=object)
 Z=F(d)*(1-F(q*(v-q),v*(v-1)))
 require(sum(P[i,i]for i in range(d))==Z,'the exact transfer trace disagrees')
 return basis,P/Z,Z

def region(v,A):
 if type(v)is not int or type(A)is not tuple or not A or any(type(i)is not int or not 0<=i<v for i in A)or tuple(sorted(set(A)))!=A or v-len(A)<2:
  raise ValueError('sorted nonempty region with at least two exterior vertices required')
 return A

def split(v,edges,A):
 graph(v,edges);region(v,A)
 aset=set(A);ea=tuple(e for e in edges if set(e)<=aset);eb=tuple(e for e in edges if not(set(e)&aset));ec=tuple(e for e in edges if len(set(e)&aset)==1)
 return ea,eb,ec

def partial_trace(v,basis,rho,A):
 if type(v)is not int or not 2<=v<=32 or type(basis)is not tuple or not basis or any(type(x)is not int or not 0<=x<(1<<v)for x in basis):
  raise ValueError('bounded exact source basis required')
 q=basis[0].bit_count()
 if basis!=masks(v,q):raise ValueError('the complete ordered source sector is required')
 region(v,A)
 if not isinstance(rho,np.ndarray)or rho.shape!=(len(basis),len(basis))or any(type(x)not in(int,F)for x in rho.flat):
  raise ValueError('a complete bounded exact rational source matrix is required')
 amask=sum(1<<i for i in A);blocks={}
 for k in range(len(A)+1):
  local=tuple(sum(1<<i for i in a)for a in combinations(A,k));where={m:i for i,m in enumerate(local)}
  R=np.zeros((len(local),len(local)),dtype=object)
  for i,x in enumerate(basis):
   if (x&amask).bit_count()!=k:continue
   for j,y in enumerate(basis):
    if x&~amask==y&~amask and (y&amask).bit_count()==k:R[where[x&amask],where[y&amask]]+=rho[i,j]
  blocks[k]=(local,R)
 return blocks

def regional_block(v,edges,q,A,k):
 if type(q)is not int or type(v)is not int or not 0<=q<=v or type(k)is not int or type(A)is not tuple or not 0<=k<=len(A):
  raise ValueError('exact supported source and regional sector indices required')
 matrix_sector(len(A),k)
 ea,eb,ec=split(v,edges,A);M=len(edges);b=v-len(A);l=q-k
 local=tuple(sum(1<<i for i in a)for a in combinations(A,k));d=len(local);where={m:i for i,m in enumerate(local)}
 if not 0<=l<=b:return local,np.zeros((d,d),dtype=object),None
 r=F(l,b);cB=F(2*len(eb)*l*(b-l),b*(b-1));beta=cB+len(ec)*r;a=1-beta/(2*M)
 D=np.zeros((d,d),dtype=object)
 for j,m in enumerate(local):
  for x,y in ea:
   if ((m>>x)^(m>>y))&1:D[j,j]+=1;D[where[m^(1<<x)^(1<<y)],j]-=1
  for e in ec:
   i=e[0] if e[0] in A else e[1]
   D[j,j]+=(1-2*r)*((m>>i)&1)
 trace_factor=1-F(q*(v-q),v*(v-1));alpha=local_tracial_eigenvalue(v,q,len(A),k)/trace_factor
 R=np.array([[alpha*(a*int(i==j)-F(D[i,j],2*M)) for j in range(d)]for i in range(d)],dtype=object)
 return local,R,{'D':D,'a':a,'r':r,'exterior_scalar':cB,'normalization_trace_factor':trace_factor,'exterior_number':l}

def modulus_bound(V,M,q,v,mA,mC,k,p):
 args=(V,M,q,v,mA,mC,k)
 if any(type(x)is not int for x in args) or not(0<=q<=V and 0<v<V-1 and 0<=k<=v and 0<=q-k<=V-v and M>0 and mA>=0 and mC>=0 and mA+mC<=M) or type(p)is not F or not 0<=p<=1:
  raise ValueError('exact supported finite graph counts and density required')
 b=V-v;l=q-k;r=F(l,b);mB=M-mA-mC
 beta=F(2*mB*l*(b-l),b*(b-1))+mC*r;a=1-beta/(2*M);ap=1-p*(1-p);L=2*mA+mC;z=F(L,2*M)/a
 if z>=1:raise ValueError('the stated logarithm-series gate is not satisfied')
 log_error=F(L*L,4*M)/(a*a*(1-z));density_error=2*mC*abs(r-p)/a+L*abs(a-ap)/(a*ap)
 return {'a':a,'target_a':ap,'r':r,'z':z,'logarithm_remainder':log_error,'density_remainder':density_error,'generator_error':log_error+density_error,'commutator_error_per_time_and_observable_norm':2*(log_error+density_error)}

def canonical_variance(v,edges,q):
 graph(v,edges)
 if v<4 or type(q)is not int or not 0<=q<=v:raise ValueError('V>=4 and exact occupation required')
 M=len(edges);degrees=Counter(i for edge in edges for i in edge);adj=sum(d*(d-1)for d in degrees.values());dis=M*(M-1)-adj
 u=F(q*(v-q),v*(v-1));w=F(q*(q-1)*(v-q)*(v-q-1),v*(v-1)*(v-2)*(v-3))
 mean=2*M*u;second=4*M*u+adj*u+4*dis*w;variance=second-mean*mean
 return {'mean':mean,'second':second,'variance':variance,'rho_variance_upper':variance/(1-u),'rho_mean':mean-variance/(2*M*(1-u))}

def cross_sector_hs_squared(v,M,q,crossings):
 if any(type(x)is not int for x in(v,M,q,crossings))or not(2<=v<=256 and M>0 and 0<=q<=v and 0<=crossings<=M):
  raise ValueError('bounded exact commutator count required')
 Z=F(choose(v,q))*(1-F(q*(v-q),v*(v-1)))
 return F(2*crossings*choose(v-2,q-1),1)/(2*M*Z)**2

def overlap_witness():
 # Relabel native triangle(0,1,5) as(0,1,2), in its one-occupation block.
 _,small=laplacian(3,((0,1),),1);_,large=laplacian(3,((0,1),(0,2),(1,2)),1)
 n=np.diag([1,0,0]).astype(object);C=(large-small)@n-n@(large-small)
 # C* C is the projection on occupations at0 and2. Hence ||C||=1 exactly.
 require(np.array_equal(C.T@C,np.diag([1,0,1])),'the exact overlap witness changed')
 return {'raw_exchange_overlap_commutator_norm':F(1),'half_density_rescaled_modular_overlap_commutator_norm':F(4,3)}


def coherent_reference_mixture(v,edges,q):
 """Exact positive mixture and retained reject law, not a compiled service.

 Uniform canonical basis preparation plus one local coherent exchange mixer
 are explicit implementation premises. The output is the accepted marginal.
 """
 basis,G=laplacian(v,edges,q);M=len(edges);d=len(basis)
 numerator=np.zeros((d,d),dtype=object);accepted=0;rejected=0
 index={x:i for i,x in enumerate(basis)}
 for i,x in enumerate(basis):
  for a,b in edges:
   x_a=(x>>a)&1;x_b=(x>>b)&1
   if x_a==x_b:
    numerator[i,i]+=1;accepted+=1
   elif x_a==1:
    y=x^(1<<a)^(1<<b);j=index[y]
    numerator[i,i]+=F(1,2);numerator[j,j]+=F(1,2)
    numerator[i,j]+=F(1,2);numerator[j,i]+=F(1,2);accepted+=1
   else:rejected+=1
 require(accepted+rejected==d*M,'the sampler lost complete outcomes')
 _,rho,Z=reference(v,edges,q)
 require(F(accepted,d*M)==1-F(q*(v-q),v*(v-1)),'the exact acceptance probability disagrees')
 require(np.array_equal(numerator/F(accepted),rho),'the sampler does not prepare its stated conditional reference')
 return {'accepted_probability':F(accepted,d*M),'rejected_probability':F(rejected,d*M),'accepted_state':numerator/F(accepted),'coherent_mixer_is_supplied':True}
