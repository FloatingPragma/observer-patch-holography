"""Fixed symmetric twin-shell signed source, all chart complements retained."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
sys.path.insert(0,str(ROOT/'reverse-engineering-reality'))
from model.observer.native import native
N=74;BITS=7;DIM=128;PAIR_DIM=DIM*DIM;M=253;DENOM=506
CELL=((0,1),(0,5),(0,7),(0,10),(0,11),(1,5),(1,7),(1,8),(1,9),(2,3),(2,4),(2,6),(2,10),(2,11),(3,4),(3,6),(3,8),(3,9),(4,5),(4,9),(4,11),(5,9),(5,11),(6,7),(6,8),(6,10),(7,8),(7,10),(8,9),(10,11))
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),(3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))
def code(p):
 m,r=divmod(int(p),37)
 return 64*m+(16*(r//12)+r%12 if r<36 else 48)
CODES=np.array([code(p) for p in range(N)],dtype=np.int64)
WHERE={int(w):p for p,w in enumerate(CODES)}
VALID=np.array([w in WHERE for w in range(DIM)])
def graph():
 edges=[]
 for m in range(2):
  o=37*m
  for shell in range(3):edges +=[(o+12*shell+a,o+12*shell+b) for a,b in CELL]
  for p in range(12):edges +=[(o+p,o+24+p),(o+12+p,o+24+p),(o+24+p,o+36)]
 edges.append((36,73));edges=tuple(sorted(edges));assert len(edges)==M and len(set(edges))==M
 return edges
SEAMS=graph()
def laplacian():
 L=np.zeros((N,N),dtype=np.int64)
 for a,b in SEAMS:L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
 return L
def ray(m=0):
 f=np.zeros(DIM,dtype=np.int64)
 for p in range(12):f[64*m+p]=1;f[64*m+16+p]=-1
 return f

def scalar(p,r,seam,coin,law='native'):
 a,b=seam;qa=int(p==a)-int(r==a);qb=int(p==b)-int(r==b)
 if law=='native':(x,y),(_,j)=native((qa,qb),(0,1),coin)
 elif law=='hold':x,y,j=qa,qb,0
 elif law=='swap':x,y,j=qb,qa,qa-qb
 else:raise ValueError(law)
 pp=p if p not in (a,b) else a if x==1 else b if y==1 else -1
 rr=r if r not in (a,b) else a if x==-1 else b if y==-1 else -1
 return pp,rr,j

def table():
 t=np.empty((2,3,3,3),dtype=np.int8)
 for c in (0,1):
  for x in (-1,0,1):
   for y in (-1,0,1):
    out,(_,j)=native((x,y),(0,1),c);t[c,x+1,y+1]=(*out,j)
 return t

def pair_event(t,seam,c):
 labels=np.arange(PAIR_DIM);wp,wr=labels//DIM,labels%DIM
 valid=(VALID[wp]&VALID[wr]&(wp!=wr))|(labels==0)
 p=np.array([WHERE.get(int(x),-2) for x in wp]);r=np.array([WHERE.get(int(x),-2) for x in wr]);p[0]=r[0]=-1
 a,b=seam;qa=(p==a).astype(np.int8)-(r==a).astype(np.int8);qb=(p==b).astype(np.int8)-(r==b).astype(np.int8)
 x,y,j=t[c,qa+1,qb+1].T
 pp=np.where((p==a)|(p==b),np.where(x==1,a,np.where(y==1,b,-1)),p)
 rr=np.where((r==a)|(r==b),np.where(x==-1,a,np.where(y==-1,b,-1)),r)
 woutp=np.where(pp>=0,CODES[np.maximum(pp,0)],0);woutr=np.where(rr>=0,CODES[np.maximum(rr,0)],0)
 out=DIM*woutp+woutr;out=np.where((pp==-1)&(rr==-1),0,out);out=np.where(valid,out,labels)
 return out.astype(np.uint16),np.where(valid,j,0).astype(np.int8),valid
