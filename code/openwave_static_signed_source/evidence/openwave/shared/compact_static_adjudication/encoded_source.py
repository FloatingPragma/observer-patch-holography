"""One declared coordinate-memory lift of literal signed endpoint repair.

No export to a second signed object occurs. These work coordinates are the
source representation. Canonical positive amplitudes are an explicit N+ law.
"""
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE_ROOT=HERE.parent/'localized_rest_join'/'standalone'
sys.path.insert(0,str(SOURCE_ROOT/'reverse-engineering-reality'))
from model.observer import native

V=240
M=630
FAIR_DENOMINATOR=2*M
PAIR_LABELS=65536

def endpoint_table():
    table=np.empty((2,3,3,3),dtype=np.int8)
    for c in (0,1):
        for x in (-1,0,1):
            for y in (-1,0,1):
                output,record=native.native((x,y),(0,1),c)
                assert record[0]==c
                table[c,x+1,y+1]=(*output,record[1])
    return table

def pair_chart():
    labels=np.arange(PAIR_LABELS,dtype=np.int64)
    p,r=labels//256,labels%256
    valid=((p<V)&(r<V)&(p!=r))|(labels==0)
    p=np.where(labels==0,-1,p)
    r=np.where(labels==0,-1,r)
    return labels,p,r,valid

def pair_event(table,seam,coin):
    labels,p,r,valid=pair_chart()
    a,b=seam
    qa=(p==a).astype(np.int8)-(r==a).astype(np.int8)
    qb=(p==b).astype(np.int8)-(r==b).astype(np.int8)
    branch=table[coin,qa+1,qb+1]
    left,right,current=branch.T
    po=np.where((p==a)|(p==b),np.where(left==1,a,np.where(right==1,b,-1)),p)
    ro=np.where((r==a)|(r==b),np.where(left==-1,a,np.where(right==-1,b,-1)),r)
    vacuum=(po==-1)&(ro==-1)
    output=np.where(vacuum,0,256*po+ro)
    # Invalid input is handled by its separate identity flag, not a body branch.
    output=np.where(valid,output,labels)
    current=np.where(valid,current,0).astype(np.int8)
    return output.astype(np.uint16),current,valid

def single_event(table,seam,coin,sign=1):
    assert sign in (-1,1)
    labels=np.arange(256,dtype=np.int64);valid=labels<V;a,b=seam
    qa=sign*(labels==a).astype(np.int8);qb=sign*(labels==b).astype(np.int8)
    left,right,current=table[coin,qa+1,qb+1].T
    output=np.where((labels==a)|(labels==b),np.where(left==sign,a,b),labels)
    return output.astype(np.uint16),np.where(valid,current,0).astype(np.int8),valid

def apply_pair_column(table,seam,coin,label):
    # Independent scalar call on the represented source. It makes the
    # mathematical chart and actual native endpoint operation explicit.
    p,r=divmod(int(label),256)
    if label==0:q={}
    elif p<V and r<V and p!=r:q={p:1,r:-1}
    else:return int(label),0,'invalid'
    a,b=seam
    (left,right),(_,j)=native.native((q.get(a,0),q.get(b,0)),(0,1),coin)
    q[a]=left;q[b]=right;q={k:v for k,v in q.items() if v}
    if not q:return 0,j,'body'
    pp=next(k for k,v in q.items() if v==1)
    rr=next(k for k,v in q.items() if v==-1)
    return 256*pp+rr,j,'body'
