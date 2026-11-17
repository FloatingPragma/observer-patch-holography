#!/usr/bin/env python3
"""Read-only reconstruction from coordinates and floor/ceiling, no producer imports."""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import argparse,copy,hashlib,json,sys,time
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent

def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):return a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1]
def dot(a,b):
    x=(0,0)
    for v,w in zip(a,b):x=add(x,mul(v,w))
    return x

def reconstruct():
    z=(0,0);one=(1,0);neg=(-1,0);phi=(0,1);minusphi=(0,-1)
    v=[(neg,phi,z),(one,phi,z),(neg,minusphi,z),(one,minusphi,z),(z,neg,phi),(z,one,phi),(z,neg,minusphi),(z,one,minusphi),(phi,z,neg),(phi,z,one),(minusphi,z,neg),(minusphi,z,one)]
    cell=[(a,b) for a in range(12) for b in range(a+1,12) if dot(v[a],v[b])==phi]
    edges=[];codes=[]
    for module in range(2):
        for shell in range(3):
            edges.extend((37*module+12*shell+a,37*module+12*shell+b) for a,b in cell)
            codes.extend(64*module+16*shell+p for p in range(12))
        for p in range(12):
            edges.extend(((37*module+p,37*module+24+p),(37*module+12+p,37*module+24+p),(37*module+24+p,37*module+36)))
        codes.append(64*module+48)
    edges.append((36,73));edges.sort();codes=np.array(codes,dtype=np.int64)
    L=np.zeros((74,74),dtype=np.int64)
    for a,b in edges:
        d=np.zeros(74,dtype=np.int64);d[a]=1;d[b]=-1;L+=np.outer(d,d)
    # Generate proper rotations from all oriented neighbor triangles, preserving their cyclic orientation via exact determinant.
    def cross(a,b):return (add(mul(a[1],b[2]),tuple(-x for x in mul(a[2],b[1]))),add(mul(a[2],b[0]),tuple(-x for x in mul(a[0],b[2]))),add(mul(a[0],b[1]),tuple(-x for x in mul(a[1],b[0]))))
    def positive(x):return float(x[0])+float(x[1])*(1+5**0.5)/2>0
    faces=[]
    for a in range(12):
        for b in range(12):
            for c in range(12):
                if a==b or b==c or a==c:continue
                if dot(v[a],v[b])==dot(v[b],v[c])==dot(v[c],v[a])==phi and positive(dot(v[a],cross(v[b],v[c]))):faces.append((a,b,c))
    base=faces[0];sig=lambda p,f:tuple(dot(v[p],v[k]) for k in f)
    perms=[]
    for face in faces:
        lookup={sig(p,face):p for p in range(12)}
        perms.append(tuple(lookup[sig(p,base)] for p in range(12)))
    perms=sorted(set(perms));assert len(perms)==60
    rotations=[]
    for g in perms:
        rotations.append(tuple(37*(p//37)+(12*((p%37)//12)+g[(p%37)%12] if p%37<36 else 36) for p in range(74)))
    return edges,codes,L,rotations

def check(r,a):
    edges,codes,L,rotations=reconstruct();assert np.array_equal(a['laplacian'],L)
    assert sorted(map(tuple,a['rotations']))==sorted(rotations)
    where=np.full(128,-2,dtype=np.int64);where[codes]=np.arange(74);valid=where>=0
    f=np.zeros(128,dtype=np.int64);f[:12]=1;f[16:28]=-1
    fidelity_num=0;escape_num=0;N=np.zeros((74,74),dtype=np.int64)
    for e,(u,v) in enumerate(edges):
        for coin in (0,1):
            output=np.arange(128);js=np.zeros(128,dtype=np.int64)
            for p in range(74):
                x=int(p==u);y=int(p==v);low=(x+y)//2;high=x+y-low
                xp,yp=(low,high) if coin==0 else (high,low)
                dest=p if p not in (u,v) else (u if xp else v)
                output[codes[p]]=codes[dest];js[codes[p]]=x-xp;N[dest,p]+=1
            assert np.array_equal(a['single_destinations'][e,coin],output)
            assert np.array_equal(a['single_currents'][e,coin],js)
            for j in (-1,0,1):
                branch=np.zeros(128,dtype=np.int64)
                for i in np.flatnonzero(valid & (js==j)):branch[output[i]]+=int(f[i])
                fidelity_num+=int(f@branch)**2;escape_num+=int(np.sum(branch[f==0]**2))
    assert np.array_equal(N,506*np.eye(74,dtype=np.int64)-L)
    assert r['unread_ray_fidelity']==str(F(fidelity_num,506*24**2))=='11867/12144'
    assert r['unread_escape_probability']==str(F(escape_num,506*24))=='1/506'
    labels=np.arange(16384);pp=where[labels//128];nn=where[labels%128]
    isvalid=(pp>=0)&(nn>=0)&(pp!=nn);isvalid[0]=True;pp[0]=nn[0]=-1
    assert int(isvalid.sum())==5403
    digest=hashlib.sha256();total=np.zeros(16384,dtype=np.int64);psi=a['pair_ray']
    expected_pair=np.zeros(16384,dtype=np.int64)
    g=np.zeros(128,dtype=np.int64);g[64:76]=1;g[80:92]=-1
    assert np.array_equal(psi,np.outer(f,g).reshape(-1))
    for u,v in edges:
        x=(pp==u).astype(np.int64)-(nn==u).astype(np.int64)
        y=(pp==v).astype(np.int64)-(nn==v).astype(np.int64)
        low=(x+y)//2;high=x+y-low
        for coin in (0,1):
            xp,yp=(low,high) if coin==0 else (high,low);j=x-xp
            po=pp.copy();no=nn.copy()
            touched=(pp==u)|(pp==v);po[touched]=-1
            po[touched & (xp==1)]=u;po[touched & (yp==1)]=v
            touched=(nn==u)|(nn==v);no[touched]=-1
            no[touched & (xp==-1)]=u;no[touched & (yp==-1)]=v
            out=labels.copy();wp=np.zeros(16384,dtype=np.int64);wn=wp.copy()
            wp[po>=0]=codes[po[po>=0]];wn[no>=0]=codes[no[no>=0]]
            out[isvalid]=(128*wp+wn)[isvalid];j[~isvalid]=0
            for current in (-1,0,1):
                ids=np.flatnonzero(isvalid & (j==current));assert len(set(map(int,out[ids])))==len(ids)
            digest.update(np.asarray(out,dtype='<u2').tobytes());digest.update(np.asarray(j,dtype=np.int8).tobytes())
            np.add.at(total,out,psi)
    assert digest.hexdigest()==r['pair_column_stream_sha256']
    assert np.array_equal(total,a['pair_transfer_sum']) and np.array_equal(total,504*psi)
    assert r['native_calls']==37462 and r['factorized_pair_columns']==8290304
    assert r['full_row_solved'] is False and r['physical_spin_generated'] is False and r['apparatus_compiler_executed'] is False
    return {'single_columns_recomputed':64768,'pair_columns_recomputed':8290304,'proper_rotations':60,'unread_ray_fidelity':'11867/12144','unread_escape':'1/506'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args();start=time.perf_counter()
    r=json.loads((HERE/'RECEIPT.json').read_text())
    with np.load(HERE/'RECORDS.npz',allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
    counts=check(r,a)
    result={'all_pass':True,'kind':'retrospective independently implemented same-group review','producer_imports':[],'counts':counts,'receipt_sha256':hashlib.sha256((HERE/'RECEIPT.json').read_bytes()).hexdigest(),'records_sha256':hashlib.sha256((HERE/'RECORDS.npz').read_bytes()).hexdigest(),'reviewer':'Codex root; producer/root identity disclosed'}
    if args.record:(HERE/'INDEPENDENT_CHECK.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    elif (HERE/'INDEPENDENT_CHECK.json').exists():assert result==json.loads((HERE/'INDEPENDENT_CHECK.json').read_text())
    print(json.dumps({'all_pass':True,'seconds':time.perf_counter()-start,**counts},indent=2))
if __name__=='__main__':main()
