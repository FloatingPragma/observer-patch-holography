import itertools
import numpy as np
from scipy.special import xlogy
from kernels import cell,s,sibling
group=sibling('all_cubic_group_for_entropy','full_cubic_sectors/group.py')

def coefficients_at(source,key,d):
    ix=d%source['L']
    c=source[key][ix[...,0],ix[...,1],ix[...,2]].copy()
    if key=='shape':
        for j in range(3):c-=source['h']*d[...,j,None]*source['gradients'][j][ix[...,0],ix[...,1],ix[...,2]]
    return c*np.exp(2j*np.pi*source['twist']*d.sum(-1)/source['L'])[...,None]

def block_at(source,key,left,right):
    d=left[:,None,:]-right[None,:,:]
    c=coefficients_at(source,key,d)
    block=np.einsum('...a,aij->...ij',c,cell.GAMMA)
    if key=='cov':block+=(np.all(d==0,axis=-1))[...,None,None]*cell.PB
    return block.transpose(0,2,1,3).reshape(4*len(left),4*len(right))

def build(source,radius):
    cells=np.array([c for c in itertools.product(range(-radius,radius+1),repeat=3) if sum(x*x for x in c)<=radius*radius])
    orbits=group.orbit_bases(cells)
    sizes={name:0 for name in ['G1','G2','H']}
    for o in orbits:
        o['columns']={}
        for name in sizes:
            n=o['sectors'][name]['multiplicity'];o['columns'][name]=np.arange(sizes[name],sizes[name]+n);sizes[name]+=n
    C={name:np.zeros((n,n),complex) for name,n in sizes.items()}
    for i,oi in enumerate(orbits):
        for j in range(i,len(orbits)):
            oj=orbits[j];full=block_at(source,'cov',oi['cells'],oj['cells'])
            for name in sizes:
                bi=oi['sectors'][name]['bases'][0];bj=oj['sectors'][name]['bases'][0]
                block=bi.conj().T@full@bj
                C[name][np.ix_(oi['columns'][name],oj['columns'][name])]=block
                if i!=j:C[name][np.ix_(oj['columns'][name],oi['columns'][name])]=block.conj().T
    K={};S=0.;floors={};eigenranges={};eigensystems={}
    for name,d in [('G1',2),('G2',2),('H',4)]:
        ev,U=np.linalg.eigh(C[name]);e=np.clip(ev,1e-13,1-1e-13)
        K[name]=(U*np.log((1-e)/e))@U.conj().T
        eigensystems[name]=(ev,U)
        physical=np.clip(ev,0,1)
        S+=d*float(np.sum(-xlogy(physical,physical)-xlogy(1-physical,1-physical)));floors[name]=int(np.sum((ev<1e-13)|(ev>1-1e-13)))
        eigenranges[name]=[float(ev[0]),float(ev[-1])]
    return dict(cells=cells,orbits=orbits,C=C,K=K,radius=radius,entropy=S,sizes=sizes,floors=floors,eigenranges=eigenranges,eigensystems=eigensystems)

def aggregates(region,source,shifts):
    r=region['radius'];side=4*r+1;offset=2*r
    out=[np.zeros((side,side,side,4,4),complex) for _ in shifts]
    orbits=region['orbits'];K=region['K']
    for i,oi in enumerate(orbits):
        for j,oj in enumerate(orbits):
            full=np.zeros((4*len(oi['cells']),4*len(oj['cells'])),complex)
            for name in K:
                block=K[name][np.ix_(oi['columns'][name],oj['columns'][name])]
                for bi,bj in zip(oi['sectors'][name]['bases'],oj['sectors'][name]['bases']):full+=bi@block@bj.conj().T
            # Trace K_ij dC_ji; d=x_j-x_i and phase belongs to x_j.
            d=oj['cells'][None,:,:]-oi['cells'][:,None,:]
            indices=(d+offset).reshape(-1,3)
            matrix=full.reshape(len(oi['cells']),4,len(oj['cells']),4).transpose(0,2,1,3)
            for acc,shift in zip(out,shifts):
                phase=np.exp(2j*np.pi*(oj['cells']@np.asarray(shift))/source['L'])
                values=(matrix*phase[None,:,None,None]).reshape(-1,4,4)
                np.add.at(acc,(indices[:,0],indices[:,1],indices[:,2]),values)
    d=np.stack(np.meshgrid(*([np.arange(-offset,offset+1)]*3),indexing='ij'),-1)
    return d,out

def trace_coefficients(acc,source,key,d):
    kernel=np.einsum('...a,aij->...ij',coefficients_at(source,key,d),cell.GAMMA)
    return complex(np.sum(acc*kernel.swapaxes(-1,-2)))

def trace_response(acc,kernel,source,d):
    ix=d%source['L'];phase=np.exp(2j*np.pi*source['twist']*d.sum(-1)/source['L'])
    values=kernel[ix[...,0],ix[...,1],ix[...,2]]*phase[...,None,None]
    return complex(np.sum(acc*values.swapaxes(-1,-2)))
