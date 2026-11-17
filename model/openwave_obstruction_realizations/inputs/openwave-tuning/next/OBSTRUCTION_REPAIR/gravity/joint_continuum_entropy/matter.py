import itertools
import numpy as np
from kernels import cell,s

def prepare_matter(source,width):
    L=source['L'];h=source['h'];R=source['R'];mass=source['mass'];twist=source['twist'];aliases=source['aliases']
    axis=2*np.pi*(np.fft.fftfreq(L)*L+twist)/(L*h)
    p0=np.stack(np.meshgrid(*([axis]*3),indexing='ij'),-1)
    G={name:np.zeros((4,4),complex) for name in ['particle','hole']}
    energy={name:np.zeros((4,4),complex) for name in G}
    fields={name:np.zeros((L,L,L,4,4),complex) for name in G}
    for alias in itertools.product(range(-aliases,aliases+1),repeat=3):
        p=p0+2*np.pi*np.array(alias)/h
        _,P,E=s.symbol(p,R,mass);w,_=cell.window(p,h)
        f=np.exp(-width*width*np.sum(p*p,axis=-1)/4)
        for name,proj in [('particle',s.I4-P),('hole',P)]:
            G[name]+=np.einsum('xyz,xyzij->ij',f*f,proj)
            energy[name]+=np.einsum('xyz,xyzij->ij',f*f*E,proj)
            fields[name]+=(w*f)[...,None,None]*proj
    records={};total_energy=0.;average=np.zeros_like(fields['particle'])
    for name in G:
        ev,U=np.linalg.eigh(G[name]);inv=(U/np.sqrt(ev))@U.conj().T
        modes=fields[name]@inv
        total_energy+=float(np.trace(inv@energy[name]@inv).real)
        average+=(1 if name=='particle' else -1)*(modes@modes.conj().swapaxes(-1,-2))
        records[name]=np.fft.ifftn(modes,axes=(0,1,2))*np.sqrt(L**3)
        records[name+'_global_gram_error']=float(np.linalg.norm(inv@G[name]@inv-s.I4))
        records[name+'_cell_captured_trace']=float(np.sum(abs(modes)**2))
    records['average_kernel']=np.fft.ifftn(average,axes=(0,1,2))
    records['width']=width;records['actual_energy']=total_energy
    return records

def response(region,source,matter,centre):
    coordinates=region['cells']+np.asarray(centre)
    ix=coordinates%source['L'];phase=np.exp(2j*np.pi*source['twist']*coordinates.sum(-1)/source['L'])
    local={name:(matter[name][ix[:,0],ix[:,1],ix[:,2]]*phase[:,None,None]).reshape(4*len(coordinates),4) for name in ['particle','hole']}
    first=0.;actual=0.;other_irrep_weight=0.
    for name,d in [('G1',2),('G2',2),('H',4)]:
        shape=(region['sizes'][name],4)
        parts={key:[np.zeros(shape,complex) for _ in range(d)] for key in local}
        for o in region['orbits']:
            for j,B in enumerate(o['sectors'][name]['bases']):
                for key,X in local.items():parts[key][j][o['columns'][name]]=B.conj().T@X[o['rows']]
        for j in range(d):
            Xp=parts['particle'][j];Xh=parts['hole'][j]
            first+=float(np.trace(Xp.conj().T@region['K'][name]@Xp-Xh.conj().T@region['K'][name]@Xh).real)
            if name!='G1':other_irrep_weight+=float(np.linalg.norm(Xp)**2+np.linalg.norm(Xh)**2)
        if np.all(np.asarray(centre)==0):
            Xp=parts['particle'][0];Xh=parts['hole'][0]
            if name=='G1':
                actual+=d*(s.entropy(region['C'][name]+Xp@Xp.conj().T-Xh@Xh.conj().T)-s.entropy(region['C'][name]))
            elif max(np.linalg.norm(Xp),np.linalg.norm(Xh))>1e-10:
                raise ValueError('Centred radial seed must transform in G1; cannot cancel a changed sector.')
    result=dict(centre=list(centre),first_derivative=first,other_irrep_mode_weight=other_irrep_weight)
    if np.all(np.asarray(centre)==0):
        result.update(actual_entropy_change=actual,relative_first_law_remainder=(first-actual)/first)
    return result
