"""Actual canonical tilted-plane entropy, with all source bands retained.

Normal direction is discretized by a complete finite Fourier CAR system.
Transverse momentum is integrated in a disk with explicit cutoff. No
negative entropy weights, transported Dirac algebra, or fitted area input.
"""
import argparse,json,time
import numpy as np
from scipy.special import xlogy
from source import tilted_covariance

parser=argparse.ArgumentParser()
parser.add_argument('--count',type=int,default=128)
parser.add_argument('--h',type=float,default=.25)
parser.add_argument('--radial',type=int,default=24)
parser.add_argument('--angular',type=int,default=12)
parser.add_argument('--cutoff',type=float,default=8.)
parser.add_argument('--speeds',type=float,nargs='+',default=[0.,.3,.6])
parser.add_argument('--orientations',nargs='+',default=['parallel','perpendicular'])
parser.add_argument('--output',default='OBSERVER_AREA.json')
args=parser.parse_args()
start=time.process_time()
count=args.count
normal=2*np.pi*np.fft.fftfreq(count,d=args.h)
region=np.arange(count//2)
delta=(region[:,None]-region[None,:])%count
nodes,weights=np.polynomial.legendre.leggauss(args.radial)
radial=(nodes+1)*args.cutoff/2
weights=weights*args.cutoff/2

def region_entropy(q):
    C,G,_=tilted_covariance(q,u)
    ker=np.fft.ifft(C,axis=0)
    CA=ker[delta].transpose(0,2,1,3).reshape(2*count,2*count)
    e=np.clip(np.linalg.eigvalsh((CA+CA.conj().T)/2),0,1)
    ent=float(np.sum(-xlogy(e,e)-xlogy(1-e,1-e)))
    return ent,float(np.max(np.abs(C@C-C)))

rows=[]
for orientation in args.orientations:
    for u in args.speeds:
        total=0.;maxprojector=0.;samples=[]
        for r,weight in zip(radial,weights):
            angles=[0.] if orientation=='parallel' else 2*np.pi*(np.arange(args.angular)+.5)/args.angular
            entropies=[]
            for angle in angles:
                q=np.zeros((count,3))
                if orientation=='parallel':
                    q[:,0]=normal;q[:,1]=r
                else:
                    q[:,0]=r*np.cos(angle);q[:,1]=r*np.sin(angle);q[:,2]=normal
                ent,proj=region_entropy(q)
                maxprojector=max(maxprojector,proj)
                entropies.append(ent)
            average=float(np.mean(entropies))
            # integral d^2q/(2pi)^2 times one-half (two interfaces).
            total+=weight*r*average/(4*np.pi)
            samples.append(dict(transverse_radius=float(r),mean_entropy=average,angular_range=float(np.ptp(entropies))))
        row=dict(orientation=orientation,speed=u,eta_per_proper_face=total,max_projector_entry_defect=maxprojector,samples=samples)
        rows.append(row)
        print(json.dumps({k:v for k,v in row.items() if k!='samples'}),flush=True)
        if time.process_time()-start>270:
            raise RuntimeError('Declared 300 CPU second resource ceiling approached; partial rows are in console, no completed claim.')
out=dict(status='ACTUAL_CANONICAL_TILTED_PLANE_ENTROPY',configuration=vars(args),physical_normal_box=count*args.h,regional_width=count*args.h/2,rows=rows,cpu_seconds=time.process_time()-start,source='Continuous positive flattened quartic projector, R=1',algebra='Full Gram-whitened sampled source fields; original finite-entropy local modes at speed zero.',full_covariant_region_net_proved=False)
with open(args.output,'w') as f:json.dump(out,f,indent=2)
