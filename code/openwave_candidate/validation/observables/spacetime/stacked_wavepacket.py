"""Declared stacked-ring extension: exact positive-band wavepacket, no derived gravity."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/openwave-matplotlib')
import sys
import numpy as np

BUNDLE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BUNDLE/'candidate'))
from oph_openwave.geometry import CARRIER_SEAMS,RING_PORTS
A=np.array([[math.sqrt(3),math.sqrt(3)/2,0],[0,1.5,0],[0,0,1.]])
D=np.array([math.sqrt(3)/2,.5,0.])
CXY2=15/352


def validate_mass(mass):
    if not np.isscalar(mass) or not np.isfinite(mass) or mass <= 0:
        raise ValueError('mass must be finite and strictly positive for this band formula')


def validate_grid(n):
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 4:
        raise ValueError('n must be an integer at least 4')


def bloch_incidence(q,vertical=True):
    q=np.asarray(q,float)
    if q.shape != (3,) or not np.all(np.isfinite(q)):
        raise ValueError('q must contain three finite momentum components')
    b=np.zeros((24,87 if vertical else 63),complex)
    for sub in (0,1):
        for e,(a,c) in enumerate(CARRIER_SEAMS):
            b[12*sub+a,30*sub+e]=1;b[12*sub+c,30*sub+e]=-1
    for e,(p,phase) in enumerate(zip(RING_PORTS,(0.,q[0],q[1]))):
        b[p,60+e]=1;b[12+p,60+e]=-np.exp(1j*phase)
    if vertical:b[:,63:]=(1-np.exp(-1j*q[2]))*np.eye(24)
    return b


def positive_mode(q,mass=.23,vertical=True):
    validate_mass(mass)
    b=bloch_incidence(q,vertical);w,u=np.linalg.eigh(b@b.conj().T);v=u[:,0]
    pivot=0 if abs(v[0])>1e-10 else int(np.argmax(abs(v)))
    v=v*np.exp(-1j*np.angle(v[pivot]));lam=max(0.,float(w[0]));energy=math.sqrt(mass*mass+lam)
    mode=np.r_[math.sqrt((energy+mass)/(2*energy))*v,b.conj().T@v/math.sqrt(2*energy*(energy+mass))]
    return energy,mode


def internal_offsets(vertical=True):
    positions=[np.zeros(3)]*12+[D]*12+[np.zeros(3)]*30+[D]*30
    positions += [(D-shift)/2 for shift in (np.zeros(3),A[:,0],A[:,1])]
    if vertical:positions += [x+np.array([0.,0.,.5]) for x in ([np.zeros(3)]*12+[D]*12)]
    return np.asarray(positions)


def band_grid(n=32,mass=.23,vertical=True):
    validate_grid(n);validate_mass(mass)
    q=2*np.pi*np.fft.fftfreq(n);count=111 if vertical else 87
    energy=np.empty((n,n,n));spinor=np.empty((n,n,n,count),complex);velocity=np.empty((n,n,n,3))
    minimum_pivot=1.;maximum_norm_error=0.
    for i,x in enumerate(q):
        for j,y in enumerate(q):
            b=bloch_incidence([x,y,0],False);values,vectors=np.linalg.eigh(b@b.conj().T)
            u=vectors[:,0];minimum_pivot=min(minimum_pivot,float(abs(u[0])))
            pivot=0 if abs(u[0])>1e-10 else int(np.argmax(abs(u)))
            u=u*np.exp(-1j*np.angle(u[pivot]));lambda_xy=max(0.,float(values[0]));beta=b.conj().T@u
            gradient=[]
            for edge,phase in ((1,x),(2,y)):
                dbeta=1j*np.exp(-1j*phase)*u[12+RING_PORTS[edge]]
                gradient.append(2*np.real(np.conj(beta[60+edge])*dbeta))
            lam=lambda_xy+(2*(1-np.cos(q)) if vertical else np.zeros(n))
            e=np.sqrt(mass*mass+lam);energy[i,j]=e
            spinor[i,j,:,:24]=np.sqrt((e+mass)/(2*e))[:,None]*u
            spinor[i,j,:,24:87]=beta/np.sqrt(2*e*(e+mass))[:,None]
            if vertical:spinor[i,j,:,87:]=((1-np.exp(1j*q))[:,None]*u)/np.sqrt(2*e*(e+mass))[:,None]
            velocity[i,j,:,0]=gradient[0]/(2*e);velocity[i,j,:,1]=gradient[1]/(2*e)
            velocity[i,j,:,2]=np.sin(q)/e if vertical else 0.
    maximum_norm_error=float(np.max(abs(np.sum(abs(spinor)**2,axis=-1)-1)))
    return {'q':q,'energy':energy,'spinor':spinor,'velocity_cell':velocity,
            'minimum_phase_pivot':minimum_pivot,'spinor_norm_max_error':maximum_norm_error}


def simulate(n=32,mass=.23,sigma=.20,q0=(.75,.45,.10),times=None,return_spatial_groups=False):
    validate_grid(n);validate_mass(mass)
    if not np.isscalar(sigma) or not np.isfinite(sigma) or sigma <= 0:
        raise ValueError('sigma must be finite and strictly positive')
    q0=np.asarray(q0,float)
    if q0.shape != (3,) or not np.all(np.isfinite(q0)):
        raise ValueError('q0 must contain three finite momentum components')
    if times is None:times=np.linspace(0,4,9)
    times=np.asarray(times,float)
    if times.ndim != 1 or not len(times) or not np.all(np.isfinite(times)):
        raise ValueError('times must be a nonempty one-dimensional finite array')
    grid=band_grid(n,mass);q=grid['q'];qq=np.stack(np.meshgrid(q,q,q,indexing='ij'),axis=-1)
    delta=(qq-np.asarray(q0)+np.pi)%(2*np.pi)-np.pi
    amplitude=np.exp(-np.sum(delta**2,axis=-1)/(4*sigma*sigma)).astype(complex)
    amplitude*=np.exp(-1j*np.sum(qq*(n//2),axis=-1));amplitude/=np.linalg.norm(amplitude.ravel())
    weight=abs(amplitude)**2;mean_velocity=A@np.sum(weight[...,None]*grid['velocity_cell'],axis=(0,1,2))
    center_energy,_=positive_mode(q0,mass)
    # Independent centered differences give the central ray; it need not equal
    # the weighted band velocity of a packet with finite momentum width.
    dq=1e-5;center_gradient=[]
    for axis in range(3):
        offset=np.eye(3)[axis]*dq
        center_gradient.append((positive_mode(np.asarray(q0)+offset,mass)[0]-positive_mode(np.asarray(q0)-offset,mass)[0])/(2*dq))
    central_velocity=A@center_gradient
    indices=np.stack(np.meshgrid(np.arange(n)-n//2,np.arange(n)-n//2,np.arange(n)-n//2,indexing='ij'),axis=-1)
    coordinates=indices@A.T;offsets=internal_offsets();edge=np.any(abs(indices)>=n//2-2,axis=-1)
    base=amplitude[...,None]*grid['spinor'];samples=[];distributions=[];spatial_groups=[]
    unique_offsets,component_bins=np.unique(offsets,axis=0,return_inverse=True)
    for t in times:
        psi=np.fft.ifftn(base*np.exp(-1j*t*grid['energy'])[...,None],axes=(0,1,2),norm='ortho')
        density=abs(psi)**2;cell=density.sum(axis=-1);components=density.sum(axis=(0,1,2))
        mean=np.sum(cell[...,None]*coordinates,axis=(0,1,2))+components@offsets
        variance=np.zeros(3)
        for axis in range(3):
            second=np.sum(cell*coordinates[...,axis]**2)+2*np.sum(density*coordinates[...,axis,None]*offsets[:,axis])+components@(offsets[:,axis]**2)
            variance[axis]=max(0,float(second-mean[axis]**2))
        samples.append({'time':float(t),'normalization':float(cell.sum()),'centroid':mean.tolist(),'position_std':np.sqrt(variance).tolist(),'edge_shell_probability':float(cell[edge].sum())})
        distributions.append(cell)
        if return_spatial_groups:spatial_groups.append(np.stack([density[...,component_bins==i].sum(axis=-1) for i in range(len(unique_offsets))],axis=-1))
    initial=np.asarray(samples[0]['centroid']);errors=[]
    for s in samples:
        ray=initial+(s['time']-times[0])*mean_velocity;s['weighted_group_velocity_ray']=ray.tolist();errors.append(float(np.linalg.norm(np.asarray(s['centroid'])-ray)))
    speed_norm=float(mean_velocity[0]**2/CXY2+mean_velocity[1]**2/CXY2+mean_velocity[2]**2)
    exact_points=[]
    for point in (q0,(.2,-.13,.31)):
        b=bloch_incidence(point);e,v=positive_mode(point,mass);h=np.block([[mass*np.eye(24),b],[b.conj().T,-mass*np.eye(87)]])
        exact_points.append({'q':list(point),'positive_energy':e,'full_H_eigen_residual':float(np.linalg.norm(h@v-e*v))})
    package=BUNDLE/'candidate/oph_openwave'
    result={'schema':'oph.explicit_stacked_ring.wavepacket.v1','scope':'Explicit 3D graph extension, uniform-record preparation, positive acoustic band and zero commit rate. Complete position distributions; centroid and reference ray, not an observed individual trajectory or derived gravity.',
      'supplied':{'cells_per_axis':n,'mass':mass,'momentum_center_cell_covectors':list(q0),'momentum_sigma':sigma,'native_times':list(map(float,times)),'vertical_seam_weight':1,'vertical_spacing':1,'commit_rate':0,'record_preparation':'uniform modular Gauss fibers, modulus6; full-sector reduction theorem applies to connected augmented graph','band_preparation':'positive acoustic eigenvector at each momentum, with a declared finite-grid component phase convention; global phase continuity is not asserted','position_apparatus':'coarse carrier locations plus whole-seam midpoints; internal ports at their carrier center','measurement_protocol':'separate identically prepared ensemble at each time; no sequential trajectory measurement','mass_or_length_empirical_anchor':None,'metric':'reference flat acoustic metric diag(-1,352/15,352/15,1), supplied interpretation of long-wavelength symbol','gravity':'none; no prescribed curvature or inferred G'},
      'geometry':{'ports_per_cell':24,'horizontal_seams':63,'vertical_seams':24,'position_components_per_cell':111,'spatial_basis_columns':A.tolist(),'internal_offsets':offsets.tolist(),'port_laplacian':'L_ring(q1,q2)+2(1-cos(q3))*I24','horizontal_acoustic_speed_squared':CXY2,'vertical_acoustic_speed_squared':1.,'speed_squared_anisotropy':1/CXY2,'all_bands':'Full generator has24 positive/negative dispersive pairs generically and63 additional -m transit modes. Only preparation occupies positive acoustic branch.'},
      'checks':{'spinor_norm_max_error':grid['spinor_norm_max_error'],'full_H_residuals':exact_points,'probability_norm_max_error':max(abs(s['normalization']-1) for s in samples),'minimum_phase_pivot':grid['minimum_phase_pivot']},
      'observables':{'samples':samples,'weighted_group_velocity':mean_velocity.tolist(),'central_group_velocity':central_velocity.tolist(),'central_energy':center_energy,'reference_ray_anchor_time':float(times[0]),'max_centroid_ray_distance':max(errors),'maximum_edge_shell_probability':max(s['edge_shell_probability'] for s in samples),'reference_metric_speed_norm':speed_norm,'reference_ray_proper_parameter_rate':math.sqrt(1-speed_norm) if speed_norm<1 else None,'geodesic_scope':'The reference metric is flat, so its connection is zero and a constant-velocity ray is geodesic. Finite-lattice wavepacket centroid comparison is reported, not asserted exact. No modular-to-proper-time conversion.'},
      'rejecting_controls':{'remove_vertical_seams':'bloch_incidence(q,vertical=False) is exactly independent of q3; three-coordinate drawing alone supplies no vertical propagation','isotropic_unit_speed_rejected':abs(1-CXY2)>.9,'fixed_band_is_not_full_theory':'other positive, negative and flat bands are present; local commits can mix bands'},
      'stored_spatial_apparatus':{'coarse_internal_offsets':unique_offsets.tolist(),'component_to_offset_bin':component_bins.tolist(),'complete_distribution':'position_distributions.npz: spatial_probability[time,cell_i,cell_j,cell_k,offset_bin]','scope':'Complete named coarse geometric-position apparatus; colocated internal port/seam labels are aggregated. Cell marginal is retained separately for plotting.'},
      'official_rows_promoted':0,'source_sha256':{str(p.relative_to(BUNDLE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),package/'geometry.py',package/'construction.py',package/'uniform_native.py')}}
    if return_spatial_groups:return result,np.asarray(distributions),coordinates,np.asarray(spatial_groups),unique_offsets
    return result,np.asarray(distributions),coordinates


def plot(result,distributions,coordinates,destination):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    samples=result['observables']['samples'];centers=np.array([s['centroid'] for s in samples]);rays=np.array([s['weighted_group_velocity_ray'] for s in samples]);times=np.array([s['time'] for s in samples])
    fig=plt.figure(figsize=(13,5),layout='constrained');ax=fig.add_subplot(121,projection='3d')
    p=distributions[-1].ravel();order=np.argsort(p)[-700:];xyz=coordinates.reshape(-1,3)[order]
    ax.scatter(*xyz.T,c=p[order],s=10,alpha=.15,cmap='Blues',label='Final cell-marginal cloud')
    ax.plot(*centers.T,'o-',color='#d65f2d',lw=2,label='Measured centroid');ax.plot(*rays.T,'--',color='#202b40',label='Weighted group-velocity ray')
    ax.set(xlabel='x',ylabel='y',zlabel='z',title='Explicit stacked ring: cell cloud and centroid');ax.legend(fontsize=8)
    bx=fig.add_subplot(122)
    for axis,color in enumerate(('#20639b','#d65f2d','#2a9d70')):
        bx.plot(times,centers[:,axis]-centers[0,axis],'-o',color=color,label='xyz'[axis]);bx.plot(times,rays[:,axis]-rays[0,axis],'--',color=color,alpha=.6)
    bx.set(xlabel='Native inverse-hop time t',ylabel='Centroid displacement',title='Three spatial components; unequal native speeds');bx.legend();bx.grid(alpha=.2)
    fig.suptitle('Positive-band wavepacket • declared vertical links • zero commit rate • no derived gravity',fontsize=12)
    fig.savefig(destination/'stacked_wavepacket.png',dpi=170);fig.savefig(destination/'stacked_wavepacket.svg');plt.close(fig)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--n',type=int,default=32);ap.add_argument('--output-dir',type=Path,default=Path(__file__).parent);args=ap.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True);result,probabilities,coordinates,spatial,offsets=simulate(args.n,return_spatial_groups=True)
    (args.output_dir/'STACKED_WAVEPACKET.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(args.output_dir/'position_distributions.npz',probability=probabilities,spatial_probability=spatial,internal_offsets=offsets,component_to_offset_bin=np.asarray(result['stored_spatial_apparatus']['component_to_offset_bin']),coordinates=coordinates,times=np.array([s['time'] for s in result['observables']['samples']]))
    plot(result,probabilities,coordinates,args.output_dir)
    print(json.dumps({'checks':result['checks'],'velocity':result['observables']['weighted_group_velocity'],'ray_error':result['observables']['max_centroid_ray_distance'],'edge_weight':result['observables']['maximum_edge_shell_probability']},indent=2))
