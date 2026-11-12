"""Complete planar wave records on a native three-dimensional lapse profile."""
from functools import lru_cache
from fractions import Fraction as F
import numpy as np
from scipy.sparse import coo_matrix,diags,bmat,eye
from scipy.sparse.linalg import expm_multiply
from scipy.fft import dct
from scipy.linalg import eigh
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import quad
from model.openwave_energy_lapse.source import prepare,BITS

@lru_cache(maxsize=2)
def profile(n=31):
 p=prepare(n,(((n//2,)*3,1),));return p,p['integer_geometry'][:,:,n//2]/2**BITS

def matrix(q,gamma=1):
 q=np.asarray(q,float);n=q.shape[0]
 if q.shape!=(n,n)or n<5 or gamma not in(-1,0,1):raise ValueError('square finite profile and gamma -1,0,1 required')
 rows=[];cols=[];data=[];weights=[];edge=0
 for axis in range(2):
  for along in range(n+1):
   for cross in range(n):
    sites=[]
    for a,sign in((along,1),(along-1,-1)):
     if 0<=a<n:
      site=(a,cross)if axis==0 else(cross,a);sites.append(site);rows.append(edge);cols.append(np.ravel_multi_index(site,(n,n)));data.append(sign)
    qedge=sum(q[site]for site in sites)/2
    weights.append(np.exp((1-gamma)*qedge));edge+=1
 D=coo_matrix((data,(rows,cols)),shape=(edge,n*n)).tocsr();S=diags(np.exp((1+3*gamma)*q.ravel()/2))
 K=(S@D.T@diags(weights)@D@S).tocsr()
 return K


def initial(n,offset=6):
 x,y=np.meshgrid(np.arange(n),np.arange(n),indexing='ij');c=n//2
 v=np.exp(-((x-(c-10))**2+(y-(c+offset))**2)/18+1j*.9*x).ravel();return v/np.linalg.norm(v)

@lru_cache(maxsize=64)
def coefficients(n,duration,kind='wave'):
 # Positive spectral enclosure valid for all stated q and gamma choices.
 low=F(7,(n+1)**2);high=F(9);center=float((low+high)/2);radius=float((high-low)/2)
 degree=4096;M=16384;theta=np.pi*(np.arange(M)+.5)/M;x=center+radius*np.cos(theta)
 f=np.exp(-1j*duration*np.sqrt(x))if kind=='wave'else np.sqrt(x)
 c=dct(f,type=2)/M;c[0]/=2
 return c[:degree+1],center,radius

def apply_function(K,v,n,duration=0,kind='wave'):
 c,center,radius=coefficients(n,duration,kind);b1=np.zeros_like(v,dtype=complex);b2=b1.copy()
 def X(w):return (K@w-center*w)/radius
 for a in c[:0:-1]:b0=2*X(b1)-b2+a*v;b2,b1=b1,b0
 return X(b1)-b2+c[0]*v

def tail_bound(n,t,kind='wave'):
 # rho=1+1/(2(n+1)) stays inside the positive-square-root analytic domain.
 import mpmath as mp
 mp.mp.dps=70;rho=1+mp.mpf(1)/(2*(n+1));B=mp.exp(mp.mpf(9)*abs(t)/20)if kind=='wave'else mp.mpf(4)
 # Truncation plus midpoint-DCT aliasing; numerical rounding is separate.
 error=2*B*rho**(-4096)/(rho-1)+8*B*rho**(-2*16384+4096)/(rho-1)
 return str(error)

def rays(q,gamma,t,offset=6):
 n=len(q);center=n//2;field=RectBivariateSpline(np.arange(n),np.arange(n),q,kx=3,ky=3);result=[]
 for name,vx in(('continuum',1.),('lattice_carrier',float(np.cos(.9/2)))):
  def acceleration(s):return -(1+gamma)*float(field.ev(center-10+vx*s,center+offset,dx=0,dy=1))
  result.append({'ray_law':name,'flat_longitudinal_speed':vx,
   'first_order_transverse_velocity_shift':float(quad(acceleration,0,t,epsabs=1e-13)[0]),
   'first_order_transverse_position_shift':float(quad(lambda s:(t-s)*acceleration(s),0,t,epsabs=1e-13)[0])})
 return result

def experiment(n=31,gamma=1,times=(8,16,24),offset=6):
 p,q=profile(n);assert np.max(abs(q))<.001
 K=matrix(q,gamma);K0=matrix(np.zeros_like(q),gamma);psi=initial(n,offset);Y=np.tile(np.arange(n),n);Hpsi=apply_function(K,psi,n,kind='sqrt');H0psi=apply_function(K0,psi,n,kind='sqrt');rows=[];B=3*(n-1)
 for t in times:
  v=apply_function(K,psi,n,t);flat=apply_function(K0,psi,n,t)
  hv=apply_function(K,v,n,kind='sqrt');hf=apply_function(K0,flat,n,kind='sqrt')
  speed=2*np.vdot(v,Y*hv).imag;speedflat=2*np.vdot(flat,Y*hf).imag
  current=[.5,float((1+speed/B)/4),float((1-speed/B)/4)]
  pos=[.5,*map(float,abs(v)**2/2)];pos0=[.5,*map(float,abs(flat)**2/2)]
  assert abs(sum(pos)-1)<2e-11 and abs(sum(current)-1)<2e-13
  # Independent local second-order wave system, with the full initial velocity.
  m=n*n;A=bmat([[None,eye(m)],[-K,None]],format='csr');start=np.r_[psi,-1j*Hpsi]
  other=expm_multiply(t*A,start)[:m];discrepancy=float(np.linalg.norm(v-other));assert discrepancy<5e-11
  shift=float((abs(v)**2-abs(flat)**2)@Y)
  rows.append({'source_time':t,'complete_vacuum_and_all_site_records':pos,'flat_control_records':pos0,
   'conditional_transverse_position_shift':shift,'conditional_transverse_velocity':float(speed),'flat_transverse_velocity':float(speedflat),
   'transverse_velocity_shift':float(speed-speedflat),'complete_vacuum_and_velocity_POVM_records':current,
   'velocity_POVM_operator_norm_bound':B,'comparison_first_order_rays':rays(q,gamma,t,offset),
   'source_energy_before':float(np.vdot(psi,Hpsi).real/2),'source_energy_after':float(np.vdot(v,hv).real/2),
   'position_measurement_work_absolute_bound':'3','velocity_measurement_work_absolute_bound':'3',
   'independent_full_second_order_wave_discrepancy':discrepancy,
   'analytic_Chebyshev_truncation_alias_bound':tail_bound(n,t),'arithmetic_roundoff_is_included_in_analytic_bound':False})
 return {'side':n,'gamma':gamma,'transverse_preparation_offset':offset,'field_carrier_dimension':n*n+1,'native_geometry_sites':n**3,
  'profile_sha256':__import__('hashlib').sha256(p['integer_geometry'].astype('<i8').tobytes()).hexdigest(),
  'complete_profile_slice':q.tolist(),'all_mode_rows':rows,'all_wave_modes_retained':True}

def run():
 rows=[experiment(31,gamma)for gamma in(1,0,-1)]+[experiment(47,1,times=(16,)),experiment(31,1,times=(16,),offset=-6)]
 primary=rows[0]['all_mode_rows'][1];reflected=rows[-1]['all_mode_rows'][0]
 assert primary['conditional_transverse_position_shift']<0 and primary['transverse_velocity_shift']<0
 assert abs(primary['conditional_transverse_position_shift']+reflected['conditional_transverse_position_shift'])<1e-10
 assert abs(primary['transverse_velocity_shift']+reflected['transverse_velocity_shift'])<1e-10
 zero=np.zeros((31,31));assert (matrix(zero,1)!=matrix(zero,-1)).nnz==0
 return {'schema':'oph.complete_native_planar_wave_deflection.v1','variants':rows,
  'finite_field_carrier':'vacuum direct-sum every planar one-quantum site; no mode projection',
  'source':'H=sqrt(S D^T A D S), S=exp((1+3gamma)q/2), A=exp((1-gamma)q_edge)',
  'preparation':'equal vacuum and normalized Gaussian one-quantum mixture; carrier momentum0.9; width3; common flat/curved state',
  'profile_origin':'full native8-level energy_lapse side31 or47 central level1 dyadic minimum; slice z=center',
  'geometry_is_static_for_this_wave_experiment':True,
  'planar_extension_is_not_claimed_as_3D_invariant_subspace':True,
  'quantum_geometry_package_source_transfer_claimed':False,
  'position_instrument':'orthogonal site/vacuum projectors with orthogonal copied records',
  'velocity_instrument':'vacuum projector and Lueders square roots of (I +/- i[H,Y]/B)/2 on excited sector, B=3(n-1)',
  'arbitrary_external_references_supported':True,'zero_profile_gamma_controls_exact':True,
  'reflection_control_position_residual':abs(primary['conditional_transverse_position_shift']+reflected['conditional_transverse_position_shift']),
  'reflection_control_velocity_residual':abs(primary['transverse_velocity_shift']+reflected['transverse_velocity_shift']),
  'ray_relation':'local principal symbol exp((1+gamma)q)*sqrt(4sin(px/2)^2+4sin(py/2)^2); explicit first-order ray comparisons, no wave-ray error theorem',
  'gamma_minus1_ray_deflection_zero_but_finite_wave_residual_retained':True,
  'physical_G_lensing_angle_or_gravity_strength_identified':False,
  'verdict':'PASS_COMPLETE_2D_NATIVE_PROFILE_WAVE_DEFLECTION'}
