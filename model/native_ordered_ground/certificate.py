"""Exact acceptance of an untrusted complete finite eigensystem proposal."""
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from . import source
from model.centered_particle_spectrum.certificate import sqrt_bounds,decimal_bounds

BITS=48


def propose():
    h,scale=source.matrix();e,v=np.linalg.eigh(h/scale)
    if v[:,0].sum()<0:v[:,0]*=-1
    return {'energies':np.rint(e*(1<<BITS)).astype(np.int64),
            'vectors':np.rint(v*(1<<BITS)).astype(np.int64)}


def verify(witness):
    h,scale=source.matrix();size=len(h)
    if set(witness)!= {'energies','vectors'}:raise ValueError('complete exact eigensystem required')
    if witness['energies'].shape!=(size,) or witness['vectors'].shape!=(size,size):raise ValueError('full matrix dimensions required')
    if any(a.dtype!=np.dtype('int64') for a in witness.values()):raise ValueError('literal int64 witness required')
    energies=witness['energies'].astype(object);v=witness['vectors'].astype(object)
    if any(a>b for a,b in zip(energies,energies[1:])):raise ValueError('ordered complete eigenvalues required')
    if any(x<=0 for x in v[:,0]):raise ValueError('strictly positive Marshall ground vector required')
    d=1<<BITS;gram=v.T@v-d*d*np.eye(size,dtype=object)
    defect=F(max(sum(abs(int(x)) for x in row) for row in gram),d*d)
    if not defect<F(1,1000000):raise ValueError('complete Gram defect too large')
    residual=h.astype(object)@v*d-scale*v*energies[None,:]
    eps=sqrt_bounds(F(sum(int(x)**2 for x in residual.flat),scale**2*d**4))[1]
    s=F(14);  # J|E|+BN = 12+2 for the exact declared example.
    if any(abs(F(int(x),d))>s for x in energies):raise ValueError('proposal outside exact Hamiltonian norm bound')
    radius=eps+2*s*defect
    e0,e1=F(int(energies[0]),d),F(int(energies[1]),d)
    if e1-e0<=2*radius:raise ValueError('ground eigenspace not separately isolated')
    norm=F(sum(int(x)**2 for x in v[:,0]),d*d)
    r0=sqrt_bounds(F(sum(int(x)**2 for x in residual[:,0]),scale**2*d**4)/norm)[1]
    error=2*r0/(e1-radius-e0)
    if error>=F(1,10**8):raise ValueError('ground vector approximation too loose')
    # In the Marshall frame the staggered X magnetization is uniform X.
    raw=sum(int(v[m,0])*int(v[m^(1<<site),0]) for m in range(size) for site in range(8))
    mag=F(raw,16*d*d)/norm
    return {'dimension':size,'vector_bits':BITS,'Gram_defect':str(defect),
      'complete_polar_Weyl_radius':str(radius),
      'ground_energy_interval':decimal_bounds(e0-r0,e0+r0,13),
      'first_excited_interval':decimal_bounds(e1-radius,e1+radius,12),
      'finite_field_gap_interval':decimal_bounds(e1-radius-e0-r0,e1+radius-e0+r0,12),
      'normalized_ground_vector_error':str(error),
      'staggered_magnetization_per_site_interval':decimal_bounds(mag-error,mag+error,12),
      'all_256_eigenvalues_certified':True,'untrusted_eigensolver_is_not_the_proof':True}
