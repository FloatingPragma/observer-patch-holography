"""Independent complete-kernel and stationary field review."""
from pathlib import Path
import hashlib,json
import numpy as np
import scipy.linalg as la
import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=ROOT/'model/openwave_intrinsic_magnetic'


def absolute(matrix):
    values,vectors=la.eigh(matrix)
    return (vectors*abs(values))@vectors.conj().T


def derivative_review():
    from model.openwave_intrinsic_magnetic.solver import dirac,native_reduction,GAMMA,BETA
    from model.native_mass_channels.algebra import folded_generators
    n=3;a=.9;N=7;r=.4;m=5.;trap=.02
    rng=np.random.default_rng(41572);field=.07*rng.normal(size=(N,N,N,3))
    d,phases=dirac(n,a,N,field,m,trap,r);d=d.toarray()
    u=s.Matrix([[s.sympify(z) for z in row] for row in native_reduction()['exact_unitary_intertwiner']])
    native,masses=folded_generators()
    # Independently recover every spatial coefficient by the Clifford trace.
    coefficients=[];nn=n**3
    tensor=d.reshape(4,nn,4,nn)
    for matrix in (*GAMMA,BETA):
        coefficients.append(np.einsum('ij,jxiy->xy',matrix,tensor)/4)
    full=sum(np.kron(np.array(g.tolist(),complex),h) for g,h in zip((*native,masses[0]),coefficients))
    unitary=np.kron(np.array(u.tolist(),complex),np.eye(nn))
    full_error=float(la.norm(unitary.conj().T@full@unitary-np.kron(np.eye(4),d)))
    assert full_error<1e-12
    full_values=la.eigvalsh(full);block_values=la.eigvalsh(d)
    assert np.max(abs(full_values-np.sort(np.tile(block_values,4))))<2e-13
    assert min(abs(full_values))>=m-1e-12
    # A random complete link direction includes both spectral signs and all
    # positive/negative cross blocks of the derivative.
    direction=rng.normal(size=field.shape);direction/=la.norm(direction.ravel())
    h=2e-5
    plus=dirac(n,a,N,field+h*direction,m,trap,r)[0].toarray()
    minus=dirac(n,a,N,field-h*direction,m,trap,r)[0].toarray()
    derivative=(plus-minus)/(2*h);kernel=absolute(d)
    x=la.solve_sylvester(kernel,kernel,d@derivative+derivative@d)
    finite=(absolute(plus)-absolute(minus))/(2*h)
    error=float(la.norm(x-finite));assert error<3e-7
    ev,v=la.eigh(d);k=v.conj().T@x@v;dd=v.conj().T@derivative@v
    ratio=(ev[:,None]+ev[None,:])/(abs(ev[:,None])+abs(ev[None,:]))
    assert la.norm(k-ratio*dd)<1e-12
    for sign in (-1,1):
        ix=np.where(sign*ev>0)[0]
        assert la.norm(k[np.ix_(ix,ix)]-sign*dd[np.ix_(ix,ix)])<1e-12
    assert la.norm(x-derivative)>1e-3  # treating both signs as positive fails
    # Full fourfold copy degeneracy may be mixed without affecting the operator.
    angle=.37;rotation=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    mix=la.block_diag(rotation,rotation);w=np.kron(mix,np.eye(nn*4))
    repeated=np.kron(np.eye(4),x)
    assert la.norm(w.conj().T@repeated@w-repeated)<1e-12
    return {'complete_native_spatial_dimension':16*nn,'exact_reduction_numerical_operator_error':full_error,
        'all_signed_spectrum_minimum_absolute_eigenvalue':float(min(abs(full_values))),
        'Sylvester_vs_absolute_kernel_finite_difference_error':error,
        'both_signs_and_copy_degeneracies_retained':True,'incorrect_positive_only_derivative_rejected':True}


def stationary_review():
    from model.openwave_intrinsic_magnetic.solver import solve,dirac
    # Direct coordinate differences, without the producer incidence helpers.
    def forward(x,j,a):return (np.roll(x,-1,axis=j)-x)/a
    def adjoint(x,j,a):return (np.roll(x,1,axis=j)-x)/a
    def curl(x,a):return np.stack((forward(x[...,1],0,a)-forward(x[...,0],1,a),
        forward(x[...,2],1,a)-forward(x[...,1],2,a),
        forward(x[...,0],2,a)-forward(x[...,2],0,a)),axis=-1)
    result=solve(n=3,spacing=1.,field_n=9,coupling=.1)
    field=result['field'];current=result['current'];a=result['spacing'];kappa=result['coupling']
    b=curl(field,a)
    curlstar=np.stack((-adjoint(b[...,0],1,a)+adjoint(b[...,2],2,a),
        adjoint(b[...,0],0,a)-adjoint(b[...,1],2,a),
        adjoint(b[...,1],1,a)-adjoint(b[...,2],0,a)),axis=-1)
    maxwell=float(la.norm((curlstar-kappa*current/a**3).ravel()))
    continuity=float(la.norm(sum(adjoint(current[...,j],j,a) for j in range(3)).ravel()))
    assert maxwell<1e-9 and continuity<1e-10
    coords=np.moveaxis(np.indices((9,9,9)),0,-1)-4
    moment=np.zeros(3)
    for j in range(3):
        positions=coords.astype(float)*a;positions[...,j]+=a/2
        edge=np.zeros_like(current);edge[...,j]=current[...,j]
        moment+=np.cross(positions,edge).sum(axis=(0,1,2))/2
    assert la.norm(moment-result['moment'])<1e-14
    energy_field=a**3*np.sum(b*b)/(2*kappa)
    assert abs(energy_field-result['field_energy'])<1e-14
    psi=result['state'];kernel=absolute(result['dirac'].toarray())
    matter_residual=float(la.norm(kernel@psi-result['eigenvalue']*psi))
    assert matter_residual<1e-10
    rng=np.random.default_rng(79183);direction=rng.normal(size=field.shape);direction/=la.norm(direction.ravel())
    # At exact continuity and zero harmonic current the projected and full
    # variations agree, so an arbitrary field direction is a valid check.
    def energy(changed):
        d=dirac(3,a,9,changed,8.,.01,.5)[0].toarray()
        return np.vdot(psi,absolute(d)@psi).real+a**3*np.sum(curl(changed,a)**2)/(2*kappa)
    step=2e-5;variation=(energy(field+step*direction)-energy(field-step*direction))/(2*step)
    assert abs(variation)<3e-9
    rest=result['eigenvalue']+energy_field
    g1=-2*8*moment[2]/result['spin'];g2=-2*rest*moment[2]/result['spin']
    assert abs(g1-result['g_supplied_mass'])<1e-13 and abs(g2-result['g_actual_rest_energy'])<1e-13
    assert rest>8 and abs(g1-g2)>1e-3
    return {'independent_small_stationary_grid':3,'current_continuity_residual':continuity,
        'direct_incidence_Maxwell_residual':maxwell,'positive_kernel_eigenstate_residual':matter_residual,
        'complete_energy_field_variation':float(variation),'moment':moment.tolist(),
        'g_supplied_mass':float(g1),'g_bound_excitation':float(g2),
        'global_minimizer_or_numerical_interval_certificate_claimed':False}


def report():
    paths=[TARGET/p for p in ('solver.py','checks.py','experiment.py','PROOF.md','results.json')]+[Path(__file__)]
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    receipt=json.loads((TARGET/'results.json').read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in receipt['input_sha256'].items())
    result={'schema':'oph.independent_intrinsic_magnetic_review.v1','input_sha256':pins,
        'positive_kernel':derivative_review(),'stationary':stationary_review(),
        'scope':'The exact finite minimizer-existence proof does not identify the displayed floating stationary branch. The architecture retains confinement, Wilson mass, transverse-only Maxwell feedback and supplied coefficients; no unconfined electron, anomaly or precision g is inferred.',
        'verdict':'PASS_INDEPENDENT_INTRINSIC_MAGNETIC_REVIEW'}
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in pins.items())
    return result


if __name__=='__main__':
    result=report();(HERE/'intrinsic_magnetic_review.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(result['verdict'])
