"""Independent full 256-state periodic-cube controls (not a thermodynamic proof)."""
import numpy as np


def matrices(field):
    if type(field) not in (int,float) or not np.isfinite(field) or field<0:
        raise ValueError('nonnegative finite real field required')
    n=8;size=1<<n; rows=np.arange(size);ops=[]
    for x in range(n):
        z=np.diag(.5-((rows>>x)&1)).astype(complex)
        sx=np.zeros((size,size),complex);sy=sx.copy()
        sx[rows^(1<<x),rows]=.5
        sy[rows^(1<<x),rows]=.5j*(1-2*((rows>>x)&1))
        ops.append((sx,sy,z))
    h=np.zeros((size,size),complex)
    for x in range(n):
        for axis in range(3):
            y=x^(1<<axis)
            # Length-two torus diagnostic: both parallel periodic bonds retained.
            unlike=((rows>>x)&1)!=((rows>>y)&1)
            h[rows,rows]+=-.5*unlike
            h[rows[unlike]^(1<<x)^(1<<y),rows[unlike]]+=.5
    order=sum(((-1)**x.bit_count()*ops[x][0] for x in range(n)),np.zeros_like(h))
    return h-field*order,ops,order


def check():
    h,ops,order=matrices(.2);vals,vecs=np.linalg.eigh(h);g=vecs[:,0];gaps=vals-vals[0]
    assert gaps[1]>1e-5
    magnetization=float(np.vdot(g,order@g).real/8)
    assert 0<magnetization<.5
    max_error=0.; checked=0
    for momentum in range(8):
        e=2*momentum.bit_count();eq=6-e
        a=sum(((-1)**(x&momentum).bit_count()*ops[x][1] for x in range(8)),np.zeros_like(h))/np.sqrt(8)
        c=sum(((-1)**(x&(momentum^7)).bit_count()*ops[x][2] for x in range(8)),np.zeros_like(h))/np.sqrt(8)
        av=a@g;cv=c@g
        ward=(-1j*np.vdot(g,(a@c-c@a)@g)).real
        max_error=max(max_error,abs(ward-magnetization))
        weights=abs(vecs.conj().T@av)**2
        assert weights[0]<1e-24
        moment=float(weights@gaps)
        assert moment<=e+.1+1e-10
        if eq:
            susceptibility=float(weights[1:]@(1/gaps[1:]))
            assert susceptibility<=1/(2*eq)+1e-10
            checked+=1
        assert abs(np.vdot(av,(h-vals[0]*np.eye(256))@av).real-moment)<1e-10
    assert max_error<1e-12
    # Full, unprojected zero-field cube verifies the biased-ground construction.
    h,ops,order=matrices(0);vals,vecs=np.linalg.eigh(h);g=vecs[:,0];m=order/8
    e=-(vals[0]+6)/8
    for momentum in range(7):
        ek=2*momentum.bit_count();eq=6-ek
        a=sum(((-1)**(x&momentum).bit_count()*ops[x][2] for x in range(8)),np.zeros_like(h))/np.sqrt(8)
        w=abs(vecs.conj().T@(a@g))**2;gap=vals-vals[0]
        assert abs(w@gap-2*e*ek/9)<1e-11
        assert w[1:]@(1/gap[1:])<=1/(4*eq)+1e-11
        assert sum(w)<=np.sqrt(e*ek/(18*eq))+1e-11
    b=float(np.vdot(m@g,m@g).real)
    psi=(g+m@g)/np.sqrt(1+b)
    assert abs(np.vdot(g,m@g))<1e-12
    assert abs(np.vdot(g,m@m@m@g))<1e-12
    assert abs(np.vdot(psi,m@psi).real-2*b/(1+b))<1e-12
    energy=float(np.vdot(psi,(h-vals[0]*np.eye(256))@psi).real)
    assert 0<energy<=3/8
    return {'full_spin_dimension':256,'periodic_bonds_counting_parallel_edges':24,
            'positive_field_unique_ground':True,'ward_momenta_checked':8,
            'nonsingular_susceptibilities_checked':checked,'biased_ground_formula_checked':True,'tight_zero_field_DLS_moments_checked':7,
            'finite_cube_proves_infinite_order':False,'native_pulses_executed':False}
