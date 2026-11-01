"""Exact quadratic kinetic stresses and quartic contact forces.

The one-particle central difference is a declared lattice momentum observable,
not a physical translation generator. Complete-Fock Clifford arithmetic binds
its finite balance to the unchanged open-box Weyl/contact action.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product

import sympy as s

from ..coarse_energy_ward import ward as w
from ..full_fock_field.field import FieldGeometry
from ..full_fock_field.local_simulation import local_tail_boundary_bound, split_bound


def rational(value):
    value=s.sympify(value)
    real,imag=value.as_real_imag()
    if not real.is_Rational or not imag.is_Rational:
        raise ValueError('Gaussian rational one-particle matrix required')
    return F(int(real.p),int(real.q)),F(int(imag.p),int(imag.q))


def lift(matrix):
    """dGamma(q)=sum_ab q_ab c_a^dagger c_b; no filled-vacuum subtraction."""
    if not isinstance(matrix,s.MatrixBase) or matrix.rows != matrix.cols:
        raise ValueError('square exact one-particle matrix required')
    ladders=[w.annihilator(p) for p in range(matrix.rows)]
    return w.total((ladders[a].adjoint()@ladders[b]).scale(rational(value))
                   for (a,b),value in matrix.todok().items())


def one_particle(geometry):
    if not isinstance(geometry,FieldGeometry) or geometry.volume > 64:
        raise ValueError('one through 64 sites in an unchanged open box required')
    if max(F(geometry.speed).denominator,F(geometry.interaction).denominator)>2**16:
        raise ValueError('binary source coefficient denominator exceeds exact domain')
    modes=geometry.modes
    h=s.MutableSparseMatrix(modes,modes,{})
    momenta=[s.MutableSparseMatrix(modes,modes,{}) for _ in range(3)]
    edges={}
    for group in geometry.groups:
        for p,q,z in group:
            z=s.Rational(F(z.real))+s.I*s.Rational(F(z.imag))
            h[p,q]+=z;h[q,p]+=s.conjugate(z)
    indices={site:x for x,site in enumerate(geometry.sites)}
    for x,site in enumerate(geometry.sites):
        for axis in range(3):
            neighbor=list(site);neighbor[axis]+=1
            y=indices.get(tuple(neighbor))
            if y is None:
                continue
            a=s.MutableSparseMatrix(modes,modes,{})
            for spin in range(2):
                p,q=2*x+spin,2*y+spin
                momenta[axis][p,q]=-s.I/2
                momenta[axis][q,p]=s.I/2
                for target in range(2):
                    a[p,2*y+target]=h[p,2*y+target]
                    a[2*y+target,p]=-h[2*y+target,p]
            edges[x,y]=s.ImmutableSparseMatrix(a)
    return s.ImmutableSparseMatrix(h),tuple(s.ImmutableSparseMatrix(q) for q in momenta),edges


def site_restriction(matrix,site):
    """The Hermitian localization {Q_x,q}/2; endpoints share each bond."""
    return s.ImmutableSparseMatrix(matrix.rows,matrix.cols,
        {(a,b):value*(int(a//2==site)+int(b//2==site))/2
         for (a,b),value in matrix.todok().items()
         if a//2==site or b//2==site})


def contact(geometry):
    u=F(geometry.interaction)
    return w.total((w.annihilator(2*x).adjoint()@w.annihilator(2*x)@
                    w.annihilator(2*x+1).adjoint()@w.annihilator(2*x+1)).scale(u)
                   for x in range(geometry.volume))


def contact_edge_force(x,y,interaction):
    """U/2 sum_s (n_x,bar(s)-n_y,bar(s)) (c_xs^dag c_ys+h.c.)."""
    u=F(interaction)
    terms=[]
    for spin in range(2):
        a,b=w.annihilator(2*x+spin),w.annihilator(2*y+spin)
        ax,bx=w.annihilator(2*x+1-spin),w.annihilator(2*y+1-spin)
        difference=ax.adjoint()@ax-bx.adjoint()@bx
        hop=a.adjoint()@b+b.adjoint()@a
        terms.append((difference@hop).scale(u/2))
    return w.total(terms)


def construct(geometry=FieldGeometry((3,3,3))):
    h,one_momenta,edges=one_particle(geometry)
    kinetic=lift(h)
    potential=contact(geometry)
    original=w.total(w.energies(geometry))
    if kinetic+potential != original:
        raise ArithmeticError('momentum balance changed the original source action')
    columns=[]
    for axis,q in enumerate(one_momenta):
        if h*q != q*h:
            raise ArithmeticError('Cartesian central difference must commute with free action')
        local=tuple(lift(site_restriction(q,x)) for x in range(geometry.volume))
        stress={edge:lift(s.I*(a*q+q*a)/2) for edge,a in edges.items()}
        stress={edge:value for edge,value in stress.items() if not value.is_zero()}
        force=tuple(w.commutator(potential,p).scale((0,1)) for p in local)
        explicit=[w.Clifford({}) for _ in geometry.sites]
        for x,y in edges:
            if geometry.sites[y][axis]-geometry.sites[x][axis] != 1:
                continue
            edge_force=contact_edge_force(x,y,geometry.interaction).scale(F(1,2))
            explicit[x]=explicit[x]+edge_force
            explicit[y]=explicit[y]+edge_force
        if force != tuple(explicit):
            raise ArithmeticError('contact force differs from the source CAR commutator')
        divergence=w.divergence(stress,geometry.volume)
        for x,p in enumerate(local):
            if p != p.adjoint() or force[x] != force[x].adjoint():
                raise ArithmeticError('momentum or force lost Hermiticity')
            if not (w.commutator(original,p).scale((0,1))+divergence[x]-force[x]).is_zero():
                raise ArithmeticError('exact forced momentum balance failed')
        if any(current != current.adjoint() for current in stress.values()):
            raise ArithmeticError('kinetic stress lost Hermiticity')
        for widths in product(*(range(1,n+1) for n in geometry.sides)):
            labels=w.block_labels(geometry.sites,widths)
            if w.divergence(w.aggregate_flows(stress,labels),max(labels)+1) != w.aggregate(divergence,labels):
                raise ArithmeticError('forced stress divergence failed to commute with blocks')
        columns.append({'local':local,'total':lift(q),'stress':stress,'force':force,
                        'total_force':w.total(force)})
    return {'H':original,'H0':kinetic,'V':potential,'h0':h,
            'one_particle_momenta':one_momenta,'columns':tuple(columns),
            'partitions':geometry.volume}


def normalized_HS_squared(operator):
    """Exact Tr(O^dag O)/dim on the entire Fock algebra, by basis orthogonality."""
    return sum((real*real+imag*imag for real,imag in operator.terms.values()),F())


def no_go_matrix_elements(sites=3):
    """N=2 doublon transitions recover every offsite quadratic coefficient.

    For each ordered pair of different sites and both spin choices, there
    is one input occupation with one fermion at each site and one output
    doublon. Only its selected c_xs^dag c_yt has a nonzero matrix element.
    Return exact signs, certifying injectivity of the contact commutator.
    """
    if type(sites) is not int or not 2 <= sites <= 64:
        raise ValueError('two through 64 sites required')
    def hop(mask,p,q):
        if not mask>>q&1:
            return None,0
        sign=(-1)**(mask&((1<<q)-1)).bit_count()
        mask ^= 1<<q
        if mask>>p&1:
            return None,0
        sign *= (-1)**(mask&((1<<p)-1)).bit_count()
        return mask^(1<<p),sign
    witnesses=[]
    for x in range(sites):
        for y in range(sites):
            if x==y:
                continue
            for spin,target in product(range(2),repeat=2):
                p,q=2*x+spin,2*y+target
                initial=(1<<(2*x+1-spin))|(1<<q)
                final=(1<<(2*x))|(1<<(2*x+1))
                found=[]
                for a in range(2*sites):
                    for b in range(2*sites):
                        output,sign=hop(initial,a,b)
                        if output==final and sign:
                            found.append((a,b,sign))
                if len(found)!=1 or found[0][:2]!=(p,q):
                    raise ArithmeticError('offsite coefficient not isolated by doublon transition')
                witnesses.append((p,q,initial,final,found[0][2]))
    return tuple(witnesses)


def conserved_quadratics(geometry=FieldGeometry((3,3,3))):
    """Exact onsite number/hidden-spin basis after the generic quadratic no-go."""
    if geometry.speed<=0 or geometry.interaction<=0:
        raise ValueError('strictly positive hopping and contact required for four-charge classification')
    h,_,_=one_particle(geometry)
    pauli=(s.eye(2),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    basis=[]
    for generator in pauli:
        matrix=s.MutableSparseMatrix(geometry.modes,geometry.modes,{})
        for x,site in enumerate(geometry.sites):
            frame=(pauli[1]**site[0])*(pauli[2]**site[1])*(pauli[3]**site[2])
            onsite=frame*generator*frame.conjugate().T
            matrix[2*x:2*x+2,2*x:2*x+2]=onsite
        matrix=s.ImmutableSparseMatrix(matrix)
        if h*matrix != matrix*h:
            raise ArithmeticError('the transported onsite generator failed free commutation')
        basis.append(matrix)
    return tuple(basis)


def cell_error(radius,duration,steps,*,samples=None,speed=F(1,8),interaction=F(1)):
    """Ambient-volume-independent error for a sampled forced momentum cell.

    The norm estimates are p<=1, stress<=2c, force<=U. Source initialization,
    local input, physical readout/hardware precision and failure energy are
    separately supplied. No source execution is asserted by this bound.
    """
    if (type(radius) is not int or not 2<radius<=4096 or type(steps) is not int or steps<1 or
            any(type(x) not in (int,F) or x<0 for x in (duration,speed,interaction)) or
            (samples is not None and (type(samples) is not int or samples<1))):
        raise ValueError('bounded radius>2, nonnegative exact data and positive integer counts required')
    t,c,u=F(duration),F(speed),F(interaction)
    epsilon=min(F(2),split_bound((2*radius+1,)*3,t,steps,c,u))
    def error(norm,support_radius,support_sites,integrated):
        tail=local_tail_boundary_bound(radius,support_radius,support_sites,t,c)
        bound=norm*min(F(2),tail+2*epsilon)
        if integrated:
            bound *= t
            if samples is not None:
                bound += t*t*norm*support_sites*(u+6*c)/(2*samples)
        return bound
    endpoint=error(F(1),1,3,False)
    stress=error(2*c,2,6,True)
    force=error(u,1,3,True)
    return {'source_unitary_error':epsilon,'momentum_endpoint_error':endpoint,
            'integrated_stress_face_error':stress,'integrated_force_error':force,
            'forced_cell_residual_bound':2*endpoint+6*stress+force}


def force_gram(columns):
    """Normalized full-Fock Hilbert--Schmidt Gram, exactly rational."""
    result=[]
    for left in columns:
        row=[]
        for right in columns:
            value=(F(),F())
            for mask,(a,b) in left['total_force'].terms.items():
                c,d=right['total_force'].terms.get(mask,(F(),F()))
                value=(value[0]+a*c+b*d,value[1]+a*d-b*c)
            if value[1]:
                raise ArithmeticError('Hermitian force Gram has an imaginary part')
            row.append(value[0])
        result.append(tuple(row))
    return tuple(result)
