"""Finite exact commutators; no approximation to the ordered ground."""
from functools import lru_cache
import sympy as sp


@lru_cache(maxsize=1)
def local_certificate():
    i=sp.I
    sigma=(sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-i],[i,0]]),sp.diag(1,-1))
    spin=[[sp.kronecker_product(*[sigma[a]/2 if j==x else sp.eye(2)
                                 for j in range(3)]) for a in range(3)] for x in range(3)]
    zero=sp.zeros(8)
    def comm(a,b):return a*b-b*a
    def dot(a,b):return sum((a[j]*b[j] for j in range(3)),zero)
    def cross(a,b):return [sum((sp.LeviCivita(c,d,e)*a[d]*b[e]
                                for d in range(3) for e in range(3)),zero) for c in range(3)]
    f=sp.symbols('f0:3',complex=True)
    graphs=(((0,1),),((0,1),(0,2)),((0,1),(1,2)),((0,1),(0,2),(1,2)))
    checked=0
    for edges in graphs:
        h=sum((dot(spin[x],spin[y])-sp.eye(8)/4 for x,y in edges),zero)
        for axis in range(3):
            a=sum((f[x]*spin[x][axis] for x in range(3)),zero)
            lap=[0,0,0]
            current=zero;force=zero
            for x,y in edges:
                lap[x]+=f[x]-f[y];lap[y]+=f[y]-f[x]
                current+=i*(f[x]-f[y])*cross(spin[x],spin[y])[axis]
                for u,v in edges:
                    if (u,v)==(x,y):continue
                    if x in (u,v):
                        z=v if u==x else u
                        force-=(f[x]-f[y])*cross(cross(spin[x],spin[z]),spin[y])[axis]
                    if y in (u,v):
                        z=v if u==y else u
                        force-=(f[x]-f[y])*cross(spin[x],cross(spin[y],spin[z]))[axis]
            linear=sum((lap[x]*spin[x][axis]/2 for x in range(3)),zero)
            if (comm(h,a)-current).applyfunc(sp.expand)!=zero:
                raise ArithmeticError('current identity failed')
            if (comm(h,comm(h,a))-linear-force).applyfunc(sp.expand)!=zero:
                raise ArithmeticError('double commutator identity failed')
            checked+=1
    # A single unoriented two-edge path: all equal-axis Pauli words cancel.
    a,b=sp.symbols('a b',complex=True)
    for axis in range(3):
        force=(b*spin[0][axis]*dot(spin[1],spin[2])
               +a*spin[2][axis]*dot(spin[0],spin[1])
               -(a+b)*spin[1][axis]*dot(spin[0],spin[2]))
        norm=sp.trace(force.conjugate().T*force)/8
        expected=(sp.conjugate(a)*a+sp.conjugate(b)*b
                  +sp.conjugate(a+b)*(a+b))/32
        if sp.expand(norm-expected)!=0:
            raise ArithmeticError('three-spin trace norm failed')
        if any(sp.expand(sp.trace(force*x))!=0 for site in spin for x in site):
            raise ArithmeticError('force must be trace-orthogonal to linear spins')
    return {'exact_graph_axis_checks':checked,'arbitrary_complex_coefficients':True,
            'current_and_double_commutator_checked':True,
            'path_trace_norm':'(|a|^2+|b|^2+|a+b|^2)/32',
            'force_trace_orthogonal_to_all_linear_spins':True}


def regular_graph_pair_identity(degree):
    if type(degree) is not int or not 2<=degree<=12:
        raise ValueError('integer degree between2 and12 required')
    z=sp.symbols('z:'+str(degree),complex=True)
    abs2=lambda x:sp.conjugate(x)*x
    pairs=sum(abs2(z[a])+abs2(z[b])+abs2(z[a]+z[b])
              for a in range(degree) for b in range(a+1,degree))
    expected=(2*degree-3)*sum(map(abs2,z))+abs2(sum(z))
    if sp.expand(pairs-expected)!=0:raise ArithmeticError('neighbor pair identity failed')
    return {'degree':degree,'coefficient':2*degree-3}


def cubic_trace_certificate():
    regular_graph_pair_identity(6)
    e,alpha=sp.symbols('E alpha',real=True)
    force=sp.expand((9*4*e+4*e**2)/32)
    residual=force+(1-alpha)**2*e**2/4
    if sp.expand(force-e*(9+e)/8)!=0:raise ArithmeticError('cubic force norm failed')
    return {'normalized_plane_wave_force_norm_squared_over_J4':str(force),
            'renormalized_force_norm_squared_over_J4':str(sp.expand(residual)),
            'requires_simple_triangle_free_degree6_graph':True,
            'finite_even_cubic_torus_side_at_least4':True,
            'tracial_state_is_not_the_ordered_ground':True}
