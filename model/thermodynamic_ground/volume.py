"""Spatial averaging, exact boundary fractions and original-mu density bounds."""
from fractions import Fraction as F

import numpy as np

from ..cube_eta_protection import cube
from ..field_symmetry import symmetry
from ..full_fock_field.field import FieldGeometry
from ..field_stress.energy import GaussianInteger as GI
from ..interacting_vacuum import ground
from .local import stationarity_error


def ground_inequality_defect(half_side, support_radius, support_sites,
                             speed=F(1, 8), interaction=F(1), observable_norm=F(1)):
    """Bound real-part negative defect AND imaginary magnitude of the average."""
    if (type(half_side) is not int or half_side < 1 or type(support_radius) is not int or
            support_radius < 0 or type(support_sites) is not int or
            not 1 <= support_sites <= (2*support_radius+1)**3 or
            any(type(x) not in (int,F) or x < 0 for x in (speed,interaction,observable_norm))):
        raise ValueError('integer local support and nonnegative exact coefficients required')
    good_side=max(0,half_side-support_radius-1)
    bad=1-F(good_side**3,half_side**3)
    return 2*F(observable_norm)**2*(F(interaction)+6*F(speed))*support_sites*bad


def averaged_dynamics_error(half_side, buffer_radius, support_radius, support_sites,
                            duration, speed=F(1, 8)):
    """Uniform state/correlation comparison for spatially averaged box grounds.

    Original box has even side2L. Good translates contain the radiusR cube.
    On good translates, compare both dynamics to that cube (2B_R); on the
    remaining fraction use2. No state-convergence rate is inferred.
    """
    if (type(half_side) is not int or type(buffer_radius) is not int or
            not 0 <= buffer_radius < half_side):
        raise ValueError('integer0<=buffer radius<half side required')
    bound = stationarity_error(buffer_radius, support_radius, support_sites, duration, speed)
    bad = 1-F((half_side-buffer_radius)**3, half_side**3)
    return min(F(2), 2*bound+2*bad)


def source_kinetic_square():
    """Exact cross-Pauli cancellation in the actual open3^3 one-body source."""
    geometry = FieldGeometry((3, 3, 3))
    pieces = []
    from scipy import sparse
    for axis in range(3):
        rows, columns, values = [], [], []
        for group in geometry.groups[2*axis:2*axis+2]:
            for p,q,z in group:
                rows.extend((p,q)); columns.extend((q,p)); values.extend((16*z,16*z.conjugate()))
        pieces.append(GI.from_exact(sparse.coo_matrix((values,(rows,columns)),
                                                     shape=(geometry.modes,geometry.modes))))
    h = GI.from_exact(16*geometry.single_h())
    square = GI.zero(geometry.modes)
    for term in pieces:
        square = square+term@term
    if not (h@h-square).is_zero():
        raise ArithmeticError('Cartesian/Pauli cancellation failed on the original source')
    for j in range(3):
        for k in range(j):
            if not (pieces[j]@pieces[k]+pieces[k]@pieces[j]).is_zero():
                raise ArithmeticError('different Cartesian generators do not anticommute')
    return {'one_body_source_dimension': geometry.modes,
            'exact_identity': 'h^2=c^2 sum_j D_j^2; Cartesian D_j commute, ||D_j||<=1',
            'all_rectangular_boxes_bound': '-sqrt3*c I<=h<=sqrt3*c I',
            'full_Fock_bound': 'H=dGamma(h)+UD>=-sqrt3*c N forU>=0'}


def particle_hole_certificate():
    """Complete two-site CAR witness of the general real-bipartite identity."""
    modes, dimension = 4, 16
    old = [symmetry._annihilator(modes,j) for j in range(modes)]
    new = [symmetry._scale(a.adjoint(),(-1)**(j//2)) for j,a in enumerate(old)]
    vacuum = np.zeros(dimension,dtype=np.int64); vacuum[-1]=1
    columns=[]
    for mask in range(dimension):
        v=vacuum.copy()
        for j in reversed(range(modes)):
            if mask>>j&1: v=new[j].adjoint().real@v
        columns.append(v)
    p=GI.from_exact(np.column_stack(columns)); identity=symmetry._identity(dimension)
    if not (p@p.adjoint()-identity).is_zero():
        raise ArithmeticError('full particle-hole map lost unitarity')
    n=[a.adjoint()@a for a in old]
    number=n[0]+n[1]+n[2]+n[3]
    contact=n[0]@n[1]+n[2]@n[3]
    kinetic=old[0].adjoint()@old[2]+old[2].adjoint()@old[0]
    kinetic=kinetic+old[1].adjoint()@old[3]+old[3].adjoint()@old[1]
    two_identity=symmetry._scale(identity,2)
    if (not (p@kinetic@p.adjoint()-kinetic).is_zero() or
            not (p@contact@p.adjoint()-contact-two_identity+number).is_zero()):
        raise ArithmeticError('original-mu particle-hole energy identity failed')
    return {'complete_Fock_dimension': dimension,
            'general_source_join': 'existing all-box onsite SU2 frame gives real spin-independent bipartite hopping without changingD orN',
            'identity': 'P H P^*=H+U(V-N)',
            'sector_minima': 'E(2V-N)=E(N)+U(V-N)',
            'consequence_for_U_positive': 'No absolute original-mu0 ground component can haveN>V; its particle-hole partner has strictly lower energy.',
            'half_filling_concluded': False}


def density_certificate():
    # The cube reference is used as a trial state in a tiling, never claimed
    # to be the large-box ground or to remain stationary across tile bonds.
    a,b=ground.reference_arrays()
    indices=np.asarray(cube.weights(4),dtype=np.int64)
    if any(int(mask)&0xAAAA for mask,aa,bb in zip(indices,a,b,strict=True) if aa or bb):
        raise ArithmeticError('tile reference has an unexpected doublon')
    h=cube.restrict(cube.operators().hubbard_H.real,4)
    if (not np.array_equal(h@a,-12*b) or not np.array_equal(h@b,-4*a) or
            int(a@a+3*b@b)!=1296 or int(a@b)!=0):
        raise ArithmeticError('the normalized source tile energy changed')
    return {'domain': 'even-sided open3D boxes, fixedc>0,U>0, all particle numbers, originalmu0',
            'tile_trial': 'normalized cube highest-weightG on every disjoint2^3 tile; fixedN4 per tile',
            'tile_energy': '-2sqrt3*c', 'tile_contact': '0',
            'cross_tile_hopping_expectation': '0 because each tile has fixed number and hence zero odd expectation',
            'energy_density_upper': '-sqrt3*c/4',
            'kinetic': source_kinetic_square(),
            'number_density_lower': '1/4',
            'particle_hole': particle_hole_certificate(),
            'number_density_upper': '1',
            'averaging_required': 'Apply spatial averaging to the finite full-ground states before taking the weak-star subsequence; a centered local limit alone need not inherit an averaged density bound.',
            'infinite_ground_density': '1/4<=omega(n_0up+n_0down)<=1',
            'density_selected_uniquely': False,
            'fixed_particle_count_continuum_theorem_applicable': False,
            'finite_state_convergence_rate_claimed': False,
            'sample_averaged_dynamic_budgets': [
                {'half_side': l, 'buffer_radius': r, 'support_radius': 0,
                 'time': '1/8', 'error_upper': str(averaged_dynamics_error(l,r,0,1,F(1,8)))}
                for l,r in ((64,16),(256,32),(1024,64),(4096,128))]}
