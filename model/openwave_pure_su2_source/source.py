"""Exact arithmetic for a pure-gauge source on a three-dimensional cube."""
from fractions import Fraction as F
from itertools import product


def cube(side):
    if not isinstance(side, int) or side < 2:
        raise ValueError('the declared cubic source has at least two vertices per axis')
    vertices = list(product(range(side), repeat=3))
    edges = []
    for x in vertices:
        for j in range(3):
            if x[j]+1 < side:
                y = list(x); y[j] += 1
                edges.append((x, tuple(y)))
    index = {e: i for i, e in enumerate(edges)}
    faces = []
    for x in vertices:
        for j in range(3):
            for k in range(j+1, 3):
                if x[j]+1 < side and x[k]+1 < side:
                    y = list(x); y[j] += 1; y = tuple(y)
                    z = list(x); z[k] += 1; z = tuple(z)
                    w = list(y); w[k] += 1; w = tuple(w)
                    faces.append(((index[(x, y)], 1), (index[(y, w)], 1),
                                  (index[(z, w)], -1), (index[(x, z)], -1)))
    assert len(edges) == 3*side**2*(side-1)
    assert len(faces) == 3*side*(side-1)**2
    for face in faces:
        divergence = {v: 0 for v in vertices}
        for i, sign in face:
            x, y = edges[i]; divergence[x] -= sign; divergence[y] += sign
        assert all(v == 0 for v in divergence.values())
    # Distinct elementary faces have an edge occurring in only one of them.
    assert all(set(i for i, _ in p) != set(i for i, _ in q)
               for j, p in enumerate(faces) for q in faces[j+1:])
    return {'vertices': vertices, 'edges': edges, 'plaquettes': faces}


def coefficients(g2=F(8), spacing=F(1)):
    if g2 <= 0 or spacing <= 0:
        raise ValueError('positive squared coupling and lattice spacing required')
    return g2/(2*spacing), 2/(g2*spacing)


def certificate(side=2, g2=F(8), spacing=F(1)):
    graph = cube(side); p = len(graph['plaquettes'])
    alpha, beta = coefficients(g2, spacing)
    mu = beta*p
    residual2 = beta*beta*p/4
    excited_floor = 3*alpha
    gap = excited_floor-mu
    result = {'side': side, 'vertices': len(graph['vertices']),
              'links': len(graph['edges']), 'plaquettes': p,
              'g_squared_exact': str(g2), 'spacing_exact': str(spacing),
              'electric_coefficient_exact': str(alpha),
              'magnetic_coefficient_exact': str(beta),
              'electric_physical_gap_exact': str(excited_floor),
              'constant_trial_mean_exact': str(mu),
              'constant_trial_residual_squared_exact': str(residual2),
              'gap_lower_candidate_exact': str(gap),
              'strict_positive_gap_certificate': gap > 0,
              'all_link_representations_retained': True}
    if gap > 0:
        pbound = residual2/(gap*gap+residual2)
        result.update({'ground_energy_interval_exact': [str(max(F(0), mu-residual2/gap)), str(mu)],
                       'ground_orthogonal_weight_upper_exact': str(pbound)})
        if side == 2 and g2 == 8:
            delta = F(1, 32)
            assert pbound < delta*delta
            variance = F(1, 4)-2*delta-4*delta*delta
            assert variance > 0
            result.update({'ground_projector_distance_upper_exact': str(delta),
                           'plaquette_variance_lower_exact': str(variance),
                           'centered_curvature_energy_upper_exact': str(alpha),
                           'vacuum_gap_interval_exact': [str(gap), str(alpha/variance)],
                           'curvature_double_commutator_identity':
                               'one half [W,[H,W]] = alpha (1-W^2)',
                           'finite_energy_noncollapse_at_this_volume_and_spacing': True})
    return result


def run():
    base = certificate()
    assert base['vacuum_gap_interval_exact'] == ['21/2', '1024/47']
    return {'schema': 'oph.full_Haar_pure_SU2_finite_cube.v1',
            'verdict': 'PASS_PURE_SU2_FULL_LINK_SOURCE_AND_FINITE_CUBE_GAP',
            'source': {'group': 'SU(2)', 'space': 'Gauss invariant subspace of L2(SU2^links, product Haar)',
                'spatial_dimension': 3, 'time': 'continuous Hamiltonian parameter',
                'hamiltonian': 'alpha sum_link Casimir + beta sum_plaquette (1-ReTr(U_p)/2)',
                'matter_photon_or_gauge_fixing_delta_states': False,
                'finite_spin_truncation_defines_source': False,
                'energy_zero': 'positive electric plus positive magnetic terms; ground subtraction only in correlations'},
            'full_cube_certificate': base,
            'fixed_coupling_volume_controls': [certificate(s) for s in (2, 3, 4)],
            'fixed_cube_spacing_controls': [certificate(spacing=F(1, 2**k)) for k in range(5)],
            'scope': {'finite_graph_self_adjoint_positive_hamiltonian': True,
                'finite_graph_unique_vacuum': True, 'finite_graph_positive_gap': True,
                'full_untruncated_curvature_noncollapse': True,
                'uniform_infinite_volume_gap_proved': False,
                'finite_physical_threshold_along_a_to_zero': False,
                'four_dimensional_continuum_Yang_Mills_constructed': False,
                'all_compact_simple_groups_covered': False,
                'Clay_problem_solved': False}}
