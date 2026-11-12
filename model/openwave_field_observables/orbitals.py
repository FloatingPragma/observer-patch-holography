"""Bound levels of the computed Dirichlet Green potential with supplied kinetics."""
from __future__ import annotations
import numpy as np
from scipy import sparse
from scipy.fft import dstn, idstn
from scipy.sparse.linalg import eigsh


def parity_frame(n, odd=False):
    if n < 3 or n % 2 != 1:
        raise ValueError('odd grid size at least three required')
    center = n//2
    coordinates = np.arange(center+int(odd), n)
    rows, cols, values = [], [], []
    for j, x in enumerate(coordinates):
        if x == center:
            rows.append(x); cols.append(j); values.append(1.)
        else:
            rows.extend([x, n-1-x]); cols.extend([j, j])
            values.extend([1/np.sqrt(2), (-1 if odd else 1)/np.sqrt(2)])
    return sparse.csc_matrix((values, (rows, cols)), shape=(n, len(coordinates))), coordinates


def dirichlet_source(n, length):
    if length <= 0:
        raise ValueError('positive physical-model box length required')
    spacing = length/(n+1)
    one = 2-2*np.cos(np.pi*np.arange(1, n+1)/(n+1))
    symbol = one[:, None, None]+one[None, :, None]+one[None, None, :]
    source = np.zeros((n, n, n)); source[n//2, n//2, n//2] = 1
    green = idstn(dstn(source, type=1, norm='ortho')/symbol, type=1, norm='ortho')
    lap = sparse.diags((-np.ones(n-1), 2*np.ones(n), -np.ones(n-1)), (-1, 0, 1), format='csr')
    # g^2=4*pi makes the continuum coefficient one; this is a declared model unit.
    potential = -4*np.pi*green/spacing
    for axis in range(3):
        assert np.max(np.abs(np.flip(potential, axis=axis)-potential)) < 1e-12
    return lap, potential, spacing


def reduced_hamiltonian(n, length=32., odd_axes=(False, False, False), potential_sign=1.):
    lap, potential, spacing = dirichlet_source(n, length)
    frames = [parity_frame(n, odd) for odd in odd_axes]
    axes = [frame[1] for frame in frames]
    kinetic = [frame[0].T@lap@frame[0] for frame in frames]
    ident = [sparse.eye(len(index), format='csr') for index in axes]
    matrix = (sparse.kron(sparse.kron(kinetic[0], ident[1]), ident[2])
              +sparse.kron(sparse.kron(ident[0], kinetic[1]), ident[2])
              +sparse.kron(sparse.kron(ident[0], ident[1]), kinetic[2]))/(2*spacing*spacing)
    return (matrix+sparse.diags(potential_sign*potential[np.ix_(*axes)].ravel())).tocsr(), spacing


def levels(n, length=32.):
    even, a = reduced_hamiltonian(n, length)
    vals, vec = eigsh(even, k=2, which='SA', tol=1e-10,
                     v0=np.random.default_rng(981).normal(size=even.shape[0]))
    odd, _ = reduced_hamiltonian(n, length, (True, False, False))
    pval, pvec = eigsh(odd, k=1, which='SA', tol=1e-10,
                      v0=np.random.default_rng(981).normal(size=odd.shape[0]))
    residuals = [np.linalg.norm(even@vec[:,j]-vals[j]*vec[:,j]) for j in range(2)]
    residuals += [np.linalg.norm(odd@pvec[:,0]-pval[0]*pvec[:,0])]
    assert max(residuals) < 1e-8
    assert vals[0] < vals[1] < 0 and vals[0] < pval[0] < 0
    # The positive discrete kinetic operator gives essential threshold zero
    # after zero-extension of the finite potential into the infinite lattice.
    # Negative finite-box Rayleigh quotients therefore witness bound states
    # for that specified extension, independently of a fitted hydrogen series.
    return {'grid': n, 'length': length, 'spacing': a,
            'full_spatial_dimension': n**3, 'even_reduced_dimension': even.shape[0],
            'ground_energy': round(float(vals[0]), 10),
            'even_excited_energy': round(float(vals[1]), 10),
            'odd_axis_energy': round(float(pval[0]), 10),
            'ground_to_odd_spacing': round(float(pval[0]-vals[0]), 10),
            'second_shell_splitting': round(float(abs(pval[0]-vals[1])), 10),
            'maximum_eigen_residual': round(float(max(residuals)), 12)}


def run():
    rows = [levels(n) for n in (41, 63, 95)]
    assert all(rows[j]['second_shell_splitting'] > rows[j+1]['second_shell_splitting'] for j in (0,1))
    # Cubic permutations identify the three odd-axis blocks exactly; check
    # their spectra independently on a small complete domain.
    p = []
    for axis in range(3):
        parity = tuple(j == axis for j in range(3))
        h, _ = reduced_hamiltonian(9, 12., parity)
        p.append(np.linalg.eigvalsh(h.toarray())[0])
    assert max(p)-min(p) < 1e-12
    wrong, _ = reduced_hamiltonian(9, 12., potential_sign=-1)
    assert np.linalg.eigvalsh(wrong.toarray())[0] > 0
    return {'task': 'Green_bound_levels', 'rows': rows,
            'checks': {'negative_ground_and_excited_states': True,
                       'exact_parity_reduction': True, 'three_cubic_odd_sectors_equal': True,
                       'source_grid_refinement_reduces_second_shell_splitting': True},
            'mutations': {'repulsive_sign_has_no_negative_level': True},
            'supplied_parameters': {'mass': 1., 'coupling_squared': '4*pi', 'box_length': 32.,
                                    'kinetic_law': 'nonrelativistic nearest-neighbour Laplacian'},
            'hydrogen_reference_spacing': .375,
            'hydrogen_precision_match_claimed': False,
            'scope': 'Conditional bound-state consumer of the computed finite Gauss Green potential. Kinetic law and mass are supplied; no physical hydrogen identification, relativistic atom, or measured atomic spectrum is claimed.'}
