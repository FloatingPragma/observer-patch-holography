"""Complete relative cubical cochains: perfect-conductor boundary, no harmonic H1."""
from itertools import combinations, product
import numpy as np
from scipy import sparse
from scipy.linalg import qr


def relative_cube(cells):
    if type(cells) is not int or cells < 2:
        raise ValueError('at least two cubic cells per axis required')
    vertices = list(product(range(1, cells), repeat=3))
    vi = {x: i for i, x in enumerate(vertices)}
    edges = []
    for axis in range(3):
        ranges = [range(1, cells)]*3
        ranges[axis] = range(cells)
        edges.extend((axis, x) for x in product(*ranges))
    ei = {x: i for i, x in enumerate(edges)}
    grad = sparse.lil_matrix((len(edges), len(vertices)), dtype=np.int64)
    for (axis, x), i in ei.items():
        y = list(x); y[axis] += 1; y = tuple(y)
        if x in vi:
            grad[i, vi[x]] = -1
        if y in vi:
            grad[i, vi[y]] = 1
    faces = []
    for a, b in combinations(range(3), 2):
        ranges = [range(1, cells)]*3
        ranges[a] = ranges[b] = range(cells)
        faces.extend((a, b, x) for x in product(*ranges))
    curl = sparse.lil_matrix((len(faces), len(edges)), dtype=np.int64)
    for row, (a, b, x) in enumerate(faces):
        xa = list(x); xa[a] += 1; xa = tuple(xa)
        xb = list(x); xb[b] += 1; xb = tuple(xb)
        for edge, sign in (((a, x), 1), ((b, xa), 1), ((a, xb), -1), ((b, x), -1)):
            if edge in ei:
                curl[row, ei[edge]] += sign
    grad, curl = grad.tocsr(), curl.tocsr()
    assert (curl@grad).nnz == 0
    assert len(vertices) == (cells-1)**3
    assert len(edges) == 3*cells*(cells-1)**2
    assert len(faces) == 3*cells*cells*(cells-1)
    positions = np.array([np.array(x, float)+np.eye(3)[a]/2 for a, x in edges])
    return {'cells': cells, 'vertices': vertices, 'edges': edges, 'faces': faces,
            'gradient': grad, 'curl': curl, 'edge_positions': positions,
            'transverse_dimension': (cells-1)**2*(2*cells+1),
            'Hodge_stiffness_lower_exact': f'8/{cells*cells}',
            'boundary': 'relative cochains; all boundary-tangential edges are fixed zero',
            'harmonic_one_forms': 0}


def transverse_frame(cells):
    data = relative_cube(cells)
    g = data['gradient'].toarray().astype(float)
    q, r = qr(g, mode='full')
    assert min(abs(np.diag(r))) > .1
    frame = q[:, g.shape[1]:]
    assert frame.shape[1] == data['transverse_dimension']
    gram = np.linalg.norm(frame.T@frame-np.eye(frame.shape[1]), 'fro')
    divergence = np.linalg.norm(g.T@frame, 'fro')
    assert gram < 1e-11 and divergence < 1e-11
    data.update({'frame': frame, 'frame_Gram_residual': gram,
                 'frame_divergence_residual': divergence})
    return data


def complete_cube(cells):
    """All vertices/links/faces; zero-charge Gauss law also holds at the boundary."""
    if type(cells) is not int or cells < 1:
        raise ValueError('positive integer number of cells required')
    vertices = list(product(range(cells+1), repeat=3))
    vi = {x: i for i, x in enumerate(vertices)}
    edges = []
    for axis in range(3):
        ranges = [range(cells+1)]*3; ranges[axis] = range(cells)
        edges.extend((axis, x) for x in product(*ranges))
    ei = {x: i for i, x in enumerate(edges)}
    grad = sparse.lil_matrix((len(edges), len(vertices)), dtype=np.int64)
    for (axis, x), i in ei.items():
        y = list(x); y[axis] += 1
        grad[i, vi[x]] = -1; grad[i, vi[tuple(y)]] = 1
    faces = []
    for a, b in combinations(range(3), 2):
        ranges = [range(cells+1)]*3; ranges[a] = ranges[b] = range(cells)
        faces.extend((a, b, x) for x in product(*ranges))
    curl = sparse.lil_matrix((len(faces), len(edges)), dtype=np.int64)
    for row, (a, b, x) in enumerate(faces):
        xa = list(x); xa[a] += 1; xa = tuple(xa)
        xb = list(x); xb[b] += 1; xb = tuple(xb)
        for edge, sign in (((a, x), 1), ((b, xa), 1), ((a, xb), -1), ((b, x), -1)):
            curl[row, ei[edge]] += sign
    grad, curl = grad.tocsr(), curl.tocsr()
    assert (curl@grad).nnz == 0
    return {'cells': cells, 'vertices': vertices, 'edges': edges, 'faces': faces,
            'gradient': grad, 'curl': curl,
            'edge_positions': np.array([np.array(x, float)+np.eye(3)[a]/2 for a, x in edges]),
            'transverse_dimension': len(edges)-len(vertices)+1,
            'Hodge_stiffness_lower_exact': f'8/{(cells+1)**2}',
            'boundary': 'complete open cubical complex; zero charge at every vertex',
            'harmonic_one_forms': 0, 'gradient_kernel_dimension': 1}


def complete_transverse_frame(cells):
    data = complete_cube(cells)
    # The omitted column is exactly redundant: every row of grad sums to zero.
    g = data['gradient'].toarray().astype(float)
    q, r = qr(g[:, 1:], mode='full')
    assert min(abs(np.diag(r))) > .01
    frame = q[:, g.shape[1]-1:]
    assert frame.shape[1] == data['transverse_dimension']
    gram = np.linalg.norm(frame.T@frame-np.eye(frame.shape[1]), 'fro')
    divergence = np.linalg.norm(g.T@frame, 'fro')
    assert gram < 1e-11 and divergence < 1e-11
    data.update({'frame': frame, 'frame_Gram_residual': gram,
                 'frame_divergence_residual': divergence})
    return data
