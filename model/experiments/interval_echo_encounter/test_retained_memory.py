"""A direct one-cell replay in the eleven-state memory and all binary/doubled two-load states."""
from itertools import combinations_with_replacement

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import carrier

from . import memory


def direct_term(matrix, edge, configurations, lookup, memories):
    """Apply each local full-basis column to the occupations with the spectators held fixed."""
    n = len(configurations)
    rows, columns, values = [], [], []
    for m in range(memories):
        for i, state in enumerate(configurations):
            incoming = tuple(state.count(x) for x in edge)
            spectators = [p for p in state if p not in edge]
            col = 9*m + law.SRC.index(incoming)
            for row in np.flatnonzero(abs(matrix[:, col]) > 1e-13):
                label, src = divmod(int(row), 9)
                x, y = law.SRC[src]
                dest = tuple(sorted(spectators + [edge[0]]*x + [edge[1]]*y))
                rows.append(label*n + lookup[dest])
                columns.append(m*n+i)
                values.append(matrix[row, col])
    return sparse.csr_matrix((values, (rows, columns)), shape=(matrix.shape[0]//9*n, memories*n))


def test_full_eleven_state_working_memory_and_doubled_loads():
    car = carrier.cell()
    states = list(combinations_with_replacement(range(car.n), 2))
    lookup = {s: i for i, s in enumerate(states)}
    dim = len(states)
    target = min(car.seams)
    u = law.echo()
    h0, h1 = law.local_term_parts()
    tagged = direct_term(h0 + 8*h1, target, states, lookup, 11)
    patch = direct_term(u, target, states, lookup, 11)
    others = sum((direct_term(law.seam_term(8), e, states, lookup, 1)
                  for e in car.seams if e != target), start=sparse.csr_matrix((dim, dim)))
    ref = np.array([8*sum((6-car.degree[i])*s.count(i)**2 for i in set(s))-96 for s in states])
    h = sparse.kron(sparse.eye(11), others+sparse.diags(ref, dtype=float), format="csr") + tagged
    free = h[:dim, :dim]
    e, v = eigh(free.toarray())
    pi = v[:, e < -5] @ v[:, e < -5].T
    ids = [i for i, s in enumerate(states) if s[0] != s[1] and s not in car.seams]
    mixture = np.eye(dim)[:, ids]
    first = np.zeros((11*dim, len(ids)), dtype=complex)
    first[:dim] = expm_multiply(-1j*free, mixture)
    written = patch @ first
    np.testing.assert_allclose(patch.conj().T @ written, first, rtol=0, atol=1e-12)
    band = sum(np.vdot(branch, pi @ branch).real
               for branch in written.reshape(11, dim, len(ids))) / len(ids)
    late = expm_multiply(-4j*h, written).reshape(11, dim, len(ids))
    occupation = np.sum(abs(late)**2, axis=(0, 2)) / len(ids)
    doubled = [lookup[i, i] for i in range(car.n)]
    assert sum(occupation[doubled]) < 1e-12
    np.testing.assert_allclose(occupation.sum(), 1., rtol=0, atol=1e-12)
    selected = occupation[[lookup[e] for e in sorted(car.seams)]]
    summary = memory.run(car)
    np.testing.assert_allclose(summary["final_distribution"]["retained"],
                               np.append(selected, 1-selected.sum()), rtol=0, atol=1e-11)
    np.testing.assert_allclose(summary["load_band_immediately_after_service"]["retained"], band,
                               rtol=0, atol=1e-11)
