"""Two sublattices on a subdivision tree with the seam-wedge hop of the scalar channel, without its level term."""

import math

import numpy as np
from scipy import sparse as sp

from model.candidate_laws.scalar_channel.carrier import require
from model.experiments.cage import construction as cell


def tree(depth, branching=4):
    """The parent-child register, in breadth-first order, with two modes per node."""
    require(type(depth) is int and depth >= 0 and branching == 4, 'invalid subdivision tree')
    nodes = (branching**(depth + 1) - 1) // (branching - 1)
    child = np.arange(1, nodes)
    shift = sp.csr_matrix((np.ones(len(child)), (child, (child - 1)//branching)), shape=(nodes, nodes))
    T = cell.HOP * (sp.eye(nodes, format='csr') - shift)
    return sp.bmat([[None, T], [T.T, None]], format='csr')


def chain(depth, length, links=True):
    """A root or sibling-difference chain, in the two-sublattice ordering."""
    require(1 <= length <= depth+1, 'invalid chain length')
    T = np.eye(length)
    if links:
        T -= 2 * np.eye(length, k=-1)
    zero = np.zeros_like(T)
    return cell.HOP * np.block([[zero, T], [T.T, zero]])


def chain_counts(depth):
    """Multiplicity of each radial or sibling-difference chain of node layers."""
    return {n: 1 if n == depth+1 else 3*4**(depth-n) for n in range(1, depth+2)}


def seed(depth, length):
    """A basis vector on the root, or a normalized difference of two siblings."""
    nodes = (4**(depth+1)-1)//3
    out = np.zeros(2*nodes)
    if length == depth+1:
        out[0] = 1
    else:
        first_at_layer = (4**(depth-length+1)-1)//3
        out[first_at_layer:first_at_layer+2] = (1/math.sqrt(2), -1/math.sqrt(2))
    return out


def outer_profile(car, name):
    """A preparation specified by seam labels and integer cell operations."""
    E, W = car['C'].shape
    if name == 'uniform':
        u = np.ones(E)
    elif name == 'port':
        require(car['ports'] == 12, 'the port profile belongs to the one-cell menu')
        Q = car['Q']
        A = Q @ Q.T - 5*np.eye(12)
        origin = np.eye(12)[:, 0]
        u = Q.T @ ((A - 5*np.eye(12)) @ ((A @ A - 5*np.eye(12)) @ origin))
    else:
        require(name == 'ring', 'unknown outer profile')
        if car['ports'] == 12:
            ring = cell.RING
        else:
            five = np.flatnonzero(car['degree'] == 5)[:12]
            ring = next(r for r in cell.rings(car, five) if all(car['degree'][i] == 5 for i in r))
        u = cell.ring_profile(car, ring).astype(float)
    u /= np.linalg.norm(u)
    image = cell.HOP * (car['C'].T @ u)
    p = float(np.linalg.norm(image))
    if name == 'uniform':
        require(p < 1e-13, 'the uniform seam mode leaks')
        return u, np.zeros(W), 0.
    partner = image / p
    require(np.max(abs(cell.HOP * (car['C'] @ partner) - p*u)) < 2e-12,
            'external profile loses a complementary mode')
    return u, partner, p


def reduced(depth, length, p, law):
    require(law in ('graded', 'no_children', 'unsigned'), 'unknown law')
    h = chain(depth, length, links=law != 'no_children')
    signs = np.ones(len(h))
    if law != 'unsigned':
        signs[length:] = -1
    return np.kron(h, np.eye(2)) + np.kron(np.diag(signs), [[0, p], [p, 0]])


def microscopic(depth, car, law):
    """Full sparse one-quantum generator including every outer wedge and every tree node."""
    require(law in ('graded', 'no_children', 'unsigned'), 'unknown law')
    B = cell.HOP * sp.csr_matrix(car['C'])
    D = sp.bmat([[None, B], [B.T, None]], format='csr')
    inner = tree(depth)
    nodes = inner.shape[0]//2
    if law == 'no_children':
        identity = cell.HOP * sp.eye(nodes)
        inner = sp.bmat([[None, identity], [identity, None]], format='csr')
    sign = np.ones(2*nodes)
    if law != 'unsigned':
        sign[nodes:] = -1
    return sp.kron(inner, sp.eye(D.shape[0]), format='csr') + sp.kron(sp.diags(sign), D, format='csr')


def inventory(depth):
    counts = chain_counts(depth)
    bands = []
    for n, copies in counts.items():
        levels = np.linalg.eigvalsh(chain(depth, n))
        positive = levels[n:] / cell.HOP
        # The light line can be far below the forward eigensolver's absolute precision.
        positive[0] = light_rate(n)
        require(0 < positive[0] <= 1+1e-13 and np.all(positive[1:] > 1), 'light-band separation failed')
        bands.append(dict(length=n, copies=copies, positive=positive.tolist()))
    nodes = (4**(depth+1)-1)//3
    require(sum(b['copies']*len(b['positive']) for b in bands) == nodes, 'incomplete signed spectrum')
    return bands


def light_rate(n):
    """Smallest singular value from the largest one of an exactly scaled inverse."""
    require(type(n) is int and n >= 1, 'a chain must be nonempty')
    i, j = np.indices((n,n))
    scaled = np.tril(np.exp2((i-j)-(n-1)))
    return 2.**(1-n)/np.linalg.norm(scaled, 2)
