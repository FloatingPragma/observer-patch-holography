"""Exact compression of the native generator to the nine admitted fixed-source colours."""

from itertools import combinations

import numpy as np

from model.executable_model import carrier


def physical_frame(band):
    """Use ordered physical port probes to remove the eigensolver's frame freedom."""
    projector = band @ band.T
    basis = []
    for v in projector.T:
        residual = v.copy()
        for u in basis:
            residual -= u * (u @ residual)
        if np.linalg.norm(residual) > 1e-9:
            basis.append(residual / np.linalg.norm(residual))
        if len(basis) == 3:
            return np.column_stack(basis)
    raise ValueError('Missing native colour axes')


def tensor(operators):
    return np.kron(np.kron(operators[0], operators[1]), operators[2])


def fixed_source_matrices(graph, ports, owners, frames, kappa, active=(True, True, True)):
    eye = np.eye(3)
    single = carrier.cell().laplacian().toarray() / 2
    separate = np.zeros((27, 27))
    contact = np.zeros_like(separate)
    for i, (owner, frame) in enumerate(zip(owners, frames, strict=True)):
        if not active[i]:
            continue
        gateway = [ports[owner, n] % 12 for n in graph[owner]]
        local = frame.T @ single @ frame + frame[gateway].T @ frame[gateway] / 2
        ops = [eye]*3
        ops[i] = local
        separate += tensor(ops)
    for i, j in combinations(range(3), 2):
        a, b = owners[i], owners[j]
        if active[i] and active[j] and graph.has_edge(a, b):
            u, v = frames[i][ports[a, b] % 12], frames[j][ports[b, a] % 12]
            ops = [eye]*3
            ops[i], ops[j] = np.outer(u, u), np.outer(v, v)
            target = separate if i == 1 else contact
            target -= (1 + 2*kappa) * tensor(ops)
    # 18*kappa is the same identity in all compared three-load operators.
    return separate, contact


def antitriplet_embedding():
    wedge = np.zeros((9, 3))
    for k, (b, c) in enumerate(combinations(range(3), 2)):
        wedge[3*b+c, k], wedge[3*c+b, k] = 1/np.sqrt(2), -1/np.sqrt(2)
    return np.kron(np.eye(3), wedge)


def spectrum(graph, ports, owners, frames, epsilon, kappa):
    sep, contact = fixed_source_matrices(graph, ports, owners, frames, kappa)
    embed = antitriplet_embedding()
    base = embed.T @ sep @ embed
    interaction = embed.T @ contact @ embed
    all_colours = np.eye(9)
    singlet = (embed.T @ epsilon.ravel()).reshape(9, 1)
    octet = np.linalg.qr(singlet, mode='complete')[0][:, 1:]
    result = {}
    for name, q in (('all', all_colours), ('singlet', singlet), ('octet', octet)):
        projected = q.T @ interaction @ q
        free = np.linalg.eigvalsh(q.T @ base @ q)
        joint = np.linalg.eigvalsh(q.T @ (base + interaction) @ q)
        result[name] = dict(operator_norm=np.linalg.norm(projected, 2),
                            ground=joint[0], separate_ground=free[0],
                            ground_shift=joint[0]-free[0])
    return result
