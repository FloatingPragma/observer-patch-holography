"""Full sparse finite carrier; numerical controls accompany an exact proof.

The 45-dimensional internal carrier and its three family copies are declared
design inputs. No change to the two-spin radial Hamiltonian is implicit.
"""
from functools import lru_cache

import numpy as np
from scipy import sparse as sp
import sympy as s

from model.common_source.chiral_candidate import action
from model.common_source import native_matter_attachment as native
from model.radial_defect_gap.finite import kernel_frame

J = np.array([[0., 1.], [-1., 0.]])


def norm(value):
    return float(sp.linalg.norm(value) if sp.issparse(value) else np.linalg.norm(value))


def native_kernel(mass_over_v, *, quark_design=True):
    """Actual neutral-Higgs tensors, with masses in the declared v(P) unit.

    Real symmetric lepton M equals Ye/sqrt(2). The optional quark matrices
    are the declared Stage5 epsilon=1/6 design, without measured mass inputs.
    The three neutrino directions are retained with exactly zero coupling.
    """
    mass = np.asarray(mass_over_v, dtype=complex)
    if mass.shape != (3, 3) or not np.isfinite(mass).all():
        raise ValueError('finite three-family matrix required')
    if norm(mass-mass.T) > 1e-13 or norm(mass.imag) > 1e-13:
        raise ValueError('this neutral-tensor assembly requires real symmetric M')
    families = {'charged_lepton': mass,
                'up': np.diag([6.**-6, 6.**-3, 1.])/np.sqrt(2),
                'down': np.diag([6.**-6, 6.**-4, 6.**-2])/np.sqrt(2)}
    if not quark_design:
        families['up'] *= 0
        families['down'] *= 0
    channels = action._representation()[4]
    result = sum((sp.kron(sp.csr_matrix(np.array(channels[name]['components'][1], complex)),
                          sp.csr_matrix(matrix), format='csr')
                  for name, matrix in families.items()), sp.csr_matrix((45, 45), dtype=complex))
    if norm(result-result.T) > 1e-13:
        raise ArithmeticError('actual all-left kernel must be symmetric')
    return result


@lru_cache(maxsize=4)
def radial(radius=1, mass=64):
    value = kernel_frame(radius, mass=mass)
    h, w = value['H'].tocsr(), sp.csr_matrix(value['frame'])
    signs = [(-1.)**sum(point) for point in value['points']]
    c = sp.kron(sp.diags(signs), sp.csr_matrix(J), format='csr')
    return h, w, c, value


def assemble(kernel, *, radius=1, mass=64, staggered=True):
    """Retain every physical annihilator and both redundant Nambu blocks."""
    k = sp.csr_matrix(kernel, dtype=complex)
    if k.shape[0] != k.shape[1] or not np.isfinite(k.data).all() or norm(k-k.T) > 1e-12:
        raise ValueError('finite symmetric internal creation kernel required')
    h, w, c, value = radial(radius, mass)
    if not staggered:
        c = sp.kron(sp.eye(len(value['points'])), sp.csr_matrix(J), format='csr')
    internal = sp.eye(k.shape[0], format='csr')
    kinetic = sp.kron(h, internal, format='csr')
    delta = sp.kron(c, k.conjugate(), format='csr')
    bdg = sp.bmat([[kinetic, delta], [delta.conj().T, -kinetic.conjugate()]], format='csr')
    e = sp.kron(w, internal, format='csr')
    inclusion = sp.block_diag((e, e.conjugate()), format='csr')
    core_delta = sp.kron(sp.csr_matrix(J), k.conjugate(), format='csr')
    core = sp.bmat([[None, core_delta], [core_delta.conj().T, None]], format='csr')
    hole_frame = sp.block_diag((sp.eye(kinetic.shape[0]), sp.kron(c.conj().T, internal)), format='csr')
    additive = sp.bmat([[kinetic, sp.kron(sp.eye(h.shape[0]), k.conjugate())],
                       [sp.kron(sp.eye(h.shape[0]), k), kinetic]], format='csr')
    return {'h': h, 'W': w, 'C': c, 'K': k, 'kinetic': kinetic,
            'delta': delta, 'bdg': bdg, 'inclusion': inclusion, 'core': core,
            'hole_frame': hole_frame, 'additive': additive, 'source': value}


def current_errors(value):
    """Anti-Hermitian generators use diag(G,conj G) in Nambu space."""
    if value['K'].shape != (45, 45):
        raise ValueError('actual native45 carrier required')
    errors = []
    for exact in action.generators():
        g = sp.csr_matrix(np.array(exact, complex))
        full = sp.kron(sp.eye(value['h'].shape[0]), g)
        small = sp.kron(sp.eye(2), g)
        full = sp.block_diag((full, full.conjugate()), format='csr')
        small = sp.block_diag((small, small.conjugate()), format='csr')
        errors.append(norm(full@value['inclusion']-value['inclusion']@small))
    p15, _ = native.action_intertwiner()
    p45 = sp.csr_matrix(np.array(s.kronecker_product(p15, s.eye(3)), complex))
    representation = native.representation()
    raw = []
    for port, exact in enumerate(representation['matter_currents']):
        source = sp.csr_matrix(np.array(s.kronecker_product(exact, s.eye(3)), complex))
        target = sp.csr_matrix(np.array(action.generators()[representation['native']['native_to_fixture_port'][port]], complex))
        raw.append(norm(target@p45-p45@source))
    return errors, raw


def residuals(value):
    h, w, c = (value[key] for key in ('h', 'W', 'C'))
    e, b, core, t = (value[key] for key in ('inclusion', 'bdg', 'core', 'hole_frame'))
    return {'radial_particle_hole': norm(h@c+c@h.conjugate()),
            'radial_kernel': norm(h@w),
            'kernel_particle_hole': norm(c@w.conjugate()-w@J),
            'symmetric_internal_kernel': norm(value['K']-value['K'].T),
            'antisymmetric_creation_pair': norm(value['delta']+value['delta'].T),
            'Hermitian_Nambu': norm(b-b.conj().T),
            'isometric_inclusion': norm(e.conj().T@e-sp.eye(e.shape[1])),
            'complete_Nambu_intertwiner': norm(b@e-e@core),
            'complete_additive_form': norm(t.conj().T@b@t-value['additive'])}
