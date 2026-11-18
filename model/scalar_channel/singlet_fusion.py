"""Executable bounded-star singlet compiler; no dense 7**r state space.

Each actual seven-state owner has vacuum0, color1:4, weak Cartesian4:7.
The neutral subrepresentation 0+1 is exactly two auxiliary spin-1/2s.
A Dyck fusion path of length2r specifies one orthonormal singlet MPS.
The first2**r paths, in down-before-up lexicographic order, are the code.
Rank0 is exactly all-owner vacuum. Arbitrary logical gates act on path ranks;
the Schur/Fourier transform is specified below by every tensor coefficient.
Unused singlets, nonsinglets and color states are retained as identity
complement by the lifted-gate definition, never silently discarded.
"""
from .carrier import require
from functools import lru_cache
from math import sqrt, comb
import json
import numpy as np
import scipy.linalg as la

def _integer(value, name, minimum=0):
    require(isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_))
            and value >= minimum, name + ' must be an integer >= ' + str(minimum))
    return int(value)


def _owners(owners, count=None):
    values = tuple(_integer(x, 'owner') for x in owners)
    require(all(x < 7 for x in values), 'owner outside seven-state register')
    require(count is None or len(values) == count, 'complete owner word required')
    return values


def _mps(state, count=None):
    require(len(state) > 0 and (count is None or len(state) == count), 'complete nonempty MPS required')
    bond = 1
    for tensor in state:
        require(isinstance(tensor, np.ndarray) and tensor.ndim == 3
                and tensor.shape[0] == 7 and tensor.shape[1] == bond
                and tensor.shape[2] > 0 and np.isfinite(tensor).all(), 'invalid physical MPS tensor')
        bond = tensor.shape[2]
    require(bond == 1, 'MPS must have closed boundary bonds')


@lru_cache(None, typed=True)
def completions(left, height):
    left = _integer(left, 'remaining steps')
    require(isinstance(height, (int, np.integer)) and not isinstance(height, (bool, np.bool_)),
            'height must be an integer')
    if left == 0:
        return int(height == 0)
    if height < 0 or height > left or (left - height) % 2:
        return 0
    return (completions(left - 1, height - 1) if height else 0) + completions(left - 1, height + 1)

def unrank(r, rank):
    r = _integer(r, 'owner count', 1)
    rank = _integer(rank, 'singlet rank')
    if not 0 <= rank < completions(2 * r, 0):
        raise ValueError('Singlet rank outside full Catalan space')
    h = 0
    path = [h]
    for left in range(2 * r, 0, -1):
        ndown = completions(left - 1, h - 1) if h else 0
        if rank < ndown:
            h -= 1
        else:
            rank -= ndown
            h += 1
        path.append(h)
    require(h == 0 and rank == 0, 'Compiler invariant: h == 0 and rank == 0')
    return tuple(path)

def rank_path(path):
    path = tuple(_integer(h, 'fusion height') for h in path)
    require(len(path) >= 3 and len(path) % 2 == 1 and path[0] == path[-1] == 0,
            'fusion path must start and finish at zero with even step count')
    r = (len(path) - 1) // 2
    rank = 0
    for k, (h, hnew) in enumerate(zip(path, path[1:])):
        h, hnew = (int(h), int(hnew))
        left = 2 * r - k
        if hnew == h + 1 and h:
            rank += completions(left - 1, h - 1)
        elif hnew != h - 1 and hnew != h + 1:
            raise ValueError('Invalid fusion step')
    return rank

def cg_tensor(h, hnew):
    """A[s, old-m-index,new-m-index], doubled spins h=2j, m=2m_j."""
    h = _integer(h, 'input spin')
    hnew = _integer(hnew, 'output spin')
    require(abs(hnew - h) == 1, 'not a spin-half fusion step')
    out = np.zeros((2, h + 1, hnew + 1), complex)
    for i, m in enumerate(range(-h, h + 1, 2)):
        for s, sigma in enumerate((-1, 1)):
            newm = m + sigma
            if abs(newm) > hnew:
                continue
            j = (newm + hnew) // 2
            if hnew == h + 1:
                z = sqrt((h + sigma * m + 2) / (2 * (h + 1)))
            elif hnew == h - 1:
                z = (-1 if sigma == 1 else 1) * sqrt((h - sigma * m) / (2 * (h + 1)))
            else:
                raise ValueError('Not a spin-half fusion step')
            out[s, i, j] = z
    return out

def owner_pair_isometry():
    """Q[down/up pair index, physical owner]. Q†Q=P_(0+weak).

    down=0,up=1. Vacuum=(up down-down up)/sqrt2. Weak xyz are
    Cartesian triplets. This is the inverse physical pairing transform.
    """
    Q = np.zeros((4, 7), complex)
    Q[2, 0] = 1 / sqrt(2)
    Q[1, 0] = -1 / sqrt(2)
    Q[0, 4] = 1 / sqrt(2)
    Q[3, 4] = -1 / sqrt(2)
    Q[0, 5] = 1j / sqrt(2)
    Q[3, 5] = 1j / sqrt(2)
    Q[1, 6] = Q[2, 6] = 1 / sqrt(2)
    return Q

@lru_cache(None, typed=True)
def mps_for_rank(r, rank):
    r = _integer(r, 'owner count', 1)
    rank = _integer(rank, 'singlet rank')
    path = unrank(r, rank)
    Q = owner_pair_isometry()
    out = []
    for k in range(r):
        A = cg_tensor(path[2 * k], path[2 * k + 1])
        B = cg_tensor(path[2 * k + 1], path[2 * k + 2])
        virtual = np.einsum('sim,tmj->stij', A, B).reshape(4, A.shape[1], B.shape[2])
        tensor = np.einsum('sa,sij->aij', Q.conj(), virtual)
        tensor.flags.writeable = False
        out.append(tensor)
    return tuple(out)

def overlap(a, b, operators=None):
    _mps(a)
    _mps(b, len(a))
    if operators is not None:
        require(isinstance(operators, dict), 'operators must map owner indices to matrices')
        for site, operator in operators.items():
            require(_integer(site, 'operator site') < len(a), 'operator site outside star')
            require(np.shape(operator) == (7, 7) and np.isfinite(operator).all(), 'invalid owner operator')
    E = np.ones((1, 1), complex)
    for k, (A, B) in enumerate(zip(a, b)):
        if operators is None or k not in operators:
            E = np.einsum('ij,aik,ajl->kl', E, A.conj(), B, optimize=True)
        else:
            E = np.einsum('ij,aik,ab,bjl->kl', E, A.conj(), operators[k], B, optimize=True)
    return E[0, 0]

def product_amplitude(r, rank, owners):
    owners = _owners(owners, r)
    v = np.ones(1, complex)
    for A, x in zip(mps_for_rank(r, rank), owners):
        v = v @ A[x]
    return v[0]

def logical_pauli(vector, xmask=0, zmask=0, phase=0):
    """Exact path-rank representation of i**phase X**xmask Z**zmask.

    Physical operation is VUV†+I-P. This only returns its code block;
    the complete complement is explicitly identity in the specification.
    """
    v = np.asarray(vector, complex)
    require(v.ndim == 1 and len(v) > 0 and len(v) & (len(v) - 1) == 0
            and np.isfinite(v).all(), 'logical vector must have finite power-of-two dimension')
    xmask = _integer(xmask, 'X mask')
    zmask = _integer(zmask, 'Z mask')
    phase = _integer(phase, 'Pauli phase')
    require(xmask < len(v) and zmask < len(v), 'Pauli mask outside logical vector')
    out = np.empty_like(v)
    for s, z in enumerate(v):
        out[s ^ xmask] = 1j ** phase * (-1) ** (s & zmask).bit_count() * z
    return out

def apply_complete_gate_to_mps_sum(r, terms, logical_action):
    """Full raw-domain action, represented as a finite linear sum of MPSs.

    Input terms are (amplitude,MPS) and may contain color or dirty spin data.
    The exact identity U_raw ψ=ψ+V(U-I)V†ψ retains the entire complement.
    This implementation never projects away, renormalizes or postselects it.
    No coefficient threshold is used. Runtime scales with2**r, a fixed local
    degree bound; it is not a complexity claim for arbitrary growing r.
    The supplied callback must be a linear unitary; output shape, finiteness
    and this input's norm are checked, without certifying arbitrary callbacks.
    """
    r = _integer(r, 'owner count', 1)
    require(completions(2 * r, 0) >= 2 ** r, 'star has insufficient singlets for r logical bits')
    terms = list(terms)
    require(len(terms) > 0 and callable(logical_action), 'nonempty state and logical action required')
    for amplitude, state in terms:
        require(np.isscalar(amplitude) and np.isfinite(amplitude), 'finite MPS amplitude required')
        _mps(state, r)
    code = [mps_for_rank(r, k) for k in range(2 ** r)]
    coeff = np.array([sum((a * overlap(v, p) for a, p in terms)) for v in code])
    acted = np.asarray(logical_action(coeff.copy()), complex)
    require(acted.shape == coeff.shape and np.isfinite(acted).all(), 'invalid logical action output')
    require(abs(np.vdot(acted, acted).real - np.vdot(coeff, coeff).real)
            <= 1e-11 * max(1.0, float(np.vdot(coeff, coeff).real)), 'logical action does not preserve input norm')
    delta = acted - coeff
    return ([*terms, *[(a, v) for a, v in zip(delta, code) if a != 0]], coeff, delta)

def product_mps(owners):
    owners = _owners(owners)
    require(len(owners) > 0, 'nonempty owner word required')
    out = []
    for x in owners:
        A = np.zeros((7, 1, 1), complex)
        A[int(x), 0, 0] = 1
        out.append(A)
    return tuple(out)

def charted_product_amplitude(r, rank, owners, chart):
    owners = _owners(owners, r)
    chart = tuple(_integer(j, 'chart index') for j in chart)
    require(sorted(chart) == list(range(r)), 'chart must be an owner permutation')
    return product_amplitude(r, rank, [owners[j] for j in chart])

def schur_orthogonality_local():
    worst = 0.0
    for h in range(31):
        cols = []
        for hp in (h - 1, h + 1) if h else (1,):
            A = cg_tensor(h, hp)
            cols.extend(A.transpose(1, 0, 2).reshape(2 * (h + 1), hp + 1).T)
        U = np.array(cols).T
        require(U.shape == (2 * (h + 1), 2 * (h + 1)), 'Compiler invariant: U.shape == (2 * (h + 1), 2 * (h + 1))')
        worst = max(worst, float(la.norm(U.conj().T @ U - np.eye(len(U)), 2)))
    return worst

def run():
    rng = np.random.default_rng(191830)
    Q = owner_pair_isometry()
    I = np.eye(7)

    def rot(axis, angle):
        x, y, z = axis / np.linalg.norm(axis)
        R = la.expm(angle * np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]]))
        return la.block_diag(np.ones((1, 1)), np.eye(3), R)
    R = rot(np.array([0.4, 0.6, -0.8]), 0.73)
    all_results = []
    for r in (10, 15):
        total = completions(2 * r, 0)
        dim = 2 ** r
        paths = np.array([unrank(r, k) for k in range(dim)], dtype=np.uint8)
        require(len({tuple(p) for p in paths}) == dim, 'Compiler invariant: len({tuple(p) for p in paths}) == dim')
        require(all((rank_path(p) == i for i, p in enumerate(paths))), 'Compiler invariant: all((rank_path(p) == i for i, p in enumerate(paths)))')
        test_indices = sorted(set([0, 1, 2, dim - 1, *rng.choice(dim, 30, replace=False).tolist()]))
        states = [mps_for_rank(r, k) for k in test_indices]
        gram = np.array([[overlap(a, b) for b in states] for a in states])
        collective = max((abs(overlap(a, a, {j: R for j in range(r)}) - 1) for a in states))
        independent = [float(np.sqrt(max(0, 2 - 2 * overlap(a, a, {j: R}).real))) for a in states for j in range(r)]
        psi = rng.normal(size=dim) + 1j * rng.normal(size=dim)
        psi /= la.norm(psi)
        pauli_error = 0.0
        for bit in range(r):
            x = logical_pauli(psi, xmask=1 << bit)
            z = logical_pauli(psi, zmask=1 << bit)
            pauli_error = max(
                pauli_error,
                float(la.norm(logical_pauli(x, xmask=1 << bit) - psi)),
                float(la.norm(logical_pauli(z, zmask=1 << bit) - psi)),
                float(la.norm(logical_pauli(x, zmask=1 << bit) + logical_pauli(z, xmask=1 << bit))),
            )
        vacuum = product_amplitude(r, 0, [0] * r)
        result = dict(
            r=r,
            code_dimension=dim,
            total_singlet_dimension=total,
            full_owner_dimension=7 ** r,
            all_code_paths_generated=True,
            all_rank_unrank_checks_passed=True,
            largest_code_MPS_bond_dimension=int(paths.max() + 1),
            all_vacuum_codeword_amplitude_real=float(vacuum.real),
            selected_code_Gram_error=float(la.norm(gram - np.eye(len(states)), 2)),
            sampled_codewords=len(states),
            collective_SO3_rotation_error=float(collective),
            largest_independent_owner_rotation_change=max(independent),
            logical_Pauli_involution_and_anticommutation_error=pauli_error,
            local_Clebsch_transform_full_unitarity_error=schur_orthogonality_local(),
            dense_7_power_r_matrix_materialized=False,
            actual_full_dense_gate_synthesis_executed=False,
            executable_gate_representation='Full Clebsch decoder, path-rank operation, inverse decoder; identity on complement.',
            raw_permutation_equivariance_claimed=False,
        )
        require(result['selected_code_Gram_error'] < 1e-12 and collective < 1e-12, "Compiler invariant: result['selected_code_Gram_error'] < 1e-12 and collective < 1e-12")
        require(abs(vacuum - 1) < 1e-12 and max(independent) > 0.1, 'Compiler invariant: abs(vacuum - 1) < 1e-12 and max(independent) > 0.1')
        require(pauli_error < 1e-12, 'Compiler invariant: pauli_error < 1e-12')
        all_results.append(result)
    Jz = np.zeros((7, 7), complex)
    Jz[4, 5] = -1j
    Jz[5, 4] = 1j
    charge_state = np.zeros(7, complex)
    charge_state[4] = 1
    charge_positive = float(np.vdot(charge_state, Jz @ Jz @ charge_state).real)
    q2_errors = []
    for r, k in ((10, 0), (10, 1023), (15, 32767)):
        a = mps_for_rank(r, k)
        val = 0j
        for j in range(r):
            val += overlap(a, a, {j: Jz @ Jz})
            for l in range(j + 1, r):
                val += 2 * overlap(a, a, {j: Jz, l: Jz})
        q2_errors.append(float(abs(val)))
    dirty = product_mps([1] + [0] * 9)
    raw = [(1 / sqrt(3), mps_for_rank(10, 1)), (1 / sqrt(3), mps_for_rank(10, 2)), (1 / sqrt(3), dirty)]
    transformed, coeff, delta = apply_complete_gate_to_mps_sum(10, raw, lambda v: logical_pauli(v, xmask=1))
    outvec = logical_pauli(coeff, xmask=1)
    require(abs(outvec[0] - 1 / sqrt(3)) < 1e-12 and abs(outvec[3] - 1 / sqrt(3)) < 1e-12, 'Compiler invariant: abs(outvec[0] - 1 / sqrt(3)) < 1e-12 and abs(outvec[3] - 1 / sqrt(3)) < 1e-12')
    dirty_amp = sum((a * overlap(dirty, p) for a, p in transformed))
    require(abs(dirty_amp - 1 / sqrt(3)) < 1e-12, 'Compiler invariant: abs(dirty_amp - 1 / sqrt(3)) < 1e-12')
    expected = np.zeros(1024, complex)
    expected[[0, 3]] = 1 / sqrt(3)
    oldowners = [4, 4] + [0] * 8
    chart = list(range(10))
    perm = [(j + 3) % 10 for j in range(10)]
    newowners = [0] * 10
    for j in range(10):
        newowners[perm[j]] = oldowners[j]
    chart_err = abs(
        charted_product_amplitude(10, 1, oldowners, chart) - charted_product_amplitude(10, 1, newowners, [perm[j] for j in chart]),
    )
    output = dict(
        schema='oph.native_singlet_fusion.v2',
        stars=all_results,
        physical_pair_transform_row_isometry_error=float(la.norm(Q @ Q.conj().T - np.eye(4), 2)),
        same_collective_Q_squared_reader=dict(
            code_errors=q2_errors,
            raw_single_weak_owner_positive_control=charge_positive,
            positive_control_is_outside_source_code=True,
            raw_single_owner_Q_squared_is_not_claimed_zero_on_entangled_code=True,
        ),
        complete_raw_domain_gate=dict(
            r=10,
            code_output_vector_error=float(la.norm(outvec - expected)),
            retained_colored_complement_amplitude_error=float(abs(dirty_amp - 1 / sqrt(3))),
            all_1024_code_overlap_coefficients_evaluated=True,
            code_probability=float(np.vdot(coeff, coeff).real),
            dirty_probability=1 / 3,
            discarded_complement=False,
            coefficient_threshold_used=False,
        ),
        passive_transported_chart_amplitude_error=float(chart_err),
        owner_initial_reset='Every owner vacuum; rank0 codeword is this product exactly.',
        gauge_condition='One collective center SO3 action on all local owners; not independent owner rotations.',
        chart_condition='Fusion ordering is declared; permutation charts require explicit code transport; raw A5 covariance not proved.',
    )
    require(max(q2_errors) < 1e-12 and charge_positive == 1, 'Compiler invariant: max(q2_errors) < 1e-12 and charge_positive == 1')
    return output
if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
