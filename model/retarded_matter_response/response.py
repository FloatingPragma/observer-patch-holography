"""Exact source CAR resolvent dilations over Q(sqrt3).

All state amplitudes are inherited from the proved original cube ground.
The small spaces below are invariant hidden-spin/number sectors of the
complete original action, not projected approximations to its evolution.
No floating-point result admits any identity in certificate().
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache

import numpy as np
from scipy import sparse

from ..cube_eta_protection import cube
from ..interacting_vacuum import charged, ground


LIMIT = np.iinfo(np.int64).max
DENOMINATOR = 36
GAP = F(7, 1024)


def _maxabs(value):
    return max(abs(int(np.min(value, initial=0))), abs(int(np.max(value, initial=0))))


def _matmul(a, b):
    """Integer product with an explicit preflight against int64 overflow."""
    a, b = np.asarray(a), np.asarray(b)
    if a.dtype != np.int64 or b.dtype != np.int64:
        raise ValueError('int64 source arrays required')
    bound = a.shape[-1]*_maxabs(a)*_maxabs(b)
    if bound > LIMIT: raise ArithmeticError('exact matrix product would overflow')
    return a@b


def gram(left, right):
    """Numerators of (A+sqrt3 B)^T(C+sqrt3 D)."""
    a, b = left; c, d = right
    ac, bd = _matmul(a.T, c), _matmul(b.T, d)
    ad, bc = _matmul(a.T, d), _matmul(b.T, c)
    if any(_maxabs(x) > LIMIT//4 for x in (ac, bd, ad, bc)):
        raise ArithmeticError('exact quadratic-field addition would overflow')
    return ac+3*bd, ad+bc


def encode(pair, denominator=1, factor=1):
    """JSON-safe exact entries [rational coefficient, sqrt3 coefficient]."""
    a, b = pair
    return [[[str(F(factor*int(x), denominator)), str(F(factor*int(y), denominator))]
             for x, y in zip(rowa, rowb, strict=True)]
            for rowa, rowb in zip(a, b, strict=True)]


@lru_cache(maxsize=None)
def sector_masks(up, down):
    if type(up) is not int or type(down) is not int or not 0 <= up <= 8 or not 0 <= down <= 8:
        raise ValueError('hidden spin occupations0..8 required')
    masks = np.array([int(mask) for mask in cube.weights(up+down)
                      if (int(mask)&0x5555).bit_count() == up], dtype=np.int64)
    masks.setflags(write=False)
    return masks


def source_block(up, down, *, interaction=True):
    """Restriction of the actual full-Fock source; U0 is only a comparator."""
    if type(interaction) is not bool: raise ValueError('Boolean interaction flag required')
    masks = sector_masks(up, down)
    source = cube.operators().hubbard_H
    if source.imag.nnz or cube.operators().scale != 16:
        raise ArithmeticError('real hidden-frame source scale changed')
    h = source.real[masks][:, masks].tocsr()
    if source.real[:, masks].nnz != h.nnz or (h-h.T).nnz:
        raise ArithmeticError('selected hidden-number space is not source invariant')
    if not interaction:
        pairs = np.array([sum(((int(mask)>>(2*x))&3) == 3 for x in range(8))
                          for mask in masks], dtype=np.int64)
        h = h-16*sparse.diags(pairs, dtype=np.int64)
    return h.tocsr()


@dataclass(frozen=True)
class Dilation:
    h16: sparse.csr_matrix
    sqrt3_shift: np.ndarray
    vectors: tuple
    norm_square: int
    reference: str
    spin: str
    interaction: bool
    addition_dimension: int
    removal_dimension: int

    def apply16(self, pair):
        """16 L applied exactly to numerator columns."""
        a, b = pair
        row_bound = int(np.max(np.asarray(abs(self.h16).sum(axis=1)), initial=0))
        magnitude = max(_maxabs(a), _maxabs(b))
        if (row_bound+3*_maxabs(self.sqrt3_shift))*magnitude > LIMIT:
            raise ArithmeticError('source polynomial would overflow')
        shift = self.sqrt3_shift[:, None]
        return self.h16@a+3*shift*b, self.h16@b+shift*a

    def moment(self, order):
        if type(order) is not int or not 0 <= order <= 4:
            raise ValueError('retained exact moment order0..4 required')
        value = self.vectors
        for _ in range(order): value = self.apply16(value)
        return encode(gram(self.vectors, value), DENOMINATOR**2*16**order, self.norm_square)


@lru_cache(maxsize=None)
def dilation(reference='ground', spin='down', *, interaction=True):
    """CAR addition plus conjugate-removal dilation, all eight site columns.

    Everything is real in the inherited hidden frame. Thus the conjugate
    removal space has the same arrays; its generator has the opposite sign.
    For the charged reference the arrays omit sqrt2 and norm_square=2
    retains that exact normalization in every Gram matrix.
    """
    if reference not in ('ground', 'charged') or spin not in ('up', 'down'):
        raise ValueError('ground/charged reference and up/down channel required')
    if type(interaction) is not bool: raise ValueError('Boolean interaction flag required')
    state = charged.reference_vector()
    up, coefficient, factor = 4, 4, 1
    if reference == 'charged':
        state = charged.car(state, 0, True)
        up, coefficient, factor = 5, 3, 2
    if charged.inner(state, state) != (1296//factor, 0):
        raise ArithmeticError('exact reference normalization changed')
    # Source stationary reference, without a chemical-potential recentering.
    h = cube.operators().hubbard_H.real
    if not (np.array_equal(h@state[0], -3*coefficient*state[1]) and
            np.array_equal(h@state[1], -coefficient*state[0])):
        raise ArithmeticError('reference is not an eigenvector of the original H')
    mode_spin = int(spin == 'down')
    add_masks = sector_masks(up+(spin == 'up'), int(spin == 'down'))
    additions = [charged.car(state, 2*x+mode_spin, True) for x in range(8)]
    add = tuple(np.stack([v[k][add_masks] for v in additions], axis=1) for k in range(2))
    add_h = source_block(up+(spin == 'up'), int(spin == 'down'), interaction=interaction)
    if spin == 'up':
        rem_masks = sector_masks(up-1, 0)
        removals = [charged.car(state, 2*x, False) for x in range(8)]
        rem = tuple(np.stack([v[k][rem_masks] for v in removals], axis=1) for k in range(2))
        values = tuple(np.vstack([a, b]) for a, b in zip(add, rem, strict=True))
        h16 = sparse.block_diag([add_h, -source_block(up-1, 0, interaction=interaction)], format='csr')
        shift = np.array([coefficient]*len(add_masks)+[-coefficient]*len(rem_masks), dtype=np.int64)
        rem_dimension = len(rem_masks)
    else:
        if any(any(np.any(a) for a in charged.car(state, 2*x+1, False)) for x in range(8)):
            raise ArithmeticError('down removal did not vanish')
        values, h16 = add, add_h
        shift = np.full(len(add_masks), coefficient, dtype=np.int64)
        rem_dimension = 0
    return Dilation(h16, shift, values, factor, reference, spin, interaction,
                    len(add_masks), rem_dimension)


def first_numerator(reference, spin, interaction=True):
    """M1 as an integer matrix divided by48, independently of its moment."""
    a = ground.band_data()[0]
    result = 3*a.copy()
    if spin == 'down' and interaction:
        density48 = [24]*8 if reference == 'ground' else [48, 32, 32, 24, 32, 24, 24, 24]
        result += np.diag(np.array(density48, dtype=np.int64))
    return result


def memory_columns(data):
    """W=L V-V M1, numerator denominator1728 (times reference sqrt factor)."""
    first = first_numerator(data.reference, data.spin, data.interaction)
    kv = data.apply16(data.vectors)
    return tuple(3*x-_matmul(v, first) for x, v in zip(kv, data.vectors, strict=True))


def channel_certificate(reference='ground', spin='down', *, interaction=True):
    data = dilation(reference, spin, interaction=interaction)
    a, b = gram(data.vectors, data.vectors)
    if not np.array_equal(data.norm_square*a, 1296*np.eye(8, dtype=np.int64)) or np.any(b):
        raise ArithmeticError('complete addition/removal CAR isometry failed')
    first = first_numerator(reference, spin, interaction)
    first_actual = gram(data.vectors, data.apply16(data.vectors))
    if (not np.array_equal(3*data.norm_square*first_actual[0], 1296*first)
            or np.any(first_actual[1])):
        raise ArithmeticError('source first retarded moment disagrees with density')
    w = memory_columns(data)
    if any(np.any(x) for x in gram(data.vectors, w)):
        raise ArithmeticError('memory columns are not orthogonal to observed columns')
    wg = gram(w, w)
    expected = ([F(0)]*8 if spin == 'up' or not interaction else
                [F(1,4)]*8 if reference == 'ground' else
                [F(0), F(2,9), F(2,9), F(1,4), F(2,9), F(1,4), F(1,4), F(1,4)])
    for i in range(8):
        for j in range(8):
            if (F(data.norm_square*int(wg[0][i,j]), 1728**2) != (expected[i] if i == j else 0)
                    or wg[1][i,j] != 0):
                raise ArithmeticError('source memory Gram matrix differs from its exact value')
    if (spin == 'up' or not interaction) and any(np.any(x) for x in w):
        raise ArithmeticError('claimed free response lacks full cyclic closure')
    return {'reference': reference, 'channel': spin, 'source_U': int(interaction),
            'addition_dimension':data.addition_dimension, 'conjugate_removal_dimension':data.removal_dimension,
            'observed_columns':8, 'CAR_isometry_exact':True,
            'reference_normalization_square_factor':data.norm_square,
            'moments':{str(n):data.moment(n) for n in range(5)},
            'memory_Gram_diagonal':[str(x) for x in expected],
            'memory_Gram_offdiagonal_exactly_zero':True,
            'memory_columns_orthogonal_to_observed_columns':True,
            'free_resolvent_by_exact_full_cyclic_closure':spin == 'up' or not interaction,
            'positive_gap_domain': '7/1024 < K < 2' if reference == 'ground' and spin == 'down' and interaction else None,
            'projection_applied_to_source_dynamics':False}


def certificate():
    ground.certificate()
    channels = [channel_certificate(reference, spin) for reference in ('ground','charged') for spin in ('up','down')]
    control = channel_certificate('ground', 'down', interaction=False)
    return {'source_geometry':[2,2,2], 'source_c':'1/8', 'source_U':'1',
            'source_full_Fock_dimension':65536, 'original_H_without_mu_or_scalar_shift':True,
            'ground_energy':'-sqrt3/4', 'charged_reference':'sqrt2 c_0up^dagger G',
            'charged_reference_energy':'-3sqrt3/16',
            'hidden_frame_unitary':'existing onsite Fock frame; physical response is its sitewise unitary conjugate',
            'retarded_definition':'-i integral_0^infinity exp(i z t)<{c_x(t),c_y^dagger}>dt; Imz>0',
            'dilation':'V^*(z-L)^(-1)V; L=(H_(N+1)-E) directsum -conjugate(H_(N-1)-E)',
            'Schur_identity':'G(z)=[zI-M1-W^*(z-Q L Q|RanQ)^(-1)W]^(-1)',
            'channels':channels, 'free_down_control':control,
            'excited_density_increment':['1/2','1/6','1/6','0','1/6','0','0','0'],
            'same_source_state_dependent_response_proved':True,
            'nonzero_dynamic_memory_proved':True,
            'common_metric_or_Einstein_premise_completed':False,
            'native_ground_preparation_or_response_measurement_executed':False,
            'physical_particle_assignments_added':0}
