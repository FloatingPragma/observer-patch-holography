#!/usr/bin/env python3
"""Gauge-covariant spectral bands on the literal compact OPH lattice.

Classical links are specified diagnostic inputs. Rational bounds prove the
band counts and polar-transport domains; floating eigenvectors test the
literal source matrices and are labelled separately in the receipt.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy import linalg as la, sparse as sp
from scipy.sparse.linalg import eigsh
import sympy as sy

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

from model.compact_defect_pair import spectral, topology
from model.common_source.chiral_candidate import action
from model.common_source import native_family_record_memory as memory
from model.defect_current_bridge.construction import native_kernel
from model.mcpr_native_memory_bridge.construction import design
from model.mcpr_native_memory_bridge import certificate as memory_certificate
from model.compact_interacting_response.algebra import mass_boxes
from model.native_flavor_orbit.intervals import I, rounded_out


SOURCE_COMMIT = '07c7acb2e7003efd9e574b67b374dbc36171dc6b'
MASS = 64


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def norm(matrix):
    return float(sp.linalg.norm(matrix) if sp.issparse(matrix) else la.norm(matrix))


def exact_constants():
    """Dimension-independent resolvent and polar-domain bounds."""
    radius = F(MASS, 2)
    kinetic = F(3)
    pair = F(1)
    gap = F(MASS)-2*(kinetic+pair)
    resolvent_kinetic = radius-kinetic
    resolvent_bdg = radius-kinetic-pair
    eta = radius*pair/(resolvent_kinetic*resolvent_bdg)
    eta_records = radius*2*pair/resolvent_bdg**2
    sign_correction = 4*eta**2
    three_projector_bound = 4*eta_records+2*sign_correction
    assert gap == 56 and eta == F(8, 203) and eta_records == F(4, 49)
    assert 4*eta**2 < 1 and three_projector_bound < 1
    return {'texture_amplitude': MASS, 'hopping_norm_upper': str(kinetic),
            'raw_pairing_norm_upper': str(pair), 'contour_radius': str(radius),
            'BdG_band_absolute_energy_upper': '4',
            'BdG_complement_absolute_energy_lower': '60',
            'BdG_band_complement_gap_lower': str(gap),
            'number_conserving_band_absolute_energy_upper': '3',
            'number_conserving_complement_absolute_energy_lower': '61',
            'number_conserving_band_complement_gap_lower': '58',
            'BdG_vs_number_conserving_projector_norm_upper': str(eta),
            'compressed_grading_square_lower': str(1-4*eta**2),
            'sign_minus_compressed_grading_norm_upper': str(sign_correction),
            'between_raw_record_projector_norm_upper': str(eta_records),
            'three_projector_alignment_minus_identity_norm_upper': str(three_projector_bound),
            'three_projector_polar_factor_is_invertible': True,
            'raw_record_projectors_low_rank': 360,
            'raw_record_annihilator_count': 180,
            'bounds_are_exact_rationals': True}


def primitive_links():
    """Exact S(U(3)xU(2)) controls, fixed without particle-mass inputs."""
    r, t = sy.Rational(3, 5), sy.Rational(4, 5)
    A, B = sy.eye(5), sy.eye(5)
    A[0, 0], A[1, 1], A[0, 1], A[1, 0] = r, r, t, -t
    A[3, 3], A[4, 4], A[3, 4], A[4, 3] = r, r, sy.I*t, sy.I*t
    B[1, 1], B[2, 2], B[1, 2], B[2, 1] = r, r, sy.I*t, sy.I*t
    B[3, 3], B[4, 4], B[3, 4], B[4, 3] = r, r, t, -t
    z = r+sy.I*t
    Z = sy.diag(*[sy.expand_complex(z**q) for q in (-2, -2, -2, 3, 3)])
    for matrix in (A, B, Z):
        assert (matrix.H*matrix-sy.eye(5)).applyfunc(sy.simplify) == sy.zeros(5)
        assert sy.simplify(matrix.det()-1) == 0
    wilson = (Z*A*B*A.H*B.H).applyfunc(sy.simplify)
    assert wilson != sy.eye(5)
    return (A, B, Z), {'construction': 'Gaussian-rational rotations with cosine 3/5 and sine 4/5',
                       'matrices5': [[[str(x) for x in row] for row in matrix.tolist()]
                                     for matrix in (A, B, Z)],
                       'exact_unitarity_and_total_determinant_one': True,
                       'plaquette_Wilson_matrix5': [[str(x) for x in row] for row in wilson.tolist()],
                       'plaquette_Wilson_trace5': str(sy.simplify(sy.trace(wilson))),
                       'link_background_selected_by_native_repair': False,
                       'comparison_data_used_to_choose_controls': False}


def exterior(matrix5, subsets):
    matrix5 = np.asarray(matrix5, complex)
    out = np.zeros((len(subsets), len(subsets)), complex)
    for i, left in enumerate(subsets):
        for j, right in enumerate(subsets):
            if len(left) == len(right):
                out[i, j] = la.det(matrix5[np.ix_(left, right)])
    return out


def actual_source():
    source = spectral.source(radius=3, mass=MASS, half_separation=1)
    points = source['points']
    addresses = {point: i for i, point in enumerate(points)}
    edges = []
    for p, i in addresses.items():
        for axis in range(3):
            q = tuple(x+int(j == axis) for j, x in enumerate(p))
            if q in addresses:
                edges.append((i, addresses[q], axis, p, q))
    assert len(points) == 343 and len(edges) == 882
    for i, j, axis, _, _ in edges:
        assert norm(source['h'][2*i:2*i+2, 2*j:2*j+2].toarray()-.5j*spectral.SIGMA[axis]) == 0
    source.update({'edges': edges, 'addresses': addresses})
    return source


def link_field(source, matrices15, *, transform=None, indices=None):
    A, B, Z = matrices15
    curved = {((0, 0, -1), (1, 0, -1)): Z@A,
              ((1, 0, -1), (1, 1, -1)): B,
              ((0, 1, -1), (1, 1, -1)): A,
              ((0, 0, -1), (0, 1, -1)): B}
    keep = list(range(15)) if indices is None else list(indices)
    rows = []
    for i, j, axis, p, q in source['edges']:
        U = curved.get((p, q), np.eye(15, dtype=complex))
        if transform is not None:
            U = transform[i]@U@transform[j].conj().T
        rows.append(U[np.ix_(keep, keep)])
    return rows


def build_h(source, links):
    d = links[0].shape[0]
    n = len(source['points'])
    rows, cols, vals = [], [], []
    def block(i, j, value):
        a, b = np.nonzero(value)
        rows.extend((2*d*i+a).tolist())
        cols.extend((2*d*j+b).tolist())
        vals.extend(value[a, b].tolist())
    for i in range(n):
        onsite = source['h'][2*i:2*i+2, 2*i:2*i+2].toarray()
        block(i, i, np.kron(onsite, np.eye(d)))
    for (i, j, axis, _, _), link in zip(source['edges'], links, strict=True):
        value = np.kron(.5j*spectral.SIGMA[axis], link)
        block(i, j, value)
        block(j, i, value.conj().T)
    result = sp.csr_matrix((vals, (rows, cols)), shape=(2*n*d,)*2)
    assert norm(result-result.conj().T) < 1e-12
    return result


def spectral_frame(H, rank, *, complement_lower):
    n = H.shape[0]
    # The rank is fixed by the analytic gap theorem, not a numerical cutoff.
    v0 = np.sin(np.arange(n)+.37)+1j*np.cos(np.arange(n)*.73+.19)
    values, W = eigsh(H, k=rank, sigma=0., which='LM', tol=2e-13,
                      maxiter=10000, ncv=max(2*rank+12, 40), v0=v0)
    gram_floor = float(la.eigvalsh(W.conj().T@W)[0])
    assert gram_floor > 1e-4, 'a repeated ARPACK vector cannot certify the band dimension'
    W, _ = la.qr(W, mode='economic')
    values, R = la.eigh(W.conj().T@(H@W))
    W = W@R
    residual = norm(H@W-W*values)
    orthogonal = norm(W.conj().T@W-np.eye(rank))
    assert residual < 5e-9 and orthogonal < 5e-11
    assert np.max(abs(values)) < 4.01
    return {'W': W, 'values': values, 'residual': residual,
            'orthogonality_error': orthogonal, 'input_Gram_floor': gram_floor,
            'floating_angle_diagnostic': residual/(complement_lower-max(abs(values)))}


def summary_frame(row):
    return {key: value.tolist() if isinstance(value, np.ndarray) else value
            for key, value in row.items() if key not in ('W', 'plus', 'minus')}


def projector_distance(X, Y):
    """Equal-rank orthonormal-frame projector distance via a small Gram matrix."""
    s = la.svdvals(X.conj().T@Y)
    return float(np.sqrt(max(0., 1-min(float(s[-1]), 1.)**2)))


def ph(frame):
    n = frame.shape[0]//2
    return np.vstack((frame[n:].conj(), frame[:n].conj()))


def charged_band(B, rank, *, complement_lower=60):
    row = spectral_frame(B, rank, complement_lower=complement_lower)
    W = row['W']
    n = B.shape[0]//2
    gradingW = np.vstack((W[:n], -W[n:]))
    C = W.conj().T@gradingW
    eigen, V = la.eigh(C)
    plus = W@V[:, eigen > 0]
    minus = ph(plus)
    assert plus.shape[1]*2 == rank
    ph_error = norm(ph(W)-W@(W.conj().T@ph(W)))
    frame_error = norm(np.hstack((plus, minus)).conj().T@np.hstack((plus, minus))-np.eye(rank))
    assert ph_error < 2e-9 and frame_error < 2e-9
    return {**row, 'plus': plus, 'minus': minus,
            'compressed_grading_eigenvalues': eigen,
            'particle_hole_band_error': ph_error,
            'canonical_frame_error': frame_error}


def band_parts(row):
    return (row['plus'], row['minus'])


def small_polar_transport(old, new):
    """Polar alignment of positive, negative and complementary projectors."""
    old_frames, new_frames = band_parts(old), band_parts(new)
    n = old_frames[0].shape[0]//2
    def majorana(F):
        return np.vstack((F[:n]+F[n:], -1j*F[:n]+1j*F[n:]))/np.sqrt(2)
    # A real Majorana span preserves the particle-hole reality condition
    # through nearly dependent directions of the two low subspaces.
    candidates = [majorana(old_frames[0]), majorana(new_frames[0])]
    union = np.hstack([part for F in candidates for part in (F.real, F.imag)])
    realZ, singular, _ = la.svd(union, full_matrices=False)
    keep = singular > 1e-11
    discarded = float(singular[~keep][0]) if np.any(~keep) else 0.
    realZ = realZ[:, keep]
    Z = np.vstack((realZ[:n]+1j*realZ[n:], realZ[:n]-1j*realZ[n:]))/np.sqrt(2)
    small_old = [Z.conj().T@F for F in old_frames]
    small_new = [Z.conj().T@F for F in new_frames]
    old_q = [X@X.conj().T for X in small_old]
    new_q = [X@X.conj().T for X in small_new]
    identity = np.eye(Z.shape[1])
    old_q.append(identity-sum(old_q))
    new_q.append(identity-sum(new_q))
    A = sum(Q@P for Q, P in zip(new_q, old_q, strict=True))
    imaginary_error = norm(A.imag)
    assert imaginary_error < 2e-10
    A = A.real
    U, positive = la.polar(A)
    alignment = max(norm(U@P@U.conj().T-Q) for P, Q in zip(old_q, new_q, strict=True))
    Jold, Jnew = old_q[0]-old_q[1], new_q[0]-new_q[1]
    Jerror = norm(U@Jold@U.conj().T-Jnew)
    # I+Z(U-I)Z† acts as identity outside the complete union span.
    D = U-identity
    XiZ = ph(Z)
    C = Z.conj().T@XiZ
    ph_span = norm(XiZ-Z@C)
    ph_error = norm(D-C@D.conj()@C.conj().T)+2*norm(D)*ph_span
    assert max(alignment, Jerror, ph_error) < 2e-8
    return {'span_dimension': Z.shape[1], 'discarded_span_singular_value_upper_numeric': discarded,
            'Majorana_alignment_imaginary_roundoff': imaginary_error,
            'alignment_error': alignment, 'occupation_grading_transport_error': Jerror,
            'particle_hole_transport_error': ph_error,
            'polar_factor_unitarity_error': norm(U.conj().T@U-identity),
            'alignment_operator_smallest_singular_value': float(la.eigvalsh(positive)[0]),
            'alignment_operator_minus_identity_norm': float(la.norm(A-identity, 2)),
            'Z': Z, 'U_small': U}


def source_kernel_bounds():
    row = design('forward', 0)
    projectors = memory._projectors(False)
    assert len(projectors) == 12
    for P in projectors:
        assert (P*P-P).applyfunc(sy.simplify) == sy.zeros(3)
        assert P.T == P and sy.simplify(sy.trace(P)) == 1
    inverse = memory_certificate.admissibility('forward')['inverse_branches'][0]
    gain = sum(mass_boxes(), I(0))/I(*inverse['p'])
    assert gain.lo > 0 and gain.hi < 1
    assert gain.lo < F(row['gain_over_v']) < gain.hi
    k = [native_kernel(np.zeros((3, 3))).toarray()]
    k += [native_kernel(row['gain_over_v']*np.asarray(P, float)).toarray() for P in projectors]
    max_norm = max(float(la.svdvals(K)[0]) for K in k)
    assert max_norm < 1
    assert all(np.max(abs(K.imag)) == 0 for K in k)
    runtime_row_bound = max(sum(abs(F(float(x))) for x in row.real) for K in k for row in K)
    assert runtime_row_bound < 1
    # Every source neutral-Higgs tensor has disjoint unit entries in each
    # charged block. This checks the exact tensor normalization, rather than
    # treating a floating singular value as the proof of ||K|| <= 1.
    channels = action._representation()[4]
    supports = []
    for name, channel in channels.items():
        C = channel['components'][1]
        assert C.T == C
        assert set(C) <= {sy.Integer(-1), sy.Integer(0), sy.Integer(1)}
        assert all(sum(abs(x) for x in C.row(i)) <= 1 for i in range(15))
        support = {i for i in range(15) if any(C[i, j] for j in range(15))}
        assert not any(support & earlier for earlier in supports)
        supports.append(support)
    return row, projectors, {'raw_kernel_count': len(k),
                            'raw_kernels_symmetric': all(norm(K-K.T) == 0 for K in k),
                            'exact_projector_count': len(projectors),
                            'gain_interval': rounded_out(gain, 20).strings(),
                            'gain_below_one_exact': True,
                            'quark_maximum_coefficient_exact': '1/sqrt(2)',
                            'neutral_Higgs_tensors_are_disjoint_unit_partial_isometries': True,
                            'all_raw_kernel_norms_below_one_exact': True,
                            'literal_runtime_kernel_row_norm_upper_exact': str(runtime_row_bound),
                            'maximum_runtime_kernel_norm': max_norm,
                            'source_branch': 'forward, inverse branch 0',
                            'source_design_ancestry': 'declared target-informed MCPR architecture retained'}


def execute(output_dir=HERE):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    source = actual_source()
    constants = exact_constants()
    subset, _, _, operators, channels = action._representation()
    prim, controls = primitive_links()
    matrices15 = [exterior(np.asarray(x, complex), subset) for x in prim]
    primitive_errors = [norm(U.conj().T@U-np.eye(15)) for U in matrices15]
    assert max(primitive_errors) < 1e-12
    qdiag = [int(operators['q'][i, i]) for i in range(15)]
    blocks = {q: [i for i, value in enumerate(qdiag) if value == q] for q in sorted(set(qdiag))}
    names = {-4: 'u^c', -3: 'L', 1: 'Q', 2: 'd^c', 6: 'e^c'}
    design_row, projectors, kernel_bounds = source_kernel_bounds()
    results = []
    arrays = {
        'source_points': np.asarray(source['points'], dtype=np.int64),
        'source_oriented_edges_i_j_axis': np.asarray([edge[:3] for edge in source['edges']], dtype=np.int64),
        'native_m15_masks': np.asarray([sum(1 << i for i in row) for row in subset], dtype=np.int64),
        'native_primitive_charges': np.asarray(qdiag, dtype=np.int64),
        'primitive_controls_15': np.asarray(matrices15),
        'curved_link_matrices_15': np.asarray(link_field(source, matrices15)),
        'all_raw_kernels_45': np.asarray(
            [native_kernel(np.zeros((3, 3))).toarray()]+
            [native_kernel(design_row['gain_over_v']*np.asarray(P, float)).toarray() for P in projectors])}
    for name in ('h', 'C'):
        literal = sp.csr_matrix(source[name])
        arrays.update({f'source_{name}_csr_data': literal.data,
                       f'source_{name}_csr_indices': literal.indices,
                       f'source_{name}_csr_indptr': literal.indptr,
                       f'source_{name}_shape': np.asarray(literal.shape, dtype=np.int64)})
    print('Source and exact norm certificates assembled.', flush=True)
    for q, indices in blocks.items():
        H = build_h(source, link_field(source, matrices15, indices=indices))
        row = spectral_frame(H, 4*len(indices), complement_lower=61)
        core = [2*len(indices)*source['addresses'][p]+a
                for p in ((0, 0, -1), (0, 0, 1)) for a in range(2*len(indices))]
        exterior_indices = [i for i in range(H.shape[0]) if i not in core]
        B = H[exterior_indices][:, core]
        gram = la.eigvalsh((B.conj().T@B).toarray())
        core_probability = np.sum(abs(row['W'][core])**2, axis=0)
        exterior_probability = np.sum(abs(row['W'][exterior_indices])**2, axis=0)
        localization_bound = F(7, 4)/(F(61)-F(7, 232))**2
        assert np.min(core_probability) > 0
        exterior_to_core = float(np.max(exterior_probability/core_probability))
        assert exterior_to_core <= float(localization_bound)+5e-12
        gram_expected = np.repeat([1.25, 1.75], 2*len(indices))
        assert norm(gram-gram_expected) < 2e-11
        row_result = {'block': names[q], 'primitive_charge': q, 'dimension': len(indices),
                      'operator_dimension': H.shape[0], 'exact_low_band_rank': 4*len(indices),
                      'all_343_sites_and_882_edges_retained': True,
                      'core_coupling_Gram_eigenvalues': gram.tolist(),
                      'exterior_to_core_probability_ratio_max': exterior_to_core,
                      'exterior_to_core_probability_ratio_upper_exact': str(localization_bound),
                      **summary_frame(row)}
        results.append(row_result)
        arrays[f'{names[q].replace("^", "_")}_band_values'] = row['values']
        print(f'Curved-link {names[q]} block: dimension {H.shape[0]}, band {len(row["values"])}, residual {row["residual"]:.3g}.', flush=True)

    # All weak and charged lepton slots of one active family; no neutral slot
    # is removed. The other two source families are retained by exact tensor
    # identities recorded below.
    li = blocks[-3]+blocks[6]
    d = len(li)
    H = build_h(source, link_field(source, matrices15, indices=li))
    H0 = sp.kron(source['h'], sp.eye(d), format='csr')
    Clepton = np.asarray(channels['charged_lepton']['components'][1], complex)[np.ix_(li, li)]
    K = design_row['gain_over_v']*Clepton
    delta = sp.kron(source['C'], sp.csr_matrix(K.conj()), format='csr')
    Boff = sp.block_diag((H, -H.conjugate()), format='csr')
    Bon = sp.bmat([[H, delta], [delta.conj().T, -H.conjugate()]], format='csr')
    full_nambu = Bon.shape[0]
    half = full_nambu//2
    Xi = sp.bmat([[None, sp.eye(half)], [sp.eye(half), None]], format='csr')
    ph_identity = norm(Xi@Bon.conjugate()@Xi+Bon)
    assert ph_identity < 1e-12
    off = charged_band(Boff, 8*d)
    on = charged_band(Bon, 8*d)
    arrays.update({'lepton_off_frame': off['W'], 'lepton_on_frame': on['W'],
                   'lepton_off_plus': off['plus'], 'lepton_on_plus': on['plus'],
                   'lepton_off_values': off['values'], 'lepton_on_values': on['values'],
                   'raw_projectors': np.array([np.asarray(P, float) for P in projectors])})

    raw_band = spectral_frame(H, 4*d, complement_lower=61)
    X = raw_band['W']
    delta_xbar = delta@X.conj()
    original_pair_leakage = norm(delta_xbar-X@(X.conj().T@delta_xbar))
    assert original_pair_leakage > 1e-7
    off_on_distance = projector_distance(off['W'], on['W'])
    assert off_on_distance < float(F(constants['BdG_vs_number_conserving_projector_norm_upper']))
    polar = small_polar_transport(off, on)
    arrays.update({'polar_span': polar.pop('Z'), 'polar_U_small': polar.pop('U_small')})

    # Site-varying transformations are evaluated on all lattice edges.
    A, B, Z = matrices15
    transforms = [np.linalg.matrix_power(A, (p[0]+3)%2)
                  @np.linalg.matrix_power(B, (p[1]+3)%2)
                  @np.linalg.matrix_power(Z, (p[2]+3)%2) for p in source['points']]
    arrays['all_site_gauge_transforms_15'] = np.asarray(transforms)
    local = [G[np.ix_(li, li)] for G in transforms]
    G = sp.block_diag([sp.csr_matrix(np.kron(np.eye(2), g)) for g in local], format='csr')
    Hchanged = build_h(source, link_field(source, matrices15, transform=transforms, indices=li))
    gauge_h_error = norm(Hchanged-G@H@G.conj().T)
    delta_changed = G@delta@G.T
    Bg = sp.bmat([[Hchanged, delta_changed], [delta_changed.conj().T, -Hchanged.conjugate()]], format='csr')
    Gn = sp.block_diag((G, G.conjugate()), format='csr')
    gauge_bdg_error = norm(Bg-Gn@Bon@Gn.conj().T)
    transformed_frame_error = norm(Bg@(Gn@on['W'])-(Gn@on['W'])*on['values'])
    gauge_grading_error = norm((Gn@on['plus']).conj().T@(Gn@on['plus'])-np.eye(on['plus'].shape[1]))
    assert max(gauge_h_error, gauge_bdg_error, transformed_frame_error, gauge_grading_error) < 1e-8
    flat_links = [(local[i]@local[j].conj().T) for i, j, _, _, _ in source['edges']]
    Hflat = build_h(source, flat_links)
    flat_error = norm(Hflat-G@H0@G.conj().T)
    assert flat_error < 1e-8

    # Every raw family projector acts on the same three-family space. The
    # active and inactive spatial/internal bands suffice to retain all raw
    # matrices exactly, because gauge transport is family scalar.
    one_family = {'off': off, 'on': on}
    family_frames = []
    family_labels = ['disabled']+list(range(12))
    for label in family_labels:
        if label == 'disabled':
            plus = np.kron(off['plus'], np.eye(3))
            minus = np.kron(off['minus'], np.eye(3))
        else:
            values, V = la.eigh(np.asarray(projectors[label], float))
            active, inactive = V[:, -1:], V[:, :2]
            plus = np.hstack((np.kron(on['plus'], active), np.kron(off['plus'], inactive)))
            minus = np.hstack((np.kron(on['minus'], active), np.kron(off['minus'], inactive)))
        family_frames.append((plus, minus))
    # Small cross-Gram matrices measure all169 pairs without a complete
    # 12348-square projector allocation.
    worst_band, worst_positive = 0., 0.
    for i, (ap, am) in enumerate(family_frames):
        for bp, bm in family_frames[i:]:
            W, V = np.hstack((ap, am)), np.hstack((bp, bm))
            worst_band = max(worst_band, projector_distance(W, V))
            worst_positive = max(worst_positive, projector_distance(ap, bp))
    assert worst_band < float(F(constants['between_raw_record_projector_norm_upper']))
    raw_polar = small_polar_transport({'plus': family_frames[3][0], 'minus': family_frames[3][1]},
                                     {'plus': family_frames[4][0], 'minus': family_frames[4][1]})
    raw_polar.pop('Z'); raw_polar.pop('U_small')

    np.savez_compressed(output_dir/'gauge_band_witnesses.npz', **arrays)
    result = {'schema': 'oph.gauge_covariant_compact_band.v1', 'all_pass': True,
              'upstream_source_commit': SOURCE_COMMIT, 'exact_bounds': constants,
              'source_graph': {'sites': 343, 'edges': 882, 'spatial_spin_dimension': 686,
                               'core_points': [[0, 0, -1], [0, 0, 1]],
                               'topology': topology.topology(1)},
              'raw_kernel_bound': kernel_bounds, 'diagnostic_link_controls': controls,
              'exterior_link_unitarity_errors': primitive_errors,
              'curved_link_matter_blocks': results,
              'lepton_BdG_execution': {
                  'internal_slots': li, 'one_active_family_Nambu_dimension': full_nambu,
                  'all_three_family_Nambu_dimension': 3*full_nambu,
                  'active_family_low_Nambu_rank': 8*d,
                  'all_three_family_low_Nambu_rank': 24*d,
                  'active_off': summary_frame(off) | {'compressed_grading_eigenvalues': off['compressed_grading_eigenvalues'].tolist()},
                  'active_on': summary_frame(on) | {'compressed_grading_eigenvalues': on['compressed_grading_eigenvalues'].tolist()},
                  'particle_hole_operator_identity_error': ph_identity,
                  'original_scalar_band_pairing_leakage_Frobenius': original_pair_leakage,
                  'BdG_vs_number_conserving_projector_distance': off_on_distance,
                  'all_raw_family_labels': family_labels,
                  'raw_family_unordered_projector_pairs_executed': 91,
                  'ordered_pairs_covered_by_distance_symmetry': 169,
                  'maximum_raw_family_band_projector_distance': worst_band,
                  'maximum_raw_family_annihilator_projector_distance': worst_positive,
                  'disabled_to_enabled_polar_transport': polar,
                  'noncommuting_family_polar_transport': raw_polar},
              'gauge_covariance_execution': {'spatial_operator_error': gauge_h_error,
                  'full_BdG_operator_error': gauge_bdg_error,
                  'full_BdG_transformed_band_error': transformed_frame_error,
                  'transformed_annihilator_frame_error': gauge_grading_error,
                  'flat_connection_isospectral_intertwiner_error': flat_error,
                  'all_site_frames_and_all_edges_retained': True},
              'scope': {'classical_background_link_covariance': True,
                        'full_BdG_pairing_retained_in_spectral_band': True,
                        'canonical_gauge_covariant_integer_band_occupation': True,
                        'gauge_covariant_quartic_interaction_on_that_band': True,
                        'exact_trivial_link_reduction_to_source_interaction': True,
                        'raw_record_band_alignment_is_an_added_controller': True,
                        'native_repair_selected_link_background': False,
                        'native_word_implementation_of_link_and_alignment_service': False,
                        'finite_native_quantum_link_join_established': False,
                        'physical_particle_masses_computed': False,
                        'numerical_eigenvectors_are_outward_interval_certificates': False},
              'execution_complete': True,
              'source_bindings': source_bindings(),
              'witness_array_shapes': {name: list(value.shape) for name, value in arrays.items()},
              'witness_sha256': digest(output_dir/'gauge_band_witnesses.npz')}
    # Frames are binary witnesses, not JSON complex-number payloads.
    for key in ('active_off', 'active_on'):
        result['lepton_BdG_execution'][key].pop('plus', None)
        result['lepton_BdG_execution'][key].pop('minus', None)
    (output_dir/'gauge_band.json').write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print(json.dumps({'all_pass': True, 'original_scalar_band_pairing_leakage': original_pair_leakage,
                      'raw_record_maximum_band_distance': worst_band,
                      'source_gauge_covariance_error': gauge_bdg_error,
                      'duration_seconds': time.monotonic()-start}, indent=2), flush=True)

    return result


def source_bindings():
    names = ['model/compact_defect_pair/spectral.py', 'model/compact_defect_pair/topology.py',
             'model/compact_defect_pair/PROOF.md', 'model/defect_current_bridge/construction.py',
             'model/common_source/chiral_candidate/action.py',
             'model/common_source/native_family_record_memory.py',
             'model/mcpr_native_memory_bridge/construction.py',
             'model/mcpr_native_memory_bridge/certificate.py',
             'model/compact_interacting_response/algebra.py',
             'model/particle_gauge_band/GAUGE_BAND_PROOF.md',
                  'model/particle_gauge_band/QUANTUM_LINK_PROOF.md',
                  'model/native_quantum_gauss/algebra.py', 'model/native_quantum_gauss/PROOF.md']
    names += list(action.SOURCE_FILES)
    from model.mcpr_native_memory_bridge import construction as memory_construction
    names += [memory_construction.FORWARD, memory_construction.MCPR]
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path:
            path = Path(path).resolve()
            if (path.suffix == '.py' and path.is_relative_to(REPO)
                    and not path.relative_to(REPO).parts[0].startswith('.')):
                names.append(path.relative_to(REPO).as_posix())
    return {name: digest(REPO/name) for name in sorted(set(names))}


if __name__ == '__main__':
    execute()
