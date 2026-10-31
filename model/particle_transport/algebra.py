"""Same-dimer charged states, graded transport and neutral pair response.

The action and candidates are exactly those of ``particle_response.dimer``:
c=1/8, U=1, K=H-N/2.  We use analytic columns, not an eigenvalue clustering
rule, to identify the two charged doublets.  ``setup`` returns the unchanged
predecessor for independent comparisons.  Floating matrix checks below are
diagnostics of the stated finite identities, not exact interval certificates.

For blocks with definite occupations n_i, the exterior lift of a block
permutation pi has sign (-1)**sum_{i<j,pi(i)>pi(j)} n_i*n_j.  This formula is
an exact integer algorithm on every occupation address.  Its composition
law is the exterior-power functor law, so arbitrary routes with the same
block permutation have the same endpoint operator.  Even ground blocks are
identity spectators; two odd candidate blocks exchange with minus the
ordinary swap of their internal doublets.  Tensoring these operator
identities by an identity proves their arbitrary-reference versions.

These are supplied permutations of decoupled identical dimers.  They do not
establish an adiabatic braid of the uniform interacting lattice, passive
protection by C/E, or electromagnetic/Lorentz particle identification.
"""
from __future__ import annotations

import ast
import math
from itertools import combinations
import numpy as np

from ..particle_response import dimer
from ..full_fock_field.field import FieldGeometry, FullFock
from ..interacting_field.lattice import SIGMA


SPEED = 1 / 8
INTERACTION = 1.
HOPPING = SPEED / 2
RADICAL = math.sqrt(INTERACTION**2 + 16*HOPPING**2)
GROUND_ENERGY = -(INTERACTION + RADICAL) / 2
LOWER_GAP = RADICAL/2 - HOPPING
LOWER_RESIDUE = (1 + 4*HOPPING/RADICAL) / 4
UPPER_RESIDUE = (1 - 4*HOPPING/RADICAL) / 4
KINDS = ("hole", "particle", "pair_hole", "pair_particle")


def _axis(axis):
    if type(axis) is not int or axis not in range(3):
        raise ValueError("integer axis in 0..2 required")
    return axis


def _kind(kind, band):
    if type(kind) is not str or kind not in KINDS:
        raise ValueError("hole, particle, pair_hole or pair_particle required")
    if type(band) is not str or band not in ("lower", "upper"):
        raise ValueError("lower or upper charged doublet required")
    if kind.startswith("pair_") and band != "lower":
        raise ValueError("the charge-two pair has only one band")


def setup(axis=0):
    axis = _axis(axis)
    sides = tuple(2 if j == axis else 1 for j in range(3))
    return dimer.setup(FieldGeometry(sides, speed=SPEED, interaction=INTERACTION))


def ground_column(axis=0):
    """Analytic ground with positive separated-singlet coefficient.

    In d_0=c_0,d_1=i sigma_axis c_1, use
    S=(d0u* d1d* - d0d* d1u*)|0>/sqrt(2),
    D=(d0u* d0d* + d1u* d1d*)|0>/sqrt(2).
    The ground is a*S+b*D, a>0,b<0, a*b=-2t/RADICAL.
    """
    axis = _axis(axis)
    space = FullFock(4)
    c = [space.annihilator(j).toarray() for j in range(4)]
    turn = 1j*SIGMA[axis]
    d = c[:2] + [sum(turn[s, j]*c[2+j] for j in range(2)) for s in range(2)]
    create = [x.conj().T for x in d]
    blank = np.eye(16, dtype=complex)[:, 0]
    separated = (create[0]@create[3] - create[1]@create[2])@blank/math.sqrt(2)
    doublon = (create[0]@create[1] + create[2]@create[3])@blank/math.sqrt(2)
    a = math.sqrt((1 + INTERACTION/RADICAL)/2)
    b = -2*HOPPING/(RADICAL*a)
    return a*separated + b*doublon


def candidate_columns(kind="hole", axis=0, band="lower"):
    """Fixed analytic charged columns; particle columns are C times holes.

    The two hole columns use the left site's ordered internal spin basis.
    No eigenvector phase or spectral visibility threshold chooses this frame.
    """
    axis = _axis(axis)
    _kind(kind, band)
    if kind.startswith("pair_"):
        out = np.zeros((16, 1), complex)
        out[0 if kind == "pair_hole" else 15, 0] = 1
        return out
    out = np.zeros((16, 2), complex)
    out[[1, 2], :] = np.eye(2)/math.sqrt(2)
    out[[4, 8], :] = (1 if band == "lower" else -1)*1j*SIGMA[axis]/math.sqrt(2)
    if kind == "particle":
        out = dimer.charge_conjugation(setup(axis))@out
    return out


def candidate_projector(kind="hole", axis=0, band="lower"):
    columns = candidate_columns(kind, axis, band)
    return columns@columns.conj().T


def candidate_gap(kind="hole", band="lower"):
    _kind(kind, band)
    if kind.startswith("pair_"):
        return (INTERACTION + RADICAL)/2
    return RADICAL/2 + (-1 if band == "lower" else 1)*HOPPING


def _permutation(permutation, block_modes):
    if (type(permutation) is not tuple or not 1 <= len(permutation) <= 12
            or any(type(j) is not int for j in permutation)
            or sorted(permutation) != list(range(len(permutation)))):
        raise ValueError("a permutation tuple of 1..12 block addresses required")
    if type(block_modes) is not int or not 1 <= block_modes <= 4:
        raise ValueError("integer block width in 1..4 required")


def block_permutation_address(address, permutation, block_modes=4):
    """Exact full-Fock address map; pi maps each old block to its new place."""
    _permutation(permutation, block_modes)
    modes = block_modes*len(permutation)
    if type(address) is not int or not 0 <= address < 1 << modes:
        raise ValueError("integer address on the entire declared Fock carrier required")
    mask = (1 << block_modes)-1
    blocks = tuple((address >> (block_modes*j)) & mask for j in range(len(permutation)))
    parities = tuple(x.bit_count() % 2 for x in blocks)
    crossings = sum(parities[i]*parities[j]
                    for i in range(len(blocks)) for j in range(i+1, len(blocks))
                    if permutation[i] > permutation[j])
    target = sum(blocks[i] << (block_modes*permutation[i]) for i in range(len(blocks)))
    return target, (-1 if crossings % 2 else 1)


def permutation_action(permutation, block_modes=4):
    """Bounded dense operator, retaining every Fock sector and scalar sign."""
    _permutation(permutation, block_modes)
    modes = block_modes*len(permutation)
    if modes > 8:
        raise ValueError("dense permutation evaluation supports at most eight modes")
    dimension = 1 << modes
    out = np.zeros((dimension, dimension), dtype=int)
    for address in range(dimension):
        target, sign = block_permutation_address(address, permutation, block_modes)
        out[target, address] = sign
    return out


def block_swap(block_modes=4):
    return permutation_action((1, 0), block_modes)


def adjacent_swap_route(permutation, block_modes=4):
    """Deterministic neighboring-mode FSWAP route for the requested blocks.

    Each returned integer is the left mode of one neighboring FSWAP.
    This is a route specification; source execution/resources are separate.
    """
    _permutation(permutation, block_modes)
    current = list(range(block_modes*len(permutation)))
    target = [None]*len(current)
    for old, new in enumerate(permutation):
        for mode in range(block_modes):
            target[new*block_modes+mode] = old*block_modes+mode
    route = []
    for place, mode in enumerate(target):
        located = current.index(mode)
        while located > place:
            left = located-1
            current[left], current[located] = current[located], current[left]
            route.append(left)
            located -= 1
    return tuple(route)


def pair_generator(axis=0):
    """Two independent identical dimers, not a uniform four-site field."""
    generator = setup(axis).generator
    return np.kron(np.eye(16), generator)+np.kron(generator, np.eye(16))


def localized_columns(kind="hole", location=0, axis=0, band="lower"):
    """One candidate and one unchanged even ground block.

    Block0 occupies the four low bits.  These isometries preserve an unknown
    internal doublet and arbitrary reference, rather than preparing a spin.
    """
    if type(location) is not int or location not in (0, 1):
        raise ValueError("integer localized block 0 or 1 required")
    columns, ground = candidate_columns(kind, axis, band), ground_column(axis)
    return np.column_stack([np.kron(ground, column) if location == 0 else np.kron(column, ground)
                            for column in columns.T])


def transfer_operator(axis=0, mode=0, quadrature="current"):
    """Neutral completed operator T=c_A,mode* c_B,mode and its quadratures.

    ``real`` is T+T*, ``current`` is -i(T-T*); ``directed`` returns T.
    No independently controlled odd field operator is needed for these even
    operators.  Axis labels the same-dimer frame, although T itself is local
    to the original mode chart and independent of the hopping axis.
    """
    _axis(axis)
    if type(mode) is not int or mode not in range(4):
        raise ValueError("integer local mode in 0..3 required")
    if type(quadrature) is not str or quadrature not in ("current", "real", "directed"):
        raise ValueError("current, real or directed transfer required")
    space = FullFock(8)
    transfer = (space.annihilator(mode).conj().T@space.annihilator(4+mode)).toarray()
    if quadrature == "directed":
        return transfer
    return transfer+transfer.conj().T if quadrature == "real" else -1j*(transfer-transfer.conj().T)


def pair_band_projector(band="lower", axis=0):
    """Both particle/hole placements; no state projection is executed."""
    _axis(axis)
    if type(band) is not str or band not in ("lower", "mixed", "upper"):
        raise ValueError("lower, mixed or upper particle-hole band required")
    choices = {"lower": (("lower", "lower"),),
               "mixed": (("lower", "upper"), ("upper", "lower")),
               "upper": (("upper", "upper"),)}[band]
    out = np.zeros((256, 256), complex)
    for a, b in choices:
        hole, particle = candidate_projector("hole", axis, a), candidate_projector("particle", axis, b)
        out += np.kron(hole, particle)+np.kron(particle, hole)
    return out


def pair_line_data(quadrature="current"):
    if type(quadrature) is not str or quadrature not in ("current", "real", "directed"):
        raise ValueError("current, real or directed transfer required")
    factor = 1 if quadrature == "directed" else 2
    return tuple((2*LOWER_GAP + j*2*HOPPING, factor*weight)
                 for j, weight in enumerate((LOWER_RESIDUE**2,
                     2*LOWER_RESIDUE*UPPER_RESIDUE, UPPER_RESIDUE**2)))


def exchange_pair_column(axis=0, mode=0, band="lower", quadrature="current"):
    """Normalized mathematical pair target, not a postselected preparation.

    The current creates an exchange-odd pair; real transfer creates an
    exchange-even pair.  Both have identical complete energy weights.
    A source preparation must implement this supplied target separately.
    """
    if quadrature not in ("current", "real") or type(quadrature) is not str:
        raise ValueError("a Hermitian current or real quadrature required")
    projector = pair_band_projector(band, axis)
    ground = ground_column(axis)
    out = projector@transfer_operator(axis, mode, quadrature)@np.kron(ground, ground)
    index = ("lower", "mixed", "upper").index(band)
    return out/math.sqrt(pair_line_data(quadrature)[index][1])


def neutral_pair_measure(axis=0, mode=0, quadrature="current"):
    """All 256 eigenstate weights, including dark and wrong-charge states.

    Tensor products of the predecessor's complete sixteen-state eigenbasis
    are used only for this independent numerical spectral check.  The three
    closed-form lines above are not inferred from frequency clustering.
    """
    operator = transfer_operator(axis, mode, quadrature)
    spectrum = setup(axis)
    ground = ground_column(axis)
    basis = np.kron(spectrum.vectors, spectrum.vectors)
    frequencies = (spectrum.energies[:, None]+spectrum.energies[None, :]).ravel()-2*GROUND_ENERGY
    weights = abs(basis.conj().T@operator@np.kron(ground, ground))**2
    numbers = (spectrum.sectors[:, None]+spectrum.sectors[None, :]).ravel()
    return {"frequencies": frequencies, "weights": weights, "numbers": numbers}


def _residual(value):
    """Frobenius residual also bounds the matrix operator norm."""
    return float(np.linalg.norm(value))


def _native_dictionary_chart():
    """Independently bind the consumed port chart and all proper rotations.

    Read the literal faces from the pinned source bytes, not the mutable loaded
    carrier. A directed edge determines the image of its oriented triangle;
    propagating across adjacent triangles determines the full port map. The
    sixty directed-edge images reconstruct the proper group without the
    carrier's graph-automorphism enumerator or cached rotation catalogue.
    """
    from ..common_source import slow_geometry
    from ..exchange_sector import statistics
    from ..observer import native

    slow_geometry.source.cache_clear()
    carrier = slow_geometry.source()
    path = slow_geometry.ARCHIVE/'oph_fpe/dynamics/self_readback_repair_closure.py'
    tree = ast.parse(path.read_text(encoding='utf-8'))
    declarations = [node.value for node in tree.body if isinstance(node, ast.Assign)
                    and any(isinstance(target, ast.Name) and target.id == 'ORIENTED_BASE_FACES'
                            for target in node.targets)]
    if len(declarations) != 1:
        raise ArithmeticError("unique pinned native face declaration required")
    faces = ast.literal_eval(declarations[0])

    def integer_rows(value, width, length):
        return (type(value) is tuple and len(value) == length
                and all(type(row) is tuple and len(row) == width
                        and all(type(item) is int for item in row) for row in value))

    if (not integer_rows(faces, 3, 20)
            or any(len(set(face)) != 3 or not set(face) <= set(range(12)) for face in faces)
            or len({frozenset(face) for face in faces}) != 20):
        raise ArithmeticError("complete pinned native face chart required")
    third = {(face[j], face[(j+1) % 3]): face[(j+2) % 3]
             for face in faces for j in range(3)}
    seams = tuple(sorted({tuple(sorted(edge)) for edge in third}))
    if (len(third) != 60 or len(seams) != 30
            or any((b, a) not in third for a, b in third)
            or any(sum(p in edge for edge in seams) != 5 for p in range(12))):
        raise ArithmeticError("pinned native oriented seam incidence changed")
    live_faces, live_seams = carrier.oriented_faces(), carrier.seams()
    if (not integer_rows(live_faces, 3, 20) or live_faces != faces
            or not integer_rows(live_seams, 2, 30) or live_seams != seams
            or not integer_rows(statistics.SEAMS, 2, 30) or statistics.SEAMS != seams):
        raise ArithmeticError("runtime native face/seam chart differs from pinned source")

    reconstructed = set()
    first, second = faces[0][:2]
    for target_edge in third:
        image = {first: target_edge[0], second: target_edge[1]}
        while True:
            before = len(image)
            for (a, b), c in third.items():
                if a not in image or b not in image:
                    continue
                target = third.get((image[a], image[b]))
                if target is None or (c in image and image[c] != target):
                    raise ArithmeticError("oriented face propagation is inconsistent")
                if c not in image:
                    if target in image.values():
                        raise ArithmeticError("oriented face propagation is not injective")
                    image[c] = target
            if len(image) == before:
                break
        if set(image) != set(range(12)) or set(image.values()) != set(range(12)):
            raise ArithmeticError("oriented faces did not determine every port image")
        reconstructed.add(tuple(image[p] for p in range(12)))
    if len(reconstructed) != 60:
        raise ArithmeticError("oriented faces did not reconstruct sixty proper rotations")
    rotations = carrier.rotations()
    if (not integer_rows(rotations, 12, 60) or len(set(rotations)) != 60
            or set(rotations) != reconstructed):
        raise ArithmeticError("complete distinct proper native rotation catalogue required")

    states = tuple(sorted(tuple(int(j in pair) for j in range(12))
                          for pair in combinations(range(12), 2)))
    pairs = tuple(tuple(j for j, occupied in enumerate(state) if occupied) for state in states)
    live_states = native.states(binary=True)
    index = {pair: j for j, pair in enumerate(pairs)}
    if (not integer_rows(live_states, 12, 66) or live_states != states
            or not integer_rows(statistics.STATES, 12, 66) or statistics.STATES != states
            or not integer_rows(statistics.PAIRS, 2, 66) or statistics.PAIRS != pairs
            or type(statistics.INDEX) is not dict or statistics.INDEX != index
            or any(type(value) is not int for value in statistics.INDEX.values())
            or not integer_rows(tuple(statistics.INDEX), 2, 66)):
        raise ArithmeticError("native binary state/pair/index chart changed")
    return rotations


def native_intertwiner_certificate():
    """Exact finite dictionary to the existing negative H66 instruments.

    At each of twelve blocks choose the same fixed normalized odd candidate
    column v, and the same even g.  W|{a,b}> has v at a,b and g elsewhere.
    W is an isometry into 48 field modes, with relative charge twice v's
    charge.  Block-number parity determines the entire permutation sign,
    so the following calculation needs only twelve parity bits.  It does
    not allocate a 48-mode state or execute a 48-mode repair instrument.

    On this embedded carrier let P_r be the unchanged native input-support
    projector.  L_r=Gamma(swap_ab)P_r/sqrt(2) for a nonzero current and
    L_r=P_r/sqrt(2) for a hold.  Then L_r W = W K_r^- for every complete raw
    flag.  Completeness is on im(W); no full-Fock extension is asserted.
    """
    from ..exchange_sector import statistics
    from ..observer.source_discrimination import circuit

    rotations = _native_dictionary_chart()
    pair_index = {pair: j for j, pair in enumerate(statistics.PAIRS)}
    if len(pair_index) != 66:
        raise ArithmeticError("native embedded pair carrier changed")
    movements, negative_movements = set(), set()
    raw_slots, support_columns = 0, 0
    for law in ("C", "E"):
        for seam in statistics.SEAMS:
            a, b = seam
            transposition = list(range(12))
            transposition[a], transposition[b] = b, a
            transposition = tuple(transposition)
            for coin in (0, 1):
                # Reconstruct all admitted records on the full native local
                # carrier, including collision inputs outside the binary image.
                # The branch builder uses circuit.record_values too; comparing
                # it only with that shared decoder would not bind abort slots.
                decoder = {}
                for left, right in ((0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0)):
                    difference, parity = left-right, (left+right) % 2
                    current = (difference-(2*coin-1)*parity)//2
                    raw = {0: 0, 1: 1, -1: 2}[current]
                    if law == "E":
                        raw += 4*parity
                    decoder[raw] = current if law == "C" else difference
                live_decoder = circuit.record_values(law, coin)
                if (type(live_decoder) is not dict or live_decoder != decoder
                        or any(type(key) is not int or type(value) is not int
                               for key, value in live_decoder.items())):
                    raise ArithmeticError("complete integer native record decoder changed")
                expected = {raw: [] for raw in range(4 if law == "C" else 8)}
                for i, pair in enumerate(statistics.PAIRS):
                    at_a, at_b = a in pair, b in pair
                    current = 0
                    target_pair, phase = pair, 1
                    if at_a != at_b:
                        destination = a if coin else b
                        origin = a if at_a else b
                        if origin != destination:
                            current = 1 if origin == a else -1
                            parity_address = sum(1 << j for j in pair)
                            target, phase = block_permutation_address(parity_address, transposition, 1)
                            target_pair = tuple(j for j in range(12) if target >> j & 1)
                            movement = (pair, target_pair)
                            movements.add(movement)
                            if phase == -1:
                                negative_movements.add(movement)
                    raw = {0: 0, 1: 1, -1: 2}[current]
                    if law == "E":
                        raw += 4*int(at_a != at_b)
                    expected[raw].append((i, pair_index[target_pair], phase))
                for raw, rows in expected.items():
                    branch = statistics.branch(law, seam, coin, raw, -1)
                    if (type(branch) is not statistics.SignedBranch
                            or type(branch.law) is not str
                            or type(branch.seam) is not tuple or len(branch.seam) != 2
                            or any(type(port) is not int for port in branch.seam)
                            or type(branch.coin) is not int or type(branch.raw) is not int
                            or type(branch.abort) is not bool
                            or any(type(column) is not tuple for column in
                                   (branch.inputs, branch.outputs, branch.signs))
                            or any(type(value) is not int or not 0 <= value < 66
                                   for column in (branch.inputs, branch.outputs) for value in column)
                            or any(type(sign) is not int or sign not in (-1, 1) for sign in branch.signs)):
                        raise ArithmeticError("strict integer native signed branch fields required")
                    actual = list(zip(branch.inputs, branch.outputs, branch.signs, strict=True))
                    if (actual != rows or branch.abort != (raw not in decoder)
                            or (branch.law, branch.seam, branch.coin, branch.raw) != (law, seam, coin, raw)):
                        raise ArithmeticError("same-dimer embedding fails a complete native signed branch")
                    raw_slots += 1
                    support_columns += len(rows)
    rotation_columns = 0
    for permutation in rotations:
        images, phases = statistics.rotation_action(permutation, -1)
        for i, pair in enumerate(statistics.PAIRS):
            target, phase = block_permutation_address(sum(1 << j for j in pair), permutation, 1)
            target_pair = tuple(j for j in range(12) if target >> j & 1)
            if (pair_index[target_pair], phase) != (images[i], phases[i]):
                raise ArithmeticError("same-dimer block permutation fails native rotation covariance")
            rotation_columns += 1
    if (len(movements), raw_slots, support_columns, len(rotations), rotation_columns) != (600, 720, 7920, 60, 3960):
        raise ArithmeticError("complete native dictionary census changed")
    return {"embedded_pair_dimension": 66, "formal_field_modes": 48,
        "candidate_fiber": "same fixed odd column in every identical dimer",
        "reference_total_CAR_number": 24,
        "relative_CAR_number_per_occupied_block": {"hole": -1, "particle": 1},
        "current_conversion": "native defect current times the candidate excess CAR number",
        "exact_directed_movement_columns": len(movements),
        "negative_directed_movement_columns": len(negative_movements),
        "complete_C_and_E_raw_slots": raw_slots,
        "complete_C_and_E_supported_columns": support_columns,
        "proper_rotations": len(rotations), "exact_rotation_columns": rotation_columns,
        "branch_scalar_amplitude": "1/sqrt(2)",
        "complete_instrument_domain": "66-dimensional image of the specified isometry",
        "every_finite_record_history_intertwines": True,
        "48_mode_source_execution_claimed": False,
        "native_memory_protection_inherited": False}


def certificate():
    swap = block_swap()
    if not np.array_equal(swap@swap, np.eye(256, dtype=int)):
        raise ArithmeticError("complete block permutation failed its exact involution")
    # Independent neighboring FSWAP action on all addresses, including |11>.
    route = adjacent_swap_route((1, 0))
    for source in range(256):
        target, sign = source, 1
        for left in route:
            low, high = (target >> left) & 1, (target >> (left+1)) & 1
            if low and high:
                sign = -sign
            if low != high:
                target ^= (1 << left) | (1 << (left+1))
        if ((target, sign) != block_permutation_address(source, (1, 0))
                or swap[target, source] != sign or np.count_nonzero(swap[:, source]) != 1):
            raise ArithmeticError("mode routing and complete graded block swap disagree")
    residuals, survival = [], []
    for axis in range(3):
        spectrum, ground = setup(axis), ground_column(axis)
        residuals += [_residual(spectrum.generator@ground-GROUND_ENERGY*ground),
                      abs(float(np.vdot(ground, ground).real)-1),
                      _residual(np.outer(ground, ground.conj())-spectrum.reference)]
        prior = dimer.candidate_projectors(spectrum)
        for kind, old in zip(KINDS, ("charge_minus_one_hole", "charge_plus_one_particle",
                                    "charge_minus_two_pair_hole", "charge_plus_two_pair_particle"), strict=True):
            columns = candidate_columns(kind, axis)
            residuals += [_residual(columns.conj().T@columns-np.eye(columns.shape[1])),
                _residual((spectrum.generator-GROUND_ENERGY*np.eye(16))@columns-candidate_gap(kind)*columns),
                _residual(columns@columns.conj().T-prior[old]),
                _residual(swap@localized_columns(kind, 0, axis)-localized_columns(kind, 1, axis))]
        generator, pair_ground = pair_generator(axis), np.kron(ground, ground)
        residuals += [_residual(swap@generator-generator@swap), _residual(swap@pair_ground-pair_ground)]
        for mode in range(4):
            for quadrature, parity in (("current", -1), ("real", 1)):
                operator = transfer_operator(axis, mode, quadrature)
                residuals.append(_residual(swap@operator@swap-parity*operator))
                measure = neutral_pair_measure(axis, mode, quadrature)
                residuals.append(float(measure["weights"][measure["numbers"] != 4].sum()))
                for band, (frequency, weight) in zip(("lower", "mixed", "upper"), pair_line_data(quadrature), strict=True):
                    projected = pair_band_projector(band, axis)@operator@pair_ground
                    residuals += [abs(float(np.vdot(projected, projected).real)-weight),
                        _residual((generator-2*GROUND_ENERGY*np.eye(256))@projected-frequency*projected)]
                for order in (0, 1, 2):
                    actual = float(np.sum(measure["weights"]*measure["frequencies"]**order))
                    residuals.append(abs(actual-sum(w*f**order for f, w in pair_line_data(quadrature))))
        # A complete allowed number read preserves charge but not this band.
        number_read = np.diag([int(address & 1 != 0) for address in range(16)])
        for kind in ("hole", "particle"):
            projector = candidate_projector(kind, axis)
            rho = projector/2
            changed = number_read@rho@number_read + (np.eye(16)-number_read)@rho@(np.eye(16)-number_read)
            kept = float(np.trace(projector@changed).real)
            survival.append(kept)
            residuals.append(abs(kept-3/4))
    if not residuals or any(not math.isfinite(value) or value > 2e-11 for value in residuals):
        raise ArithmeticError("same-candidate transport algebra failed its numerical identity checks")
    return {"scope": "SUPPLIED_DECOUPLED_IDENTICAL_DIMER_TRANSPORT",
        "full_Fock_dimension": 256, "axes_checked": 3,
        "same_lower_charged_gap": LOWER_GAP,
        "same_lower_local_CAR_residue": LOWER_RESIDUE,
        "complete_exact_permutation_columns": 256,
        "adjacent_mode_FSWAPs_per_block_exchange": len(route),
        "unknown_spin_and_arbitrary_reference_operator_identity": True,
        "even_ground_block_exchange_sign": 1,
        "two_odd_candidate_block_sign_relative_to_ordinary_swap": -1,
        "neutral_current_pair_exchange_parity": -1,
        "neutral_real_transfer_pair_exchange_parity": 1,
        "neutral_current_complete_lines": [{"frequency": f, "weight": w} for f, w in pair_line_data()],
        "spectral_states_per_setting": 256,
        "spectral_settings_checked": 24,
        "number_read_lower_doublet_survival": 3/4,
        "native_negative_lift_intertwiner": native_intertwiner_certificate(),
        "numerical_identity_admission": "Frobenius residual <= 2e-11",
        "fresh_numerical_identity_checks_passed": True,
        "passive_topological_protection": False,
        "uniform_four_site_field_transport_claimed": False,
        "physical_particle_identification": False}
