"""The retained elementary carrier, joined through its declared seam routing.

Graph geometry and gluing are supplied by the existing public source. Native
binary repair, coherent outcome sectors, conservation and continuing records
are preserved. No local fixed-particle sector is assumed invariant across a
federation seam.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib
from itertools import combinations, product
import math
from pathlib import Path
import sys

import numpy as np

from ..observer import native
from ..observer.scalar_response.experiment import canonical
from ..primitive_observer import runtime, transfer

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE_NAME = 'evidence/observer_dynamics_20260925/oph-physics-sim'
ARCHIVE = ROOT/ARCHIVE_NAME
# Literal source hashes are checked before importing the retained package.
# Package initializers load more modules than the active geometry path uses;
# those modules are also bound, without adopting their scientific claims.
SOURCE_PINS = {
    'oph_exact/__init__.py': '3ceeb20c7231ce9f8d3002529c5416523986dd9c79ea963b2e314bd13364726e',
    'oph_exact/carrier.py': '83354f120a0d40d9fd3c04721bbbb3c6140cc3c5bf10655d899e9aa4325804bb',
    'oph_exact/federation.py': '29b0bb7a1e4d8a2077f16dbb508b7e08e794a7191a4a9f51e6c1b471e726ca1d',
    'oph_fpe/__init__.py': '7f72d89ca6a367344b39089fbf24163e3eb6c742476ee5db3c813b12ff73617e',
    'oph_fpe/claims.py': '962a829a454f8909773d818b96b656fc971aa9d1579cc4c53bbe6cde3cabfe14',
    'oph_fpe/constants/__init__.py': 'bc732866538c149253d0e72723ef8891e5ed68fd676b46dfaec237875540f6aa',
    'oph_fpe/constants/oph_pixel.py': '9ab43098d6648010b6bd17cbfa7382d4355619a92fe2d32d3339d79af615ec9f',
    'oph_fpe/core/__init__.py': '56d3d5e684d53b61f38b24f32111728866fc5dc0151495dece8702d3c591d161',
    'oph_fpe/core/array_geometry.py': 'e9b0d0c20c8598cb6e49a4965d892c8323b513b36a47dbb53b0cd49a6f873cb8',
    'oph_fpe/core/echosahedral_dynamics.py': 'a20b928552bb51257bcc7b976f513fe227f5f36158b6e3b4786c315477c862da',
    'oph_fpe/core/echosahedral_federation.py': 'ead1a96812ac1b243bbe2c1d175a368cd23d5c7b29a179e8718a82baf2edeec1',
    'oph_fpe/core/graph.py': '533bb11fc2f53600d4f3bfec25c5781ddd2d1b31f41e3fbadbb73a6a219f8667',
    'oph_fpe/core/icosahedral.py': '89fac4f2d24bb4d85c4a84624b3c826e99cb717df68060dd23a4b26c6cfe9adc',
    'oph_fpe/core/patch_state.py': '85b7c4b9fe8b983571ca1944b5871c959e198ac4970c2412e85b195d4f34221f',
    'oph_fpe/core/patchnet.py': '4664fe5717aac6e1433ab8ec618786cf7fbbd381b421632153ebf9ae78918dc7',
    'oph_fpe/core/pixel_scale.py': '0c6ce1b055df4334abba2a8b150905b4bb551e8b7f18fb2674cb6bcbb06963e5',
    'oph_fpe/core/screen_microphysics.py': '8fe76b1b9bde24060afa6adcbb51e925865799866d1884490c26e4b172d59cf9',
    'oph_fpe/core/screen_ports.py': '37b79a8a949fad2d18f17e850533c6ef33581f64c8f1a72dd94e45da93cef68e',
    'oph_fpe/dynamics/__init__.py': '0f8ea883cedb9b67cb2c6af038ebdcb9fc332b0a0d62a5ea1c5cc87d0b22794c',
    'oph_fpe/dynamics/anneal.py': 'ff008ea37005cf667825a8cd8c7b2c5c6a2110523722da5ab089f8422dd366eb',
    'oph_fpe/dynamics/consensus_certificate.py': 'e6d282be435e41269ac8d581997508e5559d558ab8ffdc9e5756e3030cd7fc0f',
    'oph_fpe/dynamics/kernel_dispatcher.py': '5ecb419adbdd7a341359a8aab13bef3c2feb4490130c20defadc7a196f01e500',
    'oph_fpe/dynamics/modular.py': '8978134a5846ad215c7a69205401c479a02019b7902e31bef283f18a6290ef9e',
    'oph_fpe/dynamics/positive_geometry.py': 'b87cfb52bdc1e1d5762fa65259a026b50072d54d4e1c68bdbcac0159b05bf435',
    'oph_fpe/dynamics/repair.py': 'fdc6326cdbac39c86e00e9ad314a2d9fd6cf06714a2e06f0470b298177d80949',
    'oph_fpe/dynamics/self_readback_repair_closure.py': '7e7cda89350be6330c7425df69bf9647549e6e1890a64eaaf28289be6c640b2f',
    'oph_fpe/evidence/__init__.py': '1875a81eb0bbac4cedfd2c3df77cabe06270b2ce8313e781b1af43174cdc04f3',
    'oph_fpe/evidence/bundle.py': 'b84548c3f62187ecc19b64fa11a5dc24915bc1b9a05066b3a6496524c5f61a71',
    'oph_fpe/evidence/cross_repo_artifacts.py': 'f026e8653752e5bcb1c703ed253861d3444d7beff2fe572863b09508d7efbff6',
    'oph_fpe/evidence/hashes.py': 'e4a7ada049fa1137e2d53c8eadd36a574f03fca25793322906e854d8ef04773f',
    'oph_fpe/evidence/particle_input_policy.py': '8d3683c8a3b6d249c484d253e08fd8425899c71a8d69ab5c78ee34780a7e06d3',
    'oph_fpe/evidence/production_envelope.py': 'f34f03314efae563a3b42cae6646164354c74e7bf062f5091790a584d08c6d1e',
    'oph_fpe/evidence/verifier.py': '4997bb27fe34846b9c2a4fff14f8a8b9a6b238f4029e2ec0273f8e1f7de60f45',
    'oph_fpe/finite_groups.py': '1329f07ef90beb9d4be1efa5cd0cb76e323014670d8e90c8f7945eba54e401f6',
    'oph_fpe/gauge/__init__.py': 'ede5295ec9dac5c42571287b5971c61a24647c3c9dee0b8b7187b7f195c9bb54',
    'oph_fpe/gauge/a5_sm_certificate.py': '8d7a99cc1ae093bc570bbb3b1677a028d23b82be23e90de4635a76380bde92f1',
    'oph_fpe/gauge/covariant_overlap.py': 'f03127e7314b0be60ba18b0a0b4cc934aa8d314785636258796c68baae374aa9',
    'oph_fpe/gauge/higgs_carrier.py': 'c150b9850364a4948995eca273e3571a9808960e840ad9fa361801c2beec5958',
    'oph_fpe/gauge/mar_sieve.py': '0043cc7d784f1be263dca5e3882b3d1698b8b91d75a8f2ba28beee516086a366',
    'oph_fpe/gauge/physical_a5_sm_requirements.py': '6e4b35bdf56d025793973d6dcc15cb34d0e60b5e7e5796926fdf96fa2cd03015',
    'oph_fpe/gauge/repair_projection.py': 'e2fdd6165ae39006e3b83c0376137162d423adeabec293ab5c06c6fe0b2b0cec',
    'oph_fpe/gauge/source_to_eft.py': '6414f00bda171fad6e8e6c77295cab106fb163eae9f00a6656e0e89db4f87562',
    'oph_fpe/gauge/yang_mills_gap.py': '29d4bb3ba2ee6e25365196c48c6e6f96a61aee00325de0acde470a577098bc4c',
    'oph_fpe/groups/__init__.py': '8e99c23175f57934eecaf9103dabd7b3d20d56ed55fe38174571290b3ef2b815',
    'oph_fpe/groups/base.py': '8c3e34417a51c60fc2b3a6c7b9611a9ea476d90fd1e7c5e7816604855a90d2b2',
    'oph_fpe/groups/clock.py': '15434d1fc409eb29d585bc4d5feefd8819279276e082d9ed711de3b598de19bc',
    'oph_fpe/groups/registry.py': '7876a6087be4ca9ef4314625799a165e396f269fffbf51b7f6984805057f434f',
    'oph_fpe/groups/s3.py': '65d34b9ec87c46a1befeb4a1945e1827fb1987a61493e890bb4b3c9ba384a675',
    'oph_fpe/groups/z2.py': 'a32eb994598592ba9f9b763eee455a1e13720b0a8311a4d46c299794466a3624',
}
SOURCE_PATHS = tuple(sorted(set((*transfer.SOURCE_PATHS,
    'model/federated_observer/source.py',
    *(f'{ARCHIVE_NAME}/{path}' for path in SOURCE_PINS)))))
LETTERS = tuple(product((0, 1), (-1, 0, 1)))
CURRENT_VALUES = (-1, 0, 1)


def verify_source_pins():
    for name, expected in SOURCE_PINS.items():
        if hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'the pinned federation source changed: {name}')


@lru_cache(maxsize=1)
def _source():
    verify_source_pins()
    location = str(ARCHIVE)
    sys.path.insert(0, location)
    try:
        result = importlib.import_module('oph_exact.federation')
    finally:
        sys.path.remove(location)
    for name in SOURCE_PINS:
        module_name = name.removesuffix('.py').replace('/', '.')
        if module_name.endswith('.__init__'):
            module_name = module_name.removesuffix('.__init__')
        module = sys.modules.get(module_name)
        if module is None or Path(module.__file__).resolve() != (ARCHIVE/name).resolve():
            raise ValueError(f'foreign or missing retained federation module: {module_name}')
    return result


@dataclass(frozen=True)
class Federation:
    level: int
    cells: int
    ports: int
    seams: tuple[tuple[int, int], ...]
    neighbours: tuple[tuple[int, ...], ...]
    inter_seams: tuple[tuple[int, int], ...]

    @property
    def max_degree(self):
        return max(map(len, self.neighbours))


@lru_cache(maxsize=3, typed=True)
def federation(level):
    if type(level) is not int or level not in (0, 1, 2):
        raise ValueError('a declared federation level zero, one or two is required')
    raw = _source().build_federation(level, 'port_pair')
    seams = tuple(sorted(tuple(sorted((int(a), int(b))))
                         for a, b in zip(raw.seam_a, raw.seam_b, strict=True)))
    inter = tuple(sorted(tuple(sorted((12*int(a)+int(p), 12*int(b)+int(q))))
                         for a, p, b, q in raw.inter_pairs))
    cells = 20*4**level
    neighbours = tuple(tuple(sorted(map(int, row))) for row in raw.graph.neighbour_lists)
    if (raw.carriers != cells or raw.ports != 12*cells or len(seams) != 63*cells//2
            or len(inter) != 3*cells//2 or len(set(seams)) != len(seams)
            or any(not 0 <= a < b < raw.ports for a, b in seams)
            or sorted(map(len, neighbours)) != [5]*(9*cells)+[6]*(3*cells)):
        raise ArithmeticError('the declared elementary-carrier federation census changed')
    local = set(transfer.native_counts().seams)
    for cell in range(cells):
        found = {(a % 12, b % 12) for a, b in seams if a//12 == b//12 == cell}
        if found != local:
            raise ArithmeticError('an elementary carrier lost its original thirty seams')
    edges_from_neighbours = {(min(a, b), max(a, b))
                            for a, row in enumerate(neighbours) for b in row}
    if edges_from_neighbours != set(seams):
        raise ArithmeticError('federation adjacency disagrees with its seam list')
    seen, pending = {0}, [0]
    while pending:
        for vertex in neighbours[pending.pop()]:
            if vertex not in seen:
                seen.add(vertex); pending.append(vertex)
    if len(seen) != raw.ports:
        raise ArithmeticError('the joined federation is not connected')
    return Federation(level, cells, int(raw.ports), seams, neighbours, inter)


def cell_parents(coarse, fine):
    """Actual retained cell ancestry, independent of any port-frame extension."""
    if (type(coarse) is not int or type(fine) is not int
            or not 0 <= coarse <= fine <= 2):
        raise ValueError('declared integer cell levels 0<=coarse<=fine<=2 required')
    tower = _source().build_geodesic_icosahedral_tower(fine)
    return tuple(map(int, tower.cell_ancestor_ids(fine_level=fine, coarse_level=coarse)))


def _occupied(occupied):
    if (type(occupied) is not tuple or any(type(p) is not int or p < 0 for p in occupied)
            or tuple(sorted(set(occupied))) != occupied):
        raise ValueError('a canonical tuple of occupied binary port indices is required')
    return occupied


def _seam(seam):
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(p) is not int for p in seam) or not 0 <= seam[0] < seam[1]):
        raise ValueError('two distinct ordered nonnegative seam indices are required')
    return seam


def binary_step(occupied, seam, coin):
    """Execute exactly the existing native law on endpoints, leaving spectators.

    The caller supplies its declared seam menu. This function does not invent
    a seam or require the particles to occupy one carrier's total-two sector.
    """
    occupied, (a, b) = _occupied(occupied), _seam(seam)
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError('one exact native fair-coin value is required')
    endpoints = (int(a in occupied), int(b in occupied))
    target, (_, current) = native.native(endpoints, (0, 1), coin)
    result = set(occupied)-{a, b}
    result.update(p for p, value in zip((a, b), target, strict=True) if value)
    result = tuple(sorted(result))
    if len(result) != len(occupied):
        raise ArithmeticError('native binary repair did not conserve total occupation')
    return result, current


def local_integer_kraus():
    """Six unscaled integer Kraus matrices; actual amplitudes divide by sqrt2."""
    basis = tuple(product((0, 1), repeat=2))
    matrices = [np.zeros((4, 4), dtype=np.int64) for _ in LETTERS]
    for column, state in enumerate(basis):
        for coin in (0, 1):
            target, letter = native.native(state, (0, 1), coin)
            matrices[LETTERS.index(letter)][basis.index(target), column] = 1
    return tuple(tuple(tuple(map(int, row)) for row in matrix) for matrix in matrices)


def local_audit():
    matrices = tuple(np.asarray(matrix, dtype=np.int64) for matrix in local_integer_kraus())
    identity = np.eye(4, dtype=np.int64)
    swap = identity[:, [0, 2, 1, 3]]
    classical = sum(matrix*matrix for matrix in matrices)
    if (not np.array_equal(sum(m.T@m for m in matrices), 2*identity)
            or not np.array_equal(sum(m@m.T for m in matrices), 2*identity)
            or not np.array_equal(classical, identity+swap)):
        raise ArithmeticError('native all-input CPTP/unital/exchange identity failed')
    number = np.diag((0, 1, 1, 2))
    if any(np.any(number@m-m@number) for m in matrices):
        raise ArithmeticError('a native outcome changes total particle number')
    records = []
    for x, y, coin in product((0, 1), repeat=3):
        occupied = tuple(p for p, v in ((0, x), (1, y), (2, 1)) if v)
        target, current = binary_step(occupied, (0, 1), coin)
        if (int(0 in target)-x != -current or int(1 in target)-y != current
                or 2 not in target):
            raise ArithmeticError('native current continuity or spectator preservation failed')
        records.append({'initial_endpoints': [x, y], 'coin': coin, 'current': current,
                        'final_endpoints': [int(0 in target), int(1 in target)]})
    return {'complete_local_cases': records,
        'unscaled_kraus': [[list(row) for row in matrix] for matrix in local_integer_kraus()],
        'sum_KdaggerK_equals_I_exact': True, 'sum_KKdagger_equals_I_exact': True,
        'each_K_commutes_with_total_number_exact': True,
        'diagonal_transition_is_half_identity_plus_swap': True,
        'coherence_within_committed_outcomes_retained': True,
        'arbitrary_untouched_spectator_extension': 'tensor identity',
        'global_sector_generator': 'G_F=sum_edges (I-SWAP_ab)',
        'global_uniform_attempt_transfer': 'K_F=I-G_F/(2*number_of_seams)',
        'generator_and_transfer_domain': 'classical configuration/history Hilbert space, not the native quantum density-channel generator',
        'stationary_states': ['uniform density on each fixed-total binary sector',
            'Bernoulli product diagonal state at every declared density 0<=rho<=1'],
        'stationarity_argument': 'every local map is unital and preserves particle number; Bernoulli weights are constant within each total sector',
        'state_density_or_total_selected_by_native_rule': False,
        'local_total_two_tensor_sector_preserved_across_join': False}


def tape_data():
    graph = federation(0)
    a, b = graph.inter_seams[0]
    c = min(v for v in graph.neighbours[a] if v//12 == a//12)
    d = min(v for v in graph.neighbours[b] if v//12 == b//12)
    vertices = tuple(sorted((a, b, c, d)))
    tape = ((a, b), tuple(sorted((a, c))), tuple(sorted((b, d))), (a, b))
    if len(set(vertices)) != 4 or any(seam not in graph.seams for seam in tape):
        raise ArithmeticError('the declared crossing probe is not on the source federation')
    basis = tuple(combinations(vertices, 2))
    return graph, vertices, tape, basis, tuple(sorted((a, c)))


def tape_program():
    graph, vertices, tape, basis, initial = tape_data()
    events = []
    for step, seam in enumerate(tape):
        matrices = [np.zeros((6, 6), complex) for _ in LETTERS]
        for column, occupied in enumerate(basis):
            for coin in (0, 1):
                target, current = binary_step(occupied, seam, coin)
                # This exact membership test is the invariant-support proof;
                # no matrix compression discards a populated global address.
                if target not in basis:
                    raise ArithmeticError('native tape leaks outside its declared global support')
                matrices[LETTERS.index((coin, current))][basis.index(target), column] = 1/math.sqrt(2)
        instrument = runtime.Instrument(kraus=tuple(tuple(matrices[LETTERS.index((coin, current))]
            for coin in (0, 1)) for current in CURRENT_VALUES), visible=True,
            resources=(('native_attempts', 1), ('retained_current_records', 1)),
            binding={'global_seam': list(seam), 'outcomes': list(CURRENT_VALUES),
                     'unobserved_coin_alternatives_are_summed_as_CP_maps': True})
        events.append(runtime.Event(f'native_federation_attempt_{step+1}', instruments=(instrument,)))
    density = np.zeros((6, 6)); density[basis.index(initial), basis.index(initial)] = 1
    return runtime.Program('four_native_federation_reads', 6, density, tuple(events),
        preparation_resources=(('disclosed_basis_preparations', 1),),
        binding={'global_ports': graph.ports, 'global_total': 2,
                 'global_basis_support': [list(pair) for pair in basis],
                 'all_other_global_configurations_have_exact_zero_amplitude': True,
                 'initial_occupied_ports_disclosed': list(initial),
                 'scheduler': 'declared four-seam tape; not the uniform global-attempt schedule'})


def tape_execution():
    graph, vertices, tape, basis, initial = tape_data()
    program = tape_program()
    history = runtime.execute(program)
    prefix = runtime.execute(program, stop_after=2)
    restored = runtime.checkpoint_from_json(program, runtime.checkpoint_to_json(prefix))
    resumed = runtime.execute(program, checkpoint=restored)
    expected = {}
    first_cell, second_cell = tape[0][0]//12, tape[0][1]//12
    occupation_counts = Counter()
    crossing_counts = Counter()
    disconnected_counts = Counter()
    for coins in product((0, 1), repeat=4):
        occupied, word, field = initial, [], {p: int(p in initial) for p in vertices}
        crossing = 0
        for seam, coin in zip(tape, coins, strict=True):
            occupied, current = binary_step(occupied, seam, coin)
            word.append(CURRENT_VALUES.index(current))
            field[seam[0]] -= current; field[seam[1]] += current
            if seam[0]//12 != seam[1]//12:
                crossing += current
            if field != {p: int(p in occupied) for p in vertices}:
                raise ArithmeticError('public incident records fail to reconstruct disclosed fields')
        expected.setdefault(tuple(word), Counter())[occupied] += 1
        occupation_counts[sum(p//12 == second_cell for p in occupied)] += 1
        crossing_counts[crossing] += 1
        disconnected = initial
        for seam, coin in zip(tape, coins, strict=True):
            if seam[0]//12 == seam[1]//12:
                disconnected, _ = binary_step(disconnected, seam, coin)
        disconnected_counts[sum(p//12 == second_cell for p in disconnected)] += 1
    if occupation_counts != crossing_counts:
        raise ArithmeticError('regional particle change differs from retained boundary currents')
    if disconnected_counts != {0: 16}:
        raise ArithmeticError('removing the declared inter-carrier seams failed its control')
    maximum = 0.
    record_rows = []
    for branch in history.branches:
        target = np.zeros((6, 6))
        if branch.record in expected:
            for occupied, count in expected[branch.record].items():
                i = basis.index(occupied); target[i, i] = count/16
        maximum = max(maximum, float(np.linalg.norm(branch.density-target, 'fro')))
        if dict(branch.resources) != {'disclosed_basis_preparations': 1, 'native_attempts': 4,
                                      'retained_current_records': 4}:
            raise ArithmeticError('native federation execution resources changed')
        record_rows.append({'word': list(branch.record),
            'exact_probability': str(F(sum(expected.get(branch.record, {}).values()), 16)),
            'executed_weight_diagnostic': branch.probability})
    distance = runtime.complete_distance(history, resumed)
    if len(history.branches) != 3**4 or max(maximum, distance) > 2e-11:
        raise ArithmeticError('complete native federation records or checkpoint changed')
    return {'global_seam_tape': [list(s) for s in tape],
        'first_and_second_carriers': [first_cell, second_cell],
        'initial_occupied_ports_disclosed': list(initial),
        'invariant_global_basis_support': [list(pair) for pair in basis],
        'full_global_fixed_total_dimension': math.comb(graph.ports, 2),
        'executed_exact_invariant_support_dimension': 6,
        'support_embedding_valid_for_every_superposition_and_reference': True,
        'all_other_global_terminal_matrix_entries_exactly_zero': True,
        'complete_record_words': len(history.branches), 'positive_words': len(expected),
        'all_native_coin_trajectories_enumerated': 16,
        'coin_is_probability_summed_inside_each_current_outcome': True,
        'complete_record_law': record_rows,
        'second_carrier_final_occupation_law': {str(k): str(F(v, 16)) for k, v in sorted(occupation_counts.items())},
        'same_law_read_from_net_boundary_current': True,
        'mean_second_carrier_occupation': str(sum((F(k*v, 16) for k, v in occupation_counts.items()), F())),
        'disconnected_control_second_carrier_occupation': '0',
        'full_embedded_terminal_density_residual_diagnostic': maximum,
        'serialized_checkpoint_complete_distance_diagnostic': distance,
        'initial_loads_disclosed_but_hidden_evolving_loads_not_observer_inputs': True,
        'floating_roundoff_or_physical_service_certified': False}


def run():
    verify_source_pins()
    graphs = []
    for level in (0, 1, 2):
        graph = federation(level)
        graphs.append({'level': level, 'elementary_carriers': graph.cells, 'ports': graph.ports,
            'intra_carrier_seams': 30*graph.cells, 'inter_carrier_seams': len(graph.inter_seams),
            'total_seams': len(graph.seams), 'max_degree': graph.max_degree,
            'degree_census': {str(k): v for k, v in sorted(Counter(map(len, graph.neighbours)).items())},
            'complete_seams_sha256': hashlib.sha256(canonical(graph.seams).encode()).hexdigest(),
            'fair_uniform_attempt_denominator': 2*len(graph.seams),
            'connected': True, 'local_carrier_seams_preserved': True,
            'port_pair_gluing_and_spherical_support_supplied': True})
    return {'graph_family': graphs, 'local_operator_attachment': local_audit(),
        'complete_joined_observer_execution': tape_execution(),
        'scope': {'primitive_observers_are_federation_members': True,
            'global_screen_and_local_twelve_port_boundary_are_distinct': True,
            'global_state_is_tensor_of_local_total_two_sectors': False,
            'gluing_schedule_and_preparation_selected_by_OPH': False,
            'continuum_spacetime_or_laboratory_match_established': False}}
