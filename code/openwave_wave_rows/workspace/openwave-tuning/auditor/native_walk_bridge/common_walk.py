"""Shared LOCAL12 native coherent walk, with explicit padded source/event state.

The geometry and local C/E dilation come from existing native sources.
Only the declared local scheduler and coherent continuation are amended.
"""
from functools import lru_cache
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RER = ROOT / 'reverse-engineering-reality'
sys.path.insert(0, str(RER))
sys.path.insert(0, str(HERE))
from model.native_spatial_refinement import source as spatial

MODEL_ID = 'native-local12-current-v0-proposed'
POSITION_DIMENSION = 128
EVENT_DIMENSION = 16
ACTIVE_POSITIONS = 96
ACTIVE_EVENTS = 12
INTERNAL_MATCHINGS = (
    ((0,1),(2,3),(4,5),(6,7),(8,9),(10,11)),
    ((0,5),(1,7),(2,10),(3,9),(4,11),(6,8)),
    ((0,7),(1,5),(2,11),(3,8),(4,9),(6,10)),
    ((0,10),(1,9),(2,4),(3,6),(5,11),(7,8)),
    ((0,11),(1,8),(2,6),(3,4),(5,9),(7,10)),
)


def encode_position(vertex):
    if type(vertex) is not int or not 0 <= vertex < 96:
        raise ValueError('size-two native vertex required')
    return 16*(vertex//12)+vertex%12


def decode_position(address):
    if type(address) is not int or not 0 <= address < 128:
        raise ValueError('seven-bit position address required')
    cell, port = divmod(address, 16)
    return None if port >= 12 else 12*cell+port


@lru_cache(maxsize=1)
def graph():
    return spatial.graph(2)


@lru_cache(maxsize=1)
def geometry_tables():
    g = graph()
    assert g.ports == 96 and len(g.seams) == 288
    internal = set(map(tuple, spatial.slow_geometry.source().seams()))
    assert {e for row in INTERNAL_MATCHINGS for e in row} == internal
    assert sum(map(len, INTERNAL_MATCHINGS)) == len(internal)
    assert all(sorted(x for edge in row for x in edge) == list(range(12)) for row in INTERNAL_MATCHINGS)
    edge_index = {tuple(edge): index for index, edge in enumerate(g.seams)}
    neighbours = np.zeros((96, 6), dtype=np.int64)
    edge_ids = np.zeros((96, 6), dtype=np.int64)
    for vertex in range(96):
        cell, port = divmod(vertex, 12)
        for direction, matching in enumerate(INTERNAL_MATCHINGS):
            edge = next(e for e in matching if port in e)
            other = edge[1] if edge[0] == port else edge[0]
            neighbours[vertex, direction] = 12*cell+other
        external = [w for w in g.neighbours[vertex] if w//12 != cell]
        assert len(external) == 1
        neighbours[vertex, 5] = external[0]
        for direction, other in enumerate(neighbours[vertex]):
            edge_ids[vertex, direction] = edge_index[tuple(sorted((vertex, int(other))))]
    for vertex in range(96):
        assert set(map(int, neighbours[vertex])) == set(g.neighbours[vertex])
        for direction in range(6):
            assert neighbours[neighbours[vertex, direction], direction] == vertex
    neighbours.flags.writeable = False
    edge_ids.flags.writeable = False
    return neighbours, edge_ids


def transfer_numerators():
    result = 6*np.eye(96, dtype=np.int64)
    for vertex, row in enumerate(geometry_tables()[0]):
        for other in row:
            result[other, vertex] += 1
    return result


def native_transfer():
    """Complete local-attempt population matrix, not repeated-W populations."""
    return transfer_numerators()/12


def expected_swap_permutation():
    """Independent geometric target; never supplied to source pulse execution."""
    answer = np.arange(128*16, dtype=np.int64)
    neighbours, _ = geometry_tables()
    for vertex in range(96):
        for direction in range(6):
            before = encode_position(vertex)*16+2*direction+1
            after = encode_position(int(neighbours[vertex, direction]))*16+2*direction+1
            answer[before] = after
    assert np.array_equal(answer[answer], np.arange(len(answer)))
    return answer


def apply_walk(state):
    """Apply the source-compiled R and compiled reversible S, without H input.

    Numerical pulse errors and preparation/reset complements are reported by
    the owning compilation receipt. This array is the prepared code sector.
    Neither source nor event padding is removed or normalized away.
    """
    from compile_walk import controls
    value = np.asarray(state, dtype=complex)
    if value.shape != (128, 16) or not np.isfinite(value).all():
        raise ValueError('finite padded128-by16 common state required')
    built = controls()
    reflected = value @ built['reflection'].T
    result = np.empty(2048, dtype=complex)
    result[built['swap_permutation']] = reflected.reshape(-1)
    return result.reshape(128, 16)


def initial_cell_state():
    """Source-word preparation from the separately reset blank; no state install."""
    from compile_walk import controls
    prepared = controls()['uniform12_preparation'][:, 0]
    position = np.zeros(128, dtype=complex)
    position[:16] = prepared
    return position[:, None]*prepared[None, :]


def complete_position_records(state):
    value = np.asarray(state, dtype=complex)
    if value.shape != (128, 16) or not np.isfinite(value).all():
        raise ValueError('finite padded128-by16 common state required')
    probabilities = abs(value)**2
    valid = [encode_position(v) for v in range(96)]
    invalid = [p for p in range(128) if decode_position(p) is None]
    return {'joint_source_event_probabilities': probabilities,
            'valid_position_probabilities': probabilities[valid, :12].sum(axis=1),
            'source_padding_probability': float(probabilities[invalid].sum()),
            'event_padding_probability_on_valid_source': float(probabilities[valid, 12:].sum()),
            'original_norm_squared': float(probabilities.sum()),
            'renormalized': False,
            'service_reset_export_failures': 'additional explicit instrument flags, never included in a normalized-success claim'}


def raw_record(final_position, event, law):
    """Recover the source native labels for the final controlled local event.

    Reused coherent event amplitudes are not a fresh uniform native attempt.
    This is the inverse record-recoding measurement convention.
    """
    if type(event) is not int or not 0 <= event < 16 or law not in ('C', 'E'):
        raise ValueError('padded event address and C/E law required')
    after = decode_position(final_position)
    if after is None or event >= 12:
        return {'admitted': False, 'reason': 'source_padding' if after is None else 'event_padding'}
    direction, sigma = divmod(event, 2)
    neighbours, edges = geometry_tables()
    before = int(neighbours[after, direction]) if sigma else after
    edge_id = int(edges[after, direction])
    a, b = graph().seams[edge_id]
    coin = sigma ^ int(before == a)
    current = int(before == a)-int(after == a)
    return {'admitted': True, 'edge_index': edge_id, 'seam': [a, b],
            'source_before': before, 'source_after': after,
            'direction': direction, 'recoded_sigma': sigma,
            'law': law, 'coin': coin, 'current': current,
            'original_E_contrast': (1 if before == a else -1) if law == 'E' else None,
            'fresh_uniform_attempt_claimed': False}


def source_descriptor():
    return {'model_id': MODEL_ID, 'development_size': 2,
            'geometry_source': 'model.native_spatial_refinement.source.graph(2)',
            'state_shape': [128, 16], 'valid_source_positions': 96, 'valid_events': 12,
            'position_encoding': '(cell_index<<4)|local_port; cell_index=4*x+2*y+z; local_port<12',
            'event_encoding': '2*direction+sigma; direction0..4 internal perfect matchings; direction5 external matching',
            'internal_matchings': [[list(edge) for edge in row] for row in INTERNAL_MATCHINGS],
            'reverse_direction': 'same direction because every colour is a perfect matching',
            'step': 'W=S R with R=I_source tensor(2|uniform12><uniform12|-I16)',
            'native_population_compression': 'A* S A = I-L/12',
            'original_global_population_relation': 'T_global=(1-6/M)I+(6/M)T_local, M=288',
            'scheduler_scope': 'Local incident seam selection is a declared amendment; global nonincident attempts remain a separately reconstructed clock/tape comparison.',
            'coherent_continuation_scope': 'Repeated W reuses uncommitted registers; it is not the repeated original unread native quantum channel.',
            'clocks': ['native global or local attempted event', 'ordered M6 source-pulse word and proper-code services', 'coherent W step'],
            'physical_time_or_speed_calibrated': False,
            'field_record_extension': 'Requires the chosen joint preparation and explicit intertwiner; tracing a definite current register need not preserve W.'}
