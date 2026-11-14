"""Attach the existing graph mixing theorem to the actual observer federation.

Counts below are discarded native attempts, not modular parameters or physical
durations. A bound selects a finite preparation procedure, not a unique total
particle number, vacuum or nontrivial modular reference.
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction as F
from itertools import combinations, product
import hashlib
import math
from pathlib import Path

import numpy as np

from ..federated_observer import source
from ..observer import native

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*source.SOURCE_PATHS,
    'model/observer/QUANTUM_MIXING.md', 'model/observer/quantum_mixing.py',
    'model/native_equilibration/construction.py',
    'model/native_equilibration/__init__.py', 'model/native_equilibration/__main__.py',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py',
    'model/requirements.txt'))))


def graph_certificate(vertices, seams):
    """Exact connectivity/radius certificate; supplies a conservative gap bound."""
    if type(vertices) is not int or vertices < 2:
        raise ValueError('at least two integer vertices required')
    if (type(seams) is not tuple or not seams
            or any(type(e) is not tuple or len(e) != 2
                   or any(type(v) is not int for v in e)
                   or not 0 <= e[0] < e[1] < vertices for e in seams)
            or len(set(seams)) != len(seams)):
        raise ValueError('distinct ordered edges of a finite simple graph required')
    neighbours = [[] for _ in range(vertices)]
    for a, b in seams:
        neighbours[a].append(b); neighbours[b].append(a)
    distances = [-1]*vertices
    distances[0] = 0
    todo = deque([0])
    while todo:
        a = todo.popleft()
        for b in neighbours[a]:
            if distances[b] < 0:
                distances[b] = distances[a]+1
                todo.append(b)
    if min(distances) < 0:
        raise ValueError('connected graph required for state selection')
    diameter_upper = min(vertices-1, 2*max(distances))
    # lambda_1 >= 4/(V diameter), by the range/variance/path argument.
    gap = F(2, vertices*diameter_upper*len(seams))
    return {'vertices': vertices, 'seams': len(seams),
            'root_eccentricity': max(distances),
            'diameter_upper_bound': diameter_upper,
            'quantum_gap_lower_per_attempt': str(gap)}


def preparation_budget(dimension, gap_lower, precision_bits):
    """Integer certificate for trace distance <=2^-precision_bits, with reference.

    If b=ceil(1/gamma), Bernoulli's inequality gives
    (1-gamma)^b <= 1/(1+b*gamma) <= 1/2. No huge rational power,
    floating logarithm or uncertified numerical eigenvalue is used.
    """
    if type(dimension) is not int or dimension < 1:
        raise ValueError('positive integer sector dimension required')
    if type(gap_lower) is not F or not 0 < gap_lower <= 1:
        raise ValueError('exact positive Fraction gap at most one required')
    if type(precision_bits) is not int or not 1 <= precision_bits <= 1000:
        raise ValueError('integer precision between one and 1000 bits required')
    if dimension == 1:
        return {'attempts': 0, 'half_contraction_block': 0, 'blocks': 0,
                'dimension': 1, 'trace_distance_upper': '0'}
    block = (gap_lower.denominator+gap_lower.numerator-1)//gap_lower.numerator
    blocks = dimension.bit_length()+precision_bits
    bound = F(dimension, 2**(blocks+1))
    if block*gap_lower < 1 or bound > F(1, 2**precision_bits):
        raise ArithmeticError('integer preparation certificate failed')
    return {'attempts': block*blocks, 'half_contraction_block': block,
            'blocks': blocks, 'dimension': dimension,
            'trace_distance_upper': str(bound)}


def exact_native_channel(vertices, seams, particles, law):
    """Small-graph oracle: source Kraus enumeration, no block/gap formula."""
    graph_certificate(vertices, seams)
    if type(particles) is not int or not 0 <= particles <= vertices:
        raise ValueError('binary integer particle count in graph range required')
    if law not in ('C', 'E'):
        raise ValueError('native candidate C or E required')
    occupied = tuple(combinations(range(vertices), particles))
    if len(occupied) > 70:
        raise ValueError('small-graph oracle requires at most 70 configurations')
    basis = tuple(tuple(int(i in state) for i in range(vertices)) for state in occupied)
    indices = {x: i for i, x in enumerate(basis)}
    d = len(basis)
    counts = np.zeros((d*d, d*d), dtype=np.int64)
    for seam in seams:
        for coin in (0, 1):
            images = []
            for x in basis:
                y, (_, current) = native.native(x, seam, coin)
                label = current if law == 'C' else x[seam[0]]-x[seam[1]]
                images.append((indices[tuple(y)], label))
            for i, (a, label_a) in enumerate(images):
                for j, (b, label_b) in enumerate(images):
                    if label_a == label_b:
                        counts[a*d+b, i*d+j] += 1
    return counts, 2*len(seams), basis


def primitive_execution():
    """Actual native Kraus application from a nonstationary coherent state.

    Applies the exact integer-enumerated superoperator, then compares the
    declared pointer law with its tracial-preparation counterpart. Arithmetic
    residuals are diagnostics, separate from the analytic error certificate.
    """
    from ..primitive_observer import transfer
    seams = tuple(transfer.native_counts().seams)
    # Existing exact primitive gap is (5-sqrt(5))/60 >= 11/240.
    gap = F(11, 240)
    # (9/4)^2 > 5 certifies the sole radical comparison above.
    assert F(9, 4)**2 > 5
    budget = preparation_budget(66, gap, 10)
    # The two source states differ coherently, not just in populations.
    # Use the source channel's native matrix-unit blocks, avoiding a dense
    # 4356-square operator or exponentiation of a supplied Hamiltonian.
    from ..observer import quantum_mixing as qm
    binary = tuple(i for i, yes in enumerate(qm.BINARY) if yes)
    lookup = {i: j for j, i in enumerate(binary)}
    ket = np.zeros(66, complex)
    ket[0], ket[1] = 1/math.sqrt(2), 1j/math.sqrt(2)
    initial = np.outer(ket, ket.conj())
    rows = []
    for law in ('C', 'E'):
        final = np.zeros((66, 66), complex)
        for sector in qm.sectors(binary=True):
            coordinates = tuple((lookup[i], lookup[j]) for i, j in sector.pairs)
            v = np.array([initial[i, j] for i, j in coordinates])
            if not np.any(v):
                continue
            p = np.asarray(qm.native_block(law, sector.difference, binary=True), float)
            # Repeated native channel applications are the execution being
            # certified. The preparation record is explicitly discarded.
            for _ in range(budget['attempts']):
                v = p@v
            for (i, j), value in zip(coordinates, v, strict=True):
                final[i, j] = value
        error = final-np.eye(66)/66
        trace_distance = float(np.sum(np.abs(np.linalg.eigvalsh(error)))/2)
        normalization = abs(float(np.trace(final).real)-1)
        if normalization > 2e-12 or trace_distance > float(F(budget['trace_distance_upper']))+2e-12:
            raise ArithmeticError('native preparation failed its analytic error allowance')
        # Joint probabilities of two consecutive source current bits.
        # Direct diagonal native enumeration is independent of spectral mixing.
        basis = tuple(tuple(qm.STATES[i]) for i in binary)
        index = {x: i for i, x in enumerate(basis)}
        marked = {bit: np.zeros((66, 66)) for bit in (0, 1)}
        for j, state in enumerate(basis):
            for seam, coin in product(seams, (0, 1)):
                target, (_, current) = native.native(state, seam, coin)
                marked[int(current != 0)][index[tuple(target)], j] += 1/60
        laws = []
        for a, b in product((0, 1), repeat=2):
            effect = np.ones(66)@marked[b]@marked[a]
            laws.append({'bits': [a, b],
                         'prepared_probability': float(effect@np.diag(final).real),
                         'stationary_probability': float(effect@np.full(66, 1/66))})
        tv = sum(abs(r['prepared_probability']-r['stationary_probability']) for r in laws)/2
        if tv > trace_distance+2e-12:
            raise ArithmeticError('future instrument violates trace-distance contraction')
        rows.append({'law': law, 'discarded_attempts': budget['attempts'],
                     'trace_distance_diagnostic': trace_distance,
                     'normalization_residual': normalization,
                     'two_current_bit_law': laws, 'future_law_TV_diagnostic': tv})
    return {'budget': budget, 'gap_lower_exact': str(gap), 'executions': rows,
            'records_during_preparation': 'discarded, not conditioned upon',
            'future_protocol': 'two complete native absolute-current bits'}


def run():
    rows = []
    for level in (0, 1, 2):
        graph = source.federation(level)
        cert = graph_certificate(graph.ports, graph.seams)
        for particles in (2, graph.ports//2):
            dimension = math.comb(graph.ports, particles)
            rows.append({'level': level, 'particles': particles, **cert,
                         'preparation': preparation_budget(dimension,
                             F(cert['quantum_gap_lower_per_attempt']), 10)})
    return {'schema': 'native-equilibration-v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                          for p in SOURCE_PATHS},
        'theorem_attachment': {
            'existing_proof': 'model/observer/QUANTUM_MIXING.md equations (1)-(4)',
            'external_input': 'Caputo-Liggett-Richthammer (2010), Section 4.1.1 equation (4.1)',
            'source_url': 'https://www.stat.berkeley.edu/~aldous/206-RWG/RWGpapers/caputo.pdf',
            'native_domain': 'connected binary fixed-total sector, IID uniform seam and fair coin, C or E',
            'stationary_state': 'I/binomial(V,k), unique within the declared fixed-k sector',
            'reference_uniform_bound': 'distance((Phi^n tensor id)(rho), I/d tensor rho_R) <= min(1,d/2*(1-gamma)^n)',
            'record_comparison': 'any common future complete instrument, including adaptive interventions',
            'conditional_postselection_bound_unchanged': False,
            'state_total_or_sector_weights_selected': False,
            'finite_exact_stationarity_claimed': False,
            'tracial_modular_flow': 'identity within each full fixed-total matrix algebra',
            'nontrivial_modular_reference_selected': False,
            'attempt_parameter_is_physical_time': False},
        'federation_preparation_certificates': rows,
        'primitive_native_preparation': primitive_execution()}
