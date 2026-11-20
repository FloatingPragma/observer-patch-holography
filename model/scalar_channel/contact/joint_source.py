"""The shared 150-mode source of the material clock and radar experiments.

Six empty-load collars are in their symmetric written eigenstates. Each
contributes +1 times its scalar seam occupation, with its constant energy
removed. A distinct collar supplies the marked contact. The two retained
arms are an invariant subspace of the complete 99-state local memory law.
The scalar modes are the complete one-particle sector, without a cage
projection or a filled negative-energy reference.
One inert occupied reference mode outside these modes supplies even total
parity for the scalar compiler's code sector.
"""

from functools import lru_cache

import numpy as np

from ...polynomial_sv import law as contrast
from ..carrier import load_source, require
from ..compiler import build
from ..cage import construction as cage
from .law import arm_energies, collar_states, interaction

ENDPOINTS = (0, 5, 11, 16, 22, 29)
MARK = 2
MODES = 150


def validate_source(protocol):
    require(all(type(seam) is int for seam in protocol['endpoints'])
            and tuple(protocol['endpoints']) == ENDPOINTS,
            'the six written collar seams differ from the declared construction')
    require(type(protocol['mark']) is int and protocol['mark'] == MARK,
            'the marked collar must be seam 2, distinct from the six stations')
    require(type(protocol['kappa']) in (int, float) and protocol['kappa'] == 8,
            'the declared native collar coupling is 8')


@lru_cache(maxsize=1)
def _mapped_source():
    """Check the independently generated carrier and the native number words."""
    car = cage.cell()
    D = cage.law(car)
    data = load_source('cell')
    graph = build()
    require(D.shape == (MODES, MODES) and graph['nsource'] == MODES,
            'the source must retain all 150 one-particle modes')
    require(np.array_equal(car['edges'], data['edges'])
            and np.array_equal(car['wedges'], data['wedges'])
            and np.array_equal(data['edges'], graph['data']['edges'])
            and np.array_equal(data['wedges'], graph['data']['wedges']),
            'seam or wedge ordering differs across the source maps')
    require(np.max(abs(D - data['D'])) < 1e-14, 'scalar source maps disagree')
    compiled = np.zeros((graph['n'], graph['n']))
    compiled[np.arange(MODES), np.arange(MODES)] = D.diagonal()
    for (a, b), value in graph['source_coeff'].items():
        compiled[a, b] = compiled[b, a] = value
    reference = graph['R'][0]
    require(reference >= MODES and not np.any(compiled[reference])
            and np.max(abs(compiled[:MODES, :MODES] - D)) < 1e-14,
            'the compiled source does not preserve its inert parity reference')
    projectors = {}
    for seam in ENDPOINTS + (MARK,):
        P = interaction(graph, car['edges'][seam], 1)
        expected = np.zeros_like(D)
        expected[seam, seam] = 1
        require(np.max(abs(P - expected)) < 1e-14,
                'the native number words do not map to the declared seam')
        projectors[seam] = P
    return D, projectors


def radar_source(protocol, coupling):
    """D plus six fixed +1 potentials and a signed marked potential.

    ``coupling`` changes only the marked collar. Setting it to zero or minus
    one does not suppress or reverse the six written endpoint contacts.
    """
    validate_source(protocol)
    require(type(coupling) in (int, float) and coupling in (0, 1, -1),
            'the marked contact must be 0, +1 or -1')
    D, projectors = _mapped_source()
    result = D.copy()
    empty = arm_energies(protocol['kappa'], 1, 0, 1)
    filled = arm_energies(protocol['kappa'], 1, 1, 1)
    shift = filled[1] - empty[1]
    require(abs(shift - 1) < 1e-13 and abs(filled[0] - empty[0]) < 1e-13,
            'the native written collar does not supply its +1 occupation term')
    for seam in ENDPOINTS:
        result += shift * projectors[seam]
    marked = arm_energies(protocol['kappa'], coupling, 1, 1) - empty
    # ``empty`` is independent of contact strength because filling is zero.
    require(abs(marked[1] - coupling) < 1e-13 and abs(marked[0]) < 1e-13,
            'the native marked collar does not supply the declared contact')
    result += marked[1] * projectors[MARK]
    return result


def arm_sources(protocol, coupling, sign):
    """Both free collar arms, including their native relative energy."""
    require(sign in (1, -1), 'collar sign must be +1 or -1')
    energies = arm_energies(protocol['kappa'], coupling, 0, sign)
    shifts = arm_energies(protocol['kappa'], coupling, 1, sign) - energies
    require(np.max(abs(shifts - [0, coupling])) < 1e-13,
            'the native contact differs from the two-arm reduction')
    eye = np.eye(MODES)
    return (radar_source(protocol, 0) + energies[0] * eye,
            radar_source(protocol, coupling) + energies[1] * eye)


def source_certificate(protocol):
    """Numerical witnesses for the native invariant subspace and its source map."""
    validate_source(protocol)
    D, projectors = _mapped_source()
    exchange, square = contrast.local_term_parts()
    reduction_error = 0.0
    energies = {}
    for sign in (1, -1):
        arms = collar_states(sign)
        energies[str(sign)] = {}
        for coupling in (0, 1, -1):
            table = []
            for filling in (0, 1):
                reduced = arm_energies(8, coupling, filling, sign)
                actual = (exchange + (8 + coupling * filling) * square) @ arms
                reduction_error = max(reduction_error,
                                      float(np.max(abs(actual - arms @ np.diag(reduced)))))
                table.append(reduced.tolist())
            energies[str(sign)][str(coupling)] = table
    base = radar_source(protocol, 0)
    expected = D.copy()
    expected[list(ENDPOINTS), list(ENDPOINTS)] += 1
    mapping_error = float(np.max(abs(base - expected)))
    for coupling in (0, 1, -1):
        mapping_error = max(mapping_error, float(np.max(abs(
            radar_source(protocol, coupling) - base - coupling * projectors[MARK]))))
    return dict(modes=MODES, scalar_sector='one particle on all seams and wedges',
                parity_reference='one occupied inert reference mode R[0] outside the active 150 modes',
                endpoint_seams=list(ENDPOINTS), marked_seam=MARK,
                removed_common_endpoint_energy=float(6 * arm_energies(8, 0, 0, 1)[1]),
                native_arm_energies=energies, invariant_subspace_max_error=reduction_error,
                scalar_mapping_max_error=mapping_error,
                pointer_service='complete orthogonal scalar-mode record and collar fuel record',
                pointer_outcomes=2 * MODES,
                controls_change='marked contact only; all six endpoint potentials stay +1')
