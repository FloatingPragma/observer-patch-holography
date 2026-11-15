"""One interior Gauss constraint on the finite joint link/wave register.

Sparse keys are (link word, port, internal level, parity, seam address).
The transform and the complete projection preserve any external common
amplitude normalization used by the caller.
"""
from fractions import Fraction as F

from model.coherent_seam_repair.protection import BoundaryFlux
from .preparation import fixed_word_counts


def gamma_amplitudes(amplitudes):
    law = BoundaryFlux()
    star = law.vertex_stars[5]
    return {
        (g ^ star, port, internal, parity, address):
            (-value if port == 5 and internal >= 3 else value)
        for (g, port, internal, parity, address), value in amplitudes.items()
        if value != 0
    }


def project_amplitudes(amplitudes, positive=True):
    """Unnormalized branch of the complete {P_plus,P_minus} instrument."""
    if type(positive) is not bool:
        raise ValueError('Gauss outcome must be a boolean')
    transformed = gamma_amplitudes(amplitudes)
    sign = 1 if positive else -1
    result = {key: F(1, 2)*(amplitudes.get(key, F(0)) + sign*transformed.get(key, F(0)))
              for key in set(amplitudes) | set(transformed)}
    return {key: value for key, value in result.items() if value != 0}


def gauss_measurement(observer, state):
    """Execute both constraint outcomes without postselection or discarded weight."""
    from .model import SparseState
    observer._validate(state)
    return tuple(({'kind': 'interior_gauss_readout', 'vertex': 5, 'outcome': sign},
                  SparseState(project_amplitudes(state.as_dict(), positive=positive),
                              state.norm_denominator))
                 for sign, positive in ((1, True), (-1, False)))


def certificate():
    law = BoundaryFlux()
    words, counts = fixed_word_counts()
    star = law.vertex_stars[5]
    assert star == 6553634
    assert all(g ^ star in words and g ^ star != g for g in words)
    assert all(law.syndrome(g ^ star) == law.syndrome(g) for g in words)
    for g in words:
        for edge in law.interior:
            assert law.transport(g ^ star, edge)[0] == law.transport(g, edge)[0] ^ star
        for edge, (a, b) in enumerate(law.edges):
            before = -1 if g >> edge & 1 else 1
            after = -1 if (g ^ star) >> edge & 1 else 1
            assert after == before * (-1 if a == 5 else 1) * (-1 if b == 5 else 1)
    matches = []
    for port in range(12):
        count = 0
        for g in words:
            result = g
            for edge in law.interior:
                if port in law.edges[edge]:
                    result = law.transport(result, edge)[0]
            count += result == g ^ star
        matches.append(count)
    assert matches == [0]*12
    # Gamma has no diagonal link words; Gamma U has no diagonal link paths.
    dimension = 36000
    trace = F(-1595, 3)*sum(counts)/2
    normalization = F(5, 4)*dimension + trace
    assert trace == -23925 and normalization == 21075
    return {
        'schema': 'oph.coherent_defect_wave.interior_gauss.v1',
        'interior_vertex': 5,
        'star_mask': star,
        'gauge_orbits': len(words)//2,
        'all_local_transport_covariance_checks': len(words)*len(law.interior),
        'all_central_link_covariance_checks': len(words)*len(law.edges),
        'gamma_U_diagonal_link_paths_by_port': matches,
        'physical_dimension': dimension,
        'physical_unitary_trace': str(trace),
        'physical_reference_normalization': str(normalization),
        'physical_reference_eigenvalue_lower': str(1/(4*normalization)),
        'herald_success_from_uniform_physical_input': str(4*normalization/(9*dimension)),
        'positive_gauss_probability_from_uniform_ambient_input': '1/2',
        'joint_success_from_uniform_ambient_input': '281/2160',
        'scope': {'one_interior_constraint_preserved': True,
                  'all_boundary_gauge_constraints_supplied': False,
                  'continuous_gauge_field_derived': False,
                  'physical_particle_or_mass_identified': False},
    }


def covariant_preparation_return(observer, state, ports, internal=3):
    """Complete conditional-frame return detector preserving the Gauss sector.

    Pick the smaller word in each two-element gauge orbit. On its partner,
    undo the port-5 central phase before the wave return test and restore it
    afterwards. This is a coherent register control, not a link measurement.
    The controlled diagonal frame and this detector are explicit apparatus.
    """
    from .model import SparseState
    observer._validate(state)
    star = observer.law.vertex_stars[5]
    ports = tuple(ports)

    def frame(source):
        return SparseState({key: (-value if key[0] > (key[0] ^ star)
                                  and key[1] == 5 and key[2] >= 3 else value)
                            for key, value in source.amplitudes}, source.norm_denominator)

    return tuple(({'kind': 'covariant_wave_preparation_return', 'ports': list(ports),
                   'internal': internal, 'outcome': record['outcome']}, frame(branch))
                 for record, branch in observer.wave_preparation_measurement(frame(state), ports, internal))
