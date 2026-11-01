"""Direct shared-domain checks, without identifying distinct physical references."""
from fractions import Fraction as F
import numpy as np

from model.modular_flow.flow import ModularFlow
from model.full_fock_energy.energy import local_energies
from model.field_protocols.protocol import Read, instrument
from model.particle_response.dimer import setup


def test_energy_and_modular_services_share_the_complete_hamiltonian():
    flow = ModularFlow()
    energy = local_energies(flow.geometry)
    np.testing.assert_array_equal(energy.total.matrix().toarray()/energy.scale,
                                  flow.thermal.target()['hamiltonian'])
    assert flow.duration(F(-1,8)) == F(1,8)


def test_modular_and_adaptive_reads_share_the_same_instrument_at_matched_time():
    flow = ModularFlow()
    direct, ideal, aborts, _ = instrument(flow.geometry, Read(F(1,8), 0, 4))
    run, _ = flow.certificate(F(-1,8), 4)
    d = len(direct[0])
    # The pointer is an explicitly compiled source operation; comparison
    # with exact projectors permits only its floating execution residual.
    from model.full_fock_field.field import FullFock
    number = FullFock(flow.geometry.modes).number(0).toarray()
    for actual, projector in zip(direct, (np.eye(d)-number, number), strict=True):
        np.testing.assert_allclose(actual, projector@run.success, atol=2e-11)
    np.testing.assert_allclose(ideal[1], number@flow.unitary(F(-1,8)), atol=1e-14)
    assert aborts


def test_pure_particle_ground_is_not_the_faithful_modular_reference():
    thermal = ModularFlow().thermal.target()['reference']
    particle = setup().reference
    assert np.linalg.eigvalsh(thermal)[0] > 0.01
    assert np.linalg.matrix_rank(particle, tol=1e-10) == 1
    assert np.linalg.norm(thermal-particle, 2) > .5
