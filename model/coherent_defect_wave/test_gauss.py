"""Complete interior constraint and exact central-gauge covariance checks."""
from fractions import Fraction as F

from .gauss import certificate, gamma_amplitudes, project_amplitudes, gauss_measurement


def norm(values):
    return sum(value.conjugate()*value for value in values.values())


def test_complete_gauss_projection_retains_both_sectors():
    # The superposition includes charged matter at the constrained vertex.
    state = {(1, 5, 3, 0, 0): F(3, 5), (1, 0, 1, 1, 20): F(4, 5)}
    assert gamma_amplitudes(gamma_amplitudes(state)) == state
    positive = project_amplitudes(state)
    negative = project_amplitudes(state, positive=False)
    assert gamma_amplitudes(positive) == positive
    assert gamma_amplitudes(negative) == {key: -value for key,value in negative.items()}
    assert norm(positive) == norm(negative) == F(1, 2)
    assert norm(positive)+norm(negative) == norm(state) == 1
    assert project_amplitudes(positive, positive=False) == {}
    assert project_amplitudes(negative) == {}


def test_physical_reference_uses_its_own_dimension_and_trace():
    result = certificate()
    assert result['physical_dimension'] == 36000
    assert result['physical_unitary_trace'] == '-23925'
    assert result['physical_reference_normalization'] == '21075'
    assert result['physical_reference_eigenvalue_lower'] == '1/84300'
    assert result['joint_success_from_uniform_ambient_input'] == '281/2160'


def test_same_executed_tick_preserves_the_constraint_and_wave_contrast():
    from .model import JointObserver, SparseState
    observer = JointObserver()

    def gamma(state):
        return SparseState(gamma_amplitudes(state.as_dict()), state.norm_denominator)

    probabilities = []
    address = observer.addresses.index((0, 5))
    for word in (1, 3):
        state = observer.embedded_ports(word, (0, 5), internal=3)
        assert observer.step(gamma(state)) == gamma(observer.step(state))
        preparation = gauss_measurement(observer, state)
        assert [record['outcome'] for record, _ in preparation] == [1, -1]
        assert sum(branch.norm_squared() for _, branch in preparation) == 1
        plus = preparation[0][1]
        evolved = observer.step(plus)
        assert gamma(evolved) == evolved
        assert plus.norm_squared() == evolved.norm_squared() == F(1, 2)
        branches = observer.wave_basis_measurement(evolved, 0, internal=3, address=address)
        assert sum(branch.norm_squared() for _, branch in branches) == F(1, 2)
        assert all(gamma(branch) == branch for _, branch in branches)
        probabilities.append(branches[0][1].norm_squared()/evolved.norm_squared())
        assert observer.inverse_step(evolved) == plus
    assert probabilities == [F(1, 120), 0]


def test_same_executed_tick_preserves_gauss_projected_defect_motion():
    from .model import JointObserver, SparseState
    observer = JointObserver()
    results = []
    target = observer.law.faces.index((4, 5, 11))
    for port in (5, 6):
        initial = observer.inverse_wave(observer.embedded_port(1, port, internal=3))
        plus = SparseState(project_amplitudes(initial.as_dict()), initial.norm_denominator)
        evolved = observer.step(plus)
        assert gamma_amplitudes(evolved.as_dict()) == evolved.as_dict()
        # Actual face instrument preserves coherence between both gauge representatives.
        branches = observer.region_face_measurement(evolved)
        assert len(branches) == 5
        assert sum(branch.norm_squared() for _, branch in branches) == evolved.norm_squared()
        assert all(gamma_amplitudes(branch.as_dict()) == branch.as_dict() for _, branch in branches)
        yes = next(branch for record, branch in branches if record['face_index'] == target)
        results.append(yes.norm_squared()/evolved.norm_squared())
    assert results == [1, 0]
