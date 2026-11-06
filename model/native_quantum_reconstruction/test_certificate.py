from fractions import Fraction as Q

import pytest

from . import certificate as c


def test_native_counts_have_exact_positive_transfer_bounds():
    matrix, denominator = c.source_counts()
    assert denominator == 60
    assert all(sum(row) == 60 for row in matrix)
    assert c.source_spectral_bounds()['lower'] == '2/3'
    assert c.source_spectral_bounds()['minimum_diagonal'] == '5/6'


def test_all_matrix_units_are_generated_from_source_words():
    assert c.occupation_projector_audit()['configuration_projectors'] == 66
    paths = c.matrix_unit_paths()
    assert len(paths) == 4356
    assert max(len(path)-1 for path,_ in paths.values()) == 4
    for (i,j),(path,coefficient) in paths.items():
        assert c.generated_word(path) == {(i,j):coefficient}


def test_direct_endpoint_two_attempt_return_count_agrees():
    # Independent native branch list; no source_counts or matrix power defines it.
    branches = []
    for state in c.STATES:
        branches.append(tuple(c.INDEX[c.native.native(state, edge, coin)[0]]
                              for edge in c.quantum_mixing.seams() for coin in (0,1)))
    once = sum(target == i for i,row in enumerate(branches) for target in row)
    twice = sum(target == i for i,row in enumerate(branches)
                for middle in row for target in branches[middle])
    assert Q(once,66*60) == c.return_probability(1) == Q(28,33)
    assert Q(twice,66*60**2) == c.return_probability(2)


def test_distinct_paths_give_same_matrix_unit_moment():
    # Inserting a native self-return changes event length/probability but cancels
    # from the reconstructed moment after the explicit known-word coefficient.
    counts,denominator = c.source_counts()
    for (i,j),(path,coefficient) in list(c.matrix_unit_paths().items())[::113]:
        padded=(path[0],*path)
        coefficient2=Q(counts[i][i],denominator)*coefficient
        value=c.closed_path_probability(2,padded)/(coefficient2*c.return_probability(2))
        assert value == c.matrix_unit_moment(2,i,j)


def test_full_count_density_has_coherence_but_native_diagonal_does_not_determine_it():
    density=c.reconstructed_density(1)
    assert sum(density[i][i] for i in range(66)) == 1
    assert any(density[i][j] for i in range(66) for j in range(66) if i != j)
    values=c.example_effects(1)
    assert Q(values['weights']['plus']) > Q(values['weights']['minus'])
    assert values['weights']['y_plus'] == values['weights']['y_minus']
    assert not values['these_effects_implemented_as_native_single_trial_readouts']


def test_complex_linear_reconstruction_is_additive_with_exact_coordinates():
    a={(0,0):(Q(1,3),0),(0,1):(Q(1,7),Q(2,9))}
    b={(1,0):(Q(1,7),Q(-2,9)),(0,0):(Q(1,6),0)}
    merged={key:tuple(sum(part.get(key,(0,0))[i] for part in (a,b)) for i in (0,1))
            for key in a.keys()|b.keys()}
    left=c.sparse_moment(2,merged)
    right=tuple(x+y for x,y in zip(c.sparse_moment(2,a),c.sparse_moment(2,b)))
    assert left == right
    assert left[1] == 0


def test_clock_and_operational_boundaries_are_explicit():
    receipt=c.strict_json.load(c.HERE/'results.json')
    assert receipt['source']['clock'] == 'every native seam attempt, including zero current'
    assert receipt['scope']['born_form_on_reconstructed_complexification_derived']
    assert not receipt['scope']['arbitrary_complex_effect_native_implementation']
    assert not receipt['scope']['native_source_fields_same_as_incident_clock_dictionary']
    assert not receipt['scope']['rho_L_prepared_on_native_quantum_working_carrier']


def test_native_census_attachment_controls_rare_word_amplification():
    bound=c.configuration_attachment()
    assert bound['native_attempts_per_experiment_max'] == 11398
    assert Q(bound['minimum_word_coefficient']) == Q(1,60**4)
    assert Q(bound['coordinate_error_upper']) == Q(634,10**11)
    assert Q(bound['all_complex_effect_statistic_error_upper']) < Q(4,10**6)
    assert not bound['finite_sampling_error_included']
    assert not bound['raw_reconstructed_matrix_certified_positive']


@pytest.mark.parametrize('bad',[0,-1,True,1.5,17])
def test_invalid_loop_lengths_are_rejected(bad):
    with pytest.raises(ValueError,match='loop length'):
        c.path_counts(bad)


def test_missing_native_connections_cannot_supply_full_reconstruction(monkeypatch):
    original=c.source_counts
    c.matrix_unit_paths.cache_clear()
    monkeypatch.setattr(c,'source_counts',lambda:(tuple(tuple(60*int(i==j) for j in range(66))
                                                       for i in range(66)),60))
    try:
        with pytest.raises(ValueError,match='disconnected'):
            c.matrix_unit_paths()
    finally:
        monkeypatch.setattr(c,'source_counts',original)
        c.matrix_unit_paths.cache_clear()


def test_complex_coordinates_must_be_exact():
    with pytest.raises(ValueError,match='rational complex'):
        c.sparse_moment(1,{(0,0):(.5,0)})
