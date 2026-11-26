"""Regression checks on the actual charged source and neutral release algebra."""

from model.candidate_laws.interval_echo.openwave.row_10_oscillations.audit import (
    environment,
    release_sectors,
    sector_from_leaf_read,
    source_reflections,
)


def test_all_sharp_tip_words_have_reflections_on_the_native_source():
    source = source_reflections()
    assert source["charged_chart_dimension"] == 3584
    assert source["exact_commuting_arm_permutations"] == 6
    assert set(source["sharp_tip_word_stabilizers"]) == set(map(str, range(8)))


def test_remote_ambiguity_uses_both_leaf_loads_and_the_retained_directions():
    leaves = (165, 158)
    launch = (160, 167, 1)
    moved = (165, 167, 0)
    assert environment(launch, leaves)[1] != environment(moved, leaves)[1]
    assert (
        sector_from_leaf_read(launch, leaves)
        == sector_from_leaf_read(moved, leaves)
        == 1
    )
    opposite_launch = (160, 167, 0)
    assert environment(launch, leaves)[0] == environment(opposite_launch, leaves)[0]
    assert sector_from_leaf_read(opposite_launch, leaves) != sector_from_leaf_read(
        launch, leaves
    )


def test_invariant_on_every_native_word_and_every_reachable_word_for_both_signs():
    release = release_sectors()
    assert (release["full_release_words"], release["directed_nonzero_exchanges"]) == (
        960,
        3472,
    )
    assert release["off_diagonal_contact_entries"] == [6, 4]
    for classes in release["launch_sectors"].values():
        assert classes == dict(
            reachable_words_by_mask=[210, 211, 211, 210], environment_keys=14
        )
