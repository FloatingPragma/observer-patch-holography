"""Exact comparison bounds for the released 44-dimer wall band."""

from fractions import Fraction
from itertools import combinations

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model.carrier import federation

from ..row_20_electric.experiment import guide


def wall_sites(loads):
    chain = (1, *loads, 1)
    return {j for j in range(len(chain) - 1) if chain[j] == chain[j + 1]}


def finite_crosscheck(private):
    """Check the counting inequalities on every word of a small prefix."""
    width = len(private)
    half = width // 2
    baseline = sum(private) // 2
    unique = set()
    for occupied in combinations(range(width), half):
        word = tuple(int(j in occupied) for j in range(width))
        sites = wall_sites(word)
        walls = len(sites)
        assert walls % 2
        if walls == 1:
            unique |= sites
        extra_private_doubles = sum(private[j] for j in occupied) - baseline
        assert 2 * extra_private_doubles <= walls - 1
        less = equal = more = 0
        for j in range(width - 1):
            if word[j] == word[j + 1]:
                continue
            swapped = (*word[:j], word[j + 1], word[j], *word[j + 2 :])
            updated = wall_sites(swapped)
            assert sites ^ updated == {j, j + 2}
            difference = len(updated) - walls
            less += difference == -2
            equal += difference == 0
            more += difference == 2
        assert equal <= 2 * walls and less <= walls and more <= width - 1
        if walls == 3:
            assert less <= 2
    assert unique == set(range(0, width + 1, 2))


def verify():
    board = federation(2)
    _, guide_ports, _ = guide(board)
    internal = guide_ports[1:89]
    assert len(internal) == len(set(internal)) == 88
    assert all(
        v in board.neighbours[u]
        for u, v in zip(guide_ports[:89], guide_ports[1:90], strict=True)
    )
    private = tuple(6 - int(board.degree[p]) for p in internal)
    assert private == (0, 0, 1, 1) * 22
    finite_crosscheck(private[:12])
    native = law.seam_term(8)
    unlike = law.SRC.index((0, 2))
    reversed_pair = law.SRC.index((2, 0))
    cost = Fraction(native[unlike, unlike])
    exchange = Fraction(-native[unlike, reversed_pair])
    assert (cost, exchange) == (Fraction(65, 2), Fraction(1, 2))
    # Each positive degree-five surplus requires both one full and one
    # empty dimer. For odd W the potential loses at least (C-16)(W-1).
    alpha = Fraction(1, 5)
    if_w3 = exchange * (6 + alpha * (len(internal) - 1))
    q_gap = 2 * (cost - 16) - if_w3
    for walls in range(5, len(internal) + 2, 2):
        radius = exchange * (2 * walls + alpha * (len(internal) - 1) + walls / alpha)
        assert (cost - 16) * (walls - 1) - radius >= q_gap
    # In degree-3 wall words the graph of pairs distance two has at most
    # two edges. At most two exchange moves can enter the one-wall space.
    coupling_squared = 2 * exchange**2 * (len(internal) - 1)
    p_min_gap = q_gap - 2 * exchange
    top_above_d0 = 2 * exchange + coupling_squared / p_min_gap
    mode_one_wall_weight = 1 / (1 + coupling_squared / p_min_gap**2)
    sharp_band_weight = p_min_gap / (q_gap + top_above_d0)
    minimum_removal_frequency = cost - top_above_d0
    assert q_gap == Fraction(213, 10) and top_above_d0 == Fraction(22, 7)
    assert mode_one_wall_weight > Fraction(9, 10)
    assert sharp_band_weight > Fraction(83, 100)
    assert minimum_removal_frequency == Fraction(411, 14)
    return (
        len(internal) // 2 + 1,
        q_gap,
        minimum_removal_frequency,
        mode_one_wall_weight,
        sharp_band_weight,
    )


if __name__ == "__main__":
    count, gap, minimum, wall, sharp = verify()
    print(f"{count} isolated band modes counting multiplicity. Q ceiling = D0 - {gap}.")
    print(f"Matched-host phase frequency is at least {minimum}.")
    print(
        f"One-wall probability per mode >= {wall}. High-band probability per sharp wall >= {sharp}."
    )
