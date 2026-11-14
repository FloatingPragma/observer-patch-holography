from __future__ import annotations

from itertools import permutations

from oph_fpe.groups.base import FiniteGroup, GroupElement

Permutation = tuple[int, int, int]


def _compose(left: GroupElement, right: GroupElement) -> GroupElement:
    """Return left after right for permutations of 0,1,2."""

    l_perm = left  # type: ignore[assignment]
    r_perm = right  # type: ignore[assignment]
    return tuple(l_perm[r_perm[index]] for index in range(3))


def _inverse(element: GroupElement) -> GroupElement:
    perm = element  # type: ignore[assignment]
    inv = [0, 0, 0]
    for index, image in enumerate(perm):
        inv[image] = index
    return tuple(inv)


def _label(element: GroupElement) -> str:
    labels = {
        (0, 1, 2): "e",
        (1, 0, 2): "(01)",
        (2, 1, 0): "(02)",
        (0, 2, 1): "(12)",
        (1, 2, 0): "(012)",
        (2, 0, 1): "(021)",
    }
    return labels[element]  # type: ignore[index]


def build_s3() -> FiniteGroup:
    elements = tuple(permutations((0, 1, 2)))
    return FiniteGroup(
        name="S3",
        elements=elements,
        identity=(0, 1, 2),
        multiply_fn=_compose,
        inverse_fn=_inverse,
        label_fn=_label,
    )
