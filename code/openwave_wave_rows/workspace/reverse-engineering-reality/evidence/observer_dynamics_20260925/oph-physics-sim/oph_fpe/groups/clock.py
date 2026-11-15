from __future__ import annotations

from oph_fpe.groups.base import FiniteGroup


def build_clock(order: int) -> FiniteGroup:
    if order < 2:
        raise ValueError("clock group order must be at least 2")
    return FiniteGroup(
        name=f"C{order}",
        elements=tuple(range(order)),
        identity=0,
        multiply_fn=lambda left, right: (int(left) + int(right)) % order,
        inverse_fn=lambda element: (-int(element)) % order,
        label_fn=lambda element: str(int(element)),
    )
