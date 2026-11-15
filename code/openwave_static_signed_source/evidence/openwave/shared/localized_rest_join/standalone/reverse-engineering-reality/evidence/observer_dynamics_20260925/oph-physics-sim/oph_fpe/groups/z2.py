from __future__ import annotations

from oph_fpe.groups.base import FiniteGroup


def build_z2() -> FiniteGroup:
    return FiniteGroup(
        name="Z2",
        elements=(0, 1),
        identity=0,
        multiply_fn=lambda left, right: int(left) ^ int(right),
        inverse_fn=lambda element: int(element),
        label_fn=lambda element: str(int(element)),
    )
