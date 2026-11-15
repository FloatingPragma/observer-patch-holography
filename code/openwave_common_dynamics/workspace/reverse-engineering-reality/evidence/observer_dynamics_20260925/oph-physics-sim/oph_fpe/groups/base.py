from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Hashable, Iterable

GroupElement = Hashable


@dataclass(frozen=True)
class FiniteGroup:
    """Small finite group interface used by patch ports and holonomies."""

    name: str
    elements: tuple[GroupElement, ...]
    identity: GroupElement
    multiply_fn: Callable[[GroupElement, GroupElement], GroupElement]
    inverse_fn: Callable[[GroupElement], GroupElement]
    label_fn: Callable[[GroupElement], str] = str

    def multiply(self, left: GroupElement, right: GroupElement) -> GroupElement:
        return self.multiply_fn(left, right)

    def inverse(self, element: GroupElement) -> GroupElement:
        return self.inverse_fn(element)

    def mismatch(self, left: GroupElement, right: GroupElement) -> float:
        return 0.0 if left == right else 1.0

    def label(self, element: GroupElement) -> str:
        return self.label_fn(element)

    def parse(self, value: object) -> GroupElement:
        if value in self.elements:
            return value  # type: ignore[return-value]
        for element in self.elements:
            if self.label(element) == value:
                return element
        raise ValueError(f"{value!r} is not an element of {self.name}")

    def labels(self, elements: Iterable[GroupElement]) -> list[str]:
        return [self.label(element) for element in elements]
