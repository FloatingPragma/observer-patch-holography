#!/usr/bin/env python3
"""Exact constraint-propagating global-section enumerator.

This module is intentionally independent of the capacity evaluator so the
finite CSP/model-counting backend can be tested without circular imports.
"""
from __future__ import annotations

from typing import Any, Mapping, Sequence


def public_global_sections_csp(
    observers: Mapping[str, Sequence[str]],
    interfaces: Sequence[Mapping[str, Any]],
) -> list[dict[str, str]]:
    """Enumerate all compatible atom sections by exact backtracking.

    The routine validates the same atom/readout contract as the historical
    Cartesian-product implementation.  It changes only search order: an atom is
    rejected as soon as one already-assigned interface disagrees, and the next
    variable is selected by assigned-neighbour count and degree.
    """
    if not observers or any(not atoms for atoms in observers.values()):
        raise ValueError("every declared observer needs a nonempty atom set")
    observer_ids = sorted(observers)
    if any(len(set(observers[observer])) != len(observers[observer]) for observer in observer_ids):
        raise ValueError("local record atoms must be unique")

    incident: dict[str, list[Mapping[str, Any]]] = {observer: [] for observer in observer_ids}
    for interface in interfaces:
        left = interface["left_observer"]
        right = interface["right_observer"]
        if left not in observers or right not in observers:
            raise ValueError("interface references an unknown observer")
        left_map = interface.get("left_readout")
        right_map = interface.get("right_readout")
        if not isinstance(left_map, Mapping) or not isinstance(right_map, Mapping):
            raise ValueError("RECORD-ATOM-RESTRICTION maps are required")
        if set(left_map) != set(observers[left]) or set(right_map) != set(observers[right]):
            raise ValueError("atom readout maps must be total on endpoint atoms")
        incident[left].append(interface)
        incident[right].append(interface)

    sections: list[dict[str, str]] = []
    assignment: dict[str, str] = {}

    def other_endpoint(interface: Mapping[str, Any], observer: str) -> str:
        return (
            interface["right_observer"]
            if interface["left_observer"] == observer
            else interface["left_observer"]
        )

    def consistent(observer: str, atom: str) -> bool:
        for interface in incident[observer]:
            left = interface["left_observer"]
            right = interface["right_observer"]
            if observer == left:
                if right in assignment and interface["left_readout"][atom] != interface["right_readout"][assignment[right]]:
                    return False
            elif left in assignment and interface["right_readout"][atom] != interface["left_readout"][assignment[left]]:
                return False
        return True

    def choose_next() -> str:
        return max(
            (observer for observer in observer_ids if observer not in assignment),
            key=lambda observer: (
                sum(other_endpoint(interface, observer) in assignment for interface in incident[observer]),
                len(incident[observer]),
                observer,
            ),
        )

    def search() -> None:
        if len(assignment) == len(observer_ids):
            sections.append({observer: assignment[observer] for observer in observer_ids})
            return
        observer = choose_next()
        for atom in observers[observer]:
            if consistent(observer, atom):
                assignment[observer] = atom
                search()
                del assignment[observer]

    search()
    return sections

