#!/usr/bin/env python3
"""Pro5 exact correctable public-record capacity evaluator.

This module evaluates a paper-defined PUBLIC_CHECKPOINT_PACKET.  It uses
record-atom global sections, endogenous reachability, a frozen publicness
policy, and source-supplied joint checkpoint kernels.  Exact capacity is the
independence number of the compound confusability graph.  It never derives a
joint kernel from local marginals and never consumes Lambda, rho_op, an
electroweak target, or the supplied capacity label as a desired output.

The evaluator is finite and proof-facing.  A physical N closure remains open
until a source-derived packet and an exact finite-size slack law are supplied.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations, product
from typing import Any, Iterable, Mapping, Sequence

from public_record_csp import public_global_sections_csp


FAIL = {
    "NO_CAPACITY_READBACK",
    "AMBIGUOUS_CAPACITY_READBACK",
    "INCOMPLETE_TERMINAL_FIBER",
    "NO_RECORD_ATOM_RESTRICTION",
    "NO_PUBLIC_RECORD_REACHABILITY",
    "NO_PUBLICNESS_POLICY",
    "NO_GLOBAL_PUBLIC_CHECKPOINT_COUPLING",
    "LOCAL_MARGINAL_MISMATCH",
    "NO_CAPACITY_CARRIER_REPRESENTATION",
    "CAPACITY_EXTENSION_CONFUSABILITY_FAILED",
    "FINITE_SIZE_SELECTOR_FAILED",
    "CIRCULAR_CAPACITY_DEFINITION",
    "TARGET_TAINTED",
}


def _fail(status: str, **details: Any) -> dict[str, Any]:
    if status not in FAIL:
        raise ValueError(f"unknown status {status}")
    return {"status": status, **details}


def section_id(section: Mapping[str, str]) -> str:
    return "|".join(f"{observer}={section[observer]}" for observer in sorted(section))


def public_global_sections(
    observers: Mapping[str, Sequence[str]],
    interfaces: Sequence[Mapping[str, Any]],
) -> list[dict[str, str]]:
    """Enumerate the finite atom-level global sections."""
    if not observers or any(not atoms for atoms in observers.values()):
        raise ValueError("every declared observer needs a nonempty atom set")
    observer_ids = sorted(observers)
    if any(len(set(observers[o])) != len(observers[o]) for o in observer_ids):
        raise ValueError("local record atoms must be unique")

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

    sections: list[dict[str, str]] = []
    for values in product(*(observers[o] for o in observer_ids)):
        candidate = dict(zip(observer_ids, values, strict=True))
        compatible = all(
            interface["left_readout"][candidate[interface["left_observer"]]]
            == interface["right_readout"][candidate[interface["right_observer"]]]
            for interface in interfaces
        )
        if compatible:
            sections.append(candidate)
    return sections


# Keep the historical Cartesian implementation as a regression oracle while
# evaluating source-derived packets with exact constraint propagation.
_public_global_sections_cartesian = public_global_sections
public_global_sections = public_global_sections_csp


def reachable_public_sections(
    public_sections: Sequence[Mapping[str, str]],
    reachability_witnesses: Mapping[str, Sequence[str]],
) -> list[str]:
    """Select sections with nonempty endogenous semantic-history witnesses."""
    valid = {section_id(section) for section in public_sections}
    unknown = set(reachability_witnesses) - valid
    if unknown:
        raise ValueError(f"reachability witness references unknown sections: {sorted(unknown)}")
    return sorted(
        sid for sid, history in reachability_witnesses.items()
        if isinstance(history, Sequence) and not isinstance(history, (str, bytes)) and len(history) > 0
    )


def _channel_rows(
    channel: Mapping[str, Any], reachable: Sequence[str]
) -> dict[str, dict[str, float]]:
    rows = channel.get("rows")
    if not isinstance(rows, Mapping) or set(rows) != set(reachable):
        raise ValueError("every joint checkpoint kernel needs one row per reachable record")
    normalized: dict[str, dict[str, float]] = {}
    for source in reachable:
        row = rows[source]
        if not isinstance(row, Mapping) or not row:
            raise ValueError("checkpoint rows must be nonempty mappings")
        values = {str(out): float(prob) for out, prob in row.items()}
        if any(prob < 0 for prob in values.values()) or abs(sum(values.values()) - 1.0) > 1e-12:
            raise ValueError("checkpoint rows must be normalized probabilities")
        normalized[source] = values
    return normalized


def compound_confusability_graph(
    reachable: Sequence[str], channels: Sequence[Mapping[str, Any]]
) -> dict[str, set[str]]:
    """Return the union of all declared channel confusability graphs."""
    if not channels:
        raise ValueError("GLOBAL-PUBLIC-CHECKPOINT-COUPLING is required")
    graph = {source: set() for source in reachable}
    for channel in channels:
        rows = _channel_rows(channel, reachable)
        support = {source: {out for out, p in row.items() if p > 0} for source, row in rows.items()}
        for i, left in enumerate(reachable):
            for right in reachable[i + 1 :]:
                if support[left] & support[right]:
                    graph[left].add(right)
                    graph[right].add(left)
    return graph


def _decoder_success(
    rows: Mapping[str, Mapping[str, float]], code: Sequence[str]
) -> float:
    """Return the optimal worst-input success over deterministic decoders.

    This exhaustive routine is deliberately receipt-scale.  It implements the
    paper's worst-input definition exactly and is not intended for large
    production alphabets.
    """
    outputs = sorted({out for source in code for out in rows[source]})
    best = 0.0
    # ``None`` leaves an output unassigned.  It is redundant mathematically,
    # but makes partial decoder witnesses explicit.
    targets: tuple[str | None, ...] = (None, *code)
    for assignment in product(targets, repeat=len(outputs)):
        decoded = dict(zip(outputs, assignment, strict=True))
        success = min(
            sum(probability for out, probability in rows[source].items() if decoded[out] == source)
            for source in code
        )
        best = max(best, success)
    return best


def approximate_public_capacity(
    reachable: Sequence[str],
    channels: Sequence[Mapping[str, Any]],
    epsilon: float,
    *,
    max_vertices: int = 12,
) -> dict[str, Any]:
    """Compute finite compound worst-input ``M_epsilon`` by exhaustive search."""
    if not 0 <= epsilon <= 1:
        raise ValueError("epsilon must lie in [0,1]")
    if len(reachable) > max_vertices:
        raise ValueError("approximate evaluator is limited to receipt-scale alphabets")
    normalized = [_channel_rows(channel, reachable) for channel in channels]
    for size in range(len(reachable), 0, -1):
        for code in combinations(sorted(reachable), size):
            successes = [_decoder_success(rows, code) for rows in normalized]
            if min(successes) + 1e-12 >= 1.0 - epsilon:
                return {
                    "capacity": size,
                    "code_witness": list(code),
                    "worst_input_success_by_channel": successes,
                }
    return {"capacity": 0, "code_witness": [], "worst_input_success_by_channel": []}


def maximum_independent_set(graph: Mapping[str, set[str]]) -> list[str]:
    """Exact branch-and-bound maximum independent set for finite receipt tests."""
    vertices = tuple(sorted(graph))
    best: tuple[str, ...] = ()

    def search(candidates: tuple[str, ...], chosen: tuple[str, ...]) -> None:
        nonlocal best
        if len(chosen) + len(candidates) <= len(best):
            return
        if not candidates:
            if len(chosen) > len(best):
                best = chosen
            return
        vertex = max(candidates, key=lambda item: len(graph[item].intersection(candidates)))
        without_vertex = tuple(x for x in candidates if x != vertex)
        compatible = tuple(x for x in without_vertex if x not in graph[vertex])
        search(compatible, chosen + (vertex,))
        search(without_vertex, chosen)

    search(vertices, ())
    return list(best)


def support_relation_semigroup(
    reachable: Sequence[str], generators: Sequence[Mapping[str, Any]]
) -> set[frozenset[tuple[str, str]]]:
    """Close same-alphabet checkpoint support relations under composition."""
    universe = set(reachable)
    relations: set[frozenset[tuple[str, str]]] = set()
    for channel in generators:
        rows = _channel_rows(channel, reachable)
        relation = frozenset(
            (source, target)
            for source, row in rows.items()
            for target, probability in row.items()
            if probability > 0 and target in universe
        )
        if any(not any(left == source for left, _ in relation) for source in reachable):
            raise ValueError("indefinite support generators must continue on the record alphabet")
        relations.add(relation)

    changed = True
    while changed:
        changed = False
        known = list(relations)
        for left in known:
            for right in known:
                composed = frozenset(
                    (source, target)
                    for source, middle in left
                    for middle2, target in right
                    if middle == middle2
                )
                if composed not in relations:
                    relations.add(composed)
                    changed = True
    return relations


def certify_capacity_carrier(
    capacity_dimension: int,
    reachable: Sequence[str],
    projection_supports: Mapping[str, Sequence[int]],
) -> dict[str, Any]:
    if not isinstance(capacity_dimension, int) or capacity_dimension < 1:
        raise ValueError("capacity dimension must be positive")
    if set(projection_supports) != set(reachable):
        return _fail("NO_CAPACITY_CARRIER_REPRESENTATION")
    used: set[int] = set()
    ranks: dict[str, int] = {}
    for record in reachable:
        support = list(projection_supports[record])
        if not support or len(set(support)) != len(support):
            return _fail("NO_CAPACITY_CARRIER_REPRESENTATION", record=record)
        if any(not isinstance(i, int) or i < 0 or i >= capacity_dimension for i in support):
            return _fail("NO_CAPACITY_CARRIER_REPRESENTATION", record=record)
        if used.intersection(support):
            return _fail("NO_CAPACITY_CARRIER_REPRESENTATION", record=record)
        used.update(support)
        ranks[record] = len(support)
    return {"status": "PASS", "projection_ranks": ranks, "rank_sum": len(used)}


def evaluate_terminal(packet: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate one exact PUBLIC_CHECKPOINT_PACKET."""
    if packet.get("self_read_predicate_injected", False):
        return _fail("CIRCULAR_CAPACITY_DEFINITION")
    if packet.get("lambda_used", False) or packet.get("ew_bridge_used", False) or packet.get("rho_used", False):
        return _fail("TARGET_TAINTED")
    try:
        sections = public_global_sections(packet["observers"], packet["interfaces"])
    except (KeyError, ValueError) as exc:
        return _fail("NO_RECORD_ATOM_RESTRICTION", reason=str(exc))
    witnesses = packet.get("reachability_witnesses")
    if not isinstance(witnesses, Mapping):
        return _fail("NO_PUBLIC_RECORD_REACHABILITY")
    reachable = reachable_public_sections(sections, witnesses)
    if not reachable:
        return _fail("NO_PUBLIC_RECORD_REACHABILITY")
    policy = packet.get("publicness_policy")
    if not isinstance(policy, Sequence) or isinstance(policy, (str, bytes)) or not policy:
        return _fail("NO_PUBLICNESS_POLICY")
    observer_ids = set(packet["observers"])
    normalized_policy: list[tuple[str, ...]] = []
    for authorized in policy:
        if (
            not isinstance(authorized, Sequence)
            or isinstance(authorized, (str, bytes))
            or not authorized
            or not set(authorized).issubset(observer_ids)
        ):
            return _fail("NO_PUBLICNESS_POLICY")
        normalized_policy.append(tuple(sorted(set(authorized))))
    channels = packet.get("global_checkpoint_kernels")
    if not isinstance(channels, Sequence) or isinstance(channels, (str, bytes)) or not channels:
        return _fail("NO_GLOBAL_PUBLIC_CHECKPOINT_COUPLING")
    if not packet.get("continuation_manifest_complete", False):
        return _fail("NO_GLOBAL_PUBLIC_CHECKPOINT_COUPLING", reason="incomplete continuation manifest")
    for channel in channels:
        authorized = channel.get("authorized_observers")
        if not isinstance(authorized, Sequence) or isinstance(authorized, (str, bytes)):
            return _fail("NO_GLOBAL_PUBLIC_CHECKPOINT_COUPLING")
        if tuple(sorted(set(authorized))) not in normalized_policy or not channel.get("continuation_id"):
            return _fail("NO_GLOBAL_PUBLIC_CHECKPOINT_COUPLING")
    if not packet.get("local_marginal_consistency_passed", False):
        return _fail("LOCAL_MARGINAL_MISMATCH")
    try:
        graph = compound_confusability_graph(reachable, channels)
    except ValueError as exc:
        return _fail("NO_GLOBAL_PUBLIC_CHECKPOINT_COUPLING", reason=str(exc))
    code = maximum_independent_set(graph)
    carrier = certify_capacity_carrier(
        packet["capacity_dimension"], reachable, packet.get("projection_supports", {})
    )
    if carrier["status"] != "PASS":
        return carrier
    capacity = len(code)
    dimension = packet["capacity_dimension"]
    return {
        "status": "PASS",
        "terminal_id": packet.get("terminal_id", "UNNAMED"),
        "public_global_section_count": len(sections),
        "reachable_public_sections": reachable,
        "publicness_policy": [list(authorized) for authorized in normalized_policy],
        "confusability_graph": {k: sorted(v) for k, v in graph.items()},
        "independent_set_witness": code,
        "independence_number": capacity,
        "exact_zero_error_capacity": capacity,
        "capacity_dimension": dimension,
        "dimension_bound_passed": capacity <= dimension,
        "saturation_passed": capacity == dimension,
        "saturation_rank_one_complete": (
            capacity == dimension
            and len(reachable) == dimension
            and all(rank == 1 for rank in carrier["projection_ranks"].values())
        ),
        "capacity_slack": f"log({dimension}/{capacity})",
    }


def evaluate_terminal_fiber(
    packets: Sequence[Mapping[str, Any]], *, manifest_complete: bool
) -> dict[str, Any]:
    if not manifest_complete:
        return _fail("INCOMPLETE_TERMINAL_FIBER")
    if not packets:
        return _fail("NO_CAPACITY_READBACK")
    results = [evaluate_terminal(packet) for packet in packets]
    if any(result["status"] != "PASS" for result in results):
        return {"status": "NO_CAPACITY_READBACK", "terminal_results": results}
    capacities = [result["exact_zero_error_capacity"] for result in results]
    support = sorted(set(capacities))
    base = {
        "terminal_results": results,
        "terminal_fiber_capacity_set": support,
        "unclosed_readback_kernel": dict(Counter(capacities)),
    }
    if len(support) != 1:
        return {"status": "AMBIGUOUS_CAPACITY_READBACK", **base}
    dimension = results[0]["capacity_dimension"]
    if any(result["capacity_dimension"] != dimension for result in results):
        raise ValueError("one terminal fiber must have one supplied dimension")
    return {
        "status": "PASS",
        **base,
        "scalar_readback_dimension": support[0],
        "robust_closure": support == [dimension],
    }


def greatest_fixed_point(capacity_map: Mapping[int, int]) -> dict[str, Any]:
    domain = sorted(capacity_map)
    if not domain or domain != list(range(1, domain[-1] + 1)):
        raise ValueError("domain must be the positive chain 1..D_max")
    if any(value not in capacity_map for value in capacity_map.values()):
        raise ValueError("map must be total and closed on the declared chain")
    if any(capacity_map[d] > d for d in domain):
        raise ValueError("active capacity map must be deflationary")
    if any(capacity_map[a] > capacity_map[b] for a, b in zip(domain, domain[1:], strict=False)):
        raise ValueError("capacity extension must be monotone")
    path = [domain[-1]]
    while capacity_map[path[-1]] != path[-1]:
        path.append(capacity_map[path[-1]])
    fixed = [d for d in domain if capacity_map[d] == d]
    return {"path": path, "fixed_points": fixed, "greatest_fixed_point": path[-1]}


def no_new_confusability(
    coarse_graph: Mapping[str, set[str]],
    fine_graph: Mapping[str, set[str]],
    embedding: Mapping[str, str],
) -> bool:
    if set(embedding) != set(coarse_graph) or len(set(embedding.values())) != len(embedding):
        return False
    if any(image not in fine_graph for image in embedding.values()):
        return False
    for left in coarse_graph:
        for right in coarse_graph:
            if left == right:
                continue
            if embedding[right] in fine_graph[embedding[left]] and right not in coarse_graph[left]:
                return False
    return True


def certify_unique_slack_zero(capacity_map: Mapping[int, int], selected: int) -> dict[str, Any]:
    if selected not in capacity_map:
        raise ValueError("selected dimension is outside the certified domain")
    zeros = sorted(d for d, value in capacity_map.items() if value == d)
    if zeros != [selected]:
        return _fail("FINITE_SIZE_SELECTOR_FAILED", fixed_points=zeros)
    return {"status": "PASS", "selected_dimension": selected, "fixed_points": zeros}


def tv_robustness_bound(epsilon: float, delta: float) -> float:
    if not 0 <= epsilon <= 1 or not 0 <= delta <= 1:
        raise ValueError("epsilon and delta must lie in [0,1]")
    return min(1.0, epsilon + delta)
