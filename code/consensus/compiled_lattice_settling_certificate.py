#!/usr/bin/env python3
"""Uniform-settling realization certificate for compiled Boolean lattices (issue #328).

The abstract acyclic Boolean circuit compiler corollary
(``extra/observable_normal_forms.tex``, ``cor:boolean-circuit-compiler``,
resting on ``cor:functional-synchronous``; Lean:
``Lean/ObservableNormalForms/ObservableNormalForms/Functional.lean``,
``RankedSynchronousSystem.synchronous_depth_settling``) proves depth settling
for the explicit ranked functional system and states that correspondence with
a separately specified update rule requires its own theorem.  This packet is
that correspondence at finite scope.  It supplies:

* explicit patch encodings of a universal gate set (NAND, wire, fan-out-two),
  each an admissible patch reduct with declared ports, local state registers,
  readback, records, and a repair/update table;
* an intertwiner between patch dynamics and the ranked functional system,
  checked exhaustively per primitive over all port values and all internal
  register states (a finite bisimulation, exact);
* a composition checker over compiled circuit DAGs (acyclicity, port arity,
  single-consumer wiring, rank ladder) realizing the structural induction;
* a size-uniform settling theorem for the realized dynamics with the explicit
  constant from the primitive verification, exact settling times for three
  test circuits, and a strictly decreasing settling potential, including a
  recorded countermodel against the naive remaining-depth weighting;
* fan-out non-interference evidence and fail-closed controls: a noisy update
  table, a cross-talking fan-out, a cyclic layout, a multiply consumed port,
  and a fan-in overload.

All arithmetic is exact integer and Boolean arithmetic.  The abstract
compiler result stays separately labelled in the manifest; continuum and
physical-hardware attachment is open and named there.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
ROOT = MODULE_DIR.parents[1]

SCHEMA = "oph.compiled_lattice_settling_certificate.v1"
DEFAULT_MANIFEST = MODULE_DIR / "manifests" / "compiled_lattice_settling_reference.json"

ABSTRACT_COMPILER_RESULT = {
    "label": (
        "separately labelled abstract result; this packet cites it and does not "
        "restate or replace its proof"
    ),
    "paper": "extra/observable_normal_forms.tex",
    "paper_results": [
        "cor:boolean-circuit-compiler (Acyclic Boolean circuit compiler)",
        "cor:functional-synchronous (Synchronous settling by dependency depth)",
    ],
    "lean_module": "Lean/ObservableNormalForms/ObservableNormalForms/Functional.lean",
    "lean_declarations": [
        "ObservableNormalForms.RankedSynchronousSystem.synchronousEvolve_agrees_through_rank",
        "ObservableNormalForms.RankedSynchronousSystem.synchronous_depth_settling",
        "ObservableNormalForms.RankedSynchronousSystem.generatedExtension_unique",
    ],
    "proof_index": "Lean/ObservableNormalForms/PROOF_INDEX.md",
    "gap_this_packet_closes_at_finite_scope": (
        "cor:functional-synchronous states that correspondence with a separately "
        "specified compiler or update rule requires its own theorem; this packet "
        "supplies the realized-dynamics correspondence for the finite software "
        "patch dynamics defined here"
    ),
}


class CertificateError(ValueError):
    """Fail-closed certificate error carrying a stable code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode(
        "utf-8"
    )


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Manifest fixtures are compared byte-for-byte across platforms.
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CertificateError("JSON_READ", f"cannot read {path}: {exc}") from exc


def bit_tuples(width: int) -> Iterator[tuple[int, ...]]:
    return product((0, 1), repeat=width)


# ---------------------------------------------------------------------------
# Primitive patches
# ---------------------------------------------------------------------------
#
# A primitive patch is the finite reduct of an observer patch specialized to
# one compiled site: the port list is its access interface, the registers are
# its local state, the readback map exposes registers on output ports, the
# record is the settled register content, and the update table is its repair
# instrument.  The update table is stored over the full domain (port values
# times internal register states) so that state dependence, the cross-talk
# failure mode, is expressible and its absence is a checked fact rather than a
# representation artifact.


@dataclass(frozen=True)
class Primitive:
    name: str
    in_ports: tuple[str, ...]
    registers: tuple[str, ...]
    out_ports: tuple[str, ...]
    readback: tuple[int, ...]  # out port index -> register index
    truth: Mapping[tuple[int, ...], tuple[int, ...]]  # port values -> out port values
    update: Mapping[
        tuple[tuple[int, ...], tuple[int, ...]], tuple[int, ...]
    ]  # (port values, register state) -> next register state


def build_primitive(
    name: str,
    in_ports: Sequence[str],
    registers: Sequence[str],
    out_ports: Sequence[str],
    readback: Sequence[int],
    truth_fn: Callable[[tuple[int, ...]], tuple[int, ...]],
    next_state_fn: Callable[[tuple[int, ...], tuple[int, ...]], tuple[int, ...]],
) -> Primitive:
    truth = {inputs: truth_fn(inputs) for inputs in bit_tuples(len(in_ports))}
    update = {
        (inputs, state): next_state_fn(inputs, state)
        for inputs in bit_tuples(len(in_ports))
        for state in bit_tuples(len(registers))
    }
    return Primitive(
        name=name,
        in_ports=tuple(in_ports),
        registers=tuple(registers),
        out_ports=tuple(out_ports),
        readback=tuple(readback),
        truth=truth,
        update=update,
    )


def reference_primitives() -> dict[str, Primitive]:
    nand = build_primitive(
        "NAND",
        in_ports=("x", "y"),
        registers=("q",),
        out_ports=("z",),
        readback=(0,),
        truth_fn=lambda i: (1 - (i[0] & i[1]),),
        next_state_fn=lambda i, s: (1 - (i[0] & i[1]),),
    )
    wire = build_primitive(
        "WIRE",
        in_ports=("x",),
        registers=("q",),
        out_ports=("z",),
        readback=(0,),
        truth_fn=lambda i: (i[0],),
        next_state_fn=lambda i, s: (i[0],),
    )
    fanout2 = build_primitive(
        "FANOUT2",
        in_ports=("x",),
        registers=("q0", "q1"),
        out_ports=("z0", "z1"),
        readback=(0, 1),
        truth_fn=lambda i: (i[0], i[0]),
        next_state_fn=lambda i, s: (i[0], i[0]),
    )
    return {p.name: p for p in (nand, wire, fanout2)}


def noisy_nand_control() -> Primitive:
    """NAND with one flipped update-table entry: constant one on port image."""

    return build_primitive(
        "NOISY_NAND",
        in_ports=("x", "y"),
        registers=("q",),
        out_ports=("z",),
        readback=(0,),
        truth_fn=lambda i: (1 - (i[0] & i[1]),),
        next_state_fn=lambda i, s: (1,) if i == (1, 1) else (1 - (i[0] & i[1]),),
    )


def crosstalk_fanout_control() -> Primitive:
    """Fan-out whose first branch reads the sibling register: cross-talk."""

    return build_primitive(
        "CROSSTALK_FANOUT",
        in_ports=("x",),
        registers=("q0", "q1"),
        out_ports=("z0", "z1"),
        readback=(0, 1),
        truth_fn=lambda i: (i[0], i[0]),
        next_state_fn=lambda i, s: (i[0] ^ s[1], i[0]),
    )


def verify_primitive(prim: Primitive) -> dict[str, Any]:
    """Exhaustive per-primitive intertwiner check.

    Over every port assignment and every internal register state the check
    verifies: (a) the next state depends on the port values alone, so the
    update is local to the declared access interface and sibling branches do
    not interact; (b) one update settles the readback to the declared truth
    function; (c) the settled state is a fixed point.  Together these give the
    exact one-step commutation ``readback(update(ports, state)) = truth(ports)``
    for every state, which is the local intertwiner between the patch dynamics
    and the single-site ranked functional system carrying ``truth`` as its
    generating function.
    """

    checked_pairs = 0
    for inputs in bit_tuples(len(prim.in_ports)):
        next_states = {prim.update[(inputs, state)] for state in bit_tuples(len(prim.registers))}
        if len(next_states) != 1:
            raise CertificateError(
                "CROSS_TALK",
                f"primitive {prim.name}: next state depends on internal register "
                f"state at ports {inputs}; state-dependent images {sorted(next_states)}",
            )
        settled = next(iter(next_states))
        readout = tuple(settled[r] for r in prim.readback)
        if readout != prim.truth[inputs]:
            raise CertificateError(
                "INTERTWINER_BROKEN",
                f"primitive {prim.name}: settled readback {readout} at ports "
                f"{inputs} differs from declared truth {prim.truth[inputs]}",
            )
        if prim.update[(inputs, settled)] != settled:
            raise CertificateError(
                "INTERTWINER_BROKEN",
                f"primitive {prim.name}: settled state {settled} is not a fixed "
                f"point at ports {inputs}",
            )
        checked_pairs += 2 ** len(prim.registers)
    return {
        "ports": {"in": list(prim.in_ports), "out": list(prim.out_ports)},
        "local_state": list(prim.registers),
        "readback": {
            prim.out_ports[k]: prim.registers[prim.readback[k]]
            for k in range(len(prim.out_ports))
        },
        "record": (
            "the settled register content; at the fixed point the register "
            "carries the declared truth of the port values and one further "
            "update leaves it unchanged"
        ),
        "repair_update_rule": "synchronous table update: registers := table(port values)",
        "truth_table": {
            "".join(str(b) for b in inputs): list(outputs)
            for inputs, outputs in sorted(prim.truth.items())
        },
        "evidence_bundle": {
            "checked_port_state_pairs": checked_pairs,
            "state_independent_update": True,
            "settles_in_rounds": 1,
            "settled_state_is_fixed_point": True,
            "one_step_commutation": (
                "readback(update(ports, state)) = truth(ports) for every "
                "internal state; encode/decode are identity on register bits"
            ),
        },
    }


# ---------------------------------------------------------------------------
# Patch netlists and the composition checker
# ---------------------------------------------------------------------------

SourceRef = tuple  # ("in", input_name) or ("out", instance_name, out_port_name)


@dataclass(frozen=True)
class Instance:
    name: str
    kind: str
    sources: tuple[SourceRef, ...]


@dataclass(frozen=True)
class Netlist:
    name: str
    inputs: tuple[str, ...]
    instances: tuple[Instance, ...]
    outputs: tuple[tuple[str, str], ...]  # (output label, wire instance name)


def check_netlist(net: Netlist, prims: Mapping[str, Primitive]) -> dict[str, Any]:
    """Composition checker: the structural-induction side of the packet.

    Rejects unknown primitives, port-arity mismatches (fan-in overload),
    dangling references, multiply consumed ports (fan-out without a fan-out
    patch), and cyclic wiring.  On acceptance it returns the rank ladder, the
    compiled depth, and the settling-potential weights, with the per-patch
    weight inequality w(p) >= 1 + sum of consumer weights checked exactly.
    """

    by_name: dict[str, Instance] = {}
    for inst in net.instances:
        if inst.name in by_name or inst.name in net.inputs:
            raise CertificateError("DUPLICATE_NAME", f"{net.name}: duplicate name {inst.name}")
        by_name[inst.name] = inst

    consumed: dict[SourceRef, int] = {}
    parents: dict[str, list[str]] = {inst.name: [] for inst in net.instances}
    for inst in net.instances:
        prim = prims.get(inst.kind)
        if prim is None:
            raise CertificateError(
                "UNKNOWN_PRIMITIVE", f"{net.name}: instance {inst.name} has kind {inst.kind}"
            )
        if len(inst.sources) != len(prim.in_ports):
            raise CertificateError(
                "PORT_ARITY",
                f"{net.name}: instance {inst.name} wires {len(inst.sources)} sources "
                f"into {len(prim.in_ports)} ports of {inst.kind}",
            )
        for ref in inst.sources:
            if ref[0] == "in":
                if ref[1] not in net.inputs:
                    raise CertificateError(
                        "UNKNOWN_SIGNAL", f"{net.name}: {inst.name} reads missing input {ref[1]}"
                    )
            elif ref[0] == "out":
                src = by_name.get(ref[1])
                if src is None:
                    raise CertificateError(
                        "UNKNOWN_SIGNAL",
                        f"{net.name}: {inst.name} reads missing instance {ref[1]}",
                    )
                if ref[2] not in prims[src.kind].out_ports:
                    raise CertificateError(
                        "UNKNOWN_SIGNAL",
                        f"{net.name}: {inst.name} reads missing port {ref[2]} of {ref[1]}",
                    )
                parents[inst.name].append(ref[1])
            else:
                raise CertificateError(
                    "UNKNOWN_SIGNAL", f"{net.name}: {inst.name} has malformed ref {ref!r}"
                )
            consumed[ref] = consumed.get(ref, 0) + 1
            if consumed[ref] > 1:
                raise CertificateError(
                    "MULTI_CONSUMER",
                    f"{net.name}: source {ref!r} feeds more than one port; "
                    "replication requires an explicit fan-out patch",
                )

    for label, wire_name in net.outputs:
        if wire_name not in by_name:
            raise CertificateError(
                "UNKNOWN_SIGNAL", f"{net.name}: output {label} names missing instance {wire_name}"
            )

    # Kahn topological sort; leftover instances witness a cycle.
    remaining = {name: len(set(p for p in parents[name])) for name in parents}
    unique_parents = {name: sorted(set(parents[name])) for name in parents}
    children: dict[str, list[str]] = {name: [] for name in parents}
    for name, ps in unique_parents.items():
        for p in ps:
            children[p].append(name)
    order = sorted(name for name, n in remaining.items() if n == 0)
    queue = list(order)
    while queue:
        current = queue.pop(0)
        for child in sorted(children[current]):
            remaining[child] -= 1
            if remaining[child] == 0:
                order.append(child)
                queue.append(child)
    if len(order) != len(net.instances):
        cyclic = sorted(name for name, n in remaining.items() if n > 0)
        raise CertificateError(
            "CYCLIC_LAYOUT",
            f"{net.name}: acyclicity fails; instances on a cycle: {cyclic}",
        )

    rank: dict[str, int] = {}
    for name in order:
        ps = unique_parents[name]
        rank[name] = 1 + max((rank[p] for p in ps), default=0)
    depth = max(rank.values()) if rank else 0

    # Settling-potential weights with fan-out multiplicity: one term per
    # consuming port use.  The exact identity w(p) = 1 + sum over port uses of
    # w(consumer) implies the lemma inequality with equality.
    consumer_uses: dict[str, list[str]] = {name: [] for name in parents}
    for inst in net.instances:
        for ref in inst.sources:
            if ref[0] == "out":
                consumer_uses[ref[1]].append(inst.name)
    weight: dict[str, int] = {}
    for name in reversed(order):
        weight[name] = 1 + sum(weight[c] for c in consumer_uses[name])
    for name in order:
        if weight[name] < 1 + sum(weight[c] for c in consumer_uses[name]):
            raise CertificateError(
                "WEIGHT_INEQUALITY", f"{net.name}: weight inequality fails at {name}"
            )

    ladder: dict[str, list[str]] = {}
    for name in order:
        ladder.setdefault(str(rank[name]), []).append(name)
    for names in ladder.values():
        names.sort()

    return {
        "order": order,
        "rank": rank,
        "depth": depth,
        "weight": weight,
        "consumer_uses": consumer_uses,
        "ladder": ladder,
    }


# ---------------------------------------------------------------------------
# Realized dynamics
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CompiledNet:
    net: Netlist
    analysis: dict
    inst_order: tuple[str, ...]  # fixed instance order = register layout order
    layouts: dict  # instance -> (offset, register count)
    tables: dict  # instance -> {port tuple: register tuple}; valid after verify_primitive
    src_slots: dict  # instance -> tuple of ("i", input index) or ("r", flat register index)
    register_count: int


def compile_net(net: Netlist, prims: Mapping[str, Primitive]) -> CompiledNet:
    analysis = check_netlist(net, prims)
    inst_order = tuple(inst.name for inst in net.instances)
    by_name = {inst.name: inst for inst in net.instances}
    layouts: dict[str, tuple[int, int]] = {}
    offset = 0
    for name in inst_order:
        n = len(prims[by_name[name].kind].registers)
        layouts[name] = (offset, n)
        offset += n
    input_index = {name: k for k, name in enumerate(net.inputs)}
    tables: dict[str, dict[tuple[int, ...], tuple[int, ...]]] = {}
    src_slots: dict[str, tuple[tuple[str, int], ...]] = {}
    for name in inst_order:
        inst = by_name[name]
        prim = prims[inst.kind]
        zero = (0,) * len(prim.registers)
        tables[name] = {
            inputs: prim.update[(inputs, zero)] for inputs in bit_tuples(len(prim.in_ports))
        }
        slots = []
        for ref in inst.sources:
            if ref[0] == "in":
                slots.append(("i", input_index[ref[1]]))
            else:
                src_inst = by_name[ref[1]]
                src_prim = prims[src_inst.kind]
                port_pos = src_prim.out_ports.index(ref[2])
                reg = src_prim.readback[port_pos]
                slots.append(("r", layouts[ref[1]][0] + reg))
        src_slots[name] = tuple(slots)
    return CompiledNet(
        net=net,
        analysis=analysis,
        inst_order=inst_order,
        layouts=layouts,
        tables=tables,
        src_slots=src_slots,
        register_count=offset,
    )


def realized_step(cn: CompiledNet, inputs_vec: tuple[int, ...], state: tuple[int, ...]) -> tuple[int, ...]:
    out: list[int] = []
    for name in cn.inst_order:
        key = tuple(
            inputs_vec[idx] if kind == "i" else state[idx] for kind, idx in cn.src_slots[name]
        )
        out.extend(cn.tables[name][key])
    return tuple(out)


def extension_state(cn: CompiledNet, inputs_vec: tuple[int, ...]) -> tuple[int, ...]:
    """The generated extension: settled register values in topological order."""

    values: list[int] = [0] * cn.register_count
    for name in cn.analysis["order"]:
        key = tuple(
            inputs_vec[idx] if kind == "i" else values[idx] for kind, idx in cn.src_slots[name]
        )
        offset, n = cn.layouts[name]
        values[offset : offset + n] = cn.tables[name][key]
    return tuple(values)


def output_values(cn: CompiledNet, state: tuple[int, ...]) -> dict[str, int]:
    result = {}
    for label, wire_name in cn.net.outputs:
        offset, _ = cn.layouts[wire_name]
        result[label] = state[offset]
    return result


def unsettled_instances(
    cn: CompiledNet, state: tuple[int, ...], extension: tuple[int, ...]
) -> list[str]:
    bad = []
    for name in cn.inst_order:
        offset, n = cn.layouts[name]
        if state[offset : offset + n] != extension[offset : offset + n]:
            bad.append(name)
    return bad


def potential(cn: CompiledNet, state: tuple[int, ...], extension: tuple[int, ...]) -> int:
    weight = cn.analysis["weight"]
    return sum(weight[name] for name in unsettled_instances(cn, state, extension))


def settle_trajectory(
    cn: CompiledNet,
    inputs_vec: tuple[int, ...],
    state: tuple[int, ...],
    extension: tuple[int, ...],
) -> tuple[int, int]:
    """Run the realized dynamics to the extension.

    Returns the exact settling time and the minimal observed potential margin
    ``potential(x) - potential(step(x)) - |unsettled(x)|`` over the run (a
    nonnegative margin certifies the strict-decrease lemma with the claimed
    slack on this trajectory).  Raises fail-closed if the depth bound is
    violated.
    """

    depth = cn.analysis["depth"]
    time = 0
    min_margin = None
    while state != extension:
        if time >= depth:
            raise CertificateError(
                "SETTLING_BOUND",
                f"{cn.net.name}: not settled after depth {depth} rounds from "
                f"inputs {inputs_vec}",
            )
        phi_before = potential(cn, state, extension)
        n_unsettled = len(unsettled_instances(cn, state, extension))
        state = realized_step(cn, inputs_vec, state)
        phi_after = potential(cn, state, extension)
        margin = phi_before - phi_after - n_unsettled
        if margin < 0:
            raise CertificateError(
                "POTENTIAL_INCREASE",
                f"{cn.net.name}: potential drop {phi_before - phi_after} below the "
                f"unsettled count {n_unsettled} at round {time}",
            )
        min_margin = margin if min_margin is None else min(min_margin, margin)
        time += 1
    fixed = realized_step(cn, inputs_vec, extension)
    if fixed != extension:
        raise CertificateError(
            "INTERTWINER_BROKEN", f"{cn.net.name}: extension is not a fixed point"
        )
    return time, (0 if min_margin is None else min_margin)


# ---------------------------------------------------------------------------
# Gate-level circuits and the gate-to-patch compiler
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GateCircuit:
    name: str
    inputs: tuple[str, ...]
    gates: tuple[tuple[str, tuple[str, str]], ...]  # (gate name, (src signal, src signal))
    outputs: tuple[str, ...]  # signal names
    semantics: Callable[[Mapping[str, int]], dict[str, int]]


def compile_gates_to_patches(circuit: GateCircuit) -> Netlist:
    """Compile a NAND netlist into a patch netlist.

    Every signal with more than one consumer is replicated through a balanced
    chain of fan-out-two patches, every gate becomes a NAND patch, and every
    circuit output is exposed through a wire patch.  The compiled netlist obeys
    the single-consumer wiring rule by construction.
    """

    uses: dict[str, int] = {}
    for _, (a, b) in circuit.gates:
        uses[a] = uses.get(a, 0) + 1
        uses[b] = uses.get(b, 0) + 1
    for sig in circuit.outputs:
        uses[sig] = uses.get(sig, 0) + 1

    gate_names = {name for name, _ in circuit.gates}
    instances: list[Instance] = []

    def source_of(signal: str) -> SourceRef:
        if signal in gate_names:
            return ("out", f"g_{signal}", "z")
        if signal in circuit.inputs:
            return ("in", signal)
        raise CertificateError("UNKNOWN_SIGNAL", f"{circuit.name}: missing signal {signal}")

    taps: dict[str, list[SourceRef]] = {}
    for signal in sorted(uses):
        need = uses[signal]
        available = [source_of(signal)]
        counter = 0
        while len(available) < need:
            ref = available.pop(0)
            fname = f"f_{signal}_{counter}"
            counter += 1
            instances.append(Instance(fname, "FANOUT2", (ref,)))
            available.append(("out", fname, "z0"))
            available.append(("out", fname, "z1"))
        taps[signal] = available

    def take(signal: str) -> SourceRef:
        return taps[signal].pop(0)

    for name, (a, b) in circuit.gates:
        instances.append(Instance(f"g_{name}", "NAND", (take(a), take(b))))
    outputs = []
    for sig in circuit.outputs:
        wname = f"w_{sig}"
        instances.append(Instance(wname, "WIRE", (take(sig),)))
        outputs.append((sig, wname))

    # Fan-out patches created before their driving gate exist in list order;
    # reorder topologically for the register layout while keeping determinism.
    net = Netlist(circuit.name, circuit.inputs, tuple(instances), tuple(outputs))
    analysis = check_netlist(net, reference_primitives())
    by_name = {inst.name: inst for inst in instances}
    ordered = tuple(by_name[n] for n in analysis["order"])
    return Netlist(circuit.name, circuit.inputs, ordered, tuple(outputs))


def xor_circuit() -> GateCircuit:
    def semantics(env: Mapping[str, int]) -> dict[str, int]:
        return {"n4": env["a"] ^ env["b"]}

    return GateCircuit(
        name="xor_from_nand",
        inputs=("a", "b"),
        gates=(
            ("n1", ("a", "b")),
            ("n2", ("a", "n1")),
            ("n3", ("b", "n1")),
            ("n4", ("n2", "n3")),
        ),
        outputs=("n4",),
        semantics=semantics,
    )


def tree_circuit() -> GateCircuit:
    def semantics(env: Mapping[str, int]) -> dict[str, int]:
        na = 1 - (env["i0"] & env["i1"])
        nb = 1 - (env["i2"] & env["i3"])
        nc = 1 - (env["i4"] & env["i5"])
        nd = 1 - (env["i6"] & env["i7"])
        ma = 1 - (na & nb)
        mb = 1 - (nc & nd)
        return {"root": 1 - (ma & mb)}

    return GateCircuit(
        name="nand_tree_depth3",
        inputs=("i0", "i1", "i2", "i3", "i4", "i5", "i6", "i7"),
        gates=(
            ("na", ("i0", "i1")),
            ("nb", ("i2", "i3")),
            ("nc", ("i4", "i5")),
            ("nd", ("i6", "i7")),
            ("ma", ("na", "nb")),
            ("mb", ("nc", "nd")),
            ("root", ("ma", "mb")),
        ),
        outputs=("root",),
        semantics=semantics,
    )


def adder_circuit() -> GateCircuit:
    def semantics(env: Mapping[str, int]) -> dict[str, int]:
        total = env["a0"] + 2 * env["a1"] + env["b0"] + 2 * env["b1"]
        return {"s0": total & 1, "s1": (total >> 1) & 1, "c2": (total >> 2) & 1}

    return GateCircuit(
        name="two_bit_adder_from_nand",
        inputs=("a0", "a1", "b0", "b1"),
        gates=(
            # bit 0 half adder: sum s0, carry c1
            ("n1", ("a0", "b0")),
            ("n2", ("a0", "n1")),
            ("n3", ("b0", "n1")),
            ("s0", ("n2", "n3")),
            ("c1", ("n1", "n1")),
            # bit 1 full adder on (a1, b1, c1): sum s1, carry c2
            ("m1", ("a1", "b1")),
            ("m2", ("a1", "m1")),
            ("m3", ("b1", "m1")),
            ("m4", ("m2", "m3")),
            ("m5", ("m4", "c1")),
            ("m6", ("m4", "m5")),
            ("m7", ("c1", "m5")),
            ("s1", ("m6", "m7")),
            ("c2", ("m5", "m1")),
        ),
        outputs=("s0", "s1", "c2"),
        semantics=semantics,
    )


# ---------------------------------------------------------------------------
# Circuit-level verification
# ---------------------------------------------------------------------------

SETTLING_CONSTANT_C = 1  # rounds per compiled rank, from the primitive verification


def verify_extension_semantics(circuit: GateCircuit, cn: CompiledNet) -> int:
    """The generated extension of the compiled net equals circuit evaluation."""

    checked = 0
    for inputs_vec in bit_tuples(len(cn.net.inputs)):
        env = dict(zip(cn.net.inputs, inputs_vec))
        expected_signals = circuit.semantics(env)
        expected = {label: expected_signals[label] for label, _ in cn.net.outputs}
        got = output_values(cn, extension_state(cn, inputs_vec))
        if got != expected:
            raise CertificateError(
                "GATE_TRUTH",
                f"{circuit.name}: extension readback {got} differs from circuit "
                f"evaluation {expected} at inputs {env}",
            )
        checked += 1
    return checked


def exhaustive_settling_report(circuit: GateCircuit, cn: CompiledNet) -> dict[str, Any]:
    """All inputs times all internal register states: exact settling data."""

    depth = cn.analysis["depth"]
    worst_time = 0
    min_margin = None
    trajectories = 0
    for inputs_vec in bit_tuples(len(cn.net.inputs)):
        extension = extension_state(cn, inputs_vec)
        for init in bit_tuples(cn.register_count):
            time, margin = settle_trajectory(cn, inputs_vec, init, extension)
            worst_time = max(worst_time, time)
            if time > 0:
                min_margin = margin if min_margin is None else min(min_margin, margin)
            trajectories += 1
    if worst_time > SETTLING_CONSTANT_C * depth:
        raise CertificateError(
            "SETTLING_BOUND",
            f"{circuit.name}: worst settling time {worst_time} exceeds "
            f"c*depth = {SETTLING_CONSTANT_C * depth}",
        )
    return {
        "coverage": "exhaustive over all inputs and all internal register states",
        "trajectories": trajectories,
        "worst_settling_time": worst_time,
        "settling_bound_c_times_depth": SETTLING_CONSTANT_C * depth,
        "bound_verified": True,
        "min_potential_margin": 0 if min_margin is None else min_margin,
    }


def declared_family(cn: CompiledNet, extension: tuple[int, ...]) -> Iterator[tuple[int, ...]]:
    zero = (0,) * cn.register_count
    one = (1,) * cn.register_count
    complement = tuple(1 - b for b in extension)
    yield extension
    yield zero
    yield one
    yield complement
    for k in range(cn.register_count):
        flipped = list(extension)
        flipped[k] = 1 - flipped[k]
        yield tuple(flipped)
        flipped_zero = list(zero)
        flipped_zero[k] = 1
        yield tuple(flipped_zero)


def family_settling_report(circuit: GateCircuit, cn: CompiledNet) -> dict[str, Any]:
    """Declared initial-state family for nets too wide to exhaust.

    The family per input assignment: the extension, all-zero, all-one, the
    bitwise complement of the extension, every single-register corruption of
    the extension, and every single-register corruption of all-zero.  The
    all-initial-state depth bound for this net is carried by the composition
    checker's rank ladder together with the exhaustive primitive verification,
    not by this enumeration.
    """

    depth = cn.analysis["depth"]
    worst_time = 0
    min_margin = None
    trajectories = 0
    for inputs_vec in bit_tuples(len(cn.net.inputs)):
        extension = extension_state(cn, inputs_vec)
        for init in declared_family(cn, extension):
            time, margin = settle_trajectory(cn, inputs_vec, init, extension)
            worst_time = max(worst_time, time)
            if time > 0:
                min_margin = margin if min_margin is None else min(min_margin, margin)
            trajectories += 1
    if worst_time > SETTLING_CONSTANT_C * depth:
        raise CertificateError(
            "SETTLING_BOUND",
            f"{circuit.name}: worst settling time {worst_time} exceeds "
            f"c*depth = {SETTLING_CONSTANT_C * depth}",
        )
    return {
        "coverage": (
            "declared finite family: extension, all-zero, all-one, complement of "
            "the extension, all single-register corruptions of the extension and "
            "of all-zero, for every input assignment; the all-initial-state bound "
            "is the composition-checker rank-ladder lemma"
        ),
        "trajectories": trajectories,
        "worst_settling_time": worst_time,
        "settling_bound_c_times_depth": SETTLING_CONSTANT_C * depth,
        "bound_verified": True,
        "min_potential_margin": 0 if min_margin is None else min_margin,
    }


def rank_ladder_certification(cn: CompiledNet, prims: Mapping[str, Primitive]) -> dict[str, Any]:
    """The checked settling lemma over the DAG structure.

    Premise A (per primitive, checked exhaustively by ``verify_primitive``):
    one update from any internal state settles a patch to the truth of its
    port values, and the update reads the declared ports alone.  Premise B
    (per instance, checked here): every source of a rank-k instance has rank
    below k.  Conclusion, by induction on the rank ladder: after round k every
    instance of rank at most k holds its extension value from any initial
    register state, so the realized dynamics settles in at most
    c * depth rounds with c = 1.
    """

    rank = cn.analysis["rank"]
    by_name = {inst.name: inst for inst in cn.net.instances}
    for name in cn.inst_order:
        inst = by_name[name]
        if inst.kind not in prims:
            raise CertificateError("UNKNOWN_PRIMITIVE", f"{name} kind {inst.kind}")
        for ref in inst.sources:
            if ref[0] == "out" and rank[ref[1]] >= rank[name]:
                raise CertificateError(
                    "RANK_LADDER", f"{name} at rank {rank[name]} reads {ref[1]} at equal or higher rank"
                )
    ladder = cn.analysis["ladder"]
    return {
        "layers": {k: len(v) for k, v in sorted(ladder.items(), key=lambda kv: int(kv[0]))},
        "premise_per_primitive": "one-round settling and port locality, exhaustively checked",
        "premise_per_instance": "every source sits at strictly lower rank",
        "conclusion": (
            "after round k every instance of rank at most k is settled from any "
            "initial register state; the realized dynamics settles in at most "
            "c * depth rounds with c = 1"
        ),
        "verified": True,
    }


# ---------------------------------------------------------------------------
# Potential lemma and the remaining-depth countermodel
# ---------------------------------------------------------------------------


def potential_lemma_statement() -> dict[str, Any]:
    return {
        "weight": (
            "w(p) = 1 + sum over consuming port uses of w(consumer), computed in "
            "reverse topological order (the number of directed paths leaving p, "
            "plus one)"
        ),
        "potential": "Phi(x) = sum of w(p) over instances p not holding their extension value",
        "lemma": (
            "Phi(step(x)) <= Phi(x) - |unsettled(x)|: an instance is unsettled "
            "after the step only if some source instance was unsettled before it "
            "(port locality plus per-primitive one-round settling), so the "
            "post-step unsettled set is covered by the consumers of the pre-step "
            "unsettled set, and the weight inequality w(p) >= 1 + sum of "
            "consumer weights telescopes the sum"
        ),
        "checked_premises": [
            "per-primitive port locality and one-round settling (exhaustive)",
            "per-instance weight inequality over the DAG (exact integer check)",
        ],
    }


def depth_weight_countermodel(prims: Mapping[str, Primitive]) -> dict[str, Any]:
    """Remaining-depth weighting is not strictly decreasing under fan-out.

    The net: input x feeds a fan-out patch F1 whose branches feed fan-out
    patches F2 and F3; their four branches feed wires W1..W4, which feed wires
    V1..V4.  Witness state: only F2 and F3 unsettled, with both registers
    complemented.  One step settles F2 and F3 and unsettles W1..W4.  With the
    remaining-depth weight D - rank + 1 the potential rises from 6 to 8; with
    the path-count weight it falls from 10 to 8, matching the lemma margin.
    """

    instances = [
        Instance("F1", "FANOUT2", (("in", "x"),)),
        Instance("F2", "FANOUT2", (("out", "F1", "z0"),)),
        Instance("F3", "FANOUT2", (("out", "F1", "z1"),)),
        Instance("W1", "WIRE", (("out", "F2", "z0"),)),
        Instance("W2", "WIRE", (("out", "F2", "z1"),)),
        Instance("W3", "WIRE", (("out", "F3", "z0"),)),
        Instance("W4", "WIRE", (("out", "F3", "z1"),)),
        Instance("V1", "WIRE", (("out", "W1", "z"),)),
        Instance("V2", "WIRE", (("out", "W2", "z"),)),
        Instance("V3", "WIRE", (("out", "W3", "z"),)),
        Instance("V4", "WIRE", (("out", "W4", "z"),)),
    ]
    net = Netlist(
        "depth_weight_countermodel",
        ("x",),
        tuple(instances),
        (("v1", "V1"), ("v2", "V2"), ("v3", "V3"), ("v4", "V4")),
    )
    cn = compile_net(net, prims)
    depth = cn.analysis["depth"]
    rank = cn.analysis["rank"]
    inputs_vec = (0,)
    extension = extension_state(cn, inputs_vec)
    witness = list(extension)
    for name in ("F2", "F3"):
        offset, n = cn.layouts[name]
        for k in range(offset, offset + n):
            witness[k] = 1 - witness[k]
    state = tuple(witness)

    def depth_phi(s: tuple[int, ...]) -> int:
        return sum(
            depth - rank[name] + 1 for name in unsettled_instances(cn, s, extension)
        )

    after = realized_step(cn, inputs_vec, state)
    record = {
        "net_depth": depth,
        "witness_unsettled_before": unsettled_instances(cn, state, extension),
        "witness_unsettled_after": unsettled_instances(cn, after, extension),
        "depth_weight_potential_before": depth_phi(state),
        "depth_weight_potential_after": depth_phi(after),
        "path_count_potential_before": potential(cn, state, extension),
        "path_count_potential_after": potential(cn, after, extension),
        "finding": (
            "the potential weighted by remaining depth increases across this "
            "step, so it is not a settling potential for fan-out circuits; the "
            "path-count weight decreases by exactly the unsettled count and is "
            "the potential this certificate uses"
        ),
    }
    if not record["depth_weight_potential_after"] > record["depth_weight_potential_before"]:
        raise CertificateError("COUNTERMODEL", "depth-weight countermodel does not witness an increase")
    if not (
        record["path_count_potential_before"] - record["path_count_potential_after"]
        >= len(record["witness_unsettled_before"])
    ):
        raise CertificateError("COUNTERMODEL", "path-count potential margin fails on the witness")
    return record


# ---------------------------------------------------------------------------
# Fail-closed controls
# ---------------------------------------------------------------------------


def expect_failure(code: str, thunk: Callable[[], Any]) -> dict[str, Any]:
    try:
        thunk()
    except CertificateError as exc:
        if exc.code != code:
            raise CertificateError(
                "CONTROL_NOT_CLOSED", f"expected failure code {code}, got {exc.code}"
            )
        return {"verdict": "fails_closed", "error_code": exc.code, "witness": exc.message}
    raise CertificateError("CONTROL_NOT_CLOSED", f"control expecting {code} did not fail")


def run_controls(prims: Mapping[str, Primitive]) -> dict[str, Any]:
    controls: dict[str, Any] = {}

    controls["noisy_update_table"] = expect_failure(
        "INTERTWINER_BROKEN", lambda: verify_primitive(noisy_nand_control())
    )
    controls["noisy_update_table"]["construction"] = (
        "NAND with the (1,1) table entry flipped to 1; the intertwiner check "
        "detects the broken truth-table commutation"
    )

    controls["cross_talk_fanout"] = expect_failure(
        "CROSS_TALK", lambda: verify_primitive(crosstalk_fanout_control())
    )
    controls["cross_talk_fanout"]["construction"] = (
        "fan-out whose first branch xors the sibling register into its update; "
        "the state-independence sweep detects the sibling dependence"
    )

    cyclic = Netlist(
        "cyclic_control",
        ("x",),
        (
            Instance("W1", "WIRE", (("out", "W2", "z"),)),
            Instance("W2", "WIRE", (("out", "W1", "z"),)),
        ),
        (),
    )
    controls["cyclic_layout"] = expect_failure(
        "CYCLIC_LAYOUT", lambda: check_netlist(cyclic, prims)
    )
    controls["cyclic_layout"]["construction"] = "two wire patches feeding each other"

    overload = Netlist(
        "fan_in_overload_control",
        ("x", "y", "z"),
        (Instance("G", "NAND", (("in", "x"), ("in", "y"), ("in", "z"))),),
        (),
    )
    controls["fan_in_overload"] = expect_failure(
        "PORT_ARITY", lambda: check_netlist(overload, prims)
    )
    controls["fan_in_overload"]["construction"] = "three sources wired into a two-port NAND"

    multi = Netlist(
        "multi_consumer_control",
        ("x",),
        (Instance("G", "NAND", (("in", "x"), ("in", "x"))),),
        (),
    )
    controls["multi_consumer_port"] = expect_failure(
        "MULTI_CONSUMER", lambda: check_netlist(multi, prims)
    )
    controls["multi_consumer_port"]["construction"] = (
        "one input feeding two ports without a fan-out patch"
    )

    return controls


# ---------------------------------------------------------------------------
# Manifest assembly
# ---------------------------------------------------------------------------


def circuit_block(circuit: GateCircuit, prims: Mapping[str, Primitive], exhaustive: bool) -> dict[str, Any]:
    net = compile_gates_to_patches(circuit)
    cn = compile_net(net, prims)
    kinds: dict[str, int] = {}
    for inst in net.instances:
        kinds[inst.kind] = kinds.get(inst.kind, 0) + 1
    semantic_checks = verify_extension_semantics(circuit, cn)
    if exhaustive:
        settling = exhaustive_settling_report(circuit, cn)
    else:
        settling = family_settling_report(circuit, cn)
    ladder = rank_ladder_certification(cn, prims)
    return {
        "gate_count": len(circuit.gates),
        "patch_counts": dict(sorted(kinds.items())),
        "patch_total": len(net.instances),
        "register_count": cn.register_count,
        "compiled_depth": cn.analysis["depth"],
        "extension_semantics_checks": semantic_checks,
        "settling": settling,
        "rank_ladder": ladder,
    }


def build_manifest() -> dict[str, Any]:
    prims = reference_primitives()
    primitive_reports = {name: verify_primitive(prim) for name, prim in sorted(prims.items())}

    circuits = {
        "xor_from_nand": circuit_block(xor_circuit(), prims, exhaustive=True),
        "nand_tree_depth3": circuit_block(tree_circuit(), prims, exhaustive=True),
        "two_bit_adder_from_nand": circuit_block(adder_circuit(), prims, exhaustive=False),
    }

    fanout_report = primitive_reports["FANOUT2"]["evidence_bundle"]
    manifest = {
        "schema": SCHEMA,
        "issue": 328,
        "arithmetic": "exact integer and Boolean arithmetic; no floating point in any check",
        "abstract_compiler_result": ABSTRACT_COMPILER_RESULT,
        "primitives": primitive_reports,
        "intertwiner": {
            "statement": (
                "decode is the register readback and encode is the identity on "
                "register bits; for every primitive, every port assignment, and "
                "every internal register state, one patch update commutes exactly "
                "with the single-site ranked-system update, checked exhaustively; "
                "for compiled nets the realized synchronous step is the product of "
                "these local updates over the single-consumer wiring, so the "
                "realized dynamics is the synchronous dynamics of the ranked "
                "functional system carried by the patch DAG"
            ),
            "kind": "finite bisimulation, exact",
        },
        "composition_checker": {
            "rules": [
                "known primitive kinds only",
                "port arity matches the primitive interface (fan-in overload rejected)",
                "every source reference resolves (dangling references rejected)",
                "every port source feeds at most one consumer (replication requires a fan-out patch)",
                "acyclic wiring via topological sort (cycles rejected)",
                "rank ladder: every source sits at strictly lower rank",
                "path-count weight inequality checked per instance",
            ],
            "structural_induction": (
                "acyclic composition of exhaustively verified primitives is "
                "verified: the rank-ladder certification instantiates the "
                "induction layer by layer on each compiled circuit"
            ),
        },
        "settling_theorem": {
            "constant_c": SETTLING_CONSTANT_C,
            "constant_source": (
                "every primitive settles in exactly one round from any internal "
                "state once its ports are settled, so c = 1 round per compiled rank"
            ),
            "bound": "T(circuit) <= c * compiled_depth, uniform in circuit width",
            "size_uniformity": (
                "the bound depends on the compiled DAG depth alone; width enters "
                "only through the patch and register counts reported per circuit"
            ),
        },
        "potential": {
            **potential_lemma_statement(),
            "depth_weight_countermodel": depth_weight_countermodel(prims),
        },
        "fan_out_control": {
            "non_interference": (
                "the FANOUT2 state-independence sweep varies both branch registers "
                "over all values at every port value and requires identical branch "
                "updates, so neither branch reads its sibling"
            ),
            "checked_port_state_pairs": fanout_report["checked_port_state_pairs"],
        },
        "controls": run_controls(prims),
        "circuits": circuits,
        "resource_bounds": {
            "statement": (
                "compiled patches = gates + one fan-out patch per extra consumer "
                "of each signal + one wire patch per circuit output; registers = "
                "gates + 2 * fan-out patches + outputs; settling time is at most "
                "c * compiled depth from any initial register state"
            ),
            "per_circuit_witnesses": {
                name: {
                    "patch_total": block["patch_total"],
                    "register_count": block["register_count"],
                    "compiled_depth": block["compiled_depth"],
                    "worst_settling_time": block["settling"]["worst_settling_time"],
                }
                for name, block in sorted(circuits.items())
            },
        },
        "scope": {
            "realized": (
                "finite software patch dynamics: the three primitive patches, the "
                "compiled test circuits, and the composition checker defined in "
                "code/consensus/compiled_lattice_settling_certificate.py"
            ),
            "open": (
                "continuum-limit and physical-hardware attachment of the patch "
                "dynamics is open; no laboratory device, photonic chamber, or "
                "continuum field realization is claimed by this packet"
            ),
        },
    }
    return manifest


def verify_manifest(path: Path) -> dict[str, Any]:
    stored = load_json(path)
    recomputed = build_manifest()
    if stored != recomputed:
        raise CertificateError(
            "MANIFEST_MISMATCH", f"stored manifest {path} differs from recomputation"
        )
    return recomputed


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command")
    emit = sub.add_parser("emit", help="recompute the certificate and write the manifest")
    emit.add_argument("--output", type=Path, default=DEFAULT_MANIFEST)
    verify = sub.add_parser("verify", help="recompute and compare against a stored manifest")
    verify.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args(argv)

    if args.command == "verify":
        manifest = verify_manifest(args.manifest)
        print(f"manifest verified: {args.manifest}")
    else:
        output = args.output if args.command == "emit" else DEFAULT_MANIFEST
        manifest = build_manifest()
        write_json(output, manifest)
        print(f"manifest written: {output}")
    print(f"manifest sha256: {sha256_json(manifest)}")
    for name, block in sorted(manifest["circuits"].items()):
        settling = block["settling"]
        print(
            f"{name}: depth {block['compiled_depth']}, "
            f"worst settling time {settling['worst_settling_time']}, "
            f"bound {settling['settling_bound_c_times_depth']}, "
            f"registers {block['register_count']}"
        )
    print(f"settling constant c = {SETTLING_CONSTANT_C}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
