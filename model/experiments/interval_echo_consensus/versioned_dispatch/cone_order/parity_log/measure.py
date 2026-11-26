"""Reconstruct forcing supports and authenticated ancestry from physical tapes."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from . import apparatus
from .apparatus import MASKS, ROOTS, SHELVES, WIDTH

HERE = Path(__file__).resolve().parent


def decode(raw):
    if len(raw) != WIDTH or any(code not in "69" for code in raw):
        return None
    return sum((c == "9") << k for k, c in enumerate(raw))


def final_values(launch):
    if launch["final"][0] != launch["final"][1]:
        raise ValueError("Final heads disagree")
    words = [decode(raw) for raw in launch["final"][0]]
    if any(word is None or word >> 1 != a + 1 for a, word in enumerate(words)):
        raise ValueError("Incomplete or incorrectly stamped nominal final array")
    return [word & 1 for word in words]


def supports(table, word, mask):
    """All measurements with the given restricted input must agree at this output."""
    target = table[word][ROOTS + mask]
    forcing = [
        u
        for u in range(MASKS)
        if all(
            row[ROOTS + mask] == target
            for other, row in enumerate(table)
            if ((other ^ word) & u) == 0
        )
    ]
    least = [u for u in forcing if all(u & v == u for v in forcing)]
    if len(least) != 1:
        raise ValueError("No unique least forcing support")
    if set(forcing) != {u for u in range(MASKS) if least[0] & u == least[0]}:
        raise ValueError("Forcing sets do not form the expected upward closure")
    return least[0], sum(1 << u for u in forcing)


def authenticate(launch, mount, least):
    if len(launch["events"]) != SHELVES:
        raise ValueError("Incomplete event history")
    parent_ids = []
    ancestors = []
    heights = []
    posts = {}
    for a, event in enumerate(launch["events"]):
        eid = event["id"]
        if event["target"] != a or eid != a + 1 or a in posts:
            raise ValueError("An event lost its own fresh target")
        publish = event["publication"]
        if publish is None:
            raise ValueError("Expected a full authenticated launch")
        if decode(publish["before"]) != 0 or len(publish["peers"]) != 2:
            raise ValueError("The event did not seal a new target")
        raw = publish["after"]
        word = decode(raw)
        if not (
            word is not None and word >> 1 == eid and publish["peers"] == [raw, raw]
        ):
            raise ValueError(
                "The claimed writer or value disagrees with peer readbacks"
            )
        if launch["final"][0][a] != raw or launch["final"][1][a] != raw:
            raise ValueError("A sealed version did not persist")
        expected_change = {
            p for k, pair in enumerate(mount["words"][a]) if word >> k & 1 for p in pair
        }
        if expected_change != set(publish["changed_ports"]):
            raise ValueError("A sealed write changed the wrong ports")

        support = 0 if a < ROOTS else least[a - ROOTS]
        actual = {}
        for source, first, second in event["reads"]:
            value = decode(first)
            if source in actual or first != second or value is None:
                raise ValueError(
                    "Source reads are not independent matching addressed words"
                )
            actual[source] = value
        if set(actual) != {i for i in range(ROOTS) if support >> i & 1}:
            raise ValueError(
                "Pivotal source set and the physically read register set disagree"
            )

        parents = []
        for reg, value in actual.items():
            writer = value >> 1
            source_event = launch["events"][writer - 1] if 1 <= writer < eid else None
            if (
                source_event is None
                or source_event["target"] != reg
                or posts.get(reg) != value
            ):
                raise ValueError("Failed writer, target, or payload authentication")
            parents.append(writer)
        parent_ids.append(parents)
        generated = set(parents)
        for parent in parents:
            generated.update(ancestors[parent - 1])
        ancestors.append(sorted(generated))
        heights.append(max((heights[p - 1] + 1 for p in parents), default=0))
        posts[a] = word
    return dict(parents=parent_ids, ancestors=ancestors, heights=heights)


def interpret(bank, mount):
    by_word = sorted(bank, key=lambda r: r["input"])
    if [e["input"] for e in by_word] != list(range(MASKS)):
        raise ValueError("Missing source completion")
    table = [final_values(row) for row in by_word]
    result = []
    for word, launch in enumerate(by_word):
        if table[word][:ROOTS] != [(word >> i) & 1 for i in range(ROOTS)]:
            raise ValueError("Enrollment was not independently recovered")
        proven = [supports(table, word, m) for m in range(MASKS)]
        least = [u for u, _ in proven]
        auth = authenticate(launch, mount, least)
        result.append(
            dict(
                input=word,
                least=least,
                forcing_hex=[f"{mask:08x}" for _, mask in proven],
                **auth,
            )
        )
    return result


def report():
    bench = apparatus.enroll()
    mount = apparatus.layout(bench)
    banks = {
        mode: [apparatus.run(bench, w, mode) for w in range(MASKS)]
        for mode in ("native", "altered")
    }
    controls = [apparatus.run(bench, 21, mode) for mode in ("shutter", "stamp")]
    decoded = {mode: interpret(records, mount) for mode, records in banks.items()}
    from ....seam_word_meter.finite_program import machine
    from ....seam_word_meter.finite_program.layout import mount as routes

    paths = [
        HERE / "PROTOCOL.md",
        Path(apparatus.__file__),
        Path(__file__),
        Path(machine.__file__),
        Path(routes.__file__),
        Path(machine.LAW.__file__),
    ]
    root = HERE.parents[5]
    return dict(
        provenance="Complete deterministic sharp simulated material pointer readings. Every input is enrolled.",
        selectors="5 unconditional roots, followed by 32 simultaneously retained output shelves in increasing mask order",
        declared_masks=list(range(MASKS)),
        sources={
            str(f.resolve().relative_to(root)): hashlib.sha256(
                f.read_bytes()
            ).hexdigest()
            for f in paths
        },
        mount=mount,
        banks=banks,
        controls=controls,
        analysis=decoded,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = (json.dumps(report(), separators=(",", ":")) + "\n").encode()
    if args.output:
        with args.output.open("xb") as handle:
            handle.write(gzip.compress(data, mtime=0))
    elif gzip.decompress((HERE / "records.json.gz").read_bytes()) != data:
        raise ValueError("Single-log replay differs from the deposited record")
