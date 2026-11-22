"""Replay the three installed federations and compare their indexed records."""

import hashlib

from .cells import cell_inputs, read_cell
from .equipment import install
from .observed import need
from .priority import authenticate, compare, observe
from .source import PIN, bind_source, board_hash, canonical, fingerprint


def replay_level(level, shapes):
    screen, native = install(level, "native")
    _, reverse = install(level, "reverse")
    _, relabeled = install(level, "relabeled")
    relations = []
    transcripts = {}
    for name, service in (("native", native), ("reverse", reverse), ("relabeled", relabeled)):
        document = observe(level, screen, service)
        relations.append(authenticate(level, screen, service, document))
        if name != "relabeled":
            transcripts[name] = fingerprint(document)
    relation, opposite, same = relations
    separated = compare(screen, (native, reverse), (relation, opposite))
    need(separated["addressed_hits"] == [[1, 0], [0, 1]], "priority reversal did not separate")
    need(compare(screen, (native, relabeled), (relation, same))["equivalent"],
         "relabeling changed the incident comparisons")
    material = hashlib.sha256()
    comparison = hashlib.sha256()
    rows = []
    board = board_hash(level)
    for cell, (local, boundary) in enumerate(cell_inputs(screen, relation.order)):
        result = read_cell(screen, native, relation, board, cell, local, boundary)
        indexed_words = []
        comparison_rows = []
        for witness in result["witnesses"]:
            mask = witness["mask"]
            word = witness["internal_word"]
            indexed_words.append(word)
            row = [cell, mask, witness["signed_marks"], witness["final_loads"]]
            material.update(canonical(row).encode("ascii") + b"\n")
            bits = tuple(witness["boundary_hits"])
            external = tuple(1 - bit for bit in bits)
            accepted = tuple(sorted(local[j] for j in word))
            comparison_rows.append((cell, external, accepted, bits))
        for comparison_row in sorted(comparison_rows):
            comparison.update(repr(comparison_row).encode("ascii"))
        ports = tuple(result["ports"])
        need(shapes.setdefault(ports, indexed_words) == indexed_words,
             "repeated port layout gave a different word")
        rows.append(dict(cell=cell, ports=result["ports"], boundary=boundary,
                         compensation=result["compensation"], internal=result["internal_addresses"],
                         certificate=result["certificate"]["order"], words=indexed_words))
    return dict(level=level, board=board, seams=len(screen.edges), wedges=len(screen.wedges),
                cells=len(rows), environments=8 * len(rows), transcripts_sha256=transcripts,
                cell_records_sha256=fingerprint(rows), material_sha256=material.hexdigest(),
                comparison_sha256=comparison.hexdigest())


def replay():
    sources, profile = bind_source()
    expected = sources["replay"]
    shapes = {}
    levels = []
    for level in range(3):
        measured = replay_level(level, shapes)
        need(measured == expected["levels"][level], f"replay differs at level {level}")
        levels.append(measured)
    result = dict(source_commit=PIN, levels=levels, cells=sum(s["cells"] for s in levels),
                  environments=sum(s["environments"] for s in levels), shapes=len(shapes))
    need(result == expected, "replay summary differs")
    return dict(result, matched_source_profile=profile)
