"""Bind wedge observations to a carrier and compare the recovered relations."""

from .equipment import apparatus, collect, tensor_wedge
from .observed import need, reconstruct, separator
from .source import PIN, board_hash, canonical


def observe(level, screen, source):
    return dict(format=1, pin=PIN, level=level, board=board_hash(level),
                rows=collect(screen, apparatus(source)))


def recover_document(level, screen, doc):
    need(isinstance(doc, dict) and set(doc) == {"format", "pin", "level", "board", "rows"} and
         canonical([doc[k] for k in ("format", "pin", "level", "board")]) ==
         canonical([1, PIN, level, board_hash(level)]), "transcript binding")
    return reconstruct(screen, doc["rows"])


def authenticate(level, screen, source, doc):
    recovered = recover_document(level, screen, doc)
    measure = apparatus(source)
    for i, j, a, b in doc["rows"]:
        need([a, b] == measure(i, j, screen.wedges[i, j]), "false addressed read")
        need(int(i in recovered.higher[j]) == a and int(j in recovered.higher[i]) == b,
             "inferred winner contradicts authenticated hit")
    return recovered


def compare(screen, services, relations):
    witness = separator(screen, *relations)
    if witness is None:
        return dict(equivalent=True, experiment=None)
    pair = tuple(witness["pair"])
    evidence = [tensor_wedge(screen, source, witness) for source in services]
    for rel, bits in zip(relations, evidence):
        expected = int(pair[0] in rel.higher[pair[1]])
        need(bits == [expected, 1 - expected], "physical separator mismatch")
    need(evidence[0][0] != evidence[1][0], "same addressed readout")
    return dict(equivalent=False, experiment=witness, addressed_hits=evidence)
