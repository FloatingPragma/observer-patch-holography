"""Generate and independently check the least predecessor-closed patch of a demand."""

from .observed import need


def certificate(relation, board, requested):
    members = set(requested)
    for address in reversed(relation.order):
        if address in members:
            members.update(relation.higher[address])
    order = [address for address in relation.order if address in members]
    result = dict(board=board, requested=sorted(requested), order=order)
    validate(relation, board, requested, result)
    return result


def validate(relation, board, requested, cert):
    need(isinstance(cert, dict) and set(cert) == {"board", "requested", "order"},
         "missing or unknown certificate field")
    need(cert["board"] == board, "board identity")
    need(isinstance(cert["requested"], list) and all(type(i) is int for i in cert["requested"])
         and cert["requested"] == sorted(requested), "requested patch")
    order = cert["order"]
    need(isinstance(order, list) and all(type(i) is int and 0 <= i < len(relation.edges)
                                         for i in order), "invalid address")
    visited = set()
    for i in order:
        need(i not in visited, "duplicate address")
        need(visited.issuperset(relation.higher[i]), "missing or unordered predecessor")
        visited.add(i)
    needed = set(requested)
    frontier = list(requested)
    while frontier:
        for parent in relation.higher[frontier.pop()]:
            if parent not in needed:
                needed.add(parent)
                frontier.append(parent)
    need(visited == needed, "incorrect closure")
