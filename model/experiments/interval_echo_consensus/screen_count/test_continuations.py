"""Audit the full closure from port readings and recover each source with a label log."""

from .test_count import receipt, word


def test_all_word_lengths_and_public_labels():
    r = receipt("continuations")
    g = r["generators"]
    for j, op in enumerate(g):
        assert len(set(op)) == 256
        delta = op[0]
        assert all(y == x ^ delta for x, y in enumerate(op))
        assert delta.bit_count() == (2 if j < 7 else 1)
    for origin, raw in enumerate(r["boundary_readings"]):
        destination, exterior = word(raw, 8)
        assert destination == g[-1][origin]
        assert exterior in (0, 1)
    for name, data in r["families"].items():
        count = 7 if name == "repeat_inside" else 8
        routes = data["routes"]
        assert len(routes) == 2**count
        reached = set()
        for route in routes:
            observed = tuple(route["terminal"])
            labels = route["labels"]
            assert set(labels) <= set(range(count)) and len(set(labels)) == len(labels)
            assert observed == tuple(x ^ observed[0] for x in range(256))
            assert observed[0] not in reached
            reached.add(observed[0])
            for target in range(256):
                start = target
                for j in reversed(labels):
                    start = g[j][start]
                assert observed[start] == target
        assert all((a ^ b) in reached for a in reached for b in reached)
        assert len(data["classes"]) == (2 if count == 7 else 1)
        assert all(
            set(group) == {group[0] ^ shift for shift in reached}
            for group in data["classes"]
        )
