"""Complete finite Kraus catalog, including the invalid chart outcome."""
import json
from model import CALLS, EDGES, eraser, fine, star


def encode(ops):
    rows = [{"matrix": [[str(x) for x in (*r, 0)] for r in a] + [["0"] * 4],
             "squared_denominator": d, "outcome": k} for k, (a, d) in enumerate(ops)]
    rows.append({"matrix": [[str(int(i == 3 and j == 3)) for j in range(4)] for i in range(4)],
                 "squared_denominator": 1, "outcome": "invalid_chart"})
    return rows


catalog = {}
for law in ("native", "swap"):
    catalog[law] = {}
    for vertex in range(3):
        branches = fine(tuple(e for e in EDGES if vertex in e), law)
        catalog[law][str(vertex)] = {"fine": encode([(a, 4) for a in branches]),
                                    "erased": encode(eraser(branches))}
print(json.dumps({"catalog": catalog, "literal_endpoint_calls": CALLS,
                  "complement_convention": "invalid fourth chart state has explicit orthogonal absorbing rejection"}, indent=2))
