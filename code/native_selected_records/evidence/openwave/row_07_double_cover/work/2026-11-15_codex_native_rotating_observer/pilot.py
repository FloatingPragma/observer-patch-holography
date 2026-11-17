"""Exploratory cases frozen before execution; prints complete fixed pilot set."""
import json
from model import CALLS, density, prep, run

cases = [run(h, turns, law, direction) for h in (1, 2, 4, 6)
         for turns in (1, 2) for law in ("native", "swap") for direction in (1, -1)]
preparation = {law: {name: str(prep(law, density(v))[0])
                    for name, v in {"basis": (1, 0, 0), "uniform": (1, 1, 1),
                                    "contrast": (1, -1, 0), "third": (0, 0, 1),
                                    "large_uniform": (10**18+1, 10**18-1, 10**18)}.items()}
               for law in ("native", "swap")}
print(json.dumps({"kind": "prospectively frozen exploratory run; not confirmation",
                  "cases": cases, "preparation": preparation,
                  "literal_endpoint_calls": CALLS}, indent=2))
