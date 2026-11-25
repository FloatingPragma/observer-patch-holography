"""Read the same head source with whole-region momentum instruments at graph radii zero and one."""

from . import mounts, native

RADII = (0, 1)
HERE = mounts.HERE
EFFECTS = HERE / "regional_effects.json"


def apparatus():
    board, route, launch, patches = mounts.apparatus()
    atlas = mounts.cuts.native.Atlas(1)
    domains = []
    for p in patches:
        for radius in RADII:
            indices = [
                i
                for i, u in enumerate(p["cells"])
                if atlas.dist[p["center"]][u] <= radius
            ]
            restricted = dict(p)
            for key in ("cells", "images", "lifts"):
                restricted[key] = [p[key][i] for i in indices]
            domains.append(restricted)
    return board, route, launch, domains


def effects():
    keys = [
        [(q, n) for q in (-1, 1) for n in mounts.REGIONAL]
        for _ in range(len(RADII) * len(mounts.MOUNTS))
    ]
    arms = native.exact_rows(keys, apparatus)
    for i, a in enumerate(arms):
        a["radius"] = RADII[i % len(RADII)]
    return arms


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    a = parser.parse_args()
    records = effects()
    if a.write:
        native.store(records, EFFECTS)
    else:
        assert records == json.loads(EFFECTS.read_text())["charts"]
    print("Verified the two nested whole-region momentum instruments on each mount.")
