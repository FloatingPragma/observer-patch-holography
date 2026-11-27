"""Send authenticated words on both eligible routes before learning the fault."""

import argparse
import gzip
from itertools import combinations
import json

from ..certificates import Rules, Signers
from ..transport import encode
from .measure import BYTE_CAP, DECODE, HERE, clock_trace
from .routes import Mount


def receive(raw, source):
    bits = [DECODE.get(tuple(row)) for row in raw]
    if None in bits:
        return None
    table = bytes(sum(bits[(x >> i) & 1] << i for i in range(8)) for x in range(256))
    try:
        return json.loads(source.translate(table))
    except (UnicodeError, ValueError):
        return None


def scheduled_copies(mount, sender, receiver, disabled, cert, rules):
    packets = []
    source = encode(cert)
    for triple in combinations(range(4), 3):
        if sender not in triple or receiver not in triple:
            continue
        col = mount.columns(triple, sender, receiver, disabled=disabled)
        # A silent destination gives no pointer or acknowledgement.
        raw = None if receiver in disabled else col["raw"]
        result = None if raw is None else receive(raw, source)
        admitted = rules.cert(result, "DC")
        packets.append(
            {
                "via": triple,
                "contact_path": col["path"],
                "raw": raw,
                "bytes": len(source),
                "slots": 8 * len(source) * col["slots"],
                "admitted_value": result[2] if admitted else None,
            }
        )
    return packets


def acquire(source):
    levels = {}
    for k, board in source["levels"].items():
        mount = Mount(int(k), board["topology"]["pentagons"][0])
        cases = {}
        for key, c in board["cases"].items():
            fault = int(key)
            rules = Rules(Signers(c["signature_bank"]), 3)
            ref = next(e["packet"] for e in c["traffic"] if e["phase"] == "decision")
            cert = c["packets"][ref]
            deliveries = {}
            for a in range(4):
                if a == fault:
                    continue
                for b in range(4):
                    if a == b:
                        continue
                    copies = scheduled_copies(mount, a, b, (fault,), cert, rules)
                    accepted = [
                        p["admitted_value"]
                        for p in copies
                        if p["admitted_value"] is not None
                    ]
                    deliveries[f"{a}/{b}"] = {
                        "copies": copies,
                        "first_value": next(iter(accepted), None),
                    }
            two = {}
            for s, r in ((0, 2), (2, 0), (1, 3), (3, 1)):
                both = tuple(i for i in range(4) if i not in (s, r))
                two[f"{s}/{r}"] = scheduled_copies(mount, s, r, both, cert, rules)
            # The prefix precedes an otherwise identical first bit. Its two
            # even-length arms deliberately end with identical clock parity.
            clock = board["clock"]
            longest = max(
                p["slots"] // (8 * p["bytes"])
                for trial in deliveries.values()
                for p in trial["copies"]
            )
            deadline = 6 * 2 * 8 * BYTE_CAP * longest
            prefixes = {
                str(n): clock_trace(clock["tick"], n) for n in (2, deadline + 2)
            }
            cases[key] = {
                "decision_ref": ref,
                "deliveries": deliveries,
                "both_interiors_disabled": two,
                "phase_bound": deadline,
                "silent_prefix_clock": prefixes,
            }
        levels[k] = cases
    return {"wire_byte_cap": BYTE_CAP, "levels": levels}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = acquire(
        json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))
    )
    raw = encode(result)
    f = HERE / "duplication.json.gz"
    if args.write:
        f.write_bytes(gzip.compress(raw, compresslevel=7, mtime=0))
    elif raw != gzip.decompress(f.read_bytes()):
        raise ValueError("Dual route records differ from replay")
