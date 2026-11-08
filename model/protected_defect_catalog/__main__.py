"""Replay or export the source-specific catalogue."""
import argparse
from pathlib import Path
from tools import strict_json
from . import packet

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--check', nargs='?', const=str(packet.HERE/'catalog.json'))
p.add_argument('--write', action='store_true')
a = p.parse_args()
if a.check and a.write:
    p.error('choose --check or --write')
if a.check:
    value = strict_json.load(Path(a.check))
    packet.validate(value)
    if Path(a.check).resolve() == packet.HERE/'catalog.json':
        if (packet.HERE/'CATALOG.md').read_text() != packet.render(value):
            raise ValueError('human catalogue differs from source export')
    print('Protected-sector catalogue replay passed; physical identifications:0')
else:
    value = packet.build()
    if a.write:
        (packet.HERE/'catalog.json').write_text(packet.canonical(value))
        (packet.HERE/'CATALOG.md').write_text(packet.render(value))
    else:
        print(packet.canonical(value), end='')
