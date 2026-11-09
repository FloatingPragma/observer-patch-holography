import argparse
from pathlib import Path
from tools import strict_json
from . import packet

parser = argparse.ArgumentParser(description='Verify configured charged-field and calibrated mass readouts.')
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write', action='store_true')
mode.add_argument('--check', nargs='?', const=str(packet.OUTPUT), type=Path)
args = parser.parse_args()
value = packet.build()
if args.write:
    packet.OUTPUT.write_text(packet.canonical(value))
    packet.MARKDOWN.write_text(packet.render(value))
    packet.CANDIDATES.write_text(packet.canonical(value['certificates']['candidate_particle_map']))
else:
    packet.validate(strict_json.load(args.check))
    if packet.MARKDOWN.read_text() != packet.render(value):
        raise ValueError('charged-field rendered results differ')
    if packet.CANDIDATES.read_text() != packet.canonical(value['certificates']['candidate_particle_map']):
        raise ValueError('charged-field candidate map differs')
print('Full-CAR and free-field certificates; calibrated lepton comparison; no physical particle assignments')
