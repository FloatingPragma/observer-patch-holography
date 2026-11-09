import argparse
from pathlib import Path
from tools import strict_json
from . import packet

parser = argparse.ArgumentParser(description='Export exact source correspondences and physical admission.')
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--write', action='store_true')
mode.add_argument('--check', nargs='?', const=str(packet.OUTPUT), type=Path)
args = parser.parse_args()
value = packet.build()
if args.write:
    packet.OUTPUT.write_text(packet.canonical(value))
    packet.MARKDOWN.write_text(packet.render(value))
    packet.MASS_MARKDOWN.write_text(packet.render_masses(value))
else:
    packet.validate(strict_json.load(args.check))
    if packet.MARKDOWN.read_text() != packet.render(value):
        raise ValueError('rendered particle map differs')
    if packet.MASS_MARKDOWN.read_text() != packet.render_masses(value):
        raise ValueError('rendered source masses differ')
print('18 target channels; 15 exact source components; 0 physical particle identifications')
