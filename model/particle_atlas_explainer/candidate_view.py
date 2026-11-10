"""Current candidate SVG: frozen numerical atlas plus authenticated source families.

The numerical comparison producer remains immutable. This leaf owns the public
candidate SVG and verifies the retained base before adding confirmation evidence.
"""
import argparse
from html import escape
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from model.particle_candidate_atlas import builder as base
from . import build, families, field_panels, joint_panels, repair_panels
from . import current_retained_source_equivalence as retained
from .canvas import BG, Canvas, MUTED, TEAL, TEXT
from .parts import CW, W, X0

ROOT = build.ROOT
HERE = Path(__file__).resolve().parent
OUTPUT = ROOT/'docs/research/OPH_PARTICLE_CANDIDATE_ATLAS.svg'
AUDIT = HERE/'candidate_view_audit.json'
SOURCES = ('candidate_view.py', 'families.py', 'repair_panels.py', 'joint_panels.py', 'field_panels.py', 'canvas.py', 'parts.py', 'test_candidate_view.py', 'retained_source_equivalence.py', 'current_retained_source_equivalence.py')


def load():
    comparison, original_svg = base.build()
    confirmation = build.receipt(build.IDENTITY_TESTS)
    retained.validate_receipt_value(str(build.IDENTITY_TESTS.relative_to(ROOT)), confirmation, root=ROOT)
    if confirmation['physical_particle_map'] or confirmation['mass_predictions_GeV']:
        raise ValueError('unproved physical candidate promotion')
    return comparison, original_svg, confirmation


def render(original_svg, confirmation, repair=None, protection=None, agreement=None, joint=None, field=None):
    root = ET.fromstring(original_svg)
    original_height, width = int(root.get('height')), int(root.get('width'))
    if width != 4000:
        raise ValueError('retained numerical atlas width changed')
    c = Canvas()
    c.caps(X0, 54, 'Finite-source candidate evidence', TEAL, 14)
    c.text(X0, 106, 'The source catalogue and the tests for a particle', 40, TEXT, 700)
    y = c.para(X0, 151, CW,
        'The numerical comparisons above retain their declared action, calibration and '
        'mass conventions. The following source catalogue and complete detector results '
        'use the candidate confirmation receipt. Neither table assigns a physical species '
        'or a source-derived GeV mass.', 21, MUTED, 1.42)+72
    data = {'identity_tests': confirmation, 'repair': repair or repair_panels.load(),
            'protection': protection or repair_panels.load_protection(),
            'agreement': agreement or repair_panels.load_agreement(),
            'joint': joint or joint_panels.load(),
            'field_records': field or field_panels.load()}
    y = repair_panels.draw(c, data, y) + 70
    y = joint_panels.draw(c, data, y) + 70
    y = field_panels.draw(c, data, y) + 70
    for draw in (families.responses, families.catalogue, families.confirmation, families.chronology):
        y = draw(c, data, y)+70
    if c.overflow:
        raise ValueError('candidate extension text overflow: '+str(c.overflow))
    scale = width/W
    total_height = original_height+round((y+30)*scale)
    start = original_svg.index('<svg ')
    end = original_svg.index('>', start)+1
    tag = original_svg[start:end]
    tag = re.sub(r'height="[^"]+"', f'height="{total_height}"', tag, count=1)
    tag = re.sub(r'viewBox="[^"]+"', f'viewBox="0 0 {width} {total_height}"', tag, count=1)
    body = original_svg[:start]+tag+original_svg[end:]
    marker = '<metadata id="candidate-confirmation-extension">'+escape(json.dumps({
        'schema': 'oph.particle_candidate_current_view.v1',
        'confirmation_result_sha256': confirmation['content_sha256'],
        'retained_numerical_svg_sha256': hashlib.sha256(original_svg.encode()).hexdigest(),
        'physical_particle_map': [], 'mass_predictions_GeV': []}, sort_keys=True))+'</metadata>'
    panel = (f'<g id="candidate-confirmation-panel" transform="translate(0 {original_height}) scale({scale})">'
             f'<rect width="{W}" height="{y+30}" fill="{BG}"/>'+''.join(c.items)+'</g>')
    return body.rsplit('</svg>', 1)[0]+marker+panel+'</svg>\n'


def produce():
    comparison, old_svg, confirmation = load()
    repair = repair_panels.load()
    protection = repair_panels.load_protection()
    agreement = repair_panels.load_agreement()
    joint = joint_panels.load()
    field = field_panels.load()
    svg = render(old_svg, confirmation, repair, protection, agreement, joint, field)
    build.check_links(svg)
    pins = {**comparison['input_sha256'], **retained.authenticate_pins(confirmation['source_sha256'], root=ROOT),
            str(build.IDENTITY_TESTS.relative_to(ROOT)): build.sha(build.IDENTITY_TESTS)}
    pins.update(repair_panels.input_pins(repair))
    pins.update(repair_panels.protection_pins(protection))
    pins.update(repair_panels.agreement_pins(agreement))
    pins.update(joint_panels.input_pins(joint))
    pins.update(field['pins'])
    audit = {'schema': 'oph.particle_candidate_current_view.v1',
             'input_sha256': pins,
             'renderer_sha256': {name: build.sha(HERE/name) for name in SOURCES},
             'retained_numerical_svg_sha256': hashlib.sha256(old_svg.encode()).hexdigest(),
             'retained_numerical_audit_sha256': hashlib.sha256(base.canonical(comparison).encode()).hexdigest(),
             'confirmation_result_sha256': confirmation['content_sha256'],
             'retained_source_navigation': retained.navigation(root=ROOT),
             'protected_catalogue_count': sum(len(g['rows']) for g in confirmation['labelled_catalogue']['protected_families']),
             'source_excitation_count': sum(len(g['rows']) for g in confirmation['labelled_catalogue']['excitation_families']),
             'physical_particle_map': [], 'mass_predictions_GeV': [],
             'coherent_repair_scope': repair['scope'],
             'boundary_flux_scope': protection['scope'],
             'joint_defect_wave_scope': joint['scope'],
             'svg_sha256': hashlib.sha256(svg.encode()).hexdigest()}
    return svg, json.dumps(audit, indent=2, sort_keys=True)+'\n'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write', action='store_true')
    g.add_argument('--check', action='store_true')
    args = p.parse_args()
    svg, audit = produce()
    for path, value in ((OUTPUT, svg), (AUDIT, audit)):
        if args.write:
            path.write_text(value)
        elif not path.is_file() or path.read_text() != value:
            raise ValueError('stale current candidate view: '+str(path))
    print('PASS: current candidate SVG, retained numerical base and source-family receipt')


if __name__ == '__main__':
    main()
