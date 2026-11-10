"""Reproducible dark-mode explainer: from a protected twist to a particle."""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
from urllib.parse import unquote
import xml.etree.ElementTree as ET
import zlib

from . import families, field_panels, joint_panels, repair_panels, scenes_a, scenes_b
from . import current_retained_source_equivalence as retained
from .canvas import AMBER, BG, Canvas, TEAL, VIOLET, n
from .parts import SX, W

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SVG = ROOT / 'docs/research/OPH_PARTICLE_ATLAS_V3.svg'
PNG = SVG.with_suffix('.png')
AUDIT = HERE / 'audit.json'
OUTPUTS = ROOT / 'model/particle_observer_outputs/results.json'
MASSES = ROOT / 'model/particle_candidate_atlas/reference_values.json'
DEPENDENCIES = ROOT / 'model/particle_joint_reconstruction/dependencies.json'
MOMENTUM = ROOT / 'model/particle_primary_momentum/PROOF.md'
CANDIDATES = ROOT / 'model/standard_model_source_map/candidates/correspondence.json'
IDENTITY_TESTS = ROOT / 'model/particle_candidate_confirmation/results.json'
GAUGE = ROOT / 'model/particle_primary_reduced_sector/results.json'
SOURCES = ('__init__.py', 'canvas.py', 'parts.py', 'scenes_a.py', 'scenes_b.py', 'families.py', 'repair_panels.py', 'joint_panels.py', 'field_panels.py', 'build.py', 'retained_source_equivalence.py', 'current_retained_source_equivalence.py')

CSS = ('@keyframes pulse{0%,100%{opacity:.5}50%{opacity:1}}'
       '@keyframes hop{0%{opacity:0}1.5%,12.5%{opacity:1}14.5%,100%{opacity:0}}'
       '@keyframes glide{0%{transform:translate(var(--dx),var(--dy))}86%,100%{transform:translate(0,0)}}'
       '.pulse{animation:pulse 3.6s ease-in-out infinite}'
       '.hop{animation:hop 5.6s linear infinite}'
       '.glide{animation:glide 5.6s linear infinite}'
       '@media (prefers-reduced-motion:reduce){.pulse,.hop,.glide{animation:none}}')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def receipt(path):
    """Read a frozen receipt and check its complete producer custody."""
    if path.is_relative_to(ROOT) and str(path.relative_to(ROOT)) in retained.RETAINED_NAMES:
        name = str(path.relative_to(ROOT))
        expected = retained.load_authority(ROOT)['retained_receipts'][name]['sha256']
        if sha(path) != expected:
            raise ValueError('stale source custody: ' + name)
        return retained.read_receipt(path, root=ROOT)
    result = json.loads(path.read_text())
    pins = result['source_sha256']
    if not pins:
        raise ValueError('missing source custody: ' + str(path))
    for relative, digest in pins.items():
        source = (ROOT / relative).resolve()
        if not source.is_relative_to(ROOT) or sha(source) != digest:
            raise ValueError('stale source custody: ' + relative)
    return result


def load():
    """Certified inputs: carrier geometry, the recorded defect run, reference masses."""
    outputs = json.loads(OUTPUTS.read_text())
    run = outputs['holonomy_transport']
    if run['physical_particle_assignments'] or run['measured_particle_inputs']:
        raise ValueError('the holonomy run changed its physical scope')
    carrier = run['carrier']
    data = {
        'vertices': carrier['vertices'],
        'faces': [tuple(f) for f in carrier['faces']],
        'edges': [tuple(sorted(e)) for e in carrier['edges']],
        'path': run['graphical_path']['stages'],
        'attachment': outputs['defect_field_attachment'],
    }
    if (len(data['vertices']), len(data['edges']), len(data['faces'])) != (12, 30, 20):
        raise ValueError('carrier is not the twelve-port icosahedral patch')
    if any(s['total_energy_in_J'] != 2 for s in data['path']):
        raise ValueError('energy ledger of the recorded run is not balanced')
    refs = json.loads(MASSES.read_text())
    if refs['reference_data_are_source_predictions']:
        raise ValueError('reference masses changed their role')
    data['masses'] = {row['target']: row for row in refs['masses']}
    deps = json.loads(DEPENDENCIES.read_text())
    data['nodes'] = [{'id': node['id'], 'status': node['status']} for node in deps['nodes']]
    candidates = receipt(CANDIDATES)
    if (candidates['physical_particle_map'] or candidates['mass_predictions_GeV']
            or candidates['fit_parameters_used']):
        raise ValueError('candidate map changed its uncalibrated hypothesis scope')
    gauge = receipt(GAUGE)['selected_member_bound']['result']
    if (gauge['physical_mass_table_GeV'] or gauge['ambient_primary_vacuum_gap_proved']
            or not gauge['finite_symmetry_resolved_gauge_transition_atom_proved']):
        raise ValueError('gauge transition changed its finite-sector scope')
    data['candidates'], data['gauge'] = candidates, gauge
    identity = receipt(IDENTITY_TESTS)
    retained.validate_receipt_value(str(IDENTITY_TESTS.relative_to(ROOT)), identity, root=ROOT)
    data['identity_tests'] = identity
    data['repair'] = repair_panels.load()
    data['protection'] = repair_panels.load_protection()
    data['agreement'] = repair_panels.load_agreement()
    data['joint'] = joint_panels.load()
    data['field_records'] = field_panels.load()
    data['inputs'] = {str(p.relative_to(ROOT)): sha(p)
                      for p in (OUTPUTS, MASSES, DEPENDENCIES, MOMENTUM, CANDIDATES, GAUGE, IDENTITY_TESTS)}
    data['inputs'].update(repair_panels.input_pins(data['repair']))
    data['inputs'].update(repair_panels.protection_pins(data['protection']))
    data['inputs'].update(repair_panels.agreement_pins(data['agreement']))
    data['inputs'].update(joint_panels.input_pins(data['joint']))
    data['inputs'].update(data['field_records']['pins'])
    data['retained_source_navigation'] = retained.navigation(root=ROOT)
    data['inputs'].update(retained.authenticate_pins(identity['source_sha256'], root=ROOT))
    return data


def render(data):
    c = Canvas()
    c.nodes = []
    y = scenes_a.hero(c, data)
    y = repair_panels.draw(c, data, y + 70)
    y = joint_panels.draw(c, data, y + 70)
    y = field_panels.draw(c, data, y + 70)
    for step in (scenes_a.step1, scenes_a.step2, scenes_a.step3, scenes_a.step4,
                 scenes_b.step5, scenes_b.step67, scenes_b.step8, scenes_b.step9,
                 scenes_b.dictionary, scenes_b.roster, families.catalogue, scenes_b.evidence, families.responses,
                 families.confirmation, families.chronology, scenes_b.standing, scenes_b.footer):
        y = step(c, data, y + 70)
    height = round(y + 56)
    if c.overflow:
        raise ValueError('text overflow: ' + '; '.join(f'{d} {s[:70]!r}' for d, s in c.overflow))

    b = Canvas()
    b.rect(0, 0, W, height, BG)
    b.define('dots', f'<pattern id="dots" width="32" height="32" patternUnits="userSpaceOnUse">'
                     f'<circle cx="2" cy="2" r="1.1" fill="#19233a"/></pattern>')
    b.rect(0, 0, W, height, 'url(#dots)')
    b.glow(250, 330, 620, TEAL, 0.16)
    b.glow(1400, 380, 640, VIOLET, 0.2)
    b.glow(820, 620, 420, AMBER, 0.07)
    # the thread that joins the numbered steps
    ys = [node[0] for node in c.nodes]
    b.define('thread', f'<linearGradient id="thread" gradientUnits="userSpaceOnUse" x1="0" x2="0" '
                       f'y1="{n(ys[0])}" y2="{n(ys[-1])}"><stop offset="0" stop-color="{TEAL}"/>'
                       f'<stop offset=".45" stop-color="{AMBER}"/><stop offset="1" stop-color="{VIOLET}"/>'
                       f'</linearGradient>')
    b.rect(SX - 1.5, ys[0], 3, ys[-1] - ys[0], 'url(#thread)', op=0.55)
    for yy, num, col in c.nodes:
        b.glow(SX, yy, 44, col, 0.35)
        b.circle(SX, yy, 23, BG, col, 2.4)
        b.text(SX, yy + 7.5, num, 21, col, 800, 'middle')

    defs = ''.join(b.defs + c.defs)
    title = 'OPH particle atlas V3: from a twist to a particle'
    desc = ('A dark-mode illustrated explainer in nine numbered steps. It follows a topological '
            'defect through finite source examples and proposed physical interpretations. '
            'Each step distinguishes computed records from the proof required to identify '
            'physical position, momentum and mass. A receipt-bound panel links the candidate '
            'map, isolated gauge transition, static loop fusion and executed contact control; '
            'A receipt-bound coherent seam-repair panel connects its declared reversible '
            'memory, stationary reference, apparatus and static gauge-flux response. '
            'A separate joint construction couples movable boundary flux to the wave, '
            'with reciprocal detectors, invariant interference readout and its own reference. '
            'No source-derived GeV masses or identified physical species are asserted.')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{height}" '
            f'viewBox="0 0 {W} {height}" role="img" aria-labelledby="title desc">'
            f'<title id="title">{title}</title><desc id="desc">{desc}</desc>'
            f'<style>{CSS}</style><defs>{defs}</defs>'
            + ''.join(b.items) + ''.join(c.items) + '</svg>\n')


def check_links(svg):
    root = ET.fromstring(svg)
    for node in root.iter('{http://www.w3.org/2000/svg}a'):
        href = node.get('href', '')
        if not href:
            raise ValueError('empty graphic link')
        if href.startswith('https://'):
            continue
        path = (SVG.parent / unquote(href.split('#')[0])).resolve()
        if not path.is_relative_to(ROOT) or not path.exists():
            raise ValueError('broken graphic link: ' + href)
    return root


def build():
    data = load()
    svg = render(data)
    root = check_links(svg)
    audit = {
        'schema': 'oph.particle_atlas_explainer.v3',
        'input_sha256': data['inputs'],
        'renderer_sha256': {name: sha(HERE / name) for name in SOURCES},
        'source_physical_identifications': [],
        'reference_masses_are_predictions': False,
        'candidate_mass_predictions_GeV': data['candidates']['mass_predictions_GeV'],
        'candidate_target_count': len(data['candidates']['target_rows']),
        'candidate_identity_test_count': len(data['identity_tests']['target_rows']),
        'retained_source_navigation': data['retained_source_navigation'],
        'protected_catalogue_count': len(data['candidates']['protected_catalogue_rows']),
        'evidence_domains': {
            'isolated_gauge_transition': 'rare empty-matter reducing sector, source units',
            'static_loop_fusion': 'same primary source, no moving collision',
            'executed_contact_echo': 'separate supplied contact action and preparation',
            'named_roster': 'measured comparison inputs and candidate assignments',
            'coherent_seam_repair': 'declared reversible memory, same-walk reference and optional apparatus',
            'coherent_static_flux': 'coupled finite gauge walk, conserved internal charge and return response',
            'joint_defect_wave': 'fixed-boundary joint unitary, invariant reciprocal readouts and same-clock reference',
        },
        'coherent_repair_scope': data['repair']['scope'],
        'boundary_flux_scope': data['protection']['scope'],
        'joint_defect_wave_scope': data['joint']['scope'],
        'dependency_status': {node['id']: node['status'] for node in data['nodes']},
        'layout': {'width': root.get('width'), 'height': root.get('height')},
        'svg_sha256': hashlib.sha256(svg.encode()).hexdigest(),
    }
    return svg, json.dumps(audit, ensure_ascii=False, indent=2, sort_keys=True) + '\n'


def png_dimensions(svg):
    """Keep the vector poster full size; fit its raster preview below Cairo's limit."""
    root = ET.fromstring(svg)
    width, height = int(root.get('width')), int(root.get('height'))
    longest = max(width, height)
    if longest <= 32000:
        return width, height
    raster_width = max(1, width * 32000 // longest)
    return raster_width, max(1, (height * raster_width + width - 1) // width)


def png_bytes(svg):
    """Rasterize the complete SVG and bind its exact bytes in a PNG text chunk."""
    width, height = png_dimensions(svg)
    raw = subprocess.run(['rsvg-convert', '--format=png', '--keep-aspect-ratio',
                          '--width', str(width)], input=svg.encode(),
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True).stdout
    if raw[:8] != b'\x89PNG\r\n\x1a\n':
        raise ValueError('renderer did not produce PNG')
    payload = b'Source SVG SHA256\0' + hashlib.sha256(svg.encode()).hexdigest().encode()
    chunk = b'tEXt' + payload
    encoded = struct.pack('>I', len(payload)) + chunk + struct.pack('>I', zlib.crc32(chunk))
    return raw[:33] + encoded + raw[33:]


def slices(svg, out, step=1000):
    """Screen-sized viewBox slices for visual review."""
    root = ET.fromstring(svg)
    height = int(root.get('height'))
    out.mkdir(parents=True, exist_ok=True)
    head = svg.split('>', 1)[0]
    for i, top in enumerate(range(0, height, step)):
        new = (head.replace(f'height="{height}"', f'height="{step}"')
               .replace(f'viewBox="0 0 {W} {height}"', f'viewBox="0 {top} {W} {step}"'))
        (out / f'slice_{i:02d}.svg').write_text(new + '>' + svg.split('>', 1)[1])
    return height


def main():
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write', action='store_true')
    g.add_argument('--check', action='store_true')
    g.add_argument('--slices', type=Path, help='write review slices to this directory')
    a = p.parse_args()
    svg, audit = build()
    if a.slices:
        print('height', slices(svg, a.slices))
        return
    for path, text in ((SVG, svg), (AUDIT, audit)):
        if a.write:
            path.write_text(text)
        elif not path.is_file() or path.read_text() != text:
            raise ValueError('stale explainer atlas: ' + str(path))
    png = png_bytes(svg)
    if a.write:
        PNG.write_bytes(png)
    elif not PNG.is_file() or PNG.read_bytes() != png:
        raise ValueError('stale explainer PNG: ' + str(PNG))
    print('PASS: explainer SVG, PNG and input audit')


if __name__ == '__main__':
    main()
