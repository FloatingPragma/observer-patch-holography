"""Receipt-bound finite defect-wave interaction and its complete readouts."""
import hashlib
import json
from pathlib import Path

from .canvas import AMBER, BLUE, EDGE, MUTED, SOFT, TEAL, TEXT, VIOLET
from .parts import CW, X0, card

ROOT = Path(__file__).resolve().parents[2]
RECEIPT = ROOT / 'model/coherent_defect_wave/results.json'


def load():
    value = json.loads(RECEIPT.read_text())
    if value['schema'] != 'oph.coherent_defect_wave.v1':
        raise ValueError('joint defect-wave schema changed')
    if not value['input_sha256']:
        raise ValueError('joint source custody missing')
    for name, digest in value['input_sha256'].items():
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('stale joint defect-wave input: ' + name)
    for key in ('continuous_dynamical_gauge_theory', 'complete_boundary_gauss_constraints',
                'complete_joint_effect_tomography', 'interacting_3plus1_continuum',
                'physical_particle_species', 'physical_energy_or_GeV_mass'):
        if value['scope'][key]:
            raise ValueError('joint physical scope promotion: ' + key)
    if value['source']['dimension'] != 72000:
        raise ValueError('joint finite source domain changed')
    for key in ('wave_to_defect', 'defect_to_wave'):
        if not value['interaction'][key]['actual_full_tick_executed']:
            raise ValueError('joint interaction is not executed on the full tick')
    coherent = value['interaction']['constrained_coherence']
    if not (coherent['return_projector_commutes_with_Gamma'] and
            coherent['all_recorded_branches_preserve_positive_constraint']):
        raise ValueError('return apparatus leaves the admitted gauge sector')
    return value


def input_pins(value):
    return {**value['input_sha256'], str(RECEIPT.relative_to(ROOT)):
            hashlib.sha256(RECEIPT.read_bytes()).hexdigest()}


def draw(c, data, y):
    value = data['joint']
    interaction, reference, gauss = value['interaction'], value['reference'], value['gauss']
    c.caps(X0, y + 10, 'A joint finite source: one protected region and one wave', TEAL, 14)
    c.text(X0, y + 62, 'Defects and waves change one another', 40, TEXT, 700, maxw=CW)
    y = c.para(X0, y + 105, CW,
        'The joint program retains the link memory in the wave’s quantum state. '
        'The links change the wave step; the wave’s port controls permitted defect moves. '
        'Both happen in one fixed reversible tick on ' + f"{value['source']['dimension']:,}" +
        ' basis states. No moved-or-held outcome is measured during idle propagation.',
        21, MUTED, 1.42) + 22

    h = 135
    card(c, X0, y, CW, h, stroke=TEAL)
    w = (CW - 120) / 3
    for i, (title, body, color) in enumerate((
            ('Retained link word g', 'Ten allowed words encode five defect positions.', AMBER),
            ('Wave step Wg', 'Coherent seam transport reads the actual link signs.', BLUE),
            ('Controlled link moves', 'The wave port selects local reversible repairs.', TEAL))):
        x = X0 + 26 + i * (w + 34)
        c.text(x, y + 43, title, 25, color, 700, maxw=w)
        c.para(x, y + 77, w - 14, body, 18, SOFT, 1.35)
        if i < 2:
            c.arrow(x + w - 5, y + 44, x + w + 20, y + 44, TEAL, 2.4, 9)
    y += h + 22

    h, gap, w = 364, 22, (CW - 44) / 3
    response = interaction['defect_to_wave']['yes_and_complement_probabilities']
    coherent = interaction['constrained_coherence']['coherent_return_and_complement'][0]
    recorded = interaction['constrained_coherence']['after_complete_face_readout_return_and_complement'][0]
    field = interaction['wave_to_defect']
    face_id = str(field['deterministic_final_face_indices'][0])
    face_probabilities = [row[face_id] for row in field['complete_gauge_invariant_face_distributions']]
    for i, (tag, title, body, number, color) in enumerate((
            ('Defect → wave', 'A different defect changes a click',
             'Keep the same two-port wave preparation. Changing the allowed link word '
             'changes the region’s defect position and the output of the same wave detector. '
             'Both detector complements are retained.',
             'Click probability: ' + response[0][0] + ' versus ' + response[1][0], BLUE),
            ('Wave → defect', 'A different wave changes the defect',
             'Prepare the inverse wave step of a port-5 or port-6 state. One complete tick '
             'then moves the defect to different faces. A detector of one regional face '
             'separates these two preparations.',
             'Face detector: ' + face_probabilities[0] + ' versus ' + face_probabilities[1], TEAL),
            ('Gauge-invariant interference', 'Recording a branch changes interference',
             'A wave superposition correlates with two defect positions. The apparatus uses '
             'the retained link frame to compare amplitudes. Undoing the controls restores '
             'the return. Reading defect position first halves it; both outcomes are retained.',
             'Return: ' + coherent + ' unread; ' + recorded + ' after readout', VIOLET))):
        x = X0 + i * (w + gap)
        card(c, x, y, w, h)
        c.caps(x + 24, y + 35, tag, color, 13)
        yy = c.para(x + 24, y + 76, w - 48, title, 24, TEXT, 1.28, 700)
        c.para(x + 24, yy + 10, w - 48, body, 18, SOFT, 1.36)
        c.text(x + 24, y + h - 24, number, 17, color, 700, maxw=w - 48)
    y += h + 22

    h = 204
    card(c, X0, y, CW, h, stroke=BLUE)
    c.text(X0 + 28, y + 43, 'A reference and a gauge constraint for this same clock', 28, TEXT, 700)
    c.para(X0 + 28, y + 83, CW - 56,
        'The complete preparation filter uses this joint U, giving success probability ' +
        reference['success_probability_from_uniform_basis_ensemble'] +
        ' from its uniform basis ensemble. The admitted interior gauge constraint leaves ' +
        f"{gauss['physical_dimension']:,}" + ' basis states; imposing it first gives joint success ' +
        gauss['joint_success_from_uniform_ambient_input'] +
        '. Every failure is retained. The conditioned reference is stationary on the '
        'reduced joint source; the full record archive is not claimed stationary.', 19, SOFT, 1.4)
    c.link(X0 + 28, y + 180, 'Same-U preparation and complete outcome proof',
           '../../model/coherent_defect_wave/PREPARATION.md', 16)
    c.link(X0 + 805, y + 180, 'Interior gauge constraint', '../../model/coherent_defect_wave/GAUSS.md', 16)
    y += h + 22
    c.para(X0 + 24, y + 8, CW - 48,
        '{a|What this establishes.} A finite, reciprocal interaction with retained boundary flux, '
        'coherent propagation and complete specified detectors. Its fixed boundary, one-defect '
        'sector, quantum registers and control order are declared. Physical particle species, '
        'spin, rest masses and an interacting relativistic continuum are not identified.',
        18, MUTED, 1.42)
    y += 112
    c.link(X0 + 24, y, 'Joint evolution and reciprocal detector proofs',
           '../../model/coherent_defect_wave/PROOF.md', 17)
    c.link(X0 + 805, y, 'Exact joint receipt', '../../model/coherent_defect_wave/results.json', 17)
    return y + 22
