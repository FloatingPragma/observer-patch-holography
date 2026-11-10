"""Receipt-bound diagrams for coherent seam correction and its finite readout."""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from .canvas import AMBER, BLUE, CARD2, EDGE, GREEN, MUTED, RED, SOFT, TEAL, TEXT, VIOLET
from .parts import CW, X0, card

ROOT = Path(__file__).resolve().parents[2]
RECEIPT = ROOT / 'model/coherent_seam_repair/results.json'
PROTECTION = ROOT / 'model/coherent_seam_repair/protection_results.json'
AGREEMENT = ROOT / 'model/coherent_seam_repair/agreement_results.json'
PROOF = '../../model/coherent_seam_repair/PROOF.md'


def load():
    value = json.loads(RECEIPT.read_text())
    if value['schema'] != 'oph.coherent_seam_repair.v1':
        raise ValueError('coherent repair receipt schema changed')
    pins = value['input_sha256']
    if not pins:
        raise ValueError('coherent repair input custody missing')
    for name, digest in pins.items():
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('stale coherent repair input: ' + name)
    for key in ('physical_particle_identification', 'continuum_limit', 'GeV_mass_prediction'):
        if value['scope'][key]:
            raise ValueError('coherent repair physical scope promotion: ' + key)
    if value['walk']['full_walk_dimension'] != 1440 or value['gauged_walk']['coupled_walk_dimension'] != 7200:
        raise ValueError('coherent repair source dimensions changed')
    if not 0 < Fraction(value['reference']['success_probability']) < 1:
        raise ValueError('complete preparation requires success and failure')
    if Fraction(value['reference']['eigenvalue_lower']) <= 0:
        raise ValueError('faithful reference certificate missing')
    return value


def input_pins(value):
    return {**value['input_sha256'], str(RECEIPT.relative_to(ROOT)):
            hashlib.sha256(RECEIPT.read_bytes()).hexdigest()}


def load_protection():
    value = json.loads(PROTECTION.read_text())
    if value['schema'] != 'oph.coherent_seam_repair.boundary_flux.v1':
        raise ValueError('boundary-flux schema changed')
    for name, digest in value['input_sha256'].items():
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('stale boundary-flux input: ' + name)
    if (any(value['scope'].values()) or value['protection']['arbitrary_link_updates_protect_initial_pair']
            or value['protection']['coupled_static_wave_register_is_this_dynamic_link_program']):
        raise ValueError('boundary-flux control changed its protection or physical scope')
    return value


def protection_pins(value):
    return {**value['input_sha256'], str(PROTECTION.relative_to(ROOT)):
            hashlib.sha256(PROTECTION.read_bytes()).hexdigest()}


def load_agreement():
    value = json.loads(AGREEMENT.read_text())
    if value['schema'] != 'oph.coherent_seam_agreement.v1':
        raise ValueError('observer-agreement schema changed')
    for name, digest in value['input_sha256'].items():
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('stale observer-agreement input: ' + name)
    for claim in ('A1_complete_response_certified', 'A3_state_selection_derived',
                  'Born_rule_derived_from_classical_repair', 'physical_clock_or_growing_federation_claim'):
        if value[claim]:
            raise ValueError('observer-agreement scope promotion: ' + claim)
    if not value['direct_recharted_link_execution_checked']:
        raise ValueError('observer-agreement execution missing')
    return value


def agreement_pins(value):
    return {**value['input_sha256'], str(AGREEMENT.relative_to(ROOT)):
            hashlib.sha256(AGREEMENT.read_bytes()).hexdigest()}


def _seam(c, x, y, label, opposite, color):
    c.line(x, y, x + 148, y, color, 3)
    for offset, sign in ((0, '+'), (148, '−' if opposite else '+')):
        c.circle(x + offset, y, 24, CARD2, color, 2)
        c.text(x + offset, y + 8, sign, 28, color, 700, 'middle')
    c.text(x + 74, y + 55, label, 18, SOFT, 600, 'middle')


def draw(c, data, y):
    """One horizontal process within the portrait explainer and the candidate leaf."""
    value = data['repair']
    reference = value['reference']
    gauge = value['gauged_walk']
    apparatus = value['apparatus']
    c.caps(X0, y + 10, 'An executable repair and observer construction', TEAL, 14)
    c.text(X0, y + 60, 'Repair can move disagreement into memory', 40, TEXT, 700, maxw=CW)
    y = c.para(X0, y + 102, CW,
        'The reversible seam rule stores the removed {a|exchange disagreement} in a coherent '
        'record. A declared address reflection then makes that information participate in '
        'a finite wave. This program has its own preparation, clock and optional detector '
        'operations. Its source records do not identify a natural particle.', 21, MUTED, 1.42) + 16

    # The local correction is one reversible arrow, not loss of the minus component.
    h = 230
    card(c, X0, y, CW, h)
    _seam(c, X0 + 58, y + 92, '|−⟩ = (|10⟩ − |01⟩)/√2', True, AMBER)
    c.rect(X0 + 266, y + 60, 150, 70, CARD2, EDGE, 1.3, 12)
    c.text(X0 + 341, y + 89, 'coherent record', 15, MUTED, 400, 'middle')
    c.text(X0 + 341, y + 116, '|0⟩', 26, SOFT, 700, 'middle')
    c.arrow(X0 + 450, y + 89, X0 + 574, y + 89, TEAL, 2.8, 11)
    c.arrow(X0 + 574, y + 108, X0 + 450, y + 108, TEAL, 2.8, 11)
    c.text(X0 + 512, y + 148, 'reversible correction', 17, TEAL, 600, 'middle')
    _seam(c, X0 + 626, y + 92, '|+⟩ = (|10⟩ + |01⟩)/√2', False, TEAL)
    c.rect(X0 + 838, y + 60, 150, 70, CARD2, TEAL, 1.3, 12)
    c.text(X0 + 913, y + 89, 'coherent record', 15, MUTED, 400, 'middle')
    c.text(X0 + 913, y + 116, '|1⟩', 26, TEAL, 700, 'middle')
    c.line(X0 + 1022, y + 28, X0 + 1022, y + h - 28, EDGE, 1.2)
    c.para(X0 + 1050, y + 58, 316,
        'The seam agrees. The record holds what changed. Reversing the operation restores '
        'the input; no information is erased.', 19, SOFT, 1.42)
    c.text(X0 + 28, y + h - 25,
        'Exchange parity here labels a coherent seam state. It does not establish fermionic spin or statistics.',
        16.5, MUTED, maxw=CW - 56)
    y += h + 20

    # These two execution choices must not be shown as the same stochastic law.
    gap, w, h = 24, (CW - 24) / 2, 172
    for i, (title, body, color) in enumerate((
        ('Commit a classical outcome at every event',
         'Measuring and replacing the record gives a dissipative pump. The phase-addressed '
         'variant has a faithful finite reference and complex tomography, but its coherence '
         'decays. Repeating it does not produce an undamped wave.', AMBER),
        ('Reuse coherent memory during propagation',
         'The seam correction Q and address reflection R define the unitary step W = QR. '
         'One first-event population agrees with the native rule; later coherent histories '
         'have different probabilities. The reflection is a declared operation.', TEAL))):
        xx = X0 + i * (w + gap)
        card(c, xx, y, w, h)
        c.text(xx + 24, y + 36, title, 21, color, 700, maxw=w - 48)
        c.para(xx + 24, y + 72, w - 48, body, 18, SOFT, 1.36)
    y += h + 26

    # The one-program preparation / evolution / apparatus pipeline.
    gap, w, h = 22, (CW - 44) / 3, 355
    cards = (
        ('1  Prepare', 'A reference made with the same W',
         'One controlled walk creates success and failure outputs. Keep both records; '
         'retry only within the stated memory budget. The successful source reference '
         'is faithful and stationary under W.',
         'Success: ' + reference['success_probability'] + ' per attempt', TEAL),
        ('2  Propagate', 'Coherent records carry the wave',
         'For ' + f"{value['walk']['full_walk_dimension']:,}" + ' source-and-memory states, '
         'the exact native discriminant is I − L/60. Its walk phases satisfy '
         'cos ω = 1 − λ(L)/60. A wave tick counts one W; it is not a second in SI units.',
         'No damping under the declared unitary W', BLUE),
        ('3  Read', 'An apparatus changes the experiment',
         'Optional quarter-turn probes and a register detector separate real and imaginary '
         'coherences. Every outcome is retained. The history effects span '
         f"{apparatus['full_complex_Hermitian_effect_span_dimension']:,}" + ' Hermitian directions. '
         'Readout can disturb the state.',
         'Complete finite tomography; no arbitrary-POVM claim', VIOLET),
    )
    for i, (tag, title, body, foot, color) in enumerate(cards):
        xx = X0 + i * (w + gap)
        card(c, xx, y, w, h)
        c.caps(xx + 24, y + 34, tag, color, 13)
        yy = c.para(xx + 24, y + 72, w - 48, title, 23, TEXT, 1.3, 700)
        c.para(xx + 24, yy + 10, w - 48, body, 18, SOFT, 1.36)
        c.text(xx + 24, y + h - 21, foot, 14, color, 600, maxw=w - 48)
    y += h + 22
    c.para(X0 + 24, y + 12, CW - 48,
        'The stationary reference is a reduced walk-register state. Conditioning on retained '
        'preparation labels changes it. Modular flow is nontrivial on readable effects, but '
        'its generator is not identified with the walk clock.', 17, MUTED, 1.4)
    y += 88

    # Concrete internal holonomy changes an observer-facing return amplitude.
    h = 380
    card(c, X0, y, CW, h, fill='#131a2b', stroke=VIOLET)
    c.caps(X0 + 28, y + 35, 'A loop changes the readout in the same finite gauge walk', VIOLET, 13)
    c.text(X0 + 28, y + 76, 'Carrier rotations → retained lift → transported state → detector',
           26, TEXT, 700, maxw=CW - 56)
    left_w = 450
    yy = c.para(X0 + 28, y + 119, left_w,
        'The controller derives {t|120 lifts} from 60 proper carrier rotations. '
        'A closed two-pulse word produces Z = diag(I₃, −I₂). Its internal sign '
        'survives even when the carrier orientation returns.', 19, SOFT, 1.42)
    c.para(X0 + 28, yy + 10, left_w,
        'The finite holonomy theorem permits this cover. It does not require a '
        'single-valued spin lift of every carrier rotation.', 16.5, MUTED, 1.4)
    c.line(X0 + 510, y + 105, X0 + 510, y + h - 26, EDGE, 1.2)
    c.para(X0 + 538, y + 119, 380,
        'One declared Z link makes {a|two triangular loops} nontrivial. Coupling '
        'the five-level code to seam transport gives ' + f"{gauge['coupled_walk_dimension']:,}" +
        ' states with one unitary clock and their own stationary preparation.',
        19, SOFT, 1.42)
    c.text(X0 + 538, y + h - 34, 'Internal labels +2 and −3 are conserved.', 17, TEAL, 600)
    c.line(X0 + 962, y + 105, X0 + 962, y + h - 26, EDGE, 1.2)
    c.para(X0 + 990, y + 119, CW - 1018,
        'The loop changes the three-step return amplitude. This is an actual '
        'finite detector distinction on the coupled source.', 19, SOFT, 1.42)
    delta = gauge['static_flux_witness']['amplitude_difference']
    c.text(X0 + 990, y + 300, 'Return-amplitude change: ' + delta, 17, VIOLET, 700)
    c.text(X0 + 990, y + 333, 'The five code states are not five species.', 16, AMBER, 600)
    y += h + 20
    if 'agreement' in data:
        value_agreement = data['agreement']
        card(c, X0, y, CW, 182, stroke=BLUE)
        c.text(X0 + 28, y + 42, 'Two observers compare the same detector', 28, TEXT, 700)
        c.para(X0 + 28, y + 81, CW - 56,
            'A change of observer frame transports the state, links and detector together. '
            'The executed detector then agrees: probability ' +
            value_agreement['correct_transported_detector_probability'] +
            '. Keeping the old detector gives ' +
            value_agreement['wrong_untransported_detector_probability'] +
            ' and describes a different experiment. This proves agreement for the accepted '
            'finite chart operations; it does not certify every possible observer response.',
            19, SOFT, 1.4)
        c.link(X0 + 28, y + 160, 'Executed frame comparison and complete outcome receipt',
               '../../model/coherent_seam_repair/agreement_results.json', 16)
        y += 202
    if 'protection' in data:
        y = boundary_panel(c, data['protection'], y, joined='joint' in data) + 20
    card(c, X0, y, CW, 110, fill='#211a17', stroke=AMBER)
    c.para(X0 + 26, y + 36, CW - 52,
        '{a|Physical boundary.} The wave experiment’s link background is supplied and static. Its flux cannot '
        'be removed by a gauge rechart, but an allowed link change can erase it. Conserved '
        'internal charge is not a protected particle or a measured electric charge. '
        'No dynamical gauge backreaction, interacting 3+1-dimensional limit, physical spin '
        'or source-derived GeV mass follows from this static-link construction.', 18, SOFT, 1.4)
    y += 146
    c.link(X0 + 24, y, 'Repair, wave, reference and detector proofs', PROOF, 17)
    c.link(X0 + 695, y, 'Executable observer paper', '../../paper/executable_observer_construction.pdf', 17)
    c.link(X0 + 1120, y, 'Exact receipt', '../../model/coherent_seam_repair/results.json', 17)
    return y + 20


def boundary_panel(c, value, y, joined=False):
    h = 382
    card(c, X0, y, CW, h)
    c.caps(X0 + 28, y + 35, 'A separate finite link-memory control', TEAL, 13)
    cx, cy, radius = X0 + 184, y + 171, 88
    vertices = [(cx + radius * math.cos(-math.pi/2 + 2*math.pi*j/5),
                 cy + radius * math.sin(-math.pi/2 + 2*math.pi*j/5)) for j in range(5)]
    for i in range(5):
        triangle = [(cx, cy), vertices[i], vertices[(i+1) % 5]]
        c.poly(triangle, '#45301a' if i == 0 else CARD2, TEAL, 1.6)
    c.poly(vertices, 'none', AMBER, 3.4)
    c.circle(cx, cy, 6, TEAL)
    fx = (cx + vertices[0][0] + vertices[1][0])/3
    fy = (cy + vertices[0][1] + vertices[1][1])/3
    c.circle(fx, fy, 12, AMBER)
    c.text(fx, fy + 5, '−', 17, CARD2, 700, 'middle')
    c.text(cx, y + h - 39, 'Boundary links held fixed', 17, AMBER, 600, 'middle')
    c.text(cx, y + h - 15, 'Interior links can move the flux', 15, MUTED, 400, 'middle')
    c.line(X0 + 360, y + 68, X0 + 360, y + h - 26, EDGE, 1.2)
    c.text(X0 + 390, y + 85, 'Odd regional flux cannot repair to vacuum', 29, TEXT, 700)
    yy = c.para(X0 + 390, y + 127, CW - 418,
        f"A declared {value['link_memory_bits']}-bit link register stores central signs on the carrier. "
        'The product of face signs inside a region equals its boundary sign. If that boundary '
        'is fixed and odd, no allowed interior repair can erase every defect.', 20, SOFT, 1.42)
    boundary_scope = ('Changing the boundary can erase the pair. This control has its own link-update clock. '
                      'The joint program below couples the defect and wave in one declared evolution.' if joined else
                      'Changing the boundary can erase the pair. This dynamic link register has no proved '
                      'joint evolution or stationary reference with the coherent wave above.')
    yy = c.para(X0 + 390, yy + 10, CW - 418,
        f"The exact local program transports one defect to all {value['single_defect_positions']} region faces: "
        f"{value['involutory_transport_move_checks']} checked moves and {value['retained_no_operation_checks']} retained holds. "
        + boundary_scope, 19, MUTED, 1.42)
    c.link(X0 + 390, y + h - 20, 'Boundary identity, transport program and finite-gap bound',
           '../../model/coherent_seam_repair/PROTECTION.md', 16)
    return y + h
