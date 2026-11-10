"""Same-source charged limits and declared collective field constructions."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from .canvas import AMBER, BLUE, MUTED, SOFT, TEAL, TEXT, VIOLET
from .parts import CW, X0, card

ROOT = Path(__file__).resolve().parents[2]
PATHS = {
    'charge': 'model/coherent_defect_wave/chiral_joint/interacting_charge_localization/results.json',
    'propagation': 'model/interacting_block_propagation/results.json',
    'field': 'model/interacting_block_quadrature/results.json',
    'collision': 'model/interacting_block_scattering/results.json',
}


def validate(data):
    field = data['field']
    if field['construction']['measured_inputs'] or not field['construction']['full_physical_complement_retained']:
        raise ValueError('Field construction inputs changed')
    theorem = field['theorems']
    if (not theorem['massive_free_KG_limit_in_increasing_capacity_subsector']
            or not theorem['three_components_are_internal_scalar_labels']
            or theorem['nontrivial_interacting_relativistic_limit_proved']
            or theorem['Lorentz_vector_or_Proca_identification_proved']):
        raise ValueError('Field interpretation changed')
    for row in data['charge']['full_source_charged_bins'].values():
        law = row['energy_filtered_localization']
        if (not law['every_normalized_state_in_bin_and_every_time']
                or law['finite_bin_complement_or_higher_states_truncated']
                or law['origin_probability_lower_inferred_by_discarding_screened_complement']):
            raise ValueError('Charged localization domain changed')
    collision = data['collision']['interval_certified_five_block_collision']
    if (collision['actual_interacting_ground_input_claim']
            or collision['three_dimensional_continuum_claim']
            or collision['exact_reducing_source_dimension'] != 1024):
        raise ValueError('Finite collision interpretation changed')


def load():
    data, pins = {}, {}
    for name, relative in PATHS.items():
        path = ROOT / relative
        data[name] = json.loads(path.read_text())
        local = data[name]['input_sha256']
        for key, expected in local.items():
            target = (ROOT / key).resolve()
            if not target.is_relative_to(ROOT) or hashlib.sha256(target.read_bytes()).hexdigest() != expected:
                raise ValueError('Field evidence changed: ' + key)
        pins.update(local)
        pins[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    validate(data)
    data['pins'] = pins
    return data


def draw(c, data, y):
    d = data['field_records']
    h = 415
    card(c, X0, y, CW, h, stroke=VIOLET)
    c.text(X0+28, y+43, 'A charged energy record is not automatically a moving particle', 27, TEXT, 700)
    c.para(X0+28, y+84, CW-56,
           'In the same eight-cell source, a static reference keeps the compensating charge. '
           'Moving one unscreened charge away requires excited links. A low-energy detector '
           'therefore sees a strong restriction on that motion.', 19, SOFT, 1.4)
    for i, (key, label) in enumerate((('ec', 'eᶜ removal'), ('L', 'L removal'))):
        bound = Fraction(d['charge']['full_source_charged_bins'][key]['energy_filtered_localization']
                         ['union_remote_single_unscreened_charge_probability_upper'])
        xx = X0+28+i*650
        c.text(xx, y+190, label, 20, VIOLET, 700)
        c.text(xx, y+226, f'remote single-charge probability ≤ {float(bound):.3e}', 22, TEAL, 600)
    c.para(X0+28, y+266, CW-56,
           'This bound holds after low-bin selection at every source time. The complete detector '
           'also keeps the origin and a multi-charge cloud; screening can live in that cloud. '
           'Positive higher-energy arrivals are a different setting. These facts do not prove '
           'that a free electron has been identified.', 19, SOFT, 1.4)
    c.link(X0+28, y+h-20, 'Exact Gauss, spectral-separation and screening proof',
           '../../'+PATHS['charge'].replace('results.json', 'PROOF.md'), 16)
    y += h+35
    h = 640
    card(c, X0, y, CW, h, stroke=TEAL)
    c.text(X0+28, y+43, 'One computed neutral gap → a controlled relativistic field description', 27, TEXT, 700)
    c.para(X0+28, y+85, CW-56,
           'Take complete copies of the interacting eight-cell source. Its native current connects '
           'the ground to three neutral excitations. Declare paired exchange between those copies '
           'and retain every higher source state. The exact finite dynamics has a controlled '
           'large-capacity comparison.', 19, SOFT, 1.4)
    for j in range(3):
        xx = X0+95+j*190
        c.rect(xx, y+198, 125, 108, '#132636', TEAL, 1.5, 10)
        for a in range(2):
            for b in range(2):
                c.circle(xx+38+48*a, y+230+42*b, 10, [TEAL, BLUE, VIOLET, AMBER][2*a+b])
        if j<2:
            c.arrow(xx+132, y+251, xx+181, y+251, TEAL, 2, 7)
    c.text(X0+332, y+339, 'N full source blocks per spatial cell', 19, SOFT, 500, 'middle')
    c.arrow(X0+685, y+252, X0+775, y+252, VIOLET, 3, 9)
    c.text(X0+812, y+220, 'ω² → Δ² + (Δ²/6)|p|²', 29, TEAL, 700)
    c.para(X0+812, y+258, CW-850,
           'A massive Klein–Gordon description. Δ is the actual native gap, '
           'not a number chosen from a particle table.', 19, SOFT, 1.4)
    c.para(X0+28, y+391, CW-56,
           'Ground-state, finite-capacity, evolution and detector errors are bounded together. '
           'A retained finite orientation reference makes the field measurements relational. '
           'The three components are internal scalar labels, not the three polarizations of '
           'a measured vector boson.', 19, SOFT, 1.4)
    c.para(X0+28, y+494, CW-56,
           'The paired exchange, replication and spacing-dependent clock are explicit design '
           'choices. The established refinement becomes free. A separate number-conserving '
           'exchange gives a Schrödinger limit. Neither construction supplies a Standard Model '
           'identity or a GeV calibration.', 19, MUTED, 1.4)
    c.link(X0+28, y+h-20, 'Complete collective field and reference proof',
           '../../'+PATHS['field'].replace('results.json', 'PROOF.md'), 16)
    c.link(X0+780, y+h-20, 'Distinct number-conserving propagation proof',
           '../../'+PATHS['propagation'].replace('results.json', 'PROOF.md'), 16)
    y += h+35
    h = 500
    card(c, X0, y, CW, h, stroke=AMBER)
    c.text(X0+28, y+43, 'Finite interactions: prepare two packets and read a complete detector', 27, TEXT, 700)
    c.para(X0+28, y+85, CW-56,
           'Capacity limits couple the neutral components. A large-capacity theorem identifies '
           'a nonzero quartic channel and bounds packet admission, preparation and higher-order '
           'errors below its amplitude. The selected hit and its entire complement are retained.', 19, SOFT, 1.4)
    frames = d['collision']['interval_certified_five_block_collision']['complete_time_series']
    for j, frame in enumerate(frames):
        xx = X0+100+j*250
        lo, hi = map(Fraction, frame['native_preparation_hit_probability_interval_exact'])
        c.text(xx, y+210, 'Δt = '+frame['clock_Delta_times_t'], 18, MUTED, 500)
        c.text(xx, y+246, f'{100*float(lo):.5f}% –', 21, AMBER, 700)
        c.text(xx, y+278, f'{100*float(hi):.5f}%', 21, AMBER, 700)
    c.para(X0+28, y+321, CW-56,
           'Shown: an independent five-block, 1,024-state diagnostic, including preparation '
           'error. Its finite projection changes the free comparison’s momentum supports, '
           'so this plot alone does not isolate a free-forbidden channel. Its reference-relative occupation '
           'records are quantum probabilities, not classical trajectories.', 19, SOFT, 1.4)
    contrast = Fraction(frames[-1]['native_minus_matched_free_hit_difference_interval_exact'][0])
    c.text(X0+28, y+449,
           f'Same prepared state, same detector: final source − free probability ≥ {100*float(contrast):.5f} percentage points.',
           19, TEAL, 600)
    c.link(X0+28, y+h-20, 'Finite collision theorem, comparison domains and exact intervals',
           '../../'+PATHS['collision'].replace('results.json', 'PROOF.md'), 16)
    return y+h
