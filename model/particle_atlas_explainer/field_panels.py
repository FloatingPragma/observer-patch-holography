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
    'charged_field': 'model/affine_gauss_full_spectrum/results.json',
    'charged_clock': 'model/affine_gauss_charged_clock/results.json',
    'vector': 'model/interacting_block_proca/results.json',
    'conjugation': 'model/coherent_defect_wave/chiral_joint/interacting_charge_conjugation/results.json',
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
    charged = data['charged_field']
    if (charged['additional_state_truncation']
            or charged['physical_one_charge_dimension'] != 70928
            or charged['scope']['mass_GeV'] is not None
            or charged['scope']['Standard_Model_identity'] is not None
            or charged['scope']['translation_invariant_mass_shell_identified']):
        raise ValueError('Charged field interpretation changed')
    clock = data['charged_clock']
    if (not clock['phase_apparatus']['common_source_evolution_retained']
            or clock['phase_apparatus']['higher_source_charge_sectors_discarded']
            or clock['scope']['exact_ground_projection']
            or clock['scope']['GeV_mass'] is not None):
        raise ValueError('Charged detector domain changed')
    vector = data['vector']
    if (vector['construction']['Standard_Model_particle_assignment'] is not None
            or vector['construction']['measured_inputs']
            or vector['construction']['strict_continuum_interaction'] != 'free; fixed-capacity interactions are controlled errors'
            or not vector['complete_instruments']['all_failure_and_complement_records_retained']):
        raise ValueError('Vector readout interpretation changed')
    if data['conjugation']['reference_boundary']['same_reference_total_energy_equality_claimed']:
        raise ValueError('Conjugate reference energy was silently changed')


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
    y = charged_field_panel(c, d, y)+35
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
           'The unfiltered components are internal scalar fields. The next panel shows the '
           'additional readout and symmetry map that give a free vector description.', 19, SOFT, 1.4)
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
    y = vector_panel(c, d, y)+35
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


def charged_field_panel(c, d, y):
    h = 800
    card(c, X0, y, CW, h, stroke=BLUE)
    c.text(X0+28, y+43, 'A charge and its dynamical field: an actual finite energy', 27, TEXT, 700)
    c.para(X0+28, y+85, CW-56,
           'A declared electric/magnetic source lets charge move while its compensating reference '
           'stays fixed. Gauss constraints determine the longitudinal flux. Finite Peierls hopping '
           'couples that charge to every retained transverse field state.', 19, SOFT, 1.4)
    labels = [('13 matter modes', 'where the charge can be'),
              ('5,456 field states', 'the field can change with it'),
              ('70,928 states', 'the full one-charge calculation')]
    for i, (label, sub) in enumerate(labels):
        xx=X0+28+i*460
        c.rect(xx, y+186, 414, 86, '#132636', BLUE, 1.5, 10)
        c.text(xx+20, y+221, label, 23, BLUE, 700)
        c.text(xx+20, y+250, sub, 17, SOFT)
        if i<2:c.arrow(xx+420,y+230,xx+451,y+230,BLUE,2,7)
    lo, hi = map(Fraction, d['charged_field']['actual_charged_ground_energy_interval'])
    scale = 10**13
    lo_n, hi_n = (lo*scale).__floor__(), (hi*scale).__ceil__()
    lo_text = f'{lo_n//scale}.{lo_n%scale:013d}'
    hi_text = f'{hi_n//scale}.{hi_n%scale:013d}'
    c.text(X0+28, y+320, f'Total charged energy: [{lo_text}, {hi_text}] source units', 24, TEAL, 700)
    c.para(X0+28, y+351, CW-56,
           'The fixed reference contributes exactly 1. The field dressing lowers the zero-field '
           'compression by about 0.0001750058. This correction is computed with the same Hamiltonian; '
           'it is unrelated to a fitted lepton formula.', 18, SOFT, 1.4)
    c.text(X0+28, y+440, 'Where would a position detector click in this charged ground?', 20, TEXT, 600)
    for i,row in enumerate(d['charged_field']['actual_ground_position_records']):
        lo,hi=map(Fraction,row['actual_ground_position_probability'])
        p=float((lo+hi)/2);xx=X0+40+i*104; bh=310*p
        color=AMBER if row['source_vertex']==1 else BLUE
        c.rect(xx,y+570-bh,62,bh,color,r=3)
        c.text(xx+31,y+595,str(row['source_vertex']),16,MUTED,500,'middle')
        c.text(xx+31,y+556-bh,f'{100*p:.1f}%',15,SOFT,500,'middle')
    c.text(X0+28,y+626,'All 13 locations are retained; vertex 1 contains the fixed compensating reference.',17,MUTED)
    clock=d['charged_clock']['finite_rational_gate_record_bounds']
    fidelity_n=(Fraction(clock['accepted_source_ground_fidelity_lower'])*10**11).__floor__()
    fidelity_text=f'{fidelity_n//10**11}.{fidelity_n%10**11:011d}'
    c.para(X0+28,y+658,CW-56,
           f'A finite clock detector retains all 14 × 2²⁰ records, including 75% preparation failure. '
           f'Its accepted source has ground fidelity above {fidelity_text}. '
           'The full energy and numerical errors are bounded. This finite boundary problem '
           'does not yet identify a relativistic particle, its spin or a GeV mass.',19,SOFT,1.4)
    c.link(X0+28,y+h-20,'Full charged Hamiltonian and independent residual certificate',
           '../../'+PATHS['charged_field'].replace('results.json','PROOF.md'),16)
    c.link(X0+810,y+h-20,'Complete finite clock detector',
           '../../'+PATHS['charged_clock'].replace('results.json','PROOF.md'),16)
    return y+h


def vector_panel(c,d,y):
    h=610
    card(c,X0,y,CW,h,stroke=VIOLET)
    c.text(X0+28,y+43,'From three canonical fields to a massive vector readout',27,TEXT,700)
    c.para(X0+28,y+85,CW-56,
           'Keep the same native blocks, actual ground, clock and orientation reference. '
           'Measure specific combinations of the fields and their conjugate momenta. '
           'Those filters supply three vector polarizations and the constrained time component.',19,SOFT,1.4)
    for i,(title,body,color) in enumerate([
            ('Canonical pairs','Φ₁, Φ₂, Φ₃ and Π₁, Π₂, Π₃',TEAL),
            ('Momentum-dependent filters','A₀ uses i p·Π / Ω',BLUE),
            ('Complete vector detector','16 records, including every flag',VIOLET)]):
        xx=X0+28+i*460
        c.rect(xx,y+194,414,104,'#132636',color,1.5,10)
        c.text(xx+20,y+231,title,21,color,700)
        c.text(xx+20,y+270,body,17,SOFT)
        if i<2:c.arrow(xx+420,y+248,xx+451,y+248,color,2,7)
    c.para(X0+28,y+342,CW-56,
           'The polarization sum, Proca equations and induced spin-one Lorentz action are proved '
           'for the free limiting description. The finite source’s Lorenz defect and complete '
           'measurement errors are bounded. The time component needs canonical momentum; '
           'three field labels alone do not establish this result.',19,SOFT,1.4)
    c.text(X0+28,y+454,'Rest frequency = the native gap Δ. Effective mass parameter = Δ / c*².',23,TEAL,700)
    c.para(X0+28,y+494,CW-56,
           'The explicit finite examples bound sixteen-outcome laws using a 19-mode Gaussian '
           'calculation and full native error estimates. Their large block counts are analytic '
           'capacity bounds. This readout does not identify W, Z or an interacting Standard Model sector.',19,MUTED,1.4)
    c.link(X0+28,y+h-20,'Same-source vector filters, boosts, full records and capacity bounds',
           '../../'+PATHS['vector'].replace('results.json','PROOF.md'),16)
    return y+h
