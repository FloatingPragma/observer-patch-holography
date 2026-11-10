"""Receipt-bound finite defect-wave interaction and its complete readouts."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from .canvas import AMBER, BLUE, EDGE, MUTED, SOFT, TEAL, TEXT, VIOLET
from .parts import CW, X0, card

ROOT = Path(__file__).resolve().parents[2]
RECEIPT = ROOT / 'model/coherent_defect_wave/results.json'
EXTENSIONS = ('background', 'spectroscopy', 'spatial', 'species',
              'cubic_readout', 'energy', 'dispersion', 'physical_attachment', 'mobile_flux',
              'native_routing_join')
FORBIDDEN = {'physical_mass_identified', 'relativistic_mass_or_GeV_value',
             'GeV_mass', 'physical_particle_identity', 'SM_species_identified',
             'rest_masses_derived', 'interacting_relativistic_continuum',
             'GeV_prediction', 'physical_rest_mass_identified',
             'energies_are_physical_rest_masses', 'source_derived_GeV_masses',
             'physical_particle_or_Standard_Model_identification',
             'SM_particle_identification', 'GeV_masses', 'interacting_Lorentz_limit',
             'dynamical_chiral_gauge_theory', 'SM_identity_or_GeV_mass_derived',
             'interacting_relativistic_limit_proved', 'SM_particle_or_GeV_mass_identified',
             'positive_energy_many_body_vacuum_derived', 'physical_interacting_relativistic_limit',
             'scattering_or_SM_identity_or_GeV_mass', 'positive_energy_many_body_vacuum',
             'physical_statistics_or_chiral_matter'}


def check_extension_scope(value, name):
    """Check nested claims too: a declared rest parameter is permitted, an SM mass is not."""
    if isinstance(value, dict):
        for key, item in value.items():
            if key in FORBIDDEN and item is not False:
                raise ValueError('joined extension physical promotion: ' + name + '/' + key)
            check_extension_scope(item, name + '/' + key)
    elif isinstance(value, list):
        for item in value:
            check_extension_scope(item, name)


def check_native_join(value):
    if not value['shared_tuple']['native_U_reads_the_same_link_register_as_routing']:
        raise ValueError('native routing join has no shared link register')
    finite=value['finite_join']
    if finite['Taylor_output_renormalized'] or not finite['complete_binary_instrument']:
        raise ValueError('native routing join substitutes an incomplete approximate readout')
    rows=finite['exact_unitary_probability_enclosures']
    if len(rows)!=2 or {row['changed_face'] for row in rows}!={False,True}:
        raise ValueError('native routing join omits a detector outcome')
    bounds={row['changed_face']:(Fraction(row['lower']),Fraction(row['upper'])) for row in rows}
    if (not all(0<=lo<=hi<=1 for lo,hi in bounds.values()) or bounds[True][0]<=0
            or bounds[False][0]+bounds[True][1]!=1 or bounds[False][1]+bounds[True][0]!=1):
        raise ValueError('native routing join probability bounds are not complete')


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
    value['extensions'] = {}
    for name in EXTENSIONS:
        path = RECEIPT.parent / name / 'results.json'
        extension = json.loads(path.read_text())
        if not extension.get('input_sha256'):
            raise ValueError('missing joined extension custody: ' + name)
        for source, digest in extension['input_sha256'].items():
            target = (ROOT/source).resolve()
            if not target.is_relative_to(ROOT) or hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError('stale joined extension input: ' + source)
        check_extension_scope(extension, name)
        if name=='native_routing_join':check_native_join(extension)
        value['extensions'][name] = extension
    return value


def input_pins(value):
    pins = {**value['input_sha256'], str(RECEIPT.relative_to(ROOT)):
            hashlib.sha256(RECEIPT.read_bytes()).hexdigest()}
    for name, extension in value['extensions'].items():
        path = RECEIPT.parent/name/'results.json'
        pins.update(extension['input_sha256'])
        pins[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


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

    extensions = value['extensions']
    background = extensions['background']
    h, gap, w = 394, 22, (CW - 44) / 3
    for i, (tag, title, body, number, proof, color) in enumerate((
            ('Prepare → disturb → read', 'A fixed background with visible excitations',
             'Uniform link and wave amplitudes form a pure state fixed by the actual tick. '
             'A complete apparatus prepares it. Reading its defect face disturbs that state; '
             'a coherent filter then removes the unchanged background and retains every outcome.',
             'Fresh preparation success: ' + background['joint_success_per_fresh_ambient_trial'],
             'background', BLUE),
            ('Cell → neighbor → detector', 'Connected propagation with a finite cone',
             'Place cells on a declared cubic grid. A partial exchange carries wave information '
             'to a neighbor while keeping its link memory local and preserving its constraint. '
             'The arriving wave changes the neighbor’s defect motion. A finite cone bounds influence.',
             'Small-momentum spectator mode: quadratic', 'spatial', TEAL),
            ('One wave → several occupations', 'Several waves share the same links',
             'The occupation extension uses antisymmetric wave states and one shared defect register. '
             'Its one-wave sector is exactly the original program. Two-wave records are computed. '
             'Fermionic hopping and tensor-wave exchange are distinct spatial laws.',
             'Internal symmetry: U(3) × U(2)', 'species', VIOLET))):
        x = X0 + i * (w + gap)
        card(c, x, y, w, h)
        c.caps(x + 24, y + 35, tag, color, 13)
        yy = c.para(x + 24, y + 76, w - 48, title, 24, TEXT, 1.28, 700)
        c.para(x + 24, yy + 10, w - 48, body, 18, SOFT, 1.36)
        c.text(x + 24, y + h - 58, number, 16, color, 700, maxw=w - 48)
        c.link(x + 24, y + h - 25, 'Exact construction and domain',
               '../../model/coherent_defect_wave/' + proof + '/PROOF.md', 16)
    y += h + 22

    cubic = extensions['cubic_readout']
    energies = extensions['dispersion']
    collision = extensions['energy']['finite_dark_collision']
    records = len(cubic['experiment']['complete_joint_records'])
    refined = len(cubic['higher_momentum_experiment']['complete_joint_records'])
    ratios = ' : '.join(f'{float(x):.6f}'.rstrip('0').rstrip('.')
                        for x in energies['three_nonzero_source_energy_ratios_approximate'])
    h, gap, w = 438, 22, (CW - 44) / 3
    for i, (tag, title, body, number, proof, color) in enumerate((
            ('Whole field → momentum record', 'Three spatial directions, one clock',
             'The cubic experiment translates each whole cell: its wave and retained link data. '
             'It records momentum and clock interference together, with every alternative kept. '
             'The larger cube adds complex translation phases; its complete probabilities are '
             'computed without materializing every branch vector. Preparation failures count too.',
             f'{records} records; {refined} on the larger cube', 'cubic_readout', BLUE),
            ('One explicit positive energy law', 'Exact energies; a failed relativistic test',
             'The declared positive clock has exact carrier energies and three-dimensional '
             'dispersion. Its three populated positive gaps have the ratios below. Their '
             'curvatures imply different limiting speeds; refining this additive branch '
             'does not give a common relativistic mass shell. These are source energies, not lepton masses.',
             ratios, 'dispersion', AMBER),
            ('Two waves → three detector outcomes', 'An exact finite collision',
             'Six actual source states form an invariant interacting sector of the same '
             'positive clock. From the specified pair preparation, two waves exchange '
             'occupation between cells. The detector '
             'keeps zero, one and two waves on the left. The source contact term changes '
             'this response compared with independent waves. The carrier flux stays in its cell.',
             'Maximum double occupation: ' +
             collision['maximum_double_occupation_probability'].replace('pi', 'π').replace('**2', '²'),
             'energy', TEAL))):
        x = X0 + i * (w + gap)
        card(c, x, y, w, h)
        c.caps(x + 24, y + 35, tag, color, 12)
        yy = c.para(x + 24, y + 76, w - 48, title, 24, TEXT, 1.28, 700)
        c.para(x + 24, yy + 10, w - 48, body, 18, SOFT, 1.36)
        c.text(x + 24, y + h - 58, number, 16, color, 700, maxw=w - 48)
        c.link(x + 24, y + h - 25, 'Source, clock and complete readout',
               '../../model/coherent_defect_wave/' + proof + '/PROOF.md', 16)
    y += h + 22

    mobile = extensions['mobile_flux']
    paired, route = mobile['paired_flux'], mobile['protected_object_route']
    mass = route['positive_branch_rest_parameter']
    h, gap, w = 470, 22, (CW - 22) / 2
    for i, (tag, title, body, proof, color) in enumerate((
            ('Boundary data travel with the defect', 'The protected pair really changes cell',
             'The enlarged link domain includes a flat vacuum. An interface gate carries '
             'the interior flux and its exterior compensator together. Complete position '
             'records give equal left/right chances after one exchange and certain transfer '
             'after two. The moved link changes the receiver’s wave response. Five internal '
             'face locations describe one protected type, not five identified species.',
             'PROOF.md', AMBER),
            ('Same source object → free Dirac window', 'A rest frequency for that same flux pair',
             'Two routing lanes carry the literal paired-flux link code. Their local memory '
             'has ' + str(route['full_two_lane_local_dimension']) + ' states, including collisions. '
             'A chosen mass coin uses the actual count of ' + str(route['source_face_count']) +
             ' flipped faces. The same object’s position and Fourier readouts have a controlled '
             'free Dirac limit. Its spinor memory and coupling law are declared architecture.',
             'ROUTING_COMPLETION.md', VIOLET))):
        x = X0 + i * (w + gap)
        card(c, x, y, w, h, stroke=color)
        c.caps(x + 24, y + 35, tag, color, 12)
        yy = c.para(x + 24, y + 76, w - 48, title, 25, TEXT, 1.28, 700)
        c.para(x + 24, yy + 12, w - 48, body, 18, SOFT, 1.36)
        if i == 0:
            box_y = y + h - 143
            for box_x in (x + 44, x + w - 199):
                c.rect(box_x, box_y, 145, 62, 'none', EDGE, 1.7, 10)
                c.circle(box_x + 55, box_y + 24, 8, AMBER)
                c.circle(box_x + 90, box_y + 24, 8, 'none', AMBER, 2)
                c.line(box_x + 63, box_y + 24, box_x + 82, box_y + 24, AMBER, 1.4)
            c.arrow(x + 204, box_y + 25, x + w - 215, box_y + 25, color, 2.4, 9)
            c.text(x + 61, box_y + 52, 'inside + outside', 14, SOFT)
            c.text(x + w - 184, box_y + 52, 'inside + outside', 14, SOFT)
            c.text(x + 24, y + h - 55,
                   'P(left/right) = ' + paired['charge_position_records_after_one_exchange']['10'] +
                   ' / ' + paired['charge_position_records_after_one_exchange']['01'], 18, color, 700)
        else:
            c.text(x + 24, y + h - 115, 'ω² = |p|² + ' + str(mass*mass) + ';  |ω(0)| = ' +
                   str(mass) + ' (declared units)', 23, color, 700, maxw=w-48)
            c.para(x + 24, y + h - 79, w - 48,
                   'Bounded-window error: t·a·(3P² + ' + str(3*abs(mass)) +
                   'P). No GeV or Standard Model identity follows.',
                   17, MUTED, 1.3)
        c.link(x + 24, y + h - 25, 'Literal source encoding and proof',
               '../../model/coherent_defect_wave/mobile_flux/' + proof, 16)
    y += h + 22

    native_join=extensions['native_routing_join']
    finite_join=native_join['finite_join']
    h = 274
    card(c, X0, y, CW, h, stroke=TEAL)
    c.text(X0 + 28, y + 43, 'One finite clock joins routed flux to native repair', 28, TEXT, 700)
    c.para(X0 + 28, y + 83, CW - 56,
        'Routing and the weak native update read the very same link word. The specified '
        'joint clock retains the gauge constraints and every collision sector. In a side-three '
        'cube, native repair gives a strictly positive chance of changing the routed pair’s '
        'internal face. Both detector outcomes have rigorous probability intervals; the Taylor '
        'vector is never renormalized. Empty native registers recover the free route exactly, '
        'so this interaction gives no vacuum mass correction there. Its interacting Lorentz '
        'limit remains unproved.', 19, SOFT, 1.4)
    c.text(X0 + 28, y + h - 60, 'Certified vector error ≤ ' +
           finite_join['Taylor_vector_error_upper'] + ' · both outcomes retained',
           18, TEAL, 700, maxw=CW-56)
    c.link(X0 + 28, y + h - 25, 'Same source, reference, clock and detector',
           '../../model/coherent_defect_wave/native_routing_join/PROOF.md', 16)
    y += h + 22

    h = 228
    card(c, X0, y, CW, h, stroke=BLUE)
    c.text(X0 + 28, y + 43, 'The remaining connection is physical identification', 28, TEXT, 700)
    c.para(X0 + 28, y + 83, CW - 56,
        'Spin and currents are proved for the declared free Dirac route. Fifteen internal '
        'weights match one chiral fermion family, and their Higgs mass channels are checked; '
        'the scalar state and coupling strengths are not selected. The positive energy clock, '
        'the ordered cubic clock and the native-coupled spinor route are explicit, distinct programs. '
        'Their receipts cannot be exchanged. A shared interacting relativistic limit, a '
        'physical vacuum, statistics, measured species and GeV masses remain unproved.',
        19, SOFT, 1.4)
    c.link(X0 + 28, y + h - 25, 'Spin, currents and mass constraints',
           '../../model/coherent_defect_wave/physical_attachment/PROOF.md', 16)
    c.link(X0 + 805, y + h - 25, 'Dependencies and exact status',
           '../../model/coherent_defect_wave/PARTICLE_PATH.md', 16)
    y += h + 35
    c.link(X0 + 24, y, 'Joint evolution and reciprocal detector proofs',
           '../../model/coherent_defect_wave/PROOF.md', 17)
    c.link(X0 + 805, y, 'Exact joint receipt', '../../model/coherent_defect_wave/results.json', 17)
    return y + 22
