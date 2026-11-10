"""Receipt-bound finite defect-wave interaction and its complete readouts."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import sympy as s

from .canvas import AMBER, BLUE, EDGE, MUTED, SOFT, TEAL, TEXT, VIOLET
from .parts import CW, X0, card

ROOT = Path(__file__).resolve().parents[2]
RECEIPT = ROOT / 'model/coherent_defect_wave/results.json'
EXTENSIONS = ('background', 'spectroscopy', 'spatial', 'species',
              'cubic_readout', 'energy', 'dispersion', 'physical_attachment', 'mobile_flux',
              'native_routing_join', 'many_body_route', 'mass_instrument',
              'field_identity_audit', 'propagating_native_join', 'chiral_source',
              'chiral_regulator', 'chiral_higgs', 'chiral_particle_map', 'chiral_joint',
              'chiral_joint/two_cell')
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
             'physical_statistics_or_chiral_matter', 'interacting_relativistic_limit_or_SM_or_GeV',
             'chiral_SM_particle_identification', 'interacting_relativistic_limit',
             'source_GeV_masses', 'physical_Standard_Model_identification',
             'interacting_relativistic_mass', 'GeV_mass_derived',
             'source_mass_renormalization_or_SM_identity_or_GeV_prediction',
             'interacting_relativistic_particle_identification',
             'background_is_a_quantum_Higgs_vacuum_theorem',
             'interacting_chiral_continuum', 'measured_species_identified',
             'onsite_gap_is_relativistic_particle_mass', 'full_reference_is_physical_Lorentz_vacuum',
             'physical_pole_mass_or_GeV_prediction', 'physical_rest_mass_or_GeV_identification',
             'observed_species_identified', 'finite_gap_is_Clay_mass_gap'}


def check_extension_scope(value, name):
    """Check nested claims too: a declared rest parameter is permitted, an SM mass is not."""
    if isinstance(value, dict):
        for key, item in value.items():
            if key=='measured_GeV_predictions' and item:
                raise ValueError('joined extension physical promotion: '+name+'/'+key)
            if key=='empirically_identified_species' and item is not None:
                raise ValueError('joined extension physical promotion: '+name+'/'+key)
            absent_mass_table=(key=='source_GeV_masses' and isinstance(item,(dict,list)) and not item)
            if key in FORBIDDEN and item is not False and not absent_mass_table:
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



def check_mass_instrument(value, free):
    if value['source'] != free['program']:
        raise ValueError('mass apparatus and free field have different source clocks')
    if (not free['reference']['unique_positive_excitation_ground']
            or not free['scope']['CAR_statistics_are_declared_architecture']
            or not free['source']['all_occupations_retained']
            or free['algebra']['Wilson_clock_is_Floquet_logarithm']
            or value['continuum_rest_energy_in_source_units']!=free['reference']['full_Fock_gap']):
        raise ValueError('mass apparatus lacks its declared positive free reference')
    if value['unit_map']['source_energy_unit_in_GeV'] is not None:
        raise ValueError('mass apparatus substitutes an unowned GeV scale')
    for experiment in value['experiments']:
        records=experiment['complete_records']
        weights=[s.sympify(r['probability']) for r in records]
        if (len(records)!=36 or s.simplify(sum(weights,s.Integer(0)))!=1
                or any(p.is_nonnegative is not True for p in weights)):
            raise ValueError('mass apparatus omits complete reservoir records')
        if not experiment['filled_negative_band_reference']:
            raise ValueError('mass apparatus loses its positive reference')


def check_identity(value):
    dark=value['native_dark_current_module']
    if (not value['scope']['actual_native_R15_current_binding']
            or dark['module_prepared_or_bound_by_native_idle_law']
            or dark['native_K_on_module']!='0' or dark['internal_dimension']!=15
            or dark['empty_native_flux_code_actual_current_representation']!='trivial'):
        raise ValueError('native current module or empty flux identity changed')
    mass=value['mass_normalization_and_chirality']
    if (mass['gauge_invariant_Dirac_mass_parameter_dimension']!=5
            or not mass['Dirac_R15_contains_both_chiralities_in_same_internal_representation']):
        raise ValueError('charged Dirac identity loses its mass or chirality boundary')


def check_two_cell(value):
    if (value['schema']!='oph.coherent_defect_wave.chiral_joint.two_cell.v1'
            or value['source']['active_dimension']!=614400):
        raise ValueError('two-cell physical source changed')
    ground=next(row for row in value['full_spectra'] if row['boundary']==0)
    if (ground['full_physical_dimension']!=1053 or not ground['ground_is_certifiably_simple']
            or ground['ground_window_candidate_labels']!=[{'L':2,'twoSz':0,'eigen_index':0}]):
        raise ValueError('two-cell complete physical ground is missing')
    law=value['original_attempt_readout']
    rows=law['records']
    if (law['probabilities_renormalized'] or law['record_count_including_preparation_failure']!=len(rows)
            or len(rows)!=461 or rows[0]['outcome']!='preparation_failure'
            or rows[0]['probability_interval']!=['1052/1053']*2
            or law['preparation']['original_ground_success_probability']!='1/1053'):
        raise ValueError('two-cell original preparation failures are missing')
    if any(not 0<=s.Rational(r['probability_interval'][0])<=s.Rational(r['probability_interval'][1])<=1
           for r in rows):
        raise ValueError('two-cell invalid probability interval')
    low=sum(s.Rational(r['probability_interval'][0]) for r in rows)
    high=sum(s.Rational(r['probability_interval'][1]) for r in rows)
    if not low<=1<=high:raise ValueError('two-cell original probability law is incomplete')
    if {r['kind'] for r in value['conditional_readout']['records']}!={'add','add_fail','remove','remove_fail'}:
        raise ValueError('two-cell charged probe drops a branch')


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
        if (name in ('many_body_route','mass_instrument','field_identity_audit','propagating_native_join')
                and extension['schema']!='oph.coherent_defect_wave.'+name+'.v1'):
            raise ValueError('joined extension schema changed: '+name)
        check_extension_scope(extension, name)
        if name=='native_routing_join':check_native_join(extension)
        value['extensions'][name] = extension
    check_identity(value['extensions']['field_identity_audit'])
    check_two_cell(value['extensions']['chiral_joint/two_cell'])
    check_mass_instrument(value['extensions']['mass_instrument'],value['extensions']['many_body_route'])
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

    free=extensions['many_body_route']
    apparatus=extensions['mass_instrument']
    identity=extensions['field_identity_audit']
    rest=apparatus['continuum_rest_energy_in_source_units']
    gap=free['reference']['relative_charge_zero_gap']
    h,gap_x,w=434,22,(CW-44)/3
    for i,(tag,title,body,number,proof,color) in enumerate((
            ('A complete graded source', 'Positive particles and holes',
             'A separate declared local law retains all occupations of '
             +str(free['source']['literal_link_slots_per_cell'])+' graded flux slots per cell. '
             'Fill its negative-energy bands to obtain the lowest-energy reference. Adding an '
             'object or removing one costs positive energy; the missing object is a hole. '
             'Their free Dirac limit is controlled. The native gauge currents act trivially on these bare flux objects.',
             'Free rest energy '+rest+'; neutral-pair threshold '+gap,
             'many_body_route',BLUE),
            ('Actual source currents', 'Fifteen charged native states',
             'A two-address native wave combination is invariant under the idle law up to phase. '
             'Its even-occupation clouds form fifteen states with actual source currents and the '
             'Standard Model weight pattern. A declared gate transports the cloud with its '
             'graded flux object while conserving charge. The idle law gives zero interaction '
             'energy on this cloud and supplies no binding.',
             str(identity['native_dark_current_module']['actual_native_current_intertwiner_columns_checked'])+
             ' actual current identities; native K = 0', 'field_identity_audit',VIOLET),
            ('Two incoming objects, one link', 'A finite encounter with controlled error',
             'Added four-component memory puts both incoming objects on one clock with a '
             'common-speed free Dirac window. Their weak shared-link update changes the flux '
             'face with positive probability; complete detector bounds are retained. '
             'For one flux and one native occupation, bounded coupling gives vanishing '
             'contact response in a fixed-volume, fixed native-momentum window.',
             'Same source and clock; interacting limit open',
             'propagating_native_join',TEAL))):
        x=X0+i*(w+gap_x)
        card(c,x,y,w,h,stroke=color)
        c.caps(x+24,y+35,tag,color,12)
        yy=c.para(x+24,y+76,w-48,title,24,TEXT,1.28,700)
        c.para(x+24,yy+10,w-48,body,18,SOFT,1.36)
        c.text(x+24,y+h-58,number,16,color,700,maxw=w-48)
        c.link(x+24,y+h-25,'Exact source, reference and scope',
               '../../model/coherent_defect_wave/'+proof+'/PROOF.md',16)
    y+=h+22

    h=366
    card(c,X0,y,CW,h,stroke=BLUE)
    c.text(X0+28,y+43,'What does a source mass measurement look like?',28,TEXT,700)
    c.para(X0+28,y+83,CW-56,
        'A finite reservoir adds or removes an object from the filled reference. A battery '
        'supplies the work; the apparatus measures energy and relative pair count using the '
        'same Hamiltonian. Failed preparations are retained and can create neutral pairs. '
        'This gives an operational source mass, rather than reading a phase and calling it energy.',
        19,SOFT,1.4)
    diagram_y=y+180
    boxes=((X0+28,250,'Reservoir + battery','Prepare and supply work',TEAL),
           (X0+395,400,'Filled source reference','Add a particle / remove a hole',BLUE),
           (X0+910,350,'Energy + charge records',
            str(len(apparatus['experiments'][0]['complete_records']))+' outcomes, including failures',VIOLET))
    for xx,ww,title,subtitle,color in boxes:
        c.rect(xx,diagram_y,ww,72,'none',color,1.7,10)
        c.text(xx+ww/2,diagram_y+29,title,19,color,700,anchor='middle',maxw=ww-20)
        c.text(xx+ww/2,diagram_y+54,subtitle,15,SOFT,anchor='middle',maxw=ww-20)
    c.arrow(X0+290,diagram_y+36,X0+380,diagram_y+36,TEAL,2.4,9)
    c.arrow(X0+810,diagram_y+36,X0+895,diagram_y+36,BLUE,2.4,9)
    c.text(X0+28,y+294,'E² = |p|² + '+str(int(rest)**2)+
           '; rest energy = '+rest+' source units in this free limit',23,BLUE,700,maxw=CW-56)
    c.text(X0+28,y+326,'GeV value = '+rest+
           ' × energy unit in GeV. A physical clock reference fixes that unit; no measured particle is assigned.',
           17,MUTED,maxw=CW-56)
    c.link(X0+28,y+h-15,'Complete reservoir, battery and mass records',
           '../../model/coherent_defect_wave/mass_instrument/PROOF.md',16)
    y+=h+22

    field_map=extensions['chiral_particle_map']
    spectrum=field_map['background_spectrum']
    finite=extensions['chiral_joint/two_cell']
    charged={r['kind']:r for r in finite['conditional_readout']['records']}
    thresholds=[]
    for kind in ('add','remove'):
        first=charged[kind]['complete_energy_windows'][0]
        assert s.Rational(first['probability_interval'][0])>0
        interval=list(map(s.Rational,first['energy_transfer_interval']))
        thresholds.append(f'{float(sum(interval)/2):.8f}')
    h=678
    card(c,X0,y,CW,h,stroke=VIOLET)
    c.text(X0+28,y+43,'A chiral field dictionary from the actual source currents',28,TEXT,700)
    c.para(X0+28,y+83,CW-56,
        'A declared cell encoding retains all 30 spin/native fermion modes and every occupation. '
        'Its spectral left-Weyl propagation has a controlled free limit; its finite position '
        'kernel is nonlocal. In a specified neutral Higgs background, the source’s pair tensors '
        'join left fields to conjugate singlets. This is an exact change of field variables.',
        19,SOFT,1.4)
    table_y=y+193
    c.rect(X0+24,table_y,CW-48,182,'#142337',EDGE,1,8)
    for xx,label in ((X0+40,'Conditional field role'),(X0+540,'Electric charge'),(X0+800,'Rest parameter')):
        c.text(xx,table_y+28,label,18,VIOLET,700)
    roles=(('Up-type quark · 3 colors','2/3','up'),
           ('Down-type quark · 3 colors','−1/3','down'),
           ('Charged lepton · 1 field','−1','charged_lepton'),
           ('Neutral left lepton · 1 field','0','neutral'))
    for i,(role,charge,key) in enumerate(roles):
        yy=table_y+62+32*i
        c.text(X0+40,yy,role,20,TEXT,700)
        c.text(X0+540,yy,charge,20,SOFT)
        mass=spectrum['background_mass_energy'][key].replace('abs(', '|').replace(')', '|').replace('*','·')
        c.text(X0+800,yy,mass+' source units',20,SOFT)
    c.para(X0+28,y+408,CW-56,
        'In this background field dictionary, all fermion fields have spin ½. '
        'Weak currents act only on left components; color and electric '
        'charge act on both halves of the massive fields. Antiparticles have conjugate charges. '
        'The seven Dirac color/lepton components plus one neutral Weyl field form ONE family. '
        'No electron/muon/tau generation is selected, and v and the three Yukawa coefficients are inputs.',
        19,SOFT,1.4)
    c.para(X0+28,y+527,CW-56,
        'In the separate finite lepton control, quarks are empty spectators. Both quantum scalars '
        'and the quantum link are retained. Its 614,400 source states yield '
        '1,053 neutral physical states and one certified ground. Preparation succeeds with probability 1/1053; '
        'failures are recorded. Visible charged transfers are '+' and '.join(thresholds)+
        ' source units. Conjugate experiments reflect the boundary. These are detector energies; '
        'relativistic dispersion, electron/muon/tau identity and GeV masses are unproved.',
        18,MUTED,1.38)
    c.link(X0+28,y+h-25,'Source IDs, all currents and conditional mass formulas',
           '../../model/coherent_defect_wave/chiral_particle_map/PROOF.md',16)
    c.link(X0+805,y+h-25,'Quantum state and complete records',
           '../../model/coherent_defect_wave/chiral_joint/two_cell/READOUT_PROOF.md',16)
    y += h + 35
    c.link(X0 + 24, y, 'Joint evolution and reciprocal detector proofs',
           '../../model/coherent_defect_wave/PROOF.md', 17)
    c.link(X0 + 805, y, 'Exact joint receipt', '../../model/coherent_defect_wave/results.json', 17)
    return y + 22
