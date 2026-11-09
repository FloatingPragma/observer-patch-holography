"""Mass-first display of the common graph, with deferred gauge/gravity branches."""

ACTIVE_TARGET = 'Source-derived particle rest energies and masses, with declared SI/GeV comparison units and explicit capacity dependence.'


NODE_IDS = (
    'fixed_program', 'physical_preparation', 'local_records',
    'joint_refinement', 'relativistic_field', 'chiral_species_attachment',
    'particle_observables', 'independent_comparison', 'pure_gauge_gap',
    'geometric_gravity',
)

# These are display summaries. Titles, dependencies, status and evidence come
# from the canonical dependency document rather than a second requirement list.
SUMMARIES = {
    'fixed_program':
        'One open Wilson program fixes the action, cutoff, state, reference and clock. Capacity N=ln D distinguishes ambient and physical carriers from execution cost. A changed action declares an extension; measured mass inputs retain their ancestry.',
    'physical_preparation':
        'The exact Gauss-polynomial state retains regional boundary charges. Wilson loops establish 11 sectors with positive weight. At degree one the complete charge detector retains all 1,331 outcomes; the enormous spatial preparation tape is not claimed as executed.',
    'local_records':
        'Eight radix-five coherent meters attach opposite Q6=±6 transfers to complete energy records in [9.625,33.875] source units, with original joint success weight >1/4. They preserve source/reference histories without reading initial absolute charges. Failures and all labels are retained. Absolute terminal charge is a separate disturbing readout. Coherent offsets give finite Fourier labels, not exact momentum projections.',
    'joint_refinement':
        'Filtered and unfiltered probes share positive time dynamics. Normalized density histories with a²S→0 lose every fixed positive-energy weight. Actual Weyl probes have S=30/a², outside this obstruction; their fixed-band survival and spatial momentum require proof. The filtered relative-observation sector has commuting H≥0 and P, HΩ=PΩ=0; its sinc² energy/momentum weight may vanish.',
    'relativistic_field':
        'Three pairs of nonempty conserved occupation sectors have full-space floors above the actual ground. Two pairs are CAR odd; one is even. Pair preparation retains the opposite sector and displaced source/reference state in idle memory. Complete clock windows above the ground have probability ≥255/512. These finite sector thresholds identify no isolated particle or GeV rest mass. A separate A1u prepared window has weight >1/4.',
    'chiral_species_attachment':
        'Combined link conjugation and CAR particle/hole symmetry preserve the action, Gauss space and polynomial state. Opposite native occupations have identical spectra and conjugate complete detector laws. Proper cube rotations commute with C; odd inversion eigenstates have opposite parity under C. These finite labels identify no continuum spin or Standard Model species. Physical currents, mirror modes and chiral interactions require attachment.',
    'particle_observables':
        'An order-32 phase scan retains all 64 spatial/internal records. Eight scalar Fourier-energy events have original union weight >0.06646; complete clock reports retain >1/16. Parseval gives final-excitation union >15/256. Neighboring phase cancellation makes the pair-creation bound scale with Σcos²(qᵢ/2). These CAR-even records identify no particle species or rest mass; finite analyzer labels establish no conserved momentum or dispersion curve.',
    'independent_comparison':
        'SI units may be declared immediately. Finite mass comparison needs a resolved branch, momentum and initial-momentum control, and dispersion/readout errors. A continuum mass additionally needs uniform fixed-band survival. The primary physical-mass table in GeV is empty; calibrated inputs keep their ancestry.',
    'pure_gauge_gap':
        'Instantiate nodes 1–5 for an actual pure compact-simple gauge source or a proved source map. Dense local Euclidean decay plus finite-energy noncollapse implies a positive finite gap. Its physical hypotheses are conditional; it is separate from the photon-containing theory.',
    'geometric_gravity':
        'A declared coframe has H(I)=H_primary exactly and six strain means. Complete origin/sign records have a capacity-uniform source-marginal relative entropy D/β with 0.04099<D/β≤14175/32. The coin and reference remain in the complete output; this entropy comparison traces out the quantum reference. β is the cooling normalization. Proper geometry, Lorentz stress, diamond/area response and Einstein limits require proof.',
}


def validate(graph):
    """Reject a changed graph or scope until its presentation is reviewed."""
    if graph.get('schema') != 'oph.joint_reconstruction_dependencies.v1':
        raise ValueError('unsupported joint dependency schema')
    if (graph.get('active_target') != ACTIVE_TARGET
            or graph.get('deferred_branches') != ['pure_gauge_gap','geometric_gravity']):
        raise ValueError('mass-first target and deferred branches require explicit review')
    nodes = graph['nodes']
    if tuple(node['id'] for node in nodes) != NODE_IDS:
        raise ValueError('the ten canonical dependencies must be retained')
    for index, node in enumerate(nodes):
        expected = [] if index == 0 else [NODE_IDS[4 if index >= 8 else index-1]]
        status = ('finite_constructed' if index < 3 else
                  'conditional_criterion' if index == 8 else 'physical_proof_required')
        if (node['number'] != index+1 or node['depends_on'] != expected
                or node['status'] != status):
            raise ValueError('dependency ordering or physical scope requires review')
        if not node['acceptance'] or not node['evidence']:
            raise ValueError('dependency acceptance and evidence are required')
    branch = nodes[8]
    if (branch['required_shared_nodes'] != list(NODE_IDS[:5])
            or branch['required_for_full_particle_spectrum']
            or branch['implies_full_theory_positive_global_gap']):
        raise ValueError('pure Yang-Mills must retain its separate source domain')
    gravity = nodes[9]
    if (gravity['required_for_full_particle_spectrum']
            or gravity['proper_geometry_or_Einstein_limit_proved']):
        raise ValueError('finite source energy cannot be promoted to Einstein dynamics')
    candidate = graph['conditional_candidate']
    if (candidate['status'] != 'conditional_extension'
            or candidate['discharges_primary_physical_dependencies']
            or candidate['calibrated_tree_coefficients_are_independent_predictions']
            or candidate['full_theory_positive_global_gap_claimed']):
        raise ValueError('conditional coefficients cannot be promoted to physical predictions')
    return nodes


def render(c, data, y):
    from .drawing import W, M, INK, TEAL, ORANGE, BLUE, LINE
    graph = data['joint_dependencies']
    nodes = validate(graph)
    numbers = {node['id']: node['number'] for node in nodes}
    c.items.append('<g id="joint-field-particle-dependencies">')
    c.text(M, y, 'MASS FIRST: ONE SOURCE, ORIGINAL WEIGHTS AND DECLARED UNITS',
           80, INK, 700, serif=True)
    c.para(M, y+115, W-2*M,
           'Nodes 1–3 are finite constructions. Nodes 4–8 contain proved finite response, spatial-action and readout results, with physical mass conditions explicit. Pure Yang–Mills and gravity are deferred branches. No Standard Model species or primary mass in GeV is identified.',
           56, INK, max_y=y+240)
    chain = ' → '.join(str(node['number']) for node in nodes[:8])
    c.text(M, y+315, 'Shared particle route:  '+chain, 64, TEAL, 700)
    c.text(M+12300, y+315, 'Deferred pure gauge:  5 → 9', 64, BLUE, 700)
    c.text(M+21000, y+315, 'Deferred gravity:  5 → 10', 64, ORANGE, 700)
    width = (W-2*M)/3
    for index, node in enumerate(nodes[:9]):
        x = M+(index % 3)*width
        top = y+400+(index//3)*760
        color = TEAL if index < 3 else BLUE if index == 8 else ORANGE
        status = ('FINITE CONSTRUCTED' if index < 3 else
                  'DEFERRED · CONDITIONAL GAP CRITERION' if index == 8 else 'PHYSICAL PROOF REQUIRED')
        after = ', '.join(str(numbers[key]) for key in node['depends_on'])
        tag = status+(' · after '+after if after else ' · declared program')
        c.rect(x, top, width-70, 720, '#fffdf7', LINE, 26)
        c.rect(x+70, top+50, width-210, 130,
               '#e3efeb' if index < 3 else '#e9edf6' if index == 8 else '#fbefdd', radius=20)
        c.text(x+120, top+145, tag, 57, color, 700)
        c.para(x+80, top+280, width-230,
               f"{node['number']:02d}  {node['title']}", 70, INK, 700, max_y=top+355)
        c.para(x+80, top+410, width-230, SUMMARIES[node['id']],
               64, color, max_y=top+580)
        c.link(x+80, top+660, 'Acceptance conditions and exact evidence',
               'model/particle_joint_reconstruction/DEPENDENCIES.md', 49)
        theorem = {'local_records': 'model/particle_primary_charged_readout/PROOF.md',
                   'joint_refinement': 'model/particle_primary_weighted_survival/PROOF.md',
                   'relativistic_field': 'model/particle_primary_sector_threshold/PROOF.md',
                   'chiral_species_attachment': 'model/particle_primary_charge_conjugation/PROOF.md',
                   'particle_observables': 'model/particle_primary_populated_momentum/PROOF.md',
                   'independent_comparison': 'model/particle_primary_mass_readout/PROOF.md'}.get(
                       node['id'], node['evidence'][0]['path'])
        c.link(x+width/2, top+660, 'Scoped construction or theorem',
               theorem, 49)
        extra = {
            'fixed_program': ('Source capacity', 'model/particle_primary_capacity/PROOF.md'),
            'local_records': ('Local influence', 'model/particle_primary_locality/PROOF.md'),
            'joint_refinement': ('Origin momentum', 'model/particle_primary_origin_momentum/PROOF.md'),
            'relativistic_field': ('Exact residual', 'model/particle_primary_residual_refinement/PROOF.md'),
            'particle_observables': ('Positive window', 'model/particle_positive_window/PROOF.md'),
        }.get(node['id'])
        if extra: c.link(x+width*.77, top+660, extra[0], extra[1], 49)
        fourth = {
            'local_records': ('Joint H/P', 'model/particle_primary_momentum_transfer/PROOF.md'),
            'joint_refinement': ('Weyl records', 'model/particle_primary_local_fluctuation/PROOF.md'),
            'relativistic_field': ('Prepared A1u', 'model/particle_primary_spectral_resolution/PROOF.md'),
            'particle_observables': ('Pair channels', 'model/particle_primary_pair_channels/PROOF.md'),
        }.get(node['id'])
        if fourth: c.link(x+width*.9, top+660, fourth[0], fourth[1], 49)
        earlier = {
            'local_records': (
                ('Finite Fourier aperture', 'model/particle_primary_momentum/PROOF.md'),
                ('Complete energy clocks', 'model/particle_primary_joint_spectroscopy/PROOF.md')),
            'joint_refinement': (
                ('Spatial action', 'model/particle_primary_spatial_algebra/PROOF.md'),
                ('Occupation bound', 'model/particle_primary_refinement_response/PROOF.md'),
                ('Weight ≥3/100 in [1/(100β),46/a]', 'model/particle_primary_energy_tightness/PROOF.md'),
                ('Filtered algebra', 'model/particle_primary_field_algebra/PROOF.md'),
                ('Common time', 'model/particle_primary_current_tail/PROOF.md')),
            'relativistic_field': (
                ('Selected occupation bound', 'model/particle_primary_occupation_bound/PROOF.md'),
                ('Distinct selected A1u response', 'model/particle_primary_odd_residual/PROOF.md'),
                ('Full-source Rayleigh bound', 'model/particle_primary_variational_readout/PROOF.md')),
            'particle_observables': (
                ('Local response window', 'model/particle_primary_response_moments/PROOF.md'),
                ('Local [4,8] response', 'model/particle_primary_branch_resolution/PROOF.md'),
                ('Charged endpoint band', 'model/particle_primary_charged_string/PROOF.md'),
                ('Checkerboard A2g band', 'model/particle_primary_staggered_response/PROOF.md')),
        }.get(node['id'], ())
        for j, (label, path) in enumerate(earlier):
            c.link(x+80+j*width/len(earlier), top+610, label, path, 45)
    gravity = nodes[9]; top = y+2720
    c.rect(M, top, W-2*M-70, 530, '#fffdf7', LINE, 26)
    c.text(M+90, top+110, 'DEFERRED · PHYSICAL PROOF REQUIRED · after 5', 57, ORANGE, 700)
    c.text(M+90, top+230, '10  '+gravity['title'], 72, INK, 700)
    c.para(M+90, top+345, W-2*M-250, SUMMARIES[gravity['id']], 64, ORANGE, max_y=top+420)
    c.link(M+90, top+480, 'Actual source energy, complete records and marginal entropy',
           'model/particle_primary_stress/PROOF.md', 49)
    c.link(M+14700, top+480, 'Conditional Einstein theorem and proper-geometry hypotheses',
           'model/matched_einstein_limit/PROOF.md', 49)
    c.link(M+23000, top+480, 'Uniform D/β and source time bounds',
           'model/particle_primary_current_tail/PROOF.md', 49)
    c.link(M+9000, top+480, 'Flat coframe and strain means',
           'model/particle_primary_open_coframe/PROOF.md', 49)
    c.link(M+12000, top+480, 'Modular clock',
           'model/particle_primary_gravity_premises/PROOF.md', 49)
    c.para(M, y+3360, W-2*M,
           'Conditional candidate: a declared chiral gauge-Higgs action supplies calibrated tree/running coefficients, with a separately calibrated neutrino extension. These coefficients are not physical pole masses, full radiative matching or independent predictions. The massless photon excludes a positive global vacuum gap for the full theory.',
           60, BLUE, max_y=y+3470)
    c.link(M, y+3540, 'Canonical dependency graph and acceptance conditions',
           'model/particle_joint_reconstruction/DEPENDENCIES.md', 51)
    c.link(M+14700, y+3540, 'Conditional candidate, spectrum conventions and calibration ancestry',
           graph['conditional_candidate']['proof'], 51)
    c.items.append('</g>')
    return y+3630
