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
        'Boundary records have uniform influence bounds. Two energy readouts and a Fourier aperture retain the initial pointer, all failures and the idle reference, with controlled clock/aperture errors. The first measurement preserves the source marginal and changes purification coherences. No conserved momentum is inferred.',
    'joint_refinement':
        'Fixed-filter words have one stationary positive-energy limit and convergent finite pointer histories. Spatial smearing supplies a continuous R3 action. The selected state is not proved spatially invariant; a same-representation spatial implementation, momentum origin and nonzero finite-energy weight require proof.',
    'relativistic_field':
        'On the eight-cell source, occupation-sensitive bounds give ground energy [−727.498195, −723.7678199]; exact local moments bound the full-space trial residual by 5.231016. No level is isolated. A resolved energy–momentum branch and controlled relativistic dispersion require proof.',
    'chiral_species_attachment':
        'Identify physical currents, spin, statistics and families. Control mirror and extra modes, the Higgs vacuum and neutrino completion. Any claimed defect origin needs an actual field attachment. The declared chiral extension is distinct from the primary 60-mode carrier.',
    'particle_observables':
        'Actual eight-cell neutral pair probes: original transfer weight ≥1/2 in [2,10], final excitation weight ≥3/8 in [2,10.5] source units, before clock errors. All outcomes survive. This broad, even, gauge-neutral response has no assigned particle name, spin or rest mass. Larger-capacity transfer requires proof.',
    'independent_comparison':
        'SI units may be declared immediately. Finite mass comparison needs a resolved branch, momentum and initial-momentum control, and dispersion/readout errors. A continuum mass additionally needs uniform fixed-band survival. The primary physical-mass table in GeV is empty; calibrated inputs keep their ancestry.',
    'pure_gauge_gap':
        'Instantiate nodes 1–5 for an actual pure compact-simple gauge source or a proved source map. Dense local Euclidean decay plus finite-energy noncollapse implies a positive finite gap. Its physical hypotheses are conditional; it is separate from the photon-containing theory.',
    'geometric_gravity':
        'The primary source has exact local energy currents and a finite regional source-marginal entropy first law. The coin and reference remain in the complete output; their entropy is a different functional. Proper geometry, Lorentz stress, modular diamond and area responses, and controlled Einstein small-ball limits require proof.',
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
        theorem = {'local_records': 'model/particle_primary_joint_spectroscopy/PROOF.md',
                   'joint_refinement': 'model/particle_primary_field_algebra/PROOF.md',
                   'relativistic_field': 'model/particle_primary_occupation_bound/PROOF.md',
                   'particle_observables': 'model/particle_primary_response_moments/PROOF.md',
                   'independent_comparison': 'model/particle_primary_mass_readout/PROOF.md'}.get(
                       node['id'], node['evidence'][0]['path'])
        c.link(x+width/2, top+660, 'Scoped construction or theorem',
               theorem, 49)
        extra = {
            'fixed_program': ('Source capacity', 'model/particle_primary_capacity/PROOF.md'),
            'local_records': ('Local influence', 'model/particle_primary_locality/PROOF.md'),
            'joint_refinement': ('Spatial action', 'model/particle_primary_spatial_algebra/PROOF.md'),
            'relativistic_field': ('Exact residual', 'model/particle_primary_residual_refinement/PROOF.md'),
            'particle_observables': ('Positive window', 'model/particle_positive_window/PROOF.md'),
        }.get(node['id'])
        if extra: c.link(x+width*.77, top+660, extra[0], extra[1], 49)
        fourth = {
            'local_records': ('Fourier', 'model/particle_primary_momentum/PROOF.md'),
            'joint_refinement': ('Filters', 'model/particle_filtered_probe/PROOF.md'),
            'relativistic_field': ('Rayleigh', 'model/particle_primary_variational_readout/PROOF.md'),
            'particle_observables': ('Channels', 'model/particle_primary_channel_labels/PROOF.md'),
        }.get(node['id'])
        if fourth: c.link(x+width*.9, top+660, fourth[0], fourth[1], 49)
    gravity = nodes[9]; top = y+2720
    c.rect(M, top, W-2*M-70, 530, '#fffdf7', LINE, 26)
    c.text(M+90, top+110, 'DEFERRED · PHYSICAL PROOF REQUIRED · after 5', 57, ORANGE, 700)
    c.text(M+90, top+230, '10  '+gravity['title'], 72, INK, 700)
    c.para(M+90, top+345, W-2*M-250, SUMMARIES[gravity['id']], 64, ORANGE, max_y=top+420)
    c.link(M+90, top+480, 'Actual source energy, complete records and marginal entropy',
           gravity['evidence'][0]['path'], 49)
    c.link(M+14700, top+480, 'Conditional Einstein theorem and proper-geometry hypotheses',
           'model/matched_einstein_limit/PROOF.md', 49)
    c.para(M, y+3360, W-2*M,
           'Conditional candidate: a declared chiral gauge-Higgs action supplies calibrated tree/running coefficients, with a separately calibrated neutrino extension. These coefficients are not physical pole masses, full radiative matching or independent predictions. The massless photon excludes a positive global vacuum gap for the full theory.',
           60, BLUE, max_y=y+3470)
    c.link(M, y+3540, 'Canonical dependency graph and acceptance conditions',
           'model/particle_joint_reconstruction/DEPENDENCIES.md', 51)
    c.link(M+14700, y+3540, 'Conditional candidate, spectrum conventions and calibration ancestry',
           graph['conditional_candidate']['proof'], 51)
    c.items.append('</g>')
    return y+3630
