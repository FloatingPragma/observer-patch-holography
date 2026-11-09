"""One dependency graph for the physical field and particle target."""


NODE_IDS = (
    'fixed_program', 'physical_preparation', 'local_records',
    'joint_refinement', 'relativistic_field', 'chiral_species_attachment',
    'particle_observables', 'independent_comparison', 'pure_gauge_gap',
)

# These are display summaries. Titles, dependencies, status and evidence come
# from the canonical dependency document rather than a second requirement list.
SUMMARIES = {
    'fixed_program':
        'One open Wilson program fixes the carrier, action, cutoff and energy zero. Every consequence uses its selected state, retained conjugate reference and Hamiltonian clock. A changed action declares an extension.',
    'physical_preparation':
        'The exact Gauss-polynomial state retains regional boundary charges. Wilson loops establish 11 sectors with positive weight. At degree one the complete charge detector retains all 1,331 outcomes; the enormous spatial preparation tape is not claimed as executed.',
    'local_records':
        'The selected physical source has complete local boundary records and an influence estimate uniform in link degree and assigned spacing. Its local output and full idle reference are retained. This controls propagation without asserting convergence or Lorentz symmetry.',
    'joint_refinement':
        'Selected cooling suppresses negative-energy tails. Fixed source-time averaging of uniformly bounded physical probes gives continuous subsequential responses, possibly zero. A common joint field limit and nonzero finite-energy weight require proof, with one scaling and renormalization prescription, observable maps and complete detector laws.',
    'relativistic_field':
        'Reconstruct three spatial dimensions and time on the same correlations and clock. Prove Lorentz covariance, positivity and locality. Establish an interaction incompatible with a free/Wick field description in a fixed field dictionary, or controlled nonfree scattering.',
    'chiral_species_attachment':
        'Identify physical currents, spin, statistics and families. Control mirror and extra modes, the Higgs vacuum and neutrino completion. Any claimed defect origin needs an actual field attachment. The declared chiral extension is distinct from the primary 60-mode carrier.',
    'particle_observables':
        'Calculate masses, mixing, production, scattering and decays with original spectral weights. Specify stable branches, charged infrared thresholds, resonances and composites. Radiative matching connects action coefficients to measured particle observables.',
    'independent_comparison':
        'Fix physical units and compare complete detector distributions with controlled errors. Keep the ancestry of each fitted mass, coupling and mixing input. Agreement with a calibration input does not provide an independent spectrum prediction.',
    'pure_gauge_gap':
        'Instantiate nodes 1–5 for an actual pure compact-simple gauge source or a proved source map. Dense local Euclidean decay plus finite-energy noncollapse implies a positive finite gap. Its physical hypotheses are conditional; it is separate from the photon-containing theory.',
}


def validate(graph):
    """Reject a changed graph or scope until its presentation is reviewed."""
    if graph.get('schema') != 'oph.joint_reconstruction_dependencies.v1':
        raise ValueError('unsupported joint dependency schema')
    nodes = graph['nodes']
    if tuple(node['id'] for node in nodes) != NODE_IDS:
        raise ValueError('the nine canonical dependencies must be retained')
    for index, node in enumerate(nodes):
        expected = [] if index == 0 else [NODE_IDS[4 if index == 8 else index-1]]
        status = ('finite_constructed' if index < 3 else
                  'conditional_criterion' if index == 8 else 'physical_proof_required')
        if (node['number'] != index+1 or node['depends_on'] != expected
                or node['status'] != status):
            raise ValueError('dependency ordering or physical scope requires review')
        if not node['acceptance'] or not node['evidence']:
            raise ValueError('dependency acceptance and evidence are required')
    branch = nodes[-1]
    if (branch['required_shared_nodes'] != list(NODE_IDS[:5])
            or branch['required_for_full_particle_spectrum']
            or branch['implies_full_theory_positive_global_gap']):
        raise ValueError('pure Yang-Mills must retain its separate source domain')
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
    c.text(M, y, 'ONE DEPENDENCY GRAPH: INTERACTING FIELDS AND THE FULL SPECTRUM',
           80, INK, 700, serif=True)
    c.para(M, y+115, W-2*M,
           'Nodes 1–3 are finite constructions on the selected program. Nodes 4–8 are physical proof conditions, assumed within the conditional candidate. Node 9 has its own pure-gauge source domain. Every calculation retains its declared state, reference and clock.',
           56, INK, max_y=y+240)
    chain = ' → '.join(str(node['number']) for node in nodes[:8])
    c.text(M, y+315, 'Shared particle route:  '+chain, 64, TEAL, 700)
    c.text(M+15300, y+315, 'Separate pure-gauge branch:  5 → 9', 64, BLUE, 700)
    width = (W-2*M)/3
    for index, node in enumerate(nodes):
        x = M+(index % 3)*width
        top = y+400+(index//3)*930
        color = TEAL if index < 3 else BLUE if index == 8 else ORANGE
        status = ('FINITE CONSTRUCTED' if index < 3 else
                  'CONDITIONAL GAP CRITERION' if index == 8 else 'PHYSICAL PROOF REQUIRED')
        after = ', '.join(str(numbers[key]) for key in node['depends_on'])
        tag = status+(' · after '+after if after else ' · declared program')
        c.rect(x, top, width-70, 870, '#fffdf7', LINE, 26)
        c.rect(x+70, top+60, width-210, 155,
               '#e3efeb' if index < 3 else '#e9edf6' if index == 8 else '#fbefdd', radius=20)
        c.text(x+120, top+165, tag, 57, color, 700)
        c.para(x+80, top+315, width-230,
               f"{node['number']:02d}  {node['title']}", 72, INK, 700, max_y=top+390)
        c.para(x+80, top+485, width-230, SUMMARIES[node['id']],
               67, color, max_y=top+700)
        c.link(x+80, top+800, 'Acceptance conditions and exact evidence',
               'model/particle_joint_reconstruction/DEPENDENCIES.md', 49)
        theorem = (next(row['path'] for row in node['evidence']
                        if row['path'].endswith('/POSITIVE_ENERGY.md'))
                   if node['id'] == 'joint_refinement' else node['evidence'][0]['path'])
        c.link(x+width/2, top+800, 'Scoped construction or theorem',
               theorem, 49)
    c.para(M, y+3260, W-2*M,
           'Conditional candidate: a declared chiral gauge-Higgs action supplies calibrated tree/running coefficients, with a separately calibrated neutrino extension. These coefficients are not physical pole masses, full radiative matching or independent predictions. The massless photon excludes a positive global vacuum gap for the full theory.',
           64, BLUE, max_y=y+3450)
    c.link(M, y+3540, 'Canonical dependency graph and acceptance conditions',
           'model/particle_joint_reconstruction/DEPENDENCIES.md', 51)
    c.link(M+14700, y+3540, 'Conditional candidate, spectrum conventions and calibration ancestry',
           graph['conditional_candidate']['proof'], 51)
    c.items.append('</g>')
    return y+3630
