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
        'Five coherent qutrit meters attach native occupation increments to G1/H records without reading initial numbers; all probe, pointer and failure labels remain. Separate radix-five coherent meters attach opposite Q6=±6 transfers with original joint success weight >1/4 without reading initial absolute charges. Absolute terminal charge is a separate disturbing readout. Coherent offsets give finite Fourier labels, not exact momentum projections.',
    'joint_refinement':
        'Filtered and unfiltered probes share positive time dynamics. The relative-observation sector has commuting H≥0 and P, HΩ=PΩ=0; sinc² energy/momentum weight may vanish. Complete selector sums preserve the original energy law. A finite blank centered at 9/a retains weight >1/80 at energy scale 1/a. Odd weight escapes every fixed upper energy band; fixed-band survival and spatial momentum require proof.',
    'relativistic_field':
        'At coefficient one, exact fourth moments bound the full sector complement resolvent. A prepared odd composite has low-energy weight >1/26; the complete clock retains >1/27. Separately, one original-state creation row followed only by the final energy clock has weight >c/864 at member zero, with ρ≥cP_G. This tiny floor can collapse with capacity. Physical memories retain displaced source/reference data. No isolated pole, uniform residue or GeV rest mass follows.',
    'chiral_species_attachment':
        'A weak doublet pair plus a singlet, or a color triple, gives six odd gauge-singlet states: G1 ⊕ H. Uniform spatial preparation gives G1g ⊕ Hg; conjugates give G1u ⊕ Hu. Native-current increments and disjoint-site fermionic exchange are exact. Opposite native occupations have identical spectra; odd inversion eigenstates have opposite parity under C. These finite cube labels establish no Lorentz spin or named species. Chiral propagation and interactions require proof.',
    'particle_observables':
        'Each complete Walsh spatial row has original positive-energy weight >1/24 in [0.25,17.75] source units; G1/H parts exceed 5/384 and 11/384. All 96 successes and 8 failures remain. A separately chosen normalized phase has weight >1/3. Four uncomputed spectral measures determine each summed multiplet’s spatial response. These finite labels are not total-spin or conserved-momentum eigenvalues. Relativistic dispersion, an isolated rest branch and a GeV mass require proof.',
    'independent_comparison':
        'SI units may be declared immediately. Finite mass comparison needs a resolved branch, momentum and initial-momentum control, and dispersion/readout errors. A continuum mass additionally needs uniform fixed-band survival. The primary physical-mass table in GeV is empty; calibrated inputs keep their ancestry.',
    'pure_gauge_gap':
        'Instantiate nodes 1–5 for an actual pure compact-simple gauge source or a proved source map. Dense local Euclidean decay plus finite-energy noncollapse implies a positive finite gap. Its physical hypotheses are conditional; it is separate from the photon-containing theory.',
    'geometric_gravity':
        "A region contains its matter and internal links; crossing links remain outside. This cut has a cube-covariant state and modular flow. Entropy from boundary representation dimensions obeys 0≤S_edge≤ℓ b log 6 for word degree ℓ and b crossing links, with a strictly positive finite Wilson-loop witness. Complete center recording preserves regional entropy and can disturb global coherences. Contained regional energy differs from H'_R by an explicit boundary operator. Proper geometry, Lorentz stress, a finite universal area coefficient and Einstein dynamics require proof.",
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
        theorem = {'local_records': 'model/particle_primary_selected_spin_response/PROOF.md',
                   'joint_refinement': 'model/particle_primary_odd_low_energy/PROOF.md',
                   'relativistic_field': 'model/particle_primary_odd_resolvent/PROOF.md',
                   'chiral_species_attachment': 'model/particle_primary_composite_spin/PROOF.md',
                   'particle_observables': 'model/particle_primary_original_branch/PROOF.md',
                   'independent_comparison': 'model/particle_primary_mass_readout/PROOF.md'}.get(
                       node['id'], node['evidence'][0]['path'])
        c.link(x+width/2, top+660, 'Scoped construction or theorem',
               theorem, 49)
        extra = {
            'fixed_program': ('Source capacity', 'model/particle_primary_capacity/PROOF.md'),
            'local_records': ('Local influence', 'model/particle_primary_locality/PROOF.md'),
            'joint_refinement': ('Origin momentum', 'model/particle_primary_origin_momentum/PROOF.md'),
            'relativistic_field': ('Exact residual', 'model/particle_primary_residual_refinement/PROOF.md'),
            'particle_observables': ('Local window', 'model/particle_primary_response_moments/PROOF.md'),
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
                ('Complete energy clocks', 'model/particle_primary_joint_spectroscopy/PROOF.md'),
                ('Local charge transfers', 'model/particle_primary_charged_readout/PROOF.md')),
            'joint_refinement': (
                ('Spatial action', 'model/particle_primary_spatial_algebra/PROOF.md'),
                ('Occupation bound', 'model/particle_primary_refinement_response/PROOF.md'),
                ('Weight ≥3/100 in [1/(100β),46/a]', 'model/particle_primary_energy_tightness/PROOF.md'),
                ('Filtered algebra', 'model/particle_primary_field_algebra/PROOF.md'),
                ('Common time', 'model/particle_primary_current_tail/PROOF.md'),
                ('Observation strength', 'model/particle_primary_weighted_survival/PROOF.md')),
            'relativistic_field': (
                ('Selected occupation bound', 'model/particle_primary_occupation_bound/PROOF.md'),
                ('Distinct selected A1u response', 'model/particle_primary_odd_residual/PROOF.md'),
                ('Full-source Rayleigh bound', 'model/particle_primary_variational_readout/PROOF.md'),
                ('Sector floors', 'model/particle_primary_sector_threshold/PROOF.md'),
                ('Actual odd transport', 'model/particle_primary_composite_transport/PROOF.md')),
            'chiral_species_attachment': (
                ('Literal conjugation', 'model/particle_primary_charge_conjugation/PROOF.md'),
                ('Populated odd response', 'model/particle_primary_selected_odd_response/PROOF.md')),
            'particle_observables': (
                ('Local [4,8] response', 'model/particle_primary_branch_resolution/PROOF.md'),
                ('Charged endpoint band', 'model/particle_primary_charged_string/PROOF.md'),
                ('Checkerboard A2g band', 'model/particle_primary_staggered_response/PROOF.md'),
                ('Populated spatial scans', 'model/particle_primary_populated_momentum/PROOF.md'),
                ('Coherent creation + four kernels', 'model/particle_primary_coherent_odd_response/PROOF.md'),
                ('Four spectral measures', 'model/particle_primary_coherent_odd_response/KERNEL_PROOF.md')),
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
    c.link(M+25000, top+480, 'Induced cut and boundary entropy',
           'model/particle_primary_geometric_cut/PROOF.md', 49)
    c.link(M+9000, top+420, 'Covariant source energy and boundary improvement',
           'model/particle_primary_covariant_stress/PROOF.md', 45)
    # Four original-state results have different observable domains. Their own
    # cards prevent a density wave or an energy moment from being read as mass.
    c.text(M, y+3350, 'EXACT FINITE SOURCE RECORDS: FOUR DISTINCT TESTS',
           76, TEAL, 700, serif=True)
    card_width=(W-2*M)/4
    cards=(
        ('Odd creation and native labels',
         '>1/24 per complete Walsh row',
         'The unchanged polynomial state suppresses reverse negative-energy records. G1/H creation parts exceed 5/384 and 11/384; all 104 outcomes remain. Total cube and native-number detectors preserve the summed energy law. Every component has more than one exact transfer frequency, with variance ≥9c/4 and c=2⁻⁵⁸²⁴⁰⁷¹³³. This tiny floor neither locates a pole nor identifies a particle.',
         'particle_primary_original_branch'),
        ('Even density: real energy contrast',
         'First moments: 0 : 1 : 2 : 3',
         'Four cube classes have these exact source-energy moment ratios. The spatially constant density field is elastic; the other classes have positive inelastic weight. Complete finite clocks retain distinct mean records. These even density waves are different observables from the odd G1/H creators; the ratios are neither masses nor a relativistic dispersion law.',
         'particle_primary_original_shell_moments'),
        ('Spatial work: a controlled 3D limit',
         'Work ratios: 4 • 2 • 3/2',
         'Relative to the unit-width, amplitude-one cube: double the amplitude, all three widths, or just one width. The selected-state coefficient cancels from these limiting work ratios. The positive gradient-energy form has explicit lattice O(a²) and phase O(a⁴) errors. A limit of energy moments does not prove a fixed-energy particle residue or a four-dimensional field theory.',
         'particle_primary_critical_refinement'),
        ('Boundary records change the source',
         'Iambient = Ialgebra + 2Sedge + Ccenter',
         'Recording the region’s boundary sector preserves its marginal entropy but removes global center coherence. The actual action proves Ccenter >0 and positive recording work. The identity separates boundary representation entropy from accessible algebra information. Bounds shrink with capacity; no geometric area law, Newton constant or Einstein dynamics is inferred.',
         'particle_primary_boundary_information'),
    )
    for i,(title,equation,body,owner) in enumerate(cards):
        x=M+i*card_width;top=y+3420
        c.rect(x,top,card_width-70,620,'#edf5f1',LINE,26)
        c.text(x+80,top+100,title,58,INK,700)
        c.text(x+80,top+210,equation,60,TEAL,700)
        c.para(x+80,top+315,card_width-230,body,55,INK,max_y=top+525)
        c.link(x+80,top+575,'Exact proof, source and complete records',
               'model/'+owner+'/PROOF.md',46)
    c.link(M+card_width*.53,y+3995,'Positive-window theorem',
           'model/particle_positive_window/PROOF.md',42)
    c.para(M, y+4160, W-2*M,
           'Conditional candidate: a declared chiral gauge-Higgs action supplies calibrated tree/running coefficients, with a separately calibrated neutrino extension. These coefficients are not physical pole masses, full radiative matching or independent predictions. The massless photon excludes a positive global vacuum gap for the full theory.',
           60, BLUE, max_y=y+4270)
    c.link(M, y+4340, 'Canonical dependency graph and acceptance conditions',
           'model/particle_joint_reconstruction/DEPENDENCIES.md', 51)
    c.link(M+14700, y+4340, 'Conditional candidate, spectrum conventions and calibration ancestry',
           graph['conditional_candidate']['proof'], 51)
    c.items.append('</g>')
    return y+4430
