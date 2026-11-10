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
        'Five coherent qutrit meters attach native increments to G1/H trions. Charged regional G1 factors retain the opposite endpoint and Wilson link outside. Radix-five charge and radix-three spin counters record Q6=±6 and source spin increments ±1/2 without initial readout; each of four complete spin rows has original energy-band weight >1/16. Terminal spin/charge readout is a separate disturbing operation. Finite Fourier labels are not exact momentum projections.',
    'joint_refinement':
        'Filtered and unfiltered probes share positive time dynamics. The relative-observation sector has commuting H≥0 and P, HΩ=PΩ=0; sinc² energy/momentum weight may vanish. Complete selector sums preserve the original energy law. A finite blank centered at 9/a retains weight >1/80 at energy scale 1/a. Odd weight escapes every fixed upper energy band; fixed-band survival and spatial momentum require proof.',
    'relativistic_field':
        'A true finite gauge transition is isolated near 63 source units in the empty-matter reducing sector. The exact A1u level has transition residue >4/49 and a positive original atom. Its rare-sector probability is retained and vanishes on the selected refinement. The global ground is below −723; this gauge-sector origin cannot replace it. Odd trion matrix units reduce forced spin degeneracy, but their isolated branch is unproved. No relativistic dispersion or GeV mass follows.',
    'chiral_species_attachment':
        'Neutral trions give G1g ⊕ Hg and conjugate G1u ⊕ Hu; disjoint-site fermionic exchange is exact. Charged G1 endpoints retain an outside compensator; the even pair has A1 ⊕ T1 coefficients. Parity pairs native Γ± with the same gauge representation and identical reflected records. An equivariant field/current dictionary needs explicit breaking or changed attachment for a one-handed sector. Native Γ5 is not established Lorentz chirality.',
    'particle_observables':
        'Each complete Walsh spatial row has original positive-energy weight >1/24 in [0.25,17.75] source units; G1/H parts exceed 5/384 and 11/384. All 96 successes and 8 failures remain. A separately chosen normalized phase has weight >1/3. Four uncomputed spectral measures determine each summed multiplet’s spatial response. These finite labels are not total-spin or conserved-momentum eigenvalues. Dispersion, a rest branch or threshold and a GeV mass require proof.',
    'independent_comparison':
        'SI units may be declared immediately. Finite mass comparison needs a resolved branch or threshold, momentum and initial-momentum control, and dispersion/readout errors. A continuum mass additionally needs uniform fixed-band survival. The primary physical-mass table in GeV is empty; calibrated inputs keep their ancestry.',
    'pure_gauge_gap':
        'Instantiate nodes 1–5 for an actual pure compact-simple gauge source or a proved source map. Dense local Euclidean decay plus finite-energy noncollapse implies a positive finite gap. Its physical hypotheses are conditional; it is separate from the photon-containing theory.',
    'geometric_gravity':
        "The actual corner energy A=h′₀ generates a finite control pulse. The net incoming energy current B=i[H,A] has a strictly positive original susceptibility and a complete binary record contrast; source-clock duration and original rare weights are retained. Controlled energy balance includes applied work. A source-relative derivative gives integrated energy, without assuming a proper metric. Boundary information obeys J+C_center≤2β‖V_cut‖+r with the finite polynomial correction r. Proper geometry, Lorentz stress, a finite universal area coefficient and Einstein dynamics require proof.",
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
                   'relativistic_field': 'model/particle_primary_reduced_sector/PROOF.md',
                   'chiral_species_attachment': 'model/particle_primary_odd_sum_rules/PROOF.md',
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
                ('Full sector resolvent', 'model/particle_primary_odd_resolvent/PROOF.md'),
                ('Actual odd transport', 'model/particle_primary_composite_transport/PROOF.md')),
            'chiral_species_attachment': (
                ('Literal conjugation', 'model/particle_primary_charge_conjugation/PROOF.md'),
                ('Actual trion multiplets', 'model/particle_primary_composite_spin/PROOF.md'),
                ('Populated odd response', 'model/particle_primary_selected_odd_response/PROOF.md'),
                ('Charged regional spinor', 'model/particle_primary_regional_fermion/PROOF.md')),
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
    c.link(M+17000, top+420, 'Boundary interaction and information bound',
           'model/particle_primary_boundary_correlation/PROOF.md', 45)
    # Four original-state results have different observable domains. Their own
    # cards prevent a density wave or an energy moment from being read as mass.
    c.text(M, y+3350, 'EXACT FINITE SOURCE RECORDS: FOUR DISTINCT TESTS',
           76, TEAL, 700, serif=True)
    card_width=(W-2*M)/4
    cards=(
        ('An isolated interacting source line',
         'ΔE ∈ [63−59/1351350, 63+1/45675]',
         'In the actual empty-matter sector, the spatial-odd A1u gauge level is simple within its symmetry sector. Its Wilson-loop transition residue exceeds 4/49. Complete clocks retain original weight >c_g/13; the rare population c_g is not divided out. This is an actual spectral atom at interaction one. Its sector weight vanishes under the selected refinement, so the line is neither a surviving particle nor a GeV mass. The global vacuum is not the empty-matter vacuum.',
         'particle_primary_reduced_sector'),
        ('A charged spinor with its compensator',
         'Q6=±6 • source spin increment ±1/2',
         'One odd G1 factor lies inside the region; the opposite endpoint and Wilson link remain outside. Coherent counters record the probe’s charge and chosen-axis spin increments without reading initial values. Each of four raw spin rows has original pair-energy weight >1/16. A separate A1/T1 coefficient readout has weights >5/64 and >15/64. The complete pair is even. Pair transfer energy cannot be divided by two to obtain an endpoint mass; no electron or Lorentz spin is identified.',
         'particle_primary_regional_fermion'),
        ('Selected energy tails grow more slowly',
         'M₂ ≤ 4,000,000/a • no cutoff-degree factor',
         'The unchanged polynomial state has electric energy ≤91V/a. Native derivatives and the signed third-moment identity bound the original second moment. Complete positive-energy windows have a floor of order a and an upper endpoint of order 1/a. This does not give a nonzero fixed-energy limit. A different complete parity apparatus on this same source exactly matches the critical mean work, yet its nonelastic probability tends to zero. Critical-phase survival needs its actual current correlations.',
         'particle_primary_weighted_spectral_tail'),
        ('A local pulse changes energy-current records',
         'A=h′₀ • B=i[H,A] • applied work retained',
         'An exact finite Cayley pulse of H+A gives nonzero original contrast in the complete two-outcome current detector. A physical one-word witness proves susceptibility on the selected state; its tiny source-density factor remains in the bound. Controlled energy conservation includes the work supplied by the pulse. The source-relative derivative gives integrated energy. Identifying the control as a proper metric or obtaining Einstein dynamics requires additional geometry.',
         'particle_primary_energy_response'),
    )
    for i,(title,equation,body,owner) in enumerate(cards):
        x=M+i*card_width;top=y+3420
        c.rect(x,top,card_width-70,620,'#edf5f1',LINE,26)
        c.text(x+80,top+100,title,58,INK,700)
        c.text(x+80,top+210,equation,60,TEAL,700)
        c.para(x+80,top+315,card_width-230,body,55,INK,max_y=top+525)
        c.link(x+80,top+575,'Exact proof, source and complete records',
               'model/'+owner+'/PROOF.md',46)
    c.link(M+card_width*.53,y+3995,'Original creation bands',
           'model/particle_primary_original_branch/PROOF.md',42)
    c.link(M+card_width*1.53,y+3995,'Even density contrast',
           'model/particle_primary_original_shell_moments/PROOF.md',42)
    c.link(M+card_width*2.53,y+3995,'Controlled gradient work',
           'model/particle_primary_critical_refinement/PROOF.md',42)
    c.link(M+card_width*3.53,y+3995,'Boundary information identity',
           'model/particle_primary_boundary_information/PROOF.md',42)
    c.link(M+card_width*.74,y+3995,'Positive window',
           'model/particle_positive_window/PROOF.md',42)
    c.link(M+card_width*.31,y+4075,'Odd branch and exact matrix units',
           'model/particle_primary_branch_certificate/PROOF.md',42)
    c.link(M+card_width*1.31,y+4075,'Native parity and charged-sector scope',
           'model/particle_primary_odd_sum_rules/PROOF.md',42)
    c.link(M+card_width*2.31,y+4075,'Polynomial fluctuation and survival boundary',
           'model/particle_primary_work_spectrum/PROOF.md',42)
    c.link(M+card_width*2.70,y+4075,'Equal work ≠ spectral survival',
           'model/particle_primary_equal_work_control/PROOF.md',42)
    c.para(M, y+4160, W-2*M,
           'Conditional candidate: a declared chiral gauge-Higgs action supplies calibrated tree/running coefficients, with a separately calibrated neutrino extension. These coefficients are not physical pole masses, full radiative matching or independent predictions. The massless photon excludes a positive global vacuum gap for the full theory.',
           60, BLUE, max_y=y+4270)
    c.link(M, y+4340, 'Canonical dependency graph and acceptance conditions',
           'model/particle_joint_reconstruction/DEPENDENCIES.md', 51)
    c.link(M+14700, y+4340, 'Conditional candidate, spectrum conventions and calibration ancestry',
           graph['conditional_candidate']['proof'], 51)
    c.items.append('</g>')
    return y+4430


def render_candidate_wings(c,data):
    """Use the outer panels' unused lower space for a scoped named dictionary."""
    from fractions import Fraction as F
    from .drawing import INK, MUTED, TEAL, ORANGE, BLUE, RED, LINE
    value=data['candidate_correspondence'];targets={r['target']:r for r in value['target_rows']}
    window=data['candidate_identity_tests']['calculated_source_results']['charged']['original_compensated_dressed_energy_window']
    dressed_weight=window['original_dressed_recorded_simple_strict_lower']
    dressed_lo,dressed_hi=window['reported_transfer_interval']
    left=360;right=38360;y=30400;w=9120
    c.items.append('<g id="candidate-particle-map">')
    for x,title,subtitle in (
        (left,'BEST-GUESS PARTICLE MAP','EXACT SOURCE LABELS → TESTABLE FAMILY HYPOTHESES'),
        (right,'WHAT THE CALCULATIONS SAY','SOURCE ENERGIES, COMPENSATORS AND COLLISIONS')):
        c.rect(x,y,w,16400,'#fffdf7',LINE,35)
        c.para(x+120,y+220,w-240,title,132,INK,700,max_y=y+480)
        c.para(x+120,y+650,w-240,subtitle,74,TEAL,700,max_y=y+890)
    c.para(left+140,y+1080,w-280,
        'Gold names are physical comparison targets. Green source facts are exact in their declared domain. Orange arrows are hypotheses: a matching charge or finite spin label does not identify a relativistic particle. No fitted GeV mass is added here.',
        76,INK,max_y=y+1590)
    groups=(
        (('electron','muon','tau'),'Electron • muon • tau',
         'Exact lepton-doublet and conjugate-singlet weights. A source Higgs-conjugate scalar composite gives a dressed charged singlet. The charged regional G1 factor retains its outside compensator. Flavor and mass are unproved.',
         'charged_lepton_left + charged_lepton_right_conjugate'),
        (('electron_neutrino','muon_neutrino','tau_neutrino'),'Electron • muon • tau neutrinos',
         'Exact neutral member of the native lepton doublet. Flavor labels are not mass eigenstates. No neutrino mass or Dirac/Majorana choice follows.',
         'neutrino_left'),
        (('up','charm','top'),'Up • charm • top quarks',
         'The color-triplet weak-doublet and conjugate weak-singlet weights fit one up-type family pattern. Three names do not supply three source flavors.',
         'up_left + up_right_conjugate'),
        (('down','strange','bottom'),'Down • strange • bottom quarks',
         'The native color triplet and conjugate component weights fit one down-type family pattern. Their propagating chiral quark-field dictionary is unproved.',
         'down_left + down_right_conjugate'),
        (('photon',),'Photon',
         'Q=T3+Y is the candidate unbroken direction under a neutral Higgs vacuum. Actual polar and axial loop channels have nonzero source transitions. Their transverse photon interpretation and masslessness require proof.',
         'Q=T3+Y'),
        (('gluon',),'Gluons',
         'Eight color-current directions have the adjoint pattern. This does not turn eight rows or a neutral Wilson loop into eight observed free particles.',
         'color_current_1..8'),
        (('W_plus','W_minus'),'W⁺ and W⁻ bosons',
         'The native weak raising and lowering currents have opposite charges. Their massive vector branches need the selected vacuum and physical readout.',
         'T_plus + T_minus'),
        (('Z',),'Z boson',
         'The neutral weak/hypercharge span supports a Z-type hypothesis. Its mixing angle, massive field and detector identity are separate calculations.',
         'span(T3,Y)'),
        (('Higgs',),'Higgs boson',
         'A same-source fermion composite has Higgs-conjugate doublet charges and scalar source spin. An A1g gauge-loop scalar also couples. Their common vacuum, radial Higgs response and mass require proof.',
         'declared_Higgs_doublet_radial_direction'),
    )
    for i,(ids,title,body,source) in enumerate(groups):
        top=y+1780+i*1170;example=targets[ids[0]]['benchmark']
        charge='±1' if len(ids)==2 else example['charge']
        c.rect(left+120,top,w-240,1080,'#f8f1dc',LINE,25)
        c.text(left+270,top+150,title,94,INK,700)
        c.text(left+270,top+300,'BENCHMARK: spin '+example['spin']+' • Q='+charge+' • color '+example['color'],70,BLUE,700)
        c.para(left+270,top+460,w-560,body,76,TEAL,max_y=top+760)
        c.para(left+270,top+870,w-560,'Source ID: '+source,60,MUTED,max_y=top+1000)
    c.para(left+150,y+12490,w-300,
        'POSSIBLE FAMILY CORRESPONDENCE, NOT IDENTIFICATION. The same source IDs recur under several generation names. Exact gauge weights are internal labels; source G1 is a finite spinor, not a proved Lorentz-spin assignment. Antiparticle and flavor readouts need the same physical field dictionary.',
        76,ORANGE,700,max_y=y+13020)
    c.text(left+150,y+13320,'EXACT INTERNAL GAUGE WEIGHTS',92,INK,700)
    c.para(left+150,y+13540,w-300,
        'All-left-Weyl comparison convention: Y=Q6/6 and Q=T3+Y. These internal labels match one fermion-family pattern. Physical chirality and Lorentz spin require their own attachment.',
        75,MUTED,max_y=y+13930)
    c.text(left+170,y+14100,'COMPONENT',70,BLUE,700)
    c.text(left+2300,y+14100,'(SU3, SU2) with Y',70,BLUE,700)
    c.text(left+5800,y+14100,'Q6',70,BLUE,700)
    c.text(left+6950,y+14100,'Q components',70,BLUE,700)
    for i,row in enumerate(value['internal_gauge_multiplets']):
        yy=y+14400+i*230
        c.text(left+170,yy,row['comparison_label'],82,INK,700)
        color='3̅' if row['SU3']=='3bar' else row['SU3']
        c.text(left+2300,yy,'('+color+', '+str(row['SU2_dimension'])+')   Y='+row['Y'],82,TEAL,700)
        c.text(left+5800,yy,row['Q6'],82,TEAL)
        c.text(left+6950,yy,', '.join(row['component_Q']),82,TEAL)
    c.para(left+170,y+15570,w-340,
        'Eleven protected-object IDs retain their original source domains. They are separate from these internal weights and actual excitation IDs; none has a proved named-particle map.',
        70,ORANGE,max_y=y+15770)
    c.link(left+170,y+15900,'Complete 18-target map, source IDs and exclusions',
           'model/standard_model_source_map/candidates/CANDIDATES.md',72)
    c.link(left+170,y+16100,'Machine-readable map and exact rational energy comparisons',
           'model/standard_model_source_map/candidates/correspondence.json',66)
    c.link(left+170,y+16300,'Candidate-by-candidate tests: source facts and physical requirements',
           'model/particle_candidate_confirmation/RESULTS.md',66)

    # Regional spinor: draw both endpoints so the compensator cannot disappear.
    top=y+1050
    c.text(right+160,top+150,'ONE CHARGED FACTOR, ONE COMPLETE PAIR',94,INK,700)
    c.rect(right+260,top+410,3700,1080,'#e7f1ec',TEAL,35)
    c.rect(right+5100,top+410,3700,1080,'#e9eef7',BLUE,35)
    c.circle(right+2100,top+850,220,'#d4ece5',TEAL,8)
    c.circle(right+6950,top+850,220,'#d9e4f6',BLUE,8)
    c.line(right+2350,top+850,right+6700,top+850,ORANGE,16)
    c.text(right+4520,top+690,'Wilson link',72,ORANGE,700,anchor='middle')
    c.text(right+2100,top+870,'+6',92,TEAL,700,anchor='middle')
    c.text(right+6950,top+870,'−6',92,BLUE,700,anchor='middle')
    c.text(right+2100,top+1160,'ΔQ6',70,TEAL,700,anchor='middle')
    c.text(right+6950,top+1160,'ΔQ6',70,BLUE,700,anchor='middle')
    c.text(right+2100,top+1270,'REGION • odd G1',84,TEAL,700,anchor='middle')
    c.text(right+6950,top+1270,'OUTSIDE • compensator',78,BLUE,700,anchor='middle')
    c.para(right+160,top+1750,w-320,
        'Coherent counters attach Q6 and source-spin increments to original energy records. The full apparatus is neutral and even, with A1/T1 pair coefficients. A single endpoint is a charged-lepton-singlet-type candidate; the pair is not automatically a photon. Its transfer energy includes the link and cannot be divided by two to obtain a lepton mass. '
        f'A scalar-dressed pair also has original probability >{dressed_weight} in [{dressed_lo},{dressed_hi}] source-energy units on the same clock. All {window["complete_success_outcomes"]} success records and failure are retained. Terminal one-particle occupation is not guaranteed.',
        82,INK,max_y=top+2420)
    c.link(right+160,top+2690,'Actual charged regional source, spin and complete energy records',
           'model/particle_primary_regional_fermion/PROOF.md',68)
    c.link(right+160,top+2860,'Source scalar dressing and the original positive pair-energy window',
           'model/particle_candidate_charged/PROOF.md',66)
    top=y+4100
    c.text(right+160,top+150,'NEUTRAL COMPOSITE HYPOTHESES',94,INK,700)
    c.para(right+160,top+420,w-320,
        'LL eᶜ: neutral gauge-singlet odd trion. Its source counter increment is +1, with −1 for its conjugate. A neutral-lepton-composite analogy is plausible; a direct unbroken weak-current map to an active neutrino fails.',
        82,TEAL,max_y=top+940)
    c.para(right+160,top+1170,w-320,
        'uᶜ(dᶜ)²: neutral color-singlet odd trion with source counter increment +3; its conjugate has −3. The even composite gives ±2. The counter reads S=N₋₄+N₂+N₆: three native occupation counts, not electric charge. It separates six interventions and failure; multiple excitations can alias. Physical spin and named identities are unproved.',
        82,TEAL,max_y=top+1810)
    c.link(right+160,top+2080,'Exact conserved-charge counter, gauge currents and selection rules',
           'model/particle_candidate_composites/PROOF.md',68)
    top=y+6500
    c.rect(right+120,top,w-240,2500,'#fbefeb',RED,30)
    c.text(right+280,top+190,'A DIRECT LEPTON-RATIO GUESS FAILS',94,RED,700)
    c.para(right+280,top+440,w-560,
        'Actual neutral-sector thresholds above one common ground, in source units. These are rigorous intervals, not point masses:',
        80,INK,max_y=top+820)
    ratios=value['finite_sector_ratio_control']['result']
    def outward(q,upper=False):
        q=F(q)*10000;n=-((-q.numerator)//q.denominator) if upper else q.numerator//q.denominator
        return f'{n/10000:.4f}'
    for i,row in enumerate(ratios['source_sector_rows']):
        lo,hi=row['sector_bottom_excitation_interval'];yy=top+1040+i*270
        label={'weak_pair_charge_six':'Weak trion','color_triple':'Color trion','color_determinant_charge_six_hole':'Even four-body composite'}[row['source_id']]
        c.text(right+280,yy,label,76,INK,700)
        c.text(right+4000,yy,'['+outward(lo)+', '+outward(hi,True)+']',80,TEAL,700)
        c.line(right+6700+float(F(lo))*100,yy-25,right+6700+float(F(hi))*100,yy-25,TEAL,20)
    c.para(right+280,top+1900,w-560,
        'Every family ratio <21.16; odd-family ratios <15.31. Measured μ/e =206.7682827(46). One common unit conversion cancels. This excludes only that member-zero assignment, not other excitations, capacities or effective dictionaries.',
        82,RED,700,max_y=top+2320)
    c.link(right+250,top+2650,'Exact shared-ground ratio bounds and CODATA comparison',
           'model/standard_model_source_map/candidates/FINITE_SECTOR_RATIOS.md',68)
    top=y+9530
    c.text(right+160,top+150,'STATIC FUSION IS COMPUTED',94,INK,700)
    c.para(right+160,top+410,w-320,
        'Adjacent outward square loops cancel their common edge: two four-edge words become one six-edge loop. Their interacting 44-level cluster lies in [1311/14,1335/14], about 93.64–95.36 source units. Opposite loops retain eight edges and a separate [124,128] window.',
        82,TEAL,max_y=top+1050)
    for j,label in enumerate(('□ + □','shared edge cancels','6-edge word')):
        xx=right+200+j*2940;c.rect(xx,top+1220,2700,550,'#e5efea',TEAL,25)
        c.text(xx+1350,top+1550,label,76,TEAL,700,anchor='middle')
        if j<2:c.line(xx+2730,top+1500,xx+2880,top+1500,TEAL,6,arrow=True)
    c.para(right+160,top+2020,w-320,
        'Each complete pair row retains original weight >c_g/40. The rare sector factor remains. No moving incoming packets or physical binding mass follow from this static energy contrast.',
        80,ORANGE,max_y=top+2400)
    c.link(right+160,top+2590,'Full source loop fusion, omitted spectrum and complete clocks',
           'model/particle_primary_loop_fusion/PROOF.md',68)
    top=y+12450
    c.text(right+160,top+150,'EXECUTED INTERACTION, SEPARATE SOURCE',88,INK,700)
    c.para(right+160,top+410,w-320,
        'The supplied Weyl/contact extension executes a four-mode echo: 22,531 native pulses at c=1/8, U=1, t=4. Its numerical signal is about 0.046915; the ideal-word bound is >0.0361998125 versus free 0. All failures and the reference remain. Floating and hardware errors are separate from the rational bound.',
        80,TEAL,max_y=top+1020)
    c.para(right+160,top+1260,w-320,
        'The moving-packet route proves a complete finite observation theorem and executes 32-row five-qubit controls. Its fine-grid incoming-packet program is not executed. Neither route identifies CERN particles or primary-source cross sections.',
        80,ORANGE,max_y=top+1760)
    c.link(right+160,top+1950,'Executed contact echo', 'model/contact_scattering/README.md',68)
    c.link(right+4150,top+1950,'Moving packets and native controls', 'model/scattering_packet_readout/PROOF.md',68)
    top=y+14880
    c.text(right+160,top+150,'WHAT A DEFENSIBLE GeV MASS REQUIRES',90,INK,700)
    c.para(right+160,top+420,w-320,
        'Attach physical currents and spin; resolve the correct charged threshold, resonance or branch; control spatial momentum and dispersion on the same state and clock; quantify finite errors and capacity dependence; then attach the declared energy unit. SI conventions do not supply the missing dimensionless spectrum.',
        80,INK,max_y=top+1050)
    c.link(right+160,top+1280,'Physical mass-readout criterion and remaining source calculations',
           'docs/research/OPH_PARTICLE_MASS_RESPONSE.md#physical-mass-readout-criterion',66)
    c.items.append('</g>')
