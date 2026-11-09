"""Six output obligations, with finite evidence and physical domains separated."""


def render(c, data, y):
    from .drawing import W, M, INK, TEAL, ORANGE, BLUE
    report = data['observer_outputs']
    attachment = report['defect_field_attachment']
    locality = report['field_causality']
    gravity = report['energy_and_gravity_response']['gravity_interface']
    if attachment['physical_particle_assignments'] or locality['physical_particle_assignments']:
        raise ValueError('physical species requires a separate atlas review')
    if gravity['Einstein_dynamics_or_physical_Newton_constant_derived']:
        raise ValueError('geometric scope requires a separate atlas review')
    c.items.append('<g id="six-output-obligations">')
    c.text(M, y, 'SIX OUTPUT OBLIGATIONS: WHAT THE PROGRAM ACTUALLY SUPPLIES', 80, INK, 700, serif=True)
    c.para(M, y+115, W-2*M,
           'An architecture may be chosen. Its observable consequences must be calculated. Green text describes established finite outputs; orange text names the extra physical statement. Results on distinct prepared programs are not combined into a fictitious common vacuum.',
           47, INK, max_y=y+280)
    cards = [
        ('01  PREPARATION AND DETECTORS', 'state → interaction → complete record',
         'One declared Wilson source has a physical Gauss preparation with its conjugate reference retained. Actual loop states prove 11 nonzero regional charge sectors. The first finite member has a complete 1,331-outcome boundary-charge detector; all outcomes are kept.',
         'Finite preparation including boundary charges is constructed. Physical target: compatible states and complete records under joint refinement. The enormous full spatial matrix and preparation tape are not claimed as executed.',
         'model/particle_physical_boundary_state/PROOF.md'),
        ('02  POSITION, MOMENTUM AND CLOCK', 'same source + same state + same clock',
         'Joint translations move matter, links and the compensating reference. Their complete Fourier instrument measures conserved crystal momentum. A position detector samples an occupied mode’s site and retains vacuum; reading it changes the momentum law.',
         'Physical target: one particle’s momentum, physical spin and Lorentz clocks for the same interacting excitations. Whole-field crystal momentum and finite uncertainty do not fix a GeV mass.',
         'model/particle_joint_translation_readout/PROOF.md'),
        ('03  PROPAGATION AND RELATIVITY', 'local intervention → bounded distant change',
         'Smooth-field local records have uniform approximation bounds. With specified coupling scaling, the physical causal bound is also uniform as spacing shrinks, including growing detector supports and every occupation. Ordinary Wilson supplies a sharper allowed example.',
         'Physical target: a convergent interacting state and observable family, Lorentz symmetry and chiral matter. A causal cone alone supplies none of these identifications. Exact-sign equivalence is optional.',
         'model/particle_physical_causality/PROOF.md'),
        ('04  MASSES AND INTERACTIONS', 'source spectrum → spectral readout → matching',
         'The Q* = 6 full ground has 12 complete charged channels. First-family doublet removal and addition begin near 0.0710 and 0.4873 source-energy units. A common sector offset cannot equalize their complete spectral shapes.',
         'Physical target: mass shells, charged dressing, scattering and decay laws on the common spatial state. These finite transition lines are not electron or positron masses.',
         'model/particle_ground_charged/PROOF.md'),
        ('05  DEFECTS BECOME DISTINGUISHABLE', 'protected Q → controlled vertex → detector',
         'Declare H_Q = (3I − Q)V / 2. At one tick the scalar occupied probability is 1/4 for Q = +1 and 3/4 for Q = −1. Actual memory or regional holonomy can control this complete field response.',
         'Physical target: propagation, spin, statistics and currents identifying a species. Strengths 1 and 2 are declared. Protected motion conserves the coupled action; birth needs controller-work accounting.',
         'model/particle_defect_attachment/PROOF.md'),
        ('06  ERRORS, ENERGY AND GRAVITY', 'H = Σhₓ → Jₓᵧ = i[hₓ,hᵧ] → exact flux',
         'Same-action currents obey energy continuity. A declared coframe supplies six spatial strain operators; the shared thermal flux control has complete work and pressure records. The ten-probe Einstein decoder keeps its geometric premises explicit.',
         'Physical target: proper geometry and clocks, Lorentz stress, regional modular/area responses and controlled same-state small-ball errors. A declared mechanical volume does not supply a cosmological equation of state.',
         'model/particle_geometric_response/PROOF.md')]
    width = (W-2*M)/3
    for index, (title, chain, result, boundary, proof) in enumerate(cards):
        x = M+(index % 3)*width
        top = y+340+(index//3)*1550
        c.rect(x, top, width-70, 1460, '#fffdf7', radius=26)
        c.para(x+80, top+130, width-230, title, 65, INK, 700, max_y=top+280)
        c.rect(x+80, top+310, width-230, 240, '#e3efeb', radius=20)
        c.para(x+140, top+410, width-350, chain, 52, TEAL, 700, max_y=top+520)
        c.para(x+80, top+660, width-230, result, 54, TEAL, max_y=top+975)
        c.para(x+80, top+1080, width-230, boundary, 49, ORANGE, max_y=top+1330)
        c.link(x+80, top+1400, 'Exact operators, complete records and proof domain', proof, 43)
    c.link(M, y+3530, 'The common dependency map and all six requirements',
           'docs/research/OPH_FIELD_RECONSTRUCTION.md#six-output-obligations-evidence-and-exact-boundaries', 46)
    c.text(M+15000, y+3530, 'Finite outputs are derived; physical identities need their stated bridges.', 45, BLUE)
    c.items.append('</g>')
    return y+3630
