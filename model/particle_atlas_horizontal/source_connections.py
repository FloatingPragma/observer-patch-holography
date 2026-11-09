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
         'The diagonal scalar source has a unique full-carrier ground with Q* = 6 and complete occupation records. A separate spatial field experiment has one prepared marginal, one action and clock, and all 61 local occupation outcomes.',
         'Physical target: one spatial interacting state and detector dictionary. Closed Q* = 0 evolution cannot prepare Q* = 6. The occupation frequency 16.89261 differs from the full isolation gap > 0.0469.',
         'model/particle_scalar_ground/PROOF.md'),
        ('02  POSITION, MOMENTUM AND CLOCK', 'same source + same state + same clock',
         'Original two-time records determine a scalar gap. With the Q* = 0 prepared-sector ground, charged line separations recover family coefficients and the configured controller phase class while cancelling a sector offset.',
         'Physical target: operational spatial translations, physical spin and Lorentz clocks for those same interacting excitations. Program energy is not automatically GeV.',
         'model/particle_observer_outputs/BARE_EFFECTIVE_BOUNDARY.md'),
        ('03  PROPAGATION AND RELATIVITY', 'local intervention → bounded distant change',
         'One mixed gauge-state family has controlled local roughness and expected low Wilson-mode weight at every time and volume. Exact matching local preparations give a limit of complete polynomial detector records at fixed duration. Its occupation dynamics is nonzero.',
         'Physical target: interacting exact-sign response, chiral relativistic refinement and a physical state. Large fixed coefficients and the mixed reference are declared choices. A Wilson regulator mode is not a particle mass.',
         'model/particle_local_history_family/PROOF.md'),
        ('04  MASSES AND INTERACTIONS', 'source spectrum → spectral readout → matching',
         'The Q* = 6 full ground has 12 complete charged channels. First-family doublet removal and addition begin near 0.0710 and 0.4873 source-energy units. A common sector offset cannot equalize their complete spectral shapes.',
         'Physical target: mass shells, charged dressing, scattering and decay laws on the common spatial state. These finite transition lines are not electron or positron masses.',
         'model/particle_ground_charged/PROOF.md'),
        ('05  DEFECTS BECOME DISTINGUISHABLE', 'protected Q → controlled vertex → detector',
         'Declare H_Q = (3I − Q)V / 2. At one tick the scalar occupied probability is 1/4 for Q = +1 and 3/4 for Q = −1. Actual memory or regional holonomy can control this complete field response.',
         'Physical target: propagation, spin, statistics and currents identifying a species. Strengths 1 and 2 are declared. Protected motion conserves the coupled action; birth needs controller-work accounting.',
         'model/particle_defect_attachment/PROOF.md'),
        ('06  ERRORS, ENERGY AND GRAVITY', 'H = Σhₓ → Jₓᵧ = i[hₓ,hᵧ] → exact flux',
         'Same-action currents obey exact energy continuity. The mixed field family bounds complete local polynomial record errors by 2⁻²³ at time one, with enormous finite resources. An exact ten-probe tensor decoder applies when its same-state regional premises hold.',
         'Physical target: uniform refinement and capacity errors; geometric stress, proper clocks, regional modular/area responses and the Einstein theorem’s same-state small-ball premises.',
         'model/particle_gravity_response/PROOF.md')]
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
