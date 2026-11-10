"""Receipt-bound source families and physical candidate acceptance steps."""
from .canvas import AMBER, CARD2, EDGE, MUTED, SOFT, TEAL, TEXT, VIOLET, Canvas
from .parts import CW, X0, card


def catalogue(c, data, y):
    value = data['identity_tests']['labelled_catalogue']
    paths = {r['source_id']: r for r in data['identity_tests']['protected_particle_paths']['rows']}
    c.caps(X0, y+10, 'The source catalogue', TEAL, 14)
    c.text(X0, y+60, 'Protected objects have different kinds of labels', 40, TEXT, 700)
    y = c.para(X0, y+100, CW,
        'These eleven source entries are grouped by what their proofs protect. '
        'Only the spatial textures and kinetic nodes carry the stated winding. '
        'An internal charge or protected memory is a different object. The bare source '
        'entry is an exclusion control, not a defect.', 21, MUTED, 1.42)+12
    for family in value['protected_families']:
        body = Canvas()
        body.caps(X0+26, y+36, family['label'], TEAL, 13)
        yy = y+72
        for row in family['rows']:
            body.text(X0+26, yy, row['label'], 21, TEXT, 700)
            yy = body.para(X0+26, yy+27, CW-52, row['description'], 17, SOFT, 1.42)
            yy = body.para(X0+26, yy+5, CW-52, 'Source ID: '+row['source_id'], 13.5, MUTED, 1.35)+25
            yy = body.para(X0+26, yy, CW-52,
                'Required attachment: '+paths[row['source_id']]['first_required'], 16, AMBER, 1.42)+25
        h = yy-y
        card(c, X0, y, CW, h)
        c.items.extend(body.items)
        c.overflow.extend(body.overflow)
        y += h+18
    c.link(X0+26, y+20, 'All source domains, proofs and excitation IDs',
           '../../model/particle_candidate_confirmation/CATALOGUE.md', 17)
    c.link(X0+720, y+20, 'Per-object paths to a particle mass',
           '../../model/particle_candidate_confirmation/MASS_PATHS.md', 17)
    y += 70
    c.caps(X0, y, 'Source excitations and proposed species', VIOLET, 13)
    y = c.para(X0, y+40, CW,
        'The separate excitation table has sixteen source IDs. Internal matter components, '
        'compensated charged probes, neutral composites and a gauge transition are different '
        'families of source objects and detector responses. No protected-object ID is attached to these '
        'excitation rows, and no natural species is identified by a label alone.', 20, MUTED, 1.42)+20
    for family in value['excitation_families']:
        body = Canvas()
        body.text(X0+24, y+34, f"{family['label']} · {len(family['rows'])} source IDs", 20, TEXT, 700)
        labels = '; '.join(row['label'] for row in family['rows'])
        bottom = body.para(X0+24, y+65, CW-48, labels, 15.5, SOFT, 1.35)
        h = bottom-y+23
        c.rect(X0, y, CW, h, CARD2, EDGE, 1.2, 15)
        c.items.extend(body.items)
        c.overflow.extend(body.overflow)
        y += h+18
    return y


def confirmation(c, data, y):
    c.caps(X0, y+10, 'Candidate confirmation', AMBER, 14)
    c.text(X0, y+60, 'Six tests connect a source response to a particle', 40, TEXT, 700)
    y = c.para(X0, y+100, CW,
        'Every named candidate uses this sequence on one source, state, reference and clock. '
        'A finite positive signal confirms a detector response. Physical identification '
        'also needs its transformation, interaction and mass law.', 21, MUTED, 1.42)+16
    labels = (
        ('1 · Create and read', 'Bind the source operator, conserved currents, conjugate channel and complete detector. Keep failures and compensating charge.'),
        ('2 · Energy and momentum', 'Resolve the populated joint response using the operational translation generator and the same source clock.'),
        ('3 · Spin and statistics', 'Prove the physical Lorentz transformation, chirality and exchange law. A finite cube label alone is insufficient.'),
        ('4 · Species and interactions', 'Test current vertices, distinct flavors or polarizations, conjugates and allowed decays. Retain competing modes.'),
        ('5 · The right mass', 'Use a charged threshold, unstable resonance or specified renormalized mass. Bound volume, cutoff and numerical errors.'),
        ('6 · Independent comparison', 'Choose one common energy unit and compare complete predictions with experimental errors and declared calibration ancestry.'),
    )
    gap, h = 22, 166
    w = (CW-gap)/2
    for i, (label, body) in enumerate(labels):
        x, yy = X0+(i%2)*(w+gap), y+(i//2)*(h+gap)
        card(c, x, yy, w, h)
        c.text(x+24, yy+40, label, 23, AMBER, 700)
        c.para(x+24, yy+76, w-48, body, 17, SOFT, 1.42)
    y += 3*(h+gap)
    c.link(X0+24, y+15, 'The exact tests for every lepton, quark, gauge boson and Higgs candidate',
           '../../model/particle_candidate_confirmation/CONFIRMATION_STEPS.md', 17)
    return y+42


def chronology(c, data, y):
    c.caps(X0, y+10, 'Which appears first?', VIOLET, 14)
    c.text(X0, y+60, 'A catalogue is not a cosmological birth order', 40, TEXT, 700)
    y = c.para(X0, y+100, CW,
        'The hot early Universe supports a temperature-dependent comparison: quarks and gluons '
        'at sufficiently high temperature; photons, electrons, positrons, neutrinos and nucleons '
        'around nucleosynthesis. That does not identify the very first particles or an order '
        'of source defects.', 21, MUTED, 1.42)+16
    h = 188
    card(c, X0, y, CW, h)
    c.text(X0+26, y+43, 'Selected source state: ρ = g(H), so [ρ, H] = 0', 25, TEAL, 700)
    c.para(X0+26, y+80, CW-52,
        'A fixed one-time source effect has time-independent expectation in this preparation. '
        'Conditional apparatus histories can change. A physical early-Universe history needs '
        'an initial preparation, production and cooling dynamics, and a map to observer time. '
        'The capacity/refinement index is not that clock; a first click is not a particle birth.',
        18, SOFT, 1.42)
    y += h+32
    c.link(X0+26, y, 'Source stationarity and chronology boundary',
           '../../model/particle_candidate_confirmation/CATALOGUE.md', 16)
    c.link(X0+690, y, 'PDG: Big-Bang cosmology',
           'https://pdg.lbl.gov/2026/reviews/rpp2026-rev-bbang-cosmology.pdf', 16)
    c.link(X0+1020, y, 'PDG: nucleosynthesis',
           'https://pdg.lbl.gov/2026/reviews/rpp2026-rev-bbang-nucleosynthesis.pdf', 16)
    return y+22


def responses(c, data, y):
    source = data['identity_tests']['calculated_source_results']
    charged = source['charged_propagation']['filtered_original_window']
    scalar = source['scalar_response']
    composite = source['composite_spectrum']
    interval = '['+', '.join(charged['reported_transfer_interval'])+']'
    c.caps(X0, y+10, 'Tests on the populated source', TEAL, 14)
    c.text(X0, y+60, 'Finite detector bounds and mass tests', 40, TEXT, 700)
    y = c.para(X0, y+100, CW,
        'These bounds use the original selected source state and its clock. They keep '
        'the complete apparatus records and quantum reference. A positive record probability '
        'does not by itself identify an observed particle.', 21, MUTED, 1.42)+16
    rows = (
        ('Compensated charged pair',
         'Probability > '+charged['original_filtered_recorded_simple_strict_lower'],
         'A time-filtered local dressed pair has recorded transfer energy in '+interval+
         ' source units at member zero. The companion endpoint is retained.',
         'This is an even pair observable. The canonical bulk spatial aperture is inactive '
         'at this member. Its local positive weight does not prove a propagating charged '
         'particle or a surviving continuum mass shell.',
         '../../model/particle_candidate_charged_propagation/PROOF.md'),
        ('Neutral composite families',
         str(composite['component_count'])+' component witnesses',
         'Each signed composite component has a nonzero neighboring graded response. '
         'A complete original spatial/final-energy experiment has a strictly positive '
         'analytic contrast bound.',
         'The bound uses an extremely small state-weight floor and is not a practical '
         'shot guarantee. Final energy is not transfer energy or rest mass. A physical '
         'odd two-point dispersion has not been computed.',
         '../../model/particle_candidate_composite_spectrum/PROOF.md'),
        ('Neutral radial scalar',
         'Total probability ≥ '+scalar['uniform_original_total_success_lower'],
         'The scalar contraction has this total detection floor uniformly in capacity; '
         'member zero exceeds '+scalar['member_zero_simple_strict_lower']+
         '. A separate two-time commutator experiment has nonzero original probability.',
         'The large total includes elastic response. The positive dynamical floor is tiny '
         'and has no uniform inelastic residue. A Higgs phase, scalar resonance and mass '
         'are separate physical identifications.',
         '../../model/particle_candidate_scalar_gauge_response/PROOF.md'),
    )
    rows += (
        ('Separated charged endpoints',
         'Lengths 1 → 2 → 3',
         'Original recorded probabilities exceed '+', '.join(
             r['original_filtered_recorded_simple_strict_lower']
             for r in source['endpoint_threshold']['separations'])+
         '. Each path transports Wilson flux between opposite charged endpoints. '
         'The source action, preparation and clock are unchanged.',
         'These positive energy bands overlap. The response includes both endpoints '
         'and their dressing. A charged rest mass needs the lower energy edge at '
         'controlled momentum and a controlled treatment of the companion charge.',
         '../../model/particle_candidate_endpoint_threshold/PROOF.md'),
        ('Neutral spectral means',
         'Within-multiplet ratio < '+source['neutral_dispersion']['channels'][0][
             'within_multiplet_centroid_ratio_strict_upper'],
         'For the weak-pair/singlet trion, each G1 or H multiplet has four spatial '
         'spectral means, which may coincide, in (16/3, 38/3) source units. Their spread is < 5/9. '
         'Occupied-mode bounds constrain the original interacting source.',
         'A spectral mean averages many transitions. This ratio is not a particle '
         'mass ratio. The signed spatial response, conserved physical momentum '
         'and an isolated branch or threshold require further determination.',
         '../../model/particle_candidate_neutral_dispersion/PROOF.md'),
        ('Resolved scalar transfer',
         'Positive source cutoff: 1/40',
         'The original scalar measure assigns weight at least c × 9/(128 a² D²) '
         'to [1/40, D]. Here ρ ≥ c I and D bounds the full source energy width. '
         'Complete finite clocks retain at least three quarters in an expanded band.',
         'The floor can vanish and D grows with capacity. Elastic weight is '
         'at least 3481/160000. No compact surviving branch or Higgs mass follows. '
         'The enormous clock resources are specified; the histogram is not executed.',
         '../../model/particle_candidate_scalar_spectrum/PROOF.md'),
    )
    gap, h = 22, 442
    w = (CW-2*gap)/3
    for i, (label, headline, body, scope, href) in enumerate(rows):
        x = X0+(i%3)*(w+gap)
        yy0 = y+(i//3)*(h+gap)
        card(c, x, yy0, w, h)
        c.text(x+22, yy0+39, label, 20, TEAL, 700, maxw=w-44)
        c.text(x+22, yy0+75, headline, 18, TEXT, 700, maxw=w-44)
        yy = c.para(x+22, yy0+112, w-44, body, 17, SOFT, 1.42)
        yy = c.para(x+22, yy+20, w-44, scope, 16, MUTED, 1.42)
        if yy > yy0+h-54:
            raise ValueError('finite response card overflow: '+label)
        c.link(x+22, yy0+h-27, 'Source proof and exact bound', href, 16)
    return y+2*h+gap
