"""Declared field dynamics and complete, reference-aware measurement diagrams."""
import math
from fractions import Fraction


def render(c, data, y):
    from .drawing import W, M, INK, MUTED, TEAL, BLUE, ORANGE, RED, GOLD, LINE
    from .observer_wings import ReadableCanvas
    c = ReadableCanvas(c)
    output = data['observer_outputs']
    smooth = output['smooth_interacting_local_records']
    joint = output['joint_position_translation_records']
    control = joint['complete_reference_control']
    native = joint['native_complete_law']
    if smooth['physical_GeV_mass_chiral_or_Lorentz_limit_claimed']:
        raise ValueError('physical smooth-field scope changed')
    if joint['joint_crystal_momentum_is_physical_particle_momentum']:
        raise ValueError('joint readout became a particle identification')
    c.items.append('<g id="source-field-readouts">')
    c.text(M, y, 'ACTUAL FIELD READOUTS: WHAT CHANGES WHEN A DETECTOR RECORDS?', 88, INK, 700, serif=True)
    c.para(M, y+150, W-2*M,
           'The program supplies an action and a preparation. An apparatus performs a specified intervention. Its complete record law is calculable. A physical particle assignment requires an additional spectral and relativistic identification on that same field.',
           62, INK, max_y=y+330)
    column = (W-2*M)/3
    x0, x1, x2 = (M+i*column for i in range(3))
    top = y+420
    for x in (x0, x1, x2):
        c.rect(x, top, column-75, 3930, '#fffdf7', LINE, 30)
    width = column-340

    c.text(x0+110, top+155, '1. Choose the action; derive its records', 83, INK, 700)
    c.para(x0+110, top+300, width,
           'An abrupt sign changes sharply near a zero Wilson mode. A fixed positive smoothing is a permitted construction choice. It keeps the electric/link interaction and gives a controlled local expansion.',
           65, INK, max_y=top+620)
    # Scalar functional-calculus picture. The solid curve is the actual
    # normalized delta=1 regulator, not an empirical interpolation.
    cx, cy, sx, sy = x0+column/2, top+1170, 1300, 300
    c.line(cx-3100, cy, cx+3100, cy, LINE, 4, arrow=True)
    c.line(cx, cy+480, cx, cy-480, LINE, 4, arrow=True)
    points=[]
    for i in range(121):
        lam=-2.2+4.4*i/120
        value=math.sqrt(2)*lam/math.sqrt(lam*lam+1)
        points.append((cx+sx*lam, cy-sy*value))
    c.path('M'+' L'.join(f'{a:.2f},{b:.2f}' for a,b in points), TEAL, 12)
    c.path(f'M{cx-2860},{cy+sy} H{cx-24} M{cx+24},{cy-sy} H{cx+2860}', RED, 7, dashed=True)
    c.text(cx+2930, cy+85, 'λ', 80, INK)
    c.text(x0+150, top+1740, 'Solid: √2 λ / √(λ² + 1)     Dashed: sign(λ)', 62, TEAL, 700)
    c.para(x0+110, top+1930, width,
           'The smooth action has exponentially decreasing interaction terms. Their tails control quantum detector records, including the electric field, without assuming a uniform Wilson gap. This is a different declared action.',
           65, INK, max_y=top+2240)
    c.rect(x0+110, top+2390, width, 440, '#e5f1ec', radius=24)
    error = smooth['total_complete_local_record_error_upper']
    c.text(x0+220, top+2545, 'COMPLETE LOCAL RECORD ERROR ≤ '+error, 65, TEAL, 700)
    c.para(x0+220, top+2680, width-220,
           'At time 1 and smoothing = spacing = 1: polynomial + box + word cutoff. Finite resources are certified; they are not materialized.',
           57, INK, max_y=top+2790)
    c.para(x0+110, top+3020, width,
           'Ordinary finite-range Wilson dynamics is another allowed action. Neither choice supplies physical chirality, particle masses or a 3+1D interacting relativistic limit just from this estimate.',
           64, ORANGE, max_y=top+3330)
    c.link(x0+110, top+3550, 'Smooth local-record theorem and declared inputs',
           'model/particle_smooth_sign_field/PROOF.md', 57)
    c.link(x0+110, top+3690, 'Why low Wilson-mode weight alone is insufficient',
           'model/particle_sign_response_boundary/PROOF.md', 57)

    c.text(x1+110, top+155, '2. Move the field and its reference', 83, INK, 700)
    c.para(x1+110, top+300, width,
           'Translation moves every matter register and every link. On a charged mixed preparation it also moves the compensating reference. Leaving that reference behind can violate the joint gauge constraint.',
           65, INK, max_y=top+620)
    for row, (title, color) in enumerate((('FIELD: matter + links', TEAL), ('CONJUGATE REFERENCE', BLUE))):
        yy=top+950+row*560
        c.text(x1+170, yy-110, title, 62, color, 700)
        for j,label in enumerate(('x', 'x+b', '…', 'periodic')):
            xx=x1+1600+j*1700
            c.circle(xx, yy+100, 110, '#e6efed', color, 7)
            c.text(xx, yy+280, label, 60, color, anchor='middle')
            if j<3:c.line(xx+140, yy+100, xx+1550, yy+100, color, 8, arrow=True)
    c.text(x1+170, top+2070, 'Tjoint = Tsource ⊗ conjugate(Tsource)', 73, TEAL, 700)
    c.para(x1+110, top+2270, width,
           'A Fourier analyzer resolves all total crystal-momentum labels. The canonical referenced input has exactly k = 0. That says nothing about the momentum of an individual excitation.',
           65, INK, max_y=top+2560)
    c.rect(x1+110, top+2710, width, 570, '#e8eef7', radius=24)
    c.text(x1+220, top+2870, 'POSITION POINTER: Fx = nx / N, plus vacuum', 64, BLUE, 700)
    c.para(x1+220, top+3010, width-220,
           'Randomly sample an occupied mode and report its site. Keep the vacuum outcome. The square-root operation specifies how that record changes the state.',
           60, INK, max_y=top+3200)
    c.para(x1+110, top+3460, width,
           'Both menus are global: random-site position uses the total N; momentum translates the whole field. Their finite proofs do not supply a thermodynamic limit of these global instruments. Local occupation pointers have the separate local bound.',
           61, ORANGE, max_y=top+3720)

    c.text(x2+110, top+155, '3. Record position; compare momentum', 83, INK, 700)
    c.para(x2+110, top+300, width,
           'For the native initial state, V cells carry 60 fermion modes each. Every occupation is retained. Let A = (1 − 2^(−60V))/V and C = E[√(nx ny)/N] for distinct sites, with zero on vacuum.',
           65, INK, max_y=top+620)
    c.text(x2+160, top+860, 'p(x, k = 0) = [A + (V − 1) C] / V', 77, TEAL, 700)
    c.text(x2+160, top+1020, 'p(x, each k ≠ 0) = (A − C) / V > 0', 77, TEAL, 700)
    c.text(x2+160, top+1180, 'p(vacuum, k) = 2^(−60V) when k = 0; else 0', 68, TEAL)
    c.para(x2+110, top+1370, width,
           f'The side-{native["periodic_side"]} native example has {native["position_labels"]} position labels and {native["momentum_labels"]} momentum labels: {native["all_joint_labels"]} original joint records. These formulas describe the initial preparation; an evolved state uses its own Born law.',
           63, INK, max_y=top+1660)
    c.text(x2+130, top+1870, 'SMALL EXACT CONTROL: 2 sites, 1 mode/site, Z₂ links', 64, BLUE, 700)
    rows=control['joint_position_momentum_records']
    if sum(Fraction(r['original_probability']) for r in rows) != 1:
        raise ValueError('reference control lost complete original probabilities')
    labels=['site 0 / k0', 'site 0 / k1', 'site 1 / k0', 'site 1 / k1', 'vac / k0', 'vac / k1']
    bx, by, bw = x2+160, top+2660, width-200
    step=bw/6
    c.line(bx, by, bx+bw, by, LINE, 5)
    for i,(row,label) in enumerate(zip(rows, labels)):
        p=Fraction(row['original_probability']); h=2400*float(p)
        xx=bx+i*step
        c.rect(xx+100, by-h, step-210, max(3,h), BLUE, radius=12)
        c.text(xx+step/2, by-h-70, str(p), 68, BLUE, 700, anchor='middle')
        c.text(xx+step/2, by+100, label, 50, MUTED, anchor='middle')
    c.para(x2+110, top+2940, width,
           f'This distinct small control illustrates complete Born weights. For its phase preparation, the recorded position variance is {control["record_coordinate_variance"]}; the first-moment operator variance is {control["X_variance"]}. The extra {control["unsharp_coordinate_noise"]} is detector noise.',
           62, INK, max_y=top+3240)
    c.para(x2+110, top+3420, width,
           'Finite Robertson uncertainty applies to the actual operators. An exact finite [X,P] = iℏI is impossible. A physical particle’s position, momentum and mass require the continuum and species bridges.',
           61, ORANGE, max_y=top+3700)
    c.link(x2+110, top+3840, 'Complete native law, reference and uncertainty proof',
           'model/particle_joint_translation_readout/PROOF.md', 56)
    c.items.append('</g>')
    return top+4110
