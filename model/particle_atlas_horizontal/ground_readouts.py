"""Low-energy preparation, complete signed records and native gauge cones."""
from fractions import Fraction as F


def validate_inputs(report):
    try:
        ground = report['spatial_low_energy_preparation']
        records = report['ground_energy_spectral_records']
        gauge = report['gauge_matter_characteristics']
        native = ground['native208']
        cooled = records['cooled_source_records']
        clock = cooled['complete_clock_apparatus']
        matched = gauge['declared_common_clock_example']
        if not (native['source_dimension'] == records['source_dimension'] == clock['source_dimension'] == 208
                and F(native['beta']) == F(cooled['beta']) == 4
                and native['Taylor_order'] == cooled['Taylor_order'] == 8
                and native['dyadic_scale'] == cooled['dyadic_scale'] == 128
                and native['same_H_and_clock'] and records['same_source_H_and_clock']
                and clock['same_H_clock'] and clock['all_quantum_references_retained']):
            raise ValueError('atlas cooled source, state or clock differs')
        if native['ground_energy_interval'] != records['global_ground_interval']:
            raise ValueError('atlas ground certificates disagree')
        if not (clock['complete_original_branch_count'] == len(clock['records']) == 16
                and F(clock['probability_sum']) == 1
                and records['instrument_effect_sum_identity']
                and all(records['charge_conjugation'].values())):
            raise ValueError('atlas complete charged records differ')
        if not (matched['same_clock_gauge_Dirac_characteristic_match']
                and F(matched['electric_scale'])*F(matched['magnetic_scale']) == 15
                and all(F(v) == 1 for v in matched['all_twelve_speed_squared'])
                and len(matched['all_twelve_speed_squared']) == 12):
            raise ValueError('atlas common gauge characteristic differs')
        false_flags = (ground['physical_Lorentz_or_Einstein_claimed'],
                       ground['thermodynamic_scope']['a_to_zero_limit_proved'],
                       native['postselection'], native['rest_mass_interpretation'],
                       records['accepted_branch_normalization'],
                       records['relativistic_mass_shell_derived'],
                       records['observed_particle_identification'],
                       clock['negative_bin_is_negative_energy_measurement'],
                       gauge['scope']['nonperturbative_interacting_relativistic_limit_proved'],
                       gauge['scope']['named_particle_mass_or_running_coupling_prediction'])
        if any(flag is not False for flag in false_flags):
            raise ValueError('atlas physical interpretation requires a separate proof')
    except (KeyError, TypeError, ZeroDivisionError) as exc:
        raise ValueError('incomplete atlas ground and characteristic contract') from exc


def render(canvas, data, y):
    from .drawing import W, M, INK, TEAL, BLUE, ORANGE, LINE
    from .observer_wings import ReadableCanvas
    c = ReadableCanvas(canvas)
    report = data['observer_outputs']; validate_inputs(report)
    ground = report['spatial_low_energy_preparation']
    records = report['ground_energy_spectral_records']
    gauge = report['gauge_matter_characteristics']
    native = ground['native208']; cooled = records['cooled_source_records']
    c.items.append('<g id="ground-records-common-cone">')
    c.text(M,y,'LOW-ENERGY PREPARATION → CHARGED RECORDS → COMMON WAVE SPEED',86,INK,700,serif=True)
    c.para(M,y+170,W-2*M,
           'An excitation is measured relative to a preparation. First control how much energy that preparation can release. Then read every charged branch on the same clock. Finally compare spatial propagation: a finite transition frequency needs this connection before it can represent a particle mass.',
           60,INK,max_y=y+440)
    top=y+540; column=(W-2*M)/3; width=column-340
    x0,x1,x2=(M+i*column for i in range(3))
    for x in (x0,x1,x2): c.rect(x,top,column-80,4020,'#fffdf7',LINE,30)

    c.text(x0+100,top+150,'1. Cool the declared field',81,INK,700)
    c.para(x0+100,top+330,width,
           'Keep every fermion occupation and a finite word space on every quantum link. The dimension D and an energy interval [L,U] are bounded from the code. No guessed ground vector or spectral gap is needed.',
           59,INK,max_y=top+730)
    for j,(title,detail) in enumerate((('TARGET','exp(−βH/2)'),('FINITE PREPARATION','positive polynomial'),('STATE + REFERENCE','all branches retained'))):
        xx=x0+170+j*2820
        c.rect(xx,top+930,2600,590,'#e6f0ed',radius=25)
        c.text(xx+100,top+1090,title,44,TEAL,700)
        c.text(xx+100,top+1320,detail,48,INK)
        if j<2: c.line(xx+2640,top+1230,xx+2790,top+1230,TEAL,7,arrow=True)
    c.text(x0+150,top+1800,'mean energy − ground energy ≤ log(D)/β + (U−L)η',61,TEAL,700)
    c.para(x0+100,top+2040,width,
           'β is the chosen cooling parameter; η is the certified preparation error. A larger β and a more accurate polynomial give any positive finite energy tolerance. The construction can be extremely expensive. It does not assume passive repair performs this cooling.',
           59,INK,max_y=top+2500)
    c.para(x0+100,top+2740,width,
           'The existing finite entropy theorem supplies the log(D)/β bound. The prepared density is stationary under its own H. A conjugate reference carries the compensating gauge charge; gauge invariance does not mean every source component is a singlet.',
           58,INK,max_y=top+3190)
    c.para(x0+100,top+3390,width,
           'At fixed spacing and link cutoff, growing boxes have a stationary spatial ground subsequence. A unique state and convergent records as physical spacing shrinks require additional control.',
           60,ORANGE,max_y=top+3700)
    c.link(x0+100,top+3890,'Preparation proof and finite entropy theorem',
           'model/particle_spatial_ground_state/PROOF.md',49)

    c.text(x1+100,top+150,'2. Pair opposite source charges',81,INK,700)
    c.para(x1+100,top+330,width,
           'On the full 208-state directional control, exact charge conjugation preserves H and reverses charge. It pairs N=1 and N=3 spectra around the N=2 ground. This symmetry is checked on every occupation and flux.',
           59,INK,max_y=top+730)
    lo,hi=(F(v) for v in records['adjacent_charge_minimum_transfer_intervals']['1'])
    threshold=float((lo+hi)/2)
    base=top+1710; level=top+1110
    c.line(x1+1600,base,x1+7500,base,TEAL,8)
    c.text(x1+4550,base+170,'N = 2 ground',57,TEAL,700,anchor='middle')
    for xx,label in ((x1+2300,'N = 1'),(x1+6800,'N = 3')):
        c.line(xx-1000,level,xx+1000,level,BLUE,8)
        c.text(xx,level-150,label,60,BLUE,700,anchor='middle')
        c.line(xx,base-60,xx,level+70,ORANGE,6,arrow=True)
    c.text(x1+4550,top+2100,f'paired sector thresholds ≈ {threshold:.9f}',65,BLUE,700,anchor='middle')
    c.text(x1+4550,top+2270,'source-energy units; no GeV or species assignment',49,ORANGE,anchor='middle')
    c.para(x1+100,top+2490,width,
           'A compensated CAR probe has eight original branches. The two-bin phase apparatus keeps all 16 records and their references. A sector minimum need not have nonzero residue for that probe. A phase bin does not directly mean negative energy.',
           58,INK,max_y=top+2930)
    c.text(x1+150,top+3160,'P(energy decrease ≥ δ) ≤ ε / δ',68,TEAL,700)
    c.para(x1+100,top+3380,width,
           'ε is the prepared ground-energy excess. This Lean-checked bound uses complete original spectral weights, with no success filtering. A spatial particle interpretation also needs dispersion, currents and physical statistics.',
           57,INK,max_y=top+3730)
    c.link(x1+100,top+3890,'Exact charge symmetry, spectral bounds and complete apparatus',
           'model/particle_ground_spectral_records/PROOF.md',47)

    c.text(x2+100,top+150,'3. Match the actual gauge cone',81,INK,700)
    c.para(x2+100,top+330,width,
           'The magnetic plaquette uses the native fifteen-state trace. Its electric term uses the inverse trace metric. These factors cancel in the quadratic wave equation, giving one speed for all twelve gauge directions.',
           59,INK,max_y=top+730)
    ox=x2+2200; oy=top+1840; span=3700; rise=900
    c.line(ox,oy,ox+span+900,oy,LINE,5,arrow=True)
    c.line(ox,oy,ox,oy-rise-220,LINE,5,arrow=True)
    c.path(f'M{ox},{oy} L{ox+span},{oy-rise}',TEAL,12)
    c.path(f'M{ox},{oy} L{ox+span},{oy-rise/105**.5}',ORANGE,8,dashed=True)
    c.text(ox+span+220,oy+130,'momentum',47,INK)
    c.text(ox-130,oy-rise-300,'frequency',47,INK)
    c.text(ox+1600,oy-rise-130,'matched principal cones',50,TEAL,700)
    c.text(ox+1450,oy+240,'small-momentum comparison',49,ORANGE,700)
    c.text(x2+150,top+2330,'c²gauge = e₀κ₀ / 15;  c²matter = 1',67,TEAL,700)
    c.para(x2+100,top+2570,width,
           'Choose e₀κ₀ = 15 for a common principal speed. For κ₀ = 1/7 this means e₀ = 105. Keeping e₀ = 1 instead gives c²gauge = 1/105. Changing the time unit cannot remove that relative speed difference.',
           59,INK,max_y=top+2990)
    c.para(x2+100,top+3200,width,
           'At nonzero momentum the lattice curl has two transverse responses and a longitudinal gauge null. The nonabelian plaquette also has a nonzero quartic interaction. A shared quadratic cone does not prove the interacting Lorentz limit or identify a physical photon or gluon.',
           58,ORANGE,max_y=top+3720)
    c.link(x2+100,top+3890,'Native trace, ordered plaquette and lattice dispersion',
           'model/particle_gauge_characteristics/PROOF.md',49)
    c.items.append('</g>')
    return top+4180
