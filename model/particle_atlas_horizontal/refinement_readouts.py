"""Physical refinement bounds and one shared thermal/strain experiment."""
from fractions import Fraction
import math


def validate_inputs(report):
    """Keep the displayed shared experiment and interpretation on one contract."""
    try:
        cone=report['refinement_uniform_causality']
        state=report['interacting_field_thermal_records']
        geometry=report['declared_spatial_geometry_response']
        prep=state['finite_state_preparation']; poly=prep['positive_polynomial']
        records=state['complete_flux_records']; work=geometry['shared_quantum_rotor_work']
        thermal=state['source_thermal_identities']; boundary=geometry['metric_identification_boundary']
        coframes=geometry['compact_coframe_causality']
        if not (prep['source']==work['source']=='unchanged field_route.exact_rotor(6)'
                and prep['dimension']==work['complete_dimension']==records['source_dimension']==208
                and work['all_CAR_occupations']==16
                and work['all_electric_fluxes']==list(range(-6,7))):
            raise ValueError('atlas thermal/strain source or complete carrier changed')
        if Fraction(work['magnetic_coefficient']) != 0:
            raise ValueError('atlas directional thermal control requires kappa=0')
        if not (Fraction(prep['declared_inverse_temperature'])
                ==Fraction(poly['beta'])==Fraction(work['inverse_temperature_target'])==Fraction(1,64)
                and poly['order']==12):
            raise ValueError('atlas thermal preparation temperature or polynomial order changed')
        eps=Fraction(poly['canonical_state_half_trace_error_upper'])
        if not (eps>0 and Fraction(work['target_state_trace_error_upper'])==eps
                and Fraction(work['target_Gibbs_pressure_error_upper'])==56*eps/3):
            raise ValueError('atlas thermal and strain preparation errors disagree')
        if not (list(map(Fraction,records['measurement_times']))==[Fraction(0),Fraction(1,64)]
                and Fraction(records['duration'])==Fraction(1,64) and records['Cayley_steps']==4):
            raise ValueError('atlas flux history changed its displayed clock')
        false_flags=(cone['Lorentz_covariance_continuum_vacuum_or_chirality_proved'],
                     cone['measured_light_speed_or_particle_mass_inferred'],
                     state['physical_vacuum_or_continuum_derived'],
                     prep['joint_purification_stationary'],
                     thermal['polynomial_preparation_exactly_Gibbs'],
                     thermal['arbitrary_reference_operations_passive'],
                     boundary['proper_geometric_stress_or_Einstein_limit_claimed'],
                     boundary['physical_particle_or_GeV_assignment'],
                     coframes['compact_family']['proper_geometry_Lorentz_covariance_or_limiting_state_proved'],
                     coframes['shared_rotor_thermal_state_transferred_to_spatial_refinement'],
                     work['global_Gibbs_entropy_slope_is_regional_geometric_entropy'],
                     work['pressure_is_a_physical_cosmological_equation_of_state'],
                     work['rotor_is_a_normalizable_reducing_subspace_of_full3D_quantum_field'])
        if any(flag is not False for flag in false_flags):
            raise ValueError('atlas state or geometric interpretation requires a separate review')
    except (KeyError,TypeError,ZeroDivisionError) as exc:
        raise ValueError('incomplete atlas thermal/strain composition contract') from exc


def render(c, data, y):
    from .drawing import W, M, INK, TEAL, BLUE, ORANGE, LINE
    from .observer_wings import ReadableCanvas
    c=ReadableCanvas(c)
    report=data['observer_outputs']
    validate_inputs(report)
    cone=report['refinement_uniform_causality']
    state=report['interacting_field_thermal_records']
    geometry=report['declared_spatial_geometry_response']
    plan=cone['Wilson_fixed_time_example']
    work=geometry['shared_quantum_rotor_work']
    records=state['complete_flux_records']
    c.items.append('<g id="refinement-thermal-geometry">')
    c.text(M,y,'REFINEMENT, PREPARATION AND MECHANICAL RESPONSE',88,INK,700,serif=True)
    c.para(M,y+160,W-2*M,
           'A causal estimate constrains how records can influence one another. A prepared state fixes their probabilities. A specified deformation of the action makes work and pressure measurable. These are distinct premises; a continuum spacetime description needs their common physical identification.',
           61,INK,max_y=y+350)
    column=(W-2*M)/3; top=y+450
    x0,x1,x2=(M+i*column for i in range(3)); width=column-350
    for x in (x0,x1,x2):c.rect(x,top,column-75,4400,'#fffdf7',LINE,30)

    c.text(x0+110,top+155,'1. Refine the lattice; retain causal control',82,INK,700)
    c.para(x0+110,top+320,width,
           'Set physical distance = spacing × graph distance. For this fixed-regulator action, kinetic terms scale as 1/a. With magnetic coefficient ≤ κ₀/a, the locality rate is at most W*/a: a conservative propagation speed stays finite.',
           63,INK,max_y=top+640)
    # Plot the actual conservative rational cone of the stated Wilson example.
    ox,oy=x0+650,top+2220; sx,sy=1550,1150000
    T=float(Fraction(plan['duration']))
    d=float(Fraction(plan['physical_separation_lower']))
    v=float(Fraction(cone['ordinary_Wilson_family']['physical_speed_rational_upper']))
    c.path(f'M{ox},{oy} L{ox},{oy-sy*T} L{ox+sx*v*T},{oy-sy*T} Z',TEAL,5,'#e2f0e9')
    c.line(ox,oy,ox+sx*4.6,oy,LINE,5,arrow=True)
    c.line(ox,oy,ox,oy-sy*T-140,LINE,5,arrow=True)
    c.text(ox+sx*4.1,oy+160,'distance',54,INK)
    c.text(ox-130,oy-sy*T-200,'time',54,INK)
    px,py=ox+sx*d,oy-sy*T
    c.line(ox,py,px,py,ORANGE,5,dashed=True)
    c.circle(px,py,75,ORANGE)
    c.text(ox+sx*.35,oy-sy*T*.6,'conservative cone',57,TEAL,700)
    c.text(px-1450,py-180,'d = 4,  T = 1/1000',58,ORANGE,700)
    c.para(x0+110,top+2550,width,
           f'The declared Wilson example certifies complete-record influence ≤ 2^(−{plan["certified_error_bits"]}) for EVERY spacing a ≤ {plan["all_smaller_spacing_upper"]}. The proof includes the growing number of cells in a fixed physical region and every apparatus outcome.',
           63,TEAL,max_y=top+2920)
    c.para(x0+110,top+3140,width,
           'This is a bound on a remote intervention, not a simulated particle trajectory. Finite-spacing tails can be nonzero. Both experiments’ implementation errors must be added; choosing each ≤ a preserves vanishing outside-cone influence.',
           61,INK,max_y=top+3470)
    c.para(x0+110,top+3680,width,
           'The same causal argument covers constant spatial deformations whose stretch and compression stay within fixed bounds: electric cross terms remain within one cell. A limiting state, Lorentz boosts, physical light speed and particle masses need further identification.',
           61,ORANGE,max_y=top+4080)
    c.link(x0+110,top+4250,'Full-Fock norm, complete-history bound and all-spacing proof',
           'model/particle_physical_causality/PROOF.md',52)

    c.text(x1+110,top+155,'2. Prepare one interacting thermal control',82,INK,700)
    c.para(x1+110,top+320,width,
           'Keep all 16 Dirac occupations and 13 electric fluxes: 208 states. This directional control has no magnetic term and uses β = 1/64. A positive polynomial approximates its Gibbs amplitude; one joint unitary prepares it without a success filter.',
           63,INK,max_y=top+650)
    for j,(title,detail) in enumerate((('TARGET','exp(−βH/2)'),('FINITE CODE','positive Q₁₂'),('PREPARED','vec(Q₁₂) / √Tr(Q₁₂²)'))):
        xx=x1+180+j*2870
        c.rect(xx,top+800,2600,440,'#e8eef7',radius=22)
        c.text(xx+130,top+950,title,56,BLUE,700)
        c.text(xx+130,top+1110,detail,53,INK)
        if j<2:c.line(xx+2640,top+1020,xx+2820,top+1020,BLUE,6,arrow=True)
    c.text(x1+150,top+1450,'ALL 169 ORIGINAL FLUX RECORDS: first at 0, second at 1/64',56,BLUE,700)
    rows=records['records']; labels=records['first_flux_labels']
    assert len(rows)==169 and sum(Fraction(r['original_probability']) for r in rows)==1
    table={(r['first_flux'],r['second_flux']):Fraction(r['original_probability']) for r in rows}
    maximum=max(table.values()); cell=130; gx=x1+3700; gy=top+1680
    for i,first in enumerate(labels):
        c.text(gx+i*cell+cell/2,gy-65,str(first),39,BLUE,anchor='middle')
        c.text(gx-130,gy+i*cell+cell*.72,str(first),39,BLUE,anchor='middle')
        for j,second in enumerate(labels):
            p=table[first,second]
            intensity=math.sqrt(float(p/maximum)) if p else 0
            color='#e7e9e7' if not p else f'rgb({int(221-165*intensity)},{int(233-116*intensity)},{int(240-58*intensity)})'
            c.rect(gx+i*cell,gy+j*cell,cell-8,cell-8,color,radius=5)
    c.text(x1+180,top+2050,'first flux →',55,BLUE,700)
    c.text(x1+180,top+2240,'second flux ↓',55,BLUE,700)
    c.para(x1+180,top+2500,3100,
           f'Grey cells are retained zeros ({records["zero_probability_labels"]}). Blue uses a square-root probability scale to reveal faint transitions. No accepted subset is renormalized.',
           49,INK,max_y=top+3020)
    c.para(x1+110,top+3570,width,
           'The prepared source density is stationary and source-passive; its purification need not be stationary. The first flux read has full backaction. Exact Gibbs identities apply to the target, with a certified preparation error. The clock has a separate Cayley error. This control is not a spatial vacuum.',
           61,ORANGE,max_y=top+4070)
    c.link(x1+110,top+4250,'Preparation, all records, thermal identities and error budgets',
           'model/particle_field_state_response/PROOF.md',52)

    c.text(x2+110,top+155,'3. Deform the action; read work and pressure',82,INK,700)
    c.para(x2+110,top+320,width,
           'Declare a spatial coframe V, metric G = VᵀV and volume v = det V. Differentiate the whole source action, including its regulator and gauge terms. Six symmetric strain operators follow on the spatial carrier.',
           63,INK,max_y=top+640)
    ox,oy=x2+4000,top+1630
    for dx,dy,label in ((2300,-100,'V e₁'),(750,-700,'V e₂'),(-1550,-500,'V e₃')):
        c.line(ox,oy,ox+dx,oy+dy,TEAL,9,arrow=True)
        c.text(ox+dx,oy+dy-95,label,59,TEAL,700,anchor='middle')
    c.circle(ox,oy,65,TEAL)
    for j,label in enumerate(('Kxx','Kyy','Kzz','Kxy','Kxz','Kyz')):
        xx=x2+1800+(j%3)*2500; yy=top+2070+(j//3)*460
        c.rect(xx-600,yy-180,1200,330,'#e5f0ed',radius=24)
        c.text(xx,yy+35,label,65,TEAL,700,anchor='middle')
    c.text(x2+140,top+2980,'SHARED 208-STATE DIRECTIONAL CONTROL',60,BLUE,700)
    pressure=float(Fraction(work['pressure_prepared_at_declared_volume1']))
    ratio=float(Fraction(work['pressure_over_energy_density_prepared']))
    c.text(x2+140,top+3170,f'At declared v = 1:  p ≈ {pressure:.6f};  p/ρ ≈ {ratio:.6f}',65,TEAL,700)
    c.text(x2+140,top+3340,'p/ρ = 1/3 + ⟨H_W⟩ / (3⟨H⟩); H_W is the Wilson regulator energy.',52,TEAL)
    c.para(x2+110,top+3535,width,
           'On identically prepared inputs, alternative complete ±28 instruments read longitudinal or isotropic work. Here ρ = ⟨H⟩/v. The state and clock are shared; this directional control is distinct from the full three-dimensional quantum-link field.',
           55,INK,max_y=top+3790)
    c.para(x2+110,top+3860,width,
           'Pressure comes from the chosen action’s volume response. It is not a cosmological equation of state. Proper geometry, Lorentz stress and the regional entropy/area premises for Einstein dynamics require their own attachment.',
           59,ORANGE,max_y=top+4130)
    c.link(x2+110,top+4300,'Six strain operators, shared-state pressure and geometry boundary',
           'model/particle_geometric_response/PROOF.md',51)
    c.items.append('</g>')
    return top+4560
