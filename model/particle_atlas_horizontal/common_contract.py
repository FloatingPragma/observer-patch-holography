"""Visual contract for the selected source, state, reference and clock."""
from model.particle_common_contract import contract


def validate(report):
    contract.validate_inventory(report)
    value=report['common_source_contract']
    member=value['fixed_primary_member']
    contract.verify_member(member)
    for receipt in value['bound_consequences'].values():
        contract.verify_result(member,receipt)
    contract.verify_result(member,report['spatial_regional_modular_response']['selected_common_state_consumer'])
    if not (value['item1_fixed_construction_specification']
            and value['item2_finite_physical_preparation_with_boundary_charges']):
        raise ValueError('common construction specification is absent')
    if value['observed_particle_or_Einstein_identification'] or value['joint_physical_refinement_limit_proved']:
        raise ValueError('physical interpretation requires a separate proof')
    from model.particle_physical_boundary_state import boundary
    boundary.verify_bound(member,report['physical_boundary_state']['selected_state_boundary'])
    return member


def render(canvas,data,y):
    from .drawing import W,M,INK,TEAL,BLUE,ORANGE,LINE
    from .observer_wings import ReadableCanvas
    c=ReadableCanvas(canvas);member=validate(data['observer_outputs'])
    c.items.append('<g id="one-source-state-reference-clock">')
    c.text(M,y,'ONE PROGRAM THROUGH EVERY CLAIMED STEP',84,INK,700,serif=True)
    c.para(M,y+145,W-2*M,
        'Dependencies 1–3: one finite program has a physical Gauss preparation with boundary charges and complete local records with uniform influence bounds. The four objects below stay together in every admitted calculation. Joint physical refinement is a separate theorem.',
        58,INK,max_y=y+350)
    top=y+420;column=(W-2*M)/4;inside=column-270
    rows=[
        ('SOURCE','Open Wilson lattice',
         'All 60 fermion modes at each cell; native quantum links; open boundaries with zero external charge. Electric coefficient 105 and magnetic coefficient 1/7 give matched quadratic gauge and matter speeds.',
         'Space, cutoff, action and energy zero are fixed for each member.',
         'model/particle_common_contract/PROOF.md'),
        ('STATE','2. Physical state + boundary charges',
         'The selected density Q² / Tr Q² lies in the Gauss kernel. Across a regional cut, inside charge pairs with opposite outside charge; neither is erased. Actual Wilson loops show 11 sectors with positive weight. At degree one, a complete detector retains all 1,331 charge outcomes.',
         'Exact algebraic preparation is specified; the enormous spatial tape is not executed.',
         'model/particle_physical_boundary_state/PROOF.md'),
        ('REFERENCE','Keep the other half',
         'The canonical purification retains a conjugate quantum register. Its idle Hamiltonian is zero. This physical preparation makes source and reference separately gauge invariant. No branch or reference is discarded.',
         'A stationary source marginal does not make its joint purification stationary.',
         'model/particle_common_contract/PROOF.md#reference-and-time'),
        ('CLOCK','Use this H throughout',
         'Evolution is exp(−itH), in source units with ℏ = 1. Its preparation, interventions and complete records must refer to this clock. No fitted GeV conversion is inserted. Finite implementations carry their clock error.',
         'Modular K = log(Tr Q²) − 2 log Q on physical support; it is not an exact Gibbs clock.',
         'model/particle_spatial_modular_response/PROOF.md')]
    for j,(tag,title,body,boundary,link) in enumerate(rows):
        x=M+j*column
        c.rect(x,top,column-70,1360,'#fffdf7',LINE,28)
        c.rect(x+75,top+75,inside,190,'#e6f0ed',radius=20)
        c.text(x+150,top+200,tag,62,TEAL,700)
        c.text(x+110,top+380,title,70,INK,700)
        c.para(x+110,top+520,inside,body,54,INK,max_y=top+960)
        c.para(x+110,top+1030,inside,boundary,51,ORANGE,max_y=top+1230)
        c.link(x+110,top+1310,'Definition, actual consumer and proof',link,43)
        if j<3:c.line(x+column-45,top+680,x+column+30,top+680,TEAL,6,arrow=True)
    c.para(M,top+1500,W-2*M,
        'Verified joins: fixed source → physical preparation → complete boundary-charge records with uniform causal bounds and the full idle reference retained. Eleven supported charge sectors are not eleven particle species. The common quadratic speed does not prove Lorentz covariance. The interacting particle theory and the pure Yang–Mills gap have separate physical proof conditions.',
        56,BLUE,max_y=top+1710)
    c.items.append('</g>')
    return top+1750
