"""Current joint native observer construction with explicit physical boundaries."""
from copy import deepcopy
from fractions import Fraction as F

NAMES=('native_observer_program','native_charged_propagation','native_electric_gauss',
       'native_gauss_walk','native_moving_clock','native_gauss_clock',
       'native_transport_spectrum','native_spatial_refinement','native_criterion_packet')
TITLES=('Complete native observer program and compiled reflection',
        'Charged transport, contact interaction and localized apparatus modes',
        'Integer boundary charge and complete electric-loop records',
        'One evolving charged source and integer electric field',
        'Transported clock on the undressed source',
        'Complete clock observations on the same charged field state',
        'Conditional transport energies and enlargement controls',
        'Declared three-dimensional routing and fixed-mode scalar dispersion',
        'Literal charge-criterion review packet')

def compact(value):
    if isinstance(value,dict):
        return {k:compact(v)for k,v in value.items()
                if k!='source_sha256'and not k.endswith('_zlib_base64')
                and k!='initial_full_region_charge_table'}
    if isinstance(value,list):return [compact(v)for v in value]
    return deepcopy(value)

def require(condition,message):
    if not condition:raise ValueError(message)

def project(data):
    require(all(n in data for n in NAMES),'missing joint source owner')
    for n in NAMES:require(str(data[n].get('verdict','')).startswith('PASS'),'unaccepted joint source: '+n)
    p=data[NAMES[0]];charged=data[NAMES[1]];electric=data[NAMES[2]];walk=data[NAMES[3]]
    bare=data[NAMES[4]];clock=data[NAMES[5]];spectrum=data[NAMES[6]];cubic=data[NAMES[7]];packet=data[NAMES[8]]
    identity=p['program_sha256'];extension=electric['extension_descriptor']
    require(identity==walk['common_program_sha256']==clock['common_program_sha256']
            ==packet['common_program_sha256']==bare['single_program_execution']['program_sha256'],
            'different common native program')
    require(extension==walk['integer_extension']==clock['integer_extension_descriptor']
            ==p['integer_extension_protocol']['extension_descriptor'],'different integer field extension')
    require(p['program']['hardware']==charged['construction']==cubic['local_hardware'],
            'different local native hardware')
    require(walk['preparation']['link_limit']==p['program']['integer_extension']['link_window']==4,
            'different integer-link cutoff')
    require(walk['preparation']['matching_labels']==[0,1]and walk['preparation']['active_character']==1
            and walk['preparation']['uniform_events']==60,'different primary source preparation')
    require(clock['scope']['same_charged_source_field_and_clock'],'clock is not on the joint field state')
    require(not clock['scope']['local_autonomous_clock_proved']
            and not clock['scope']['physical_rest_mass_identified']
            and clock['scope']['official_OpenWave_awards']==0,'tracked code promoted to physical clock/mass')
    require(not bare['scope']['literal_OpenWave_row3_complete'],'bare tracked clock promoted')
    require(not spectrum['scope']['spectrum_of_flux_dressed_joint_electric_program']
            and not spectrum['unit_contract']['physical_rest_mass_assignment']
            and not spectrum['unit_contract']['source_clock_and_walk_have_proved_common_stationary_generator'],
            'undressed conditional energies promoted to joint rest masses')
    require(not cubic['scope']['whole_interacting_quantum_field_theory_constructed']
            and not cubic['scope']['four_dimensional_Yang_Mills_mass_gap']
            and not cubic['scope']['Einstein_gravity_derived'],'scalar routing limit promoted')
    require(not packet['criterion1']['physical_acceptance_complete'],'candidate charge packet promoted')
    require(not electric['integer_cut_charge']['row1_review_scope']['official_row_promotion']
            and walk['scope']['OpenWave_awards']==0
            and p['physical_boundaries']['official_OpenWave_awards']==0,'finite native source awarded physics')
    require(not walk['scope']['reset_in_inverse_test']and not walk['scope']['undressed_source_spectrum_transferred'],
            'reset or unproved spectrum transfer admitted')
    require(F(walk['observed_electric_phase_response']['two_step_source_pair_total_variation'])==0
            and F(walk['observed_electric_phase_response']['three_step_response_lower_with_native_event_word_bound'])>0,
            'electric phase response scope changed')
    authority=packet['official_authority']
    require(authority['counts']=={'negative':0,'partial':0,'scored':0,'untested':31,'validated':0}
            and authority['rank']==7 and not authority['local_receipts_change_authority'],
            'official status differs from the pinned authority')
    require(len(authority['all31_criteria'])==31
            and {r['id']for r in authority['all31_criteria']}==set(range(1,32))
            and all(r['upstream_status']['M10']=='🚧'for r in authority['all31_criteria']),
            'literal upstream row status changed')
    require(len(clock['exact_transported_matrix_units'])==3
            and all(r['complete_matrix_units']==144 and r['exact_complete_Gram_identity']
                    and r['full_inverse_basis_checks']==12
                    for r in clock['exact_transported_matrix_units']),
            'complete transported matrix-unit contract changed')
    selected={n:compact(data[n])for n in NAMES}
    return {'source_packages':selected,'common_program_sha256':identity,
      'shared_integer_extension':deepcopy(extension),
      'proved_finite_interfaces':{
        'native_event_reflection_has_source_word_bounds':True,
        'integer_charge_has_complete_boundary_readout':True,
        'charged_source_and_electric_links_evolve_in_one_state':True,
        'tracked_clock_has_complete_same_joint_state_records':True,
        'full_inverse_keeps_all_source_and_field_information':True,
        'declared_cubic_routing_has_fixed_Fourier_scalar_dispersion_limit':True},
      'domains':{
        'joint_local_program':'Two-load native source, active memory, signed tag, integer links, coherent event and declared tracking/phase apparatus.',
        'conditional_transport_energies':'Undressed source spectrum, electron anchor and an extra declared event-time map; not joint electric-field rest masses.',
        'cubic_routing':'Same local native law on a declared three-dimensional routing family; the controlled scalar dispersion limit holds for any finite fixed set of Fourier modes. Whole joint interacting field/clock transfer is unproved.'},
      'official_status':{'validated':0,'untested':31,'rank':7,'authority':deepcopy(packet['official_authority'])},
      'physical_premises_closed':[],'physical_identification_ids_closed':[],
      'independent_Nature_matches_added':0,'new_independent_GeV_rest_masses':[],
      'complete_physical_mission_proved':False}

def decimal_lower(value,places=12):
    q=F(value);scale=10**places;n=q.numerator*scale//q.denominator
    return f'{n//scale}.{n%scale:0{places}d}'

def render(value):
    sources=value['source_packages'];walk=sources['native_gauss_walk'];electric=sources['native_electric_gauss']
    spectrum=sources['native_transport_spectrum'];cubic=sources['native_spatial_refinement']
    tv=walk['observed_electric_phase_response']['three_step_source_pair_total_variation']
    lower=decimal_lower(walk['observed_electric_phase_response']['three_step_response_lower_with_native_event_word_bound'])
    lines=['## One native observer construction with charged field and clock records','',
      'One extended native construction evolves the charged source, active memory, transported charge tag, integer electric links and event register together. Complete clock observations use an explicit reversible tracking apparatus on this same joint state. Architecture choices are declared; uniqueness is not a premise.','',
      '- [Joint charged-field evolution](../model/native_gauss_walk/PROOF.md). Every event preserves Gauss charge. Two steps contain 7,320 nonzero joint arcs, and the exact full inverse restores every source, tag, link and event amplitude without resetting a subsystem. An electric phase on these same links leaves position records equal for two steps; the third-step position distributions differ by total variation '+tv+'. The native event-word bound retains a strictly positive lower limit '+lower+'.',
      '- [Same-state tracked clock](../model/native_gauss_clock/PROOF.md). Complete transported matrix units and two-time outcomes include the canonical-dressing complement and all declared apparatus resources. The inverse full walk and coherent dressing decoder are part of the measurement; a local autonomous physical clock is not inferred.',
      '- [Integer boundary charge](../model/native_electric_gauss/PROOF.md). All 4,096 regional cuts measure integer electric flux. The controlled-hop experiment has regional expectations one half while its complete sharp records are (1,0) and (0,1); an equal half-charge split has zero probability. The integer extension and unit are declared, and physical electron identification is unproved.',
      '- [Compiled event reflection and complete records](../model/native_observer_program/PROOF.md). The common program retains source, active-memory and event padding, processor failures and original C/E record commitment. The recoded identity/swap event is distinct from the native final-endpoint coin.',
      '- [Charged transport and exact contact interaction](../model/native_charged_propagation/PROOF.md). The complete two-load source retains every charged sector and fixed event mode. Compact persistent event modes require their own apparatus preparation and have no physical particle identity.',
      '- [Conditional transport-energy intervals](../model/native_transport_spectrum/PROOF.md). The six active characters have exact finite phase spectra. Applying the electron anchor and a declared event-time map gives source-energy intervals; the clock and walk have no proved common stationary generator, and these are not flux-dressed particle rest masses.',
      '- [Three-dimensional scalar dispersion](../model/native_spatial_refinement/PROOF.md). A declared cubic routing uses the same local native rule, with controlled lowest scalar branch |k|^2/8 and quartic remainder. The dispersion limit concerns any finite fixed set of Fourier modes. The original spherical federation is unchanged. The whole charged-field and clock construction has no proved four-dimensional continuum or gravity transfer.',
      '- [Charge-criterion review packet](../model/native_criterion_packet/PROOF.md). Literal charge/flux records and indivisibility controls support a finite review candidate. Official OpenWave status is 0 validated and 31 untested; no local calculation changes an upstream icon.','']
    return lines
