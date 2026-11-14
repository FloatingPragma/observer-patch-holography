"""Literal size-four native spatial source and finite precision refinement."""
from connection import *
from itertools import product

def save(name,obj):
    (HERE/name).write_text(json.dumps(obj,indent=2,default=str)+'\n');print(name,'PASS',flush=True)

def native_plane():
    width=2;size=4;mode=(1,0,0);h=F(1,4);info=spatial_layout(width)
    prep=[]
    for axis in range(3):
        for j in range(width):
            wire=(2-axis)*width+j
            prep+=list(words.remap(native.primitive_templates()['hadamard'],(wire,)))
            prep.append(('Z',(wire,),native.angle(-F(mode[axis]*(1<<j),size))))
    plane_events=events(prep);wave,failed,gram,prep_bound=execute(plane_events,6,np.eye(64)[:,0]);wave=wave[:,0]
    # Track the preparation's exact scalar phase as part of the vector.
    phase=np.exp(-1j*math.pi*sum(F(mode[a]*(size-1),size) for a in range(3)))
    idealwave=np.array([phase*np.exp(2j*np.pi*(mode[0]*(c//16)+mode[1]*((c//4)%4)+mode[2]*(c%4))/size)/8 for c in range(64)])
    assert np.linalg.norm(wave-idealwave)<3e-12
    # Flavor, cell and spin occupy source bits9..0. All sixteen spin/flavor
    # columns include the retained mass-padding block.
    embedding=np.zeros((1024,16),complex)
    for flavor in range(4):
        for spin in range(4):
            for cell,value in enumerate(wave):embedding[spin+4*cell+256*flavor,spin+4*flavor]=value
    tape=events(spatial_matter_word(width,h));out,fail,complete_gram,bound=execute(tape,10,embedding)
    block=np.zeros((16,16),complex)
    for flavor,mass in enumerate((*MASSES,F())):block[4*flavor:4*flavor+4,4*flavor:4*flavor+4]=bloch(2*np.pi*np.array(mode)/size,h,mass)
    error=float(np.linalg.norm(out-embedding@block,2));assert error<3e-10
    compressed=embedding.conj().T@out;assert np.linalg.norm(compressed-block,2)<3e-10
    save('SPATIAL_NATIVE_PROGRAM.json',{'plane_preparation':plane_events,'whole_spatial_macro':tape})
    return {'L':4,'mode_integer':mode,'h':str(h),'native_source_bits':10,'all_spin_flavor_columns':16,
        'physical_source_positions':64,'retained_flavor_padding_columns':4,'source_native_pulses':len(tape),
        'plane_native_pulses':len(plane_events),'full_spatial_to_Fourier_intertwiner_residual':error,
        'complete_Gram_residual':complete_gram,'first_export_failure_Gram_trace':float(np.trace(fail).real),
        'source_native_error_exact':str(bound),'native_plane_error_exact':str(prep_bound),
        'phase_alignment_applied':False,'all_first_source_failures_retained':True}

def refinement_contract():
    # A literal finite native program is defined for every dyadic L; only
    # small members are executed. This prices its controlled apparatus too.
    size=1<<18;width=18;h=F(1,size);T=F(1,4);steps=int(T/h);mode=(1,2,3)
    p_upper=[F(44,7)*abs(x) for x in mode];m=max(MASSES)
    values=[m]+p_upper;S=sum((values[i]*values[j] for i in range(4) for j in range(i+1,4)),F())
    source_bound=T*h*S
    mass_boxes=read(SPECIES/'REGISTER.json')['effect_intervals']
    mass_radius=max((F(hi)-F(lo))/2 for lo,hi in mass_boxes)
    coefficient_bound=T*mass_radius
    primitive=native.primitive_templates();logical_counts={}
    info=spatial_layout(width)
    for axis in range(3):
        for name,targets in signed_shift_gates(width,info['coordinates'][axis],1,info['scratch']):logical_counts[name]=logical_counts.get(name,0)+1
    mass_pulses=len(list(species.flavor_mass_word(MASSES,h,info['flavor'],0)))
    rotations=2*sum(len(list(words.spin_diagonalizer(a))) for a in range(3))
    per_macro=mass_pulses+rotations+sum(count*len(primitive[name]) for name,count in logical_counts.items())
    # Every primitive-controlled template in connection.controlled has <=69
    # elementary M6 pulses (J is the largest eight-Pauli expansion). Check the
    # finite primitive kinds instead of assuming a gate-library constant.
    factors={kind:len(list(controlled([(kind,tuple(range(2 if kind=='J'else 1)),native.angle(0,0,F(1,7)))],3))) for kind in ('I','X','Y','Z','H','J')}
    control_factor=max(factors.values())
    # Native Fourier preparation uses H plus Z per coordinate bit, spin
    # tomography and flavor flips at most24 pulses; clock/probe apparatus and
    # both pointers at most40. Upper bound includes all declared source wires.
    # Split each integer mode phase into |nu_axis| equal Z rotations.
    # Every q_pi then has magnitude <=1/2, inside the universal pulse bound.
    prep_pulses=width*(6+sum(abs(x) for x in mode))+24
    total=steps*(control_factor*per_macro+1)+prep_pulses+40
    precision=128
    angle=max(F(4),h*m)
    while total*words.universal_pulse_error(precision,angle)>F(1,10**10):precision+=32
    native_bound=total*words.universal_pulse_error(precision,angle)
    reset=native.reset_repetition_plan(info['bits']+4,F(1,10**9))
    bound=min(F(1),source_bound+coefficient_bound+native_bound+F(reset['complete_timeout_union_bound_exact']))
    return {'L':size,'h':str(h),'T':str(T),'macro_count':steps,'mode_integer':mode,
        'uniform_physical_flavor_source_error_exact':str(source_bound),'source_error':float(source_bound),
        'mass_effect_coefficient_radius_exact':str(mass_radius),'coefficient_transfer_error_exact':str(coefficient_bound),
        'native_pulses_upper':total,'control_template_max_factor':control_factor,'per_primitive_control_pulses':factors,
        'precision_bits':precision,'reduced_source_memory_bits':info['bits'],'reset_registers_including_pointer_work':info['bits']+4,
        'native_operator_error_upper_exact':str(native_bound),'reset_contract':reset,
        'complete_observed_record_TV_bound_exact':str(bound),'complete_observed_record_TV_bound':float(bound),
        'huge_native_word_executed':False,'full_field_state_allocated':False,
        'spatial_family':'Actual root-port routes with dyadic signed binary translations; all four spin and flavor code blocks, including padding, retained.',
        'field_domain':'Observed free-factor instrument; arbitrary autonomous spectator reference. Full coupled-field failure tapes are distinct.'}

def full_field_terminal_contract():
    # Native complete field holding at the same finite size2 remains priced
    # even when no field measurement is requested. The static contract keeps
    # all zero-point phases and normalized harmonic packets.
    contract=words.static_holding_contract(F(4),F(1,10**8))
    native_error=F(contract['native_export_operator_error_exact'])
    reset=F(contract['preparation']['native_reset']['complete_timeout_union_bound_exact'])
    # Cutoff and splitting change only this independent spectator unitary in
    # a reduced matter experiment. Native failure/control error is retained.
    reduced=native_error+reset
    return {'source_member_L':2,'lattice_holding_time':4,'physical_holding_time':'4h',
        'field_contract':contract,'additional_reduced_record_TV_bound_exact':str(reduced),
        'all_field_native_words_executed':False,
        'record_scope':'Fresh coalesced complete terminal implementation. Every field first-export failure and preparation timeout is retained; no per-macro failure-history identity is asserted.',
        'full_state_boundary':'Field full-state or electric-observable comparison requires the stated harmonic/spectral transfers. Cancellation applies only to matter/clock marginal laws.',
        'larger_L_field_bound_transferred':False}

def main():
    pins=verify();save('SPATIAL_NATIVE.json',native_plane())
    save('REFINEMENT_CONTRACT.json',refinement_contract())
    save('FULL_FIELD_TERMINAL.json',full_field_terminal_contract())
    assert verify()==pins;save('SPATIAL_CHECKS.json',{'all_checks_pass':True,'producer_pins_unchanged':True})
if __name__=='__main__':main()
