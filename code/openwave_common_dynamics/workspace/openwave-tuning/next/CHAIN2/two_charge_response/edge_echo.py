"""Two-charge connected phase from literal shared-current native subwords.

This is a declared inverse-transport interferometer, not an uninterrupted
full-V trajectory. All four coherent branches contain exactly two charges.
"""
from fractions import Fraction as F
import hashlib
import json
import math

import numpy as np
import two_charge as t
import native_contract as nc


def run():
    mb=t.af.mb;compiler=mb.compiler
    field_program=json.loads((nc.BASE/'vector_field/FIELD_PROGRAM.json').read_text())
    block=t.af.mode_field.field_block(16)
    primitives=mb.compiled.controls()['primitive_words']
    # Compress the four supported global modes0,1,16,17 into0,1,2,3.
    # The intervening global spectator modes are exactly empty. Local
    # mode2 is a genuine occupied spectator in the exchange phase.
    mapping={0:4,1:5,2:6,3:7,4:1,5:3,6:8,7:9,8:10}
    transport=[]
    for kind,targets in block['program']['gates']:
        transport+=compiler.remap(primitives[kind],[mapping[q] for q in targets])
    transport+=mb.fock_swap_word(1,3)
    electric=compiler.remap(field_program['electric_word'],[4,5,6,7])
    alpha=.17;F0=1
    preparation=[]
    for bit in (0,2,4):preparation+=compiler.remap(primitives['flip'],[bit])
    preparation+=mb.givens_word(0,1,math.pi/4)+mb.givens_word(2,3,math.pi/4)
    # Phase differences relative to E(F0) are local, independently fixed
    # by the one-moving-charge branches. The connected phase survives.
    rephase=(mb.occupation_phase_word(1,alpha*(F0+.5))
             +mb.occupation_phase_word(3,alpha*(-F0+.5)))
    analysis=mb.givens_word(0,1,-math.pi/4)+mb.givens_word(2,3,-math.pi/4)
    source=np.zeros(1<<11,complex);source[0]=1
    prepared=compiler.execute(preparation,11,source)
    echo=transport+electric+compiler.inverse(transport)
    final=compiler.execute(echo+rephase+analysis,11,prepared)
    phases=[];truth=[];columns=np.zeros((1<<11,4),complex)
    index=0
    for left in (0,1):
        for right in (0,1):
            bitmask=(1<<left)|(1<<(2+right))|(F0<<4)
            columns[bitmask,index]=1
            current=left-right
            energy=(F0+current)**2/2
            phases.append(np.exp(-1j*alpha*energy))
            truth.append({'left_moving':left,'right_moving':right,'charges':2,
                          'current_on_shared_edge':current,'shared_electric_energy':energy,
                          'initial_and_restored_field':F0})
            index+=1
    actual=compiler.execute(echo,11,columns)
    wanted=columns*np.asarray(phases)[None,:]
    echo_error=float(np.linalg.norm(actual-wanted));assert echo_error<3e-10
    connected=phases[0]*phases[3]/(phases[1]*phases[2])
    assert abs(connected-np.exp(1j*alpha))<1e-14
    probabilities=[]
    for row in truth:
        left,right=row['left_moving'],row['right_moving']
        mask=(1<<left)|(1<<(2+right))|(F0<<4)
        probabilities.append(float(abs(final[mask])**2))
    expected=[(10+6*math.cos(alpha))/16]+[(2-2*math.cos(alpha))/16]*3
    residual=float(np.max(abs(np.array(probabilities)-expected)))
    assert residual<2e-11 and abs(np.vdot(final,final).real-1)<2e-10
    # Additive phase control preserves both one-charge phase calibrations.
    additive=phases.copy();additive[3]=phases[1]*phases[2]/phases[0]
    additive_state=np.sum(columns*np.asarray(additive)[None,:],axis=1)/2
    control=compiler.execute(rephase+analysis,11,additive_state)
    control_probability=float(abs(control[(1<<0)|(1<<2)|(F0<<4)])**2)
    assert abs(control_probability-1)<2e-11
    # Independently split fields make cross energy vanish. The same four
    # inputs give E(F0+j_left)+E(j_right), a separable phase truth table.
    independent=[np.exp(-1j*alpha*((F0+l)**2+r*r)/2) for l in (0,1) for r in (0,1)]
    independent_connected=independent[0]*independent[3]/(independent[1]*independent[2])
    assert abs(independent_connected-1)<1e-14
    bound=nc.massive.mass_coin.source_error.exact_template_bound
    base=json.loads((nc.BASE/'vector_field/COMPOSITION.json').read_text())
    root=json.loads((t.INPUT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    field_errors=json.loads((nc.BASE/'vector_field/NATIVE_FIELD.json').read_text())
    flip=F(root['primitive_operator_errors_exact']['flip'])
    prep_error=3*flip
    analysis_error=F()
    for a,b in ((0,1),(2,3)):
        prep_error+=bound(mb.givens_word(a,b,math.pi/4),nc.composition.givens_template(a,b,nc.composition.cert._pi(1,4)))
        analysis_error+=bound(mb.givens_word(a,b,-math.pi/4),nc.composition.givens_template(a,b,nc.composition.cert._pi(-1,4)))
    for mode,angle in ((1,F(alpha)*F(3,2)),(3,-F(alpha)/2)):
        analysis_error+=bound(mb.occupation_phase_word(mode,float(angle)),nc.composition.phase_template(mode,angle))
    # Conservative full-carrier embedding: all transport edges, all electric
    # registers and all source/reflection factors are charged by the supplied
    # full-matter/current bound, despite invariant zero spectators.
    source_error=2*F(base['complete_matter_current_layer_error_exact'])+288*F(field_errors['electric_per_edge_error_exact'])
    apparatus=nc.errors(1,True,0,0)
    total=F(apparatus['complete_native_TV_upper_exact'])-F(apparatus['preparation_operator_error_exact'])+prep_error+source_error+analysis_error
    tv=1-probabilities[0]
    result={'schema':'oph.two_charge.shared_field_edge_echo.v1','candidate_id':t.MODEL,
            'status':'executed_native_subword_interferometer','alpha':alpha,'mass_coin_role':'Its inverse is an explicit preparation interface if entering through a full macrostep; this local echo directly executes the frozen current/transport and electric subwords.',
            'initial_state':'One charge per endpoint, each in equal held/moving orbital; common F=P1 has edge0=1 and all others0, fixed G=2delta0. All four branches have N=2.',
            'source_words':'Frozen MODE_CURRENT_PROGRAMS16 field block; frozen full-CAR fock_swap_word; frozen FIELD_PROGRAM electric_word; literal inverse transport; explicit local native one-body rephasing and analysis.',
            'global_supported_modes':[0,1,16,17],'local_active_wires':11,
            'complete_source_wires':3203,'complete_readout_bits':3200,
            'truth_table':truth,'connected_phase':[float(connected.real),float(connected.imag)],
            'connected_phase_expected':'exp(i alpha)','connected_phase_angle':alpha,
            'local_native_echo_residual':echo_error,'native_interference_probability_residual':residual,
            'all_four_N2_event_probabilities':probabilities,
            'complete_law':'These four occupation pairs and restored field have the displayed ideal probabilities; all other ideal outcomes are zero. Native reset, export, padding and pointer complements have the complete bound below. No conditioning.',
            'unconditional_executed_norm_squared':float(np.vdot(final,final).real),
            'additive_one_charge_phase_control_probability_00':control_probability,
            'independent_field_connected_phase':[float(independent_connected.real),float(independent_connected.imag)],
            'event_law_TV_vs_additive_control':tv,
            'native_complete_TV_upper_exact':str(total),'native_complete_TV_upper':math.nextafter(float(total),math.inf),
            'native_pulses_executed_local':len(preparation)+len(echo)+len(rephase)+len(analysis),
            'source_preparation_word':preparation,'transport_word':transport,'electric_word':electric,
            'inverse_transport_word':compiler.inverse(transport),'rephase_word':rephase,'analysis_word':analysis,
            'full_embedding_cost':'Conservatively charge two full frozen matter/current layers, all288electric words, full3200-bit native readout and all3203blank/3200ready resets. The phase-readout interferometer is an additional intervention, not a free projector.',
            'scope':'A concrete connected two-charge quantum interaction in native shared-field subwords with a complete inverse/readout intervention. It is not the uninterrupted V(.17,.41) trajectory, a first-particle displacement measurement, reciprocal force, Coulomb law, physical mass or continuum result.',
            'source_snapshot_sha256':hashlib.sha256((t.HERE/'INPUT_SNAPSHOT.json').read_bytes()).hexdigest()}
    (t.HERE/'EDGE_ECHO.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('local_native_echo_residual','event_law_TV_vs_additive_control','native_complete_TV_upper','native_pulses_executed_local')}),flush=True)


if __name__=='__main__':run()
