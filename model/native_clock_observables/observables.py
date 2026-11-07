"""Complete local Pauli instruments using the admitted native clock services.

The fixed biased H66 carrier, source chart and six-level primitives are reused.
The feedback correction is an explicit observer protocol, not passive repair.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib,json
import numpy as np
import sympy as s

from model.local_repair_family import local_bias as local, mixing
from model.local_repair_family import modular_clock as clock, learned_flow
from model.native_program import interpreter
from model.native_modular_reference import cap_local_clock as cap, phase_readout as phase

DIMENSION=66
PROBE=(1,7)
AXES=('X','Y','Z')


def read_program():
    if clock.PROBE != PROBE or cap.SEAMS[0] != PROBE:
        raise ArithmeticError('the existing native clock probe changed')
    command={'op':'repair','selector':{'kind':'fixed','seam':list(PROBE)},'law':'E','repeat':1}
    result=interpreter.program(commands=[command])
    if tuple(result['parameters']['port_weights']) != (25,25)+(8,)*10:
        raise ArithmeticError('the fixed biased native source changed')
    return result


@dataclass(frozen=True)
class Branch:
    kind:str
    law:str|None
    coin:int|None
    raw:int|None
    contrast:int|None
    current:int|None
    weight:F
    pairs:tuple[tuple[int,int],...]


def branches():
    program=read_program()
    marker={(i,i,i,i):1 for i in range(66)}
    rows=[]
    for out in interpreter._outcomes(program,program['commands'][0],marker):
        rec=out.record
        weight=local.coin_probability(clock.parameters(),PROBE,rec.coin) if rec.law=='E' else F()
        pairs=[]
        for (j,i,k,l),coefficient in out.density.items():
            if j!=k or i!=l or coefficient!=weight:
                raise ArithmeticError('the literal native branch marker changed')
            pairs.append((i,j))
        rows.append(Branch(rec.kind,rec.law,rec.coin,rec.raw,rec.label,rec.current,weight,tuple(sorted(pairs))))
    if len(rows)!=25:
        raise ArithmeticError('all native hold, C and E slots must be retained')
    return tuple(rows)


def data():
    states=local.basis(); pairrows=cap.pairs(0)
    X=s.zeros(66);Y=s.zeros(66);Z=s.zeros(66);Q=s.zeros(66);J=s.zeros(66)
    for i,j,code_i,code_j in pairrows:
        if (states[i][1],states[i][7],states[j][1],states[j][7])!=(1,0,0,1):
            raise ArithmeticError('the complete source pair orientation changed')
        if any(states[i][k]!=states[j][k] for k in range(12) if k not in PROBE):
            raise ArithmeticError('a source analyzer changes the spectator')
        X[i,j]=X[j,i]=1;Y[i,j]=-s.I;Y[j,i]=s.I
        Z[i,i]=1;Z[j,j]=-1;Q[i,i]=Q[j,j]=1;J[i,j]=1;J[j,i]=-1
    if len(pairrows)!=10 or Q.trace()!=20:
        raise ArithmeticError('all ten spectator pairs are required')
    zero=s.eye(66)-Q;c=s.sqrt(2)/2
    rotations={'Z':s.eye(66),'X':zero+c*(Q+J),'Y':zero+c*(Q-s.I*X)}
    P={'Z':{-1:(Q-Z)/2,0:zero,1:(Q+Z)/2},
       'X':{-1:(Q-X)/2,0:zero,1:(Q+X)/2},
       'Y':{-1:(Q-Y)/2,0:zero,1:(Q+Y)/2}}
    return {'states':states,'pairs':pairrows,'branches':branches(),
            'X':s.ImmutableSparseMatrix(X),'Y':s.ImmutableSparseMatrix(Y),
            'Z':s.ImmutableSparseMatrix(Z),'Q':s.ImmutableSparseMatrix(Q),
            'zero':s.ImmutableSparseMatrix(zero),'swap':s.ImmutableSparseMatrix(zero+J),
            'rotations':{a:s.ImmutableSparseMatrix(U) for a,U in rotations.items()},
            'projectors':{a:{v:s.ImmutableSparseMatrix(M) for v,M in row.items()} for a,row in P.items()},
            'pi':s.ImmutableSparseMatrix(s.diag(*local.reference(clock.parameters())))}


def axis(value):
    if type(value) is not str or value not in AXES:
        raise ValueError('an axis X, Y or Z is required')
    return value


def branch_operator(d, row, selected, *, corrected=False):
    """Integer/algebraic branch factor; squared amplitude row.weight is separate."""
    selected=axis(selected)
    if type(corrected) is not bool:
        raise ValueError('an original Boolean correction selector is required')
    V=s.ImmutableSparseMatrix(66,66,{(j,i):1 for i,j in row.pairs})
    U=d['rotations'][selected]
    if not corrected:
        return V*U
    correction=d['swap'] if row.current else s.eye(66)
    return s.ImmutableSparseMatrix(U.H*correction*V*U)


def instrument_certificate(d=None):
    d=data() if d is None else d
    if d['branches']!=branches():
        raise ArithmeticError('the complete native branch labels, weights or maps changed')
    I=s.eye(66);zero=s.zeros(66)
    if any(M.H!=M for M in (d['X'],d['Y'],d['Z'])):
        raise ArithmeticError('an actual local observable is not Hermitian')
    if any(M*M!=d['Q'] for M in (d['X'],d['Y'],d['Z'])):
        raise ArithmeticError('the even endpoint sector was not retained')
    if d['X']*d['Y']-d['Y']*d['X']!=2*s.I*d['Z']:
        raise ArithmeticError('the Pauli commutator orientation changed')
    if d['rotations']['X']**2!=d['swap']:
        raise ArithmeticError('the recorded displacement correction is not two existing rotations')
    rows=[]
    for selected in AXES:
        U=d['rotations'][selected]
        if U.H*U!=I or U.H*d['Z']*U!=d[selected]:
            raise ArithmeticError('the complete analyzer effect failed')
        effects={v:s.zeros(66) for v in (-1,0,1)}
        normal=s.zeros(66)
        for row in d['branches']:
            K=branch_operator(d,row,selected)
            normal+=row.weight*K.H*K
            if row.pairs:
                if row.law!='E' or row.contrast not in (-1,0,1):
                    raise ArithmeticError('a nonzero source outcome has an invalid contrast')
                effects[row.contrast]+=row.weight*K.H*K
                restored=branch_operator(d,row,selected,corrected=True)
                sign=-1 if row.current==-1 else 1
                if restored!=sign*d['projectors'][selected][row.contrast]:
                    raise ArithmeticError('a fine retained branch is not its full-carrier Lueders map')
                rows.append((selected,row.coin,row.raw,row.contrast,row.current,sign))
            elif K!=zero:
                raise ArithmeticError('an empty native address has a nonzero operator')
        if normal!=I or effects!=d['projectors'][selected]:
            raise ArithmeticError('complete native effects or normalization failed')
    return {'source_dimension':66,'raw_slots_per_axis':25,
            'zero_native_slots_per_axis':sum(not row.pairs for row in d['branches']),
            'active_dimension':20,'even_zero_dimension':46,'spectator_pairs':10,
            'exact_full_carrier_fine_branch_identities':len(rows),
            'matrix_units_per_identity':4356,'fine_branch_sign_table':[list(row) for row in rows],
            'raw_effects':'P_axis,delta with delta=-1,0,+1; P_0 retains both even sectors',
            'raw_update':'destructive biased placement; not a Lueders instrument',
            'correction':'undo recorded displacement by two balanced cap rotations when current != 0, then undo analyzer',
            'corrected_fine_map':'weight_coin * P_axis,delta rho P_axis,delta; scalar sign cancels',
            'source_and_reference_domain':'every operator on H66 tensor any untouched finite reference',
            'noncommuting_effect_norm':'||[P_X,+,P_Y,+]||=1/2',
            'new_primitive_gate_or_action':False,'observer_protocol_design_is_declared':True}


def uncertainty_certificate(d=None):
    d=data() if d is None else d
    rho=d['pi'];P=d['Q'];X,Y,Z=(d[k] for k in AXES)
    mean={k:s.trace(rho*d[k]) for k in AXES}
    variance={k:s.trace(rho*d[k]*d[k])-mean[k]**2 for k in AXES}
    comm=s.trace(rho*(X*Y-Y*X))
    if s.trace(rho)!=1 or any(v<0 for v in rho.diagonal()):
        raise ArithmeticError('the native reference is not a density matrix')
    active=s.trace(rho*P)
    if active!=F(3201,7505) or mean['Z']!=F(1649,7505) or comm!=2*s.I*mean['Z']:
        raise ArithmeticError('the source-derived uncertainty witness changed')
    lower=s.simplify(s.conjugate(comm)*comm/4)
    margin=variance['X']*variance['Y']-lower
    if margin<0:
        raise ArithmeticError('finite state Robertson inequality failed')
    return {'reference':'actual biased canonical rho_*; no odd-sector conditioning',
            'mean':{k:str(v) for k,v in mean.items()},'variance':{k:str(v) for k,v in variance.items()},
            'active_probability':str(active),'commutator_expectation':str(comm),
            'Robertson_lower':str(lower),'Robertson_margin':str(margin),
            'finite_theorem':'EventAlgebra.Robertson.finite_state_robertson_commutator',
            'position_momentum_CCR_identified':False,'universal_operational_Born_effect_quotient_supplied':False}


def native_output(row, density):
    """Arbitrary source/reference block after one exact native raw outcome."""
    density=np.asarray(density,complex)
    if density.ndim!=2 or not density.shape[0] or density.shape[0]%66:
        raise ValueError('a nonempty H66 tensor reference operator is required')
    ref=density.shape[0]//66
    if density.shape!=(66*ref,66*ref):
        raise ValueError('a square H66 tensor reference operator is required')
    out=np.zeros_like(density)
    if row.pairs:
        inputs=np.array([i*ref+r for i,_ in row.pairs for r in range(ref)])
        outputs=np.array([j*ref+r for _,j in row.pairs for r in range(ref)])
        out[np.ix_(outputs,outputs)]=float(row.weight)*density[np.ix_(inputs,inputs)]
    return out


def ideal_complete_output(d, density, selected, *, corrected=False):
    selected=axis(selected);ref=len(density)//66;U=np.kron(np.asarray(d['rotations'][selected],complex),np.eye(ref))
    prepared=U@density@U.conj().T;result={}
    for k,row in enumerate(d['branches']):
        out=native_output(row,prepared)
        if corrected:
            correction=np.kron(np.asarray(d['swap'] if row.current else s.eye(66),complex),np.eye(ref))
            out=U.conj().T@correction@out@correction.conj().T@U
        result[k]=out
    return result


def words():
    return {'X':cap.rotation_word(0,1),'X_inverse':cap.rotation_word(0,-1),
            'Y':phase.quadrature_word(0),'Y_inverse':cap.compiler.inverse(phase.quadrature_word(0))}


def _inverse_y():
    result=np.eye(512,dtype=complex)
    for _,_,a,b in reversed(cap.pairs(0)):
        tape=cap.compiler.inverse(phase.quadrature_pair_word(a,b))
        result=cap.compiler.execute(tape,9)@result
    return result


def source_arrays():
    return {'X':cap.rotation(0,1),'X_inverse':cap.rotation(0,-1),
            'Y':phase.quadrature_rotation(0),'Y_inverse':_inverse_y()}


def multiplication_error(pulses, blocks=10):
    """Outward bound from returned products to exact products of stored pulses.

    IEEE binary64 round-to-nearest complex products: local dots have at most
    four entries and pair-block composition dots have 512. Gamma_8d covers
    their real operations; a separate absolute allowance covers subnormals.
    """
    if type(pulses) is not int or not 0 < pulses <= 12000 or blocks != 10:
        raise ValueError('the certified ten-pair word resource envelope is required')
    unit=F(1,2**53)
    small=32*unit/(1-32*unit);large=4096*unit/(1-4096*unit)
    gate=F(2,10**11)
    envelope=pulses*(gate+64*small*(1+gate))+blocks*1024*large+F(1,2**900)
    if envelope>=F(1,4):
        raise ArithmeticError('the product norm envelope exceeds its proof domain')
    error=128*pulses*small*(1+gate)+2048*blocks*large+F(1,2**900)
    return {'operator_error_upper':str(error),'pulse_count':pulses,
            'pair_block_products':blocks,'floating_unit_roundoff':str(unit),
            'local_complex_dot_gamma':str(small),'full512_complex_dot_gamma':str(large),
            'intermediate_Frobenius_bound':'32','propagation_gain_bound':'2',
            'IEEE_binary64_round_to_nearest_arithmetic_assumed':True,
            'matrix_exponential_vs_exact_analytic_exponential_included':False}


def compile_certificate(arrays=None):
    arrays=source_arrays() if arrays is None else arrays
    tape=words();result={}
    for name,word in tape.items():
        admitted=interpreter.source_word(word,9)
        error=(cap.certify_rotation(0,1 if name=='X' else -1,arrays[name]) if name.startswith('X')
               else phase.certify_quadrature(0,arrays[name] if name=='Y' else arrays[name].conj().T))
        if error>cap.ROTATION_GATE:
            raise ArithmeticError('the literal source analyzer failed its full512 certificate')
        gate,failure=clock._processor_failure_bound(len(word))
        menu={(p['kind'],p['angle'],len(p['targets'])) for p in word}
        for kind,angle,bits in menu:
            leakage,gain=clock._processor_bounds(cap.compiler.native_pulse(kind,angle),1<<bits)
            if max(leakage,gain)>gate:
                raise ArithmeticError('a complete processor export failed its per-pulse certificate')
        returned_error=cap.ROTATION_GATE
        arithmetic=multiplication_error(len(word))
        error=returned_error+F(arithmetic['operator_error_upper'])
        delta=(error*(2+error)+failure)/2
        result[name]={'operator_error_upper':str(error),'returned_array_operator_error_upper':str(returned_error),
                      'product_arithmetic':arithmetic,'first_processor_failure_mass_upper':str(failure),
                      'complete_halfdiamond_upper':str(delta),'channel_norm_upper':str((1+error)**2+failure),
                      'word_digest':interpreter.digest(admitted),'resources':admitted['declared_resources'],
                      'unique_processor_pulses_checked':len(menu)}
    eX=F(result['X']['operator_error_upper']);eY=F(result['Y']['operator_error_upper'])
    a=eX*(2+eX);b=eY*(2+eY)
    lower=F(1,2)-2*(a+b+a*b)
    if lower<=0:
        raise ArithmeticError('compiled effects have no certified incompatibility margin')
    return {'words':result,'full512_addresses_checked':512,'quantum_work_qubits':9,
            'raw_noncommuting_effect_commutator_norm_lower':str(lower),
            'all_code_complement_first_abort_labels_retained_in_bound':True,
            'every_long_word_abort_matrix_materialized':False,
            'physical_hardware_error_included':False}


def composed_error(cert, stages):
    error=F();gain=F(1);pulses=0
    for stage in stages:
        row=cert['words'][stage];g=F(row['channel_norm_upper'])
        error=g*error+F(row['complete_halfdiamond_upper']);gain*=g
        pulses+=row['resources']['native_pulses']
    return {'complete_halfdiamond_upper':str(error),'channel_norm_upper':str(gain),
            'source_pulses_upper':pulses,'proper_code_load_export_calls_upper':2*pulses}


def protocol_certificate(cert):
    result={}
    for selected in AXES:
        stages=[] if selected=='Z' else [selected]
        raw=composed_error(cert,stages)
        restored=composed_error(cert,stages+['X','X']+([] if selected=='Z' else [selected+'_inverse']))
        result[selected]={'raw':raw,'corrected_worst_branch':restored,
                          'native_detector_attempts':1,'native_raw_slots':25,
                          'current_feedback_uses_only_retained_record':True}
    return result


def compiled_output(d, density, selected, arrays, *, corrected=False):
    """Literal successful-code arrays with every address-export abort retained.

    Proper-code primitive first failures have separate analytic complete bounds.
    No surviving branch or address export is normalized.
    """
    selected=axis(selected)
    if type(corrected) is not bool:
        raise ValueError('an original Boolean correction selector is required')
    density=np.asarray(density,complex)
    if density.ndim!=2 or density.shape[0]==0 or density.shape[0]%66 or density.shape[1]!=density.shape[0]:
        raise ValueError('a square H66 tensor reference operator is required')
    result={}
    prepared=density
    if selected!='Z':
        prepared,abort=clock._export(arrays[selected],prepared)
        result[('analyzer_export_abort',)]=abort
    for k,row in enumerate(d['branches']):
        out=native_output(row,prepared)
        if corrected:
            stages=(['X','X'] if row.current else [])+([] if selected=='Z' else [selected+'_inverse'])
            for index,name in enumerate(stages):
                out,abort=clock._export(arrays[name],out)
                result[('correction_export_abort',k,index,name)]=abort
        result[('native',k)]=out
    return result


def ideal_output_with_same_exports(d, density, selected, *, corrected=False):
    ref=len(density)//66
    blocks=ideal_complete_output(d,density,selected,corrected=corrected)
    result={('native',k):block for k,block in blocks.items()}
    if selected!='Z':
        result[('analyzer_export_abort',)]=np.zeros((ref,ref),complex)
    if corrected:
        for k,row in enumerate(d['branches']):
            stages=(['X','X'] if row.current else [])+([] if selected=='Z' else [selected+'_inverse'])
            for index,name in enumerate(stages):
                result[('correction_export_abort',k,index,name)]=np.zeros((ref,ref),complex)
    return result


def modular_clock_certificate(d=None):
    d=data() if d is None else d
    weights=local.reference(clock.parameters())
    for i,j,_,_ in d['pairs']:
        if weights[i]/weights[j]!=F(25,8):
            raise ArithmeticError('the accessible pair has a different native modular gap')
    # Polynomial identities in cos(t), sin(t), valid for every real t.
    c,t=s.symbols('c t',real=True)
    ux=s.Matrix([[c-s.I*t,0],[0,1]])
    xx=s.Matrix([[0,1],[1,0]]);yy=s.Matrix([[0,-s.I],[s.I,0]])
    if s.expand(ux.H*xx*ux-(c*xx-t*yy))!=s.zeros(2):
        raise ArithmeticError('the source modular Heisenberg sign changed')
    # U_s has diagonal ratio exp(-i*s*log(25/8)); its Heisenberg action
    # is cos(s*Delta) X - sin(s*Delta) Y.
    return {'generator':'K=-log(rho_*); U_s=rho_*^(-i*s)=exp(i*s*K)',
            'pair_population_ratio':'25/8','Delta':'log(25/8)',
            'Heisenberg_X':'cos(s*Delta)*X - sin(s*Delta)*Y',
            'Heisenberg_Y':'sin(s*Delta)*X + cos(s*Delta)*Y',
            'source_program':'model.local_repair_family.learned_flow.source_program',
            'certified_existing_clock_parameters':['1/2','1'],
            'arbitrary_real_parameter_claims_only_the_ideal_algebra':True,
            'working_modular_K_identified_with_record_transfer_minus_log_S':False,
            'physical_time_identified':False}


def preparation_certificate():
    result=mixing.mixing_bound_summary(clock.parameters(),1600)
    if result['half_diamond_upper']!='542159/1000000000000':
        raise ArithmeticError('the source preparation bound changed')
    return {**result,'joint_preparation_archive_included':False,
            'future_compiled_protocol_halftrace_bound':'channel_norm_upper*preparation_halfdiamond + protocol_complete_halfdiamond',
            'arbitrary_prior_reference_retained':True,
            'protocol_adaptation_to_preparation_archive_covered':False}
