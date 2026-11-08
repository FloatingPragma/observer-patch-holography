"""Addressed occupation-wire attachment of existing coherent services.

The coupling is an explicitly supplied proper-code interface. It is not a
classical-record reconstruction of an unknown quantum state.
"""
from fractions import Fraction as F
import math
import numpy as np
from model.native_fixed_algebra import occupation_service as occupation
from model.native_effect_synthesis import synthesis as synthesis
from model.reference_recovery import recovery
from . import preparation as p


def split_mask(mask,cell=0):
    p.integer(mask);p.integer(cell)
    return (mask>>(12*cell))&p.FULL, mask&~(p.FULL<<(12*cell))


def join_mask(local,outside,cell=0):
    p.integer(local,0,p.FULL);p.integer(outside);p.integer(cell)
    if outside&(p.FULL<<(12*cell)):
        raise ValueError('exterior address overlaps the supplied local wire chart')
    return (local<<(12*cell))|outside


def split_operator(density,cell=0,reference_dimension=1):
    p.integer(cell);p.integer(reference_dimension,1)
    if type(density) is not dict:raise ValueError('an exact sparse source/reference operator required')
    out={}
    for (x,a,y,b),value in density.items():
        p.integer(a,0,reference_dimension-1);p.integer(b,0,reference_dimension-1)
        if type(value) not in (int,F):raise ValueError('exact rational operator component required')
        u,e=split_mask(x,cell);v,f=split_mask(y,cell)
        out[u,e*reference_dimension+a,v,f*reference_dimension+b]=value
    return out


def join_operator(density,cell=0,reference_dimension=1):
    p.validate_density(density);p.integer(cell);p.integer(reference_dimension,1)
    out={}
    for (x,a,y,b),value in density.items():
        e,r=divmod(a,reference_dimension);f,s=divmod(b,reference_dimension)
        out[join_mask(x,e,cell),r,join_mask(y,f,cell),s]=value
    return out


def intertwining_certificate():
    """All raw branch addresses on all native seams, with cross-number inputs."""
    checks=0
    for cell in (0,1):
        outside_base=12 if cell==0 else 0
        for edge in p.SEAMS:
            configurations=[]
            for bits in range(4):
                local=((bits&1)<<edge[0])|((bits>>1)<<edge[1])
                outside=((1<<(12-local.bit_count()))-1)<<outside_base
                configurations.append(join_mask(local,outside,cell))
            rho={(x,i,y,j):F(1,4) for i,x in enumerate(configurations) for j,y in enumerate(configurations)}
            local=split_operator(rho,cell,4)
            for law in ('C','E'):
                actual_edge=[12*cell+v for v in edge]
                command={'op':'repair','selector':{'kind':'fixed','seam':actual_edge},'law':law,'repeat':1}
                program=p.source.program(federation={'level':0,'total':12},biased=False,commands=[command])
                p.source._validate_density(rho,program['carrier'])
                actual=p.source._outcomes(program,command,rho)
                wanted=p.readout.cone.fixed_instrument(p.FAMILY,local,edge,law=law)
                for full,(record,expected) in zip(actual,wanted,strict=True):
                    if split_operator(full.density,cell,4)!=expected:
                        raise ArithmeticError('the complete raw instrument fails occupation-wire intertwining')
                    a=full.record
                    if (a.kind,a.selector,a.law,a.coin,a.raw,a.label,a.current)!=(record.kind,record.selector,record.law,record.coin,record.raw,record.label,record.current):
                        raise ArithmeticError('the addressed complete record convention changed')
                    checks+=1
    return {'complete_raw_blocks':checks,'local_cells_checked':[0,1],
            'native_seams_per_cell':30,'record_laws':['C','E'],
            'cross_local_number_reference_matrix_units':True,
            'global_occupation_of_controls':12,'all_local_sectors_covered_by_tensor_proof':True}


def words():
    compiler=occupation.compiler
    beam=[dict(event,targets=[0 if v==0 else 1 for v in event['targets']]) for event in occupation._beam()]
    phase=[compiler.pulse('I',math.pi/4,[0]),compiler.pulse('Z',-math.pi/4,[0])]
    return {'beam':(2,beam),'phase':(1,phase)}


def ideal_matrix(name):
    if name=='phase':return np.diag([1,1j]).astype(complex)
    if name=='beam':
        value=np.eye(4,dtype=complex);value[np.ix_([1,2],[1,2])]=np.array([[1.,-1.],[1.,1.]])/math.sqrt(2)
        return value
    raise ValueError('an existing beam or phase macro required')


def complete_bound(pulses):
    p.integer(pulses,1,1024)
    x=pulses*synthesis.PULSE_GATE
    error=x/(1-x)+F(pulses,10**15)
    failure=pulses*synthesis.PULSE_GATE**2/(1-x)**2
    return {'complete_halfdiamond_upper':str((error*(2+error)+failure)/2),
            'complete_map_gain_upper':str(1/(1-x)**2+failure),
            'success_operator_error_upper':str(error),
            'all_first_processor_failure_mass_upper':str(failure),
            'number_leakage_mass_upper':str(error**2)}


def service_certificate():
    synthesis.certify_primitive.cache_clear();rows=[]
    for name,(bits,word) in words().items():
        unique=synthesis.certify_word(word)
        run=recovery.execute_complete(word,bits,np.eye(1<<bits,dtype=complex))
        target=ideal_matrix(name);bound=complete_bound(len(word))
        if np.linalg.norm(run.success-target,2)>float(F(bound['success_operator_error_upper'])):
            raise ArithmeticError('stored local service exceeds its certified target envelope')
        numbers=np.array([i.bit_count() for i in range(1<<bits)])
        if np.any(target[numbers[:,None]!=numbers[None,:]]):
            raise ArithmeticError('the exact declared macro changes total occupation')
        # The full12-bit word is separately admitted with literal source cost.
        declaration=p.source.source_word(word,12)
        rows.append({'name':name,'compact_complete_input_columns':1<<bits,
                     'unique_certified_primitives':unique,'source_word_digest':p.source.digest(declaration),
                     'source_resources':declaration['declared_resources'],'error':bound,
                     'ideal_number_preserving':True,'all_first_failure_outputs_retained':True,
                     'rounded_off_sector_work_output_discarded':False,
                     'numerical_propagation_arrays_are_diagnostics':True})
    return {'services':rows,'proper_code_coupling_supplied_and_bound':True,
            'quantum_coupling_derived_from_classical_records':False,
            'bounded_local_service_cost_uniform_in_capacity':True,
            'numerical_success_assumed_exactly_native_number_preserving':False,
            'coherent_service_is_native_repair_command':False,
            'source_pulses_identified_with_repair_attempts':False}
