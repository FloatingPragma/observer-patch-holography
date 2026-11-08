"""Versioned dyadic source control with an exact occupation-sector return.

This is a separately declared M1 coherent-control service. It does not
amend the frozen native_program source_word alphabet or its V2 verdict.
"""
from fractions import Fraction as F
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.certified_reference import source
from model.native_effect_synthesis import intervals as iv
from model.native_program import interpreter as native
from model.observer.source_scalar import compiler
from . import bell

SCHEMA='oph.derived_occupation_exchange.v1'
ROOT=Path(compiler.ROOT)
GATE=F(1,10**13)


def exact_identity():
    spec=source.specification()
    coefficients=(F(0),)*4+(F(-1,2),F(-1,2),F(1,2),F(1,2))+(F(0),)*4
    result=[]
    for i in range(6):
        row=[]
        for j in range(6):
            cell=tuple(tuple(sum(c*g[i][j][part][irr] for c,g in zip(coefficients,spec,strict=True))
                             for irr in (0,1))for part in (0,1))
            wanted=1 if (i,j) in ((1,2),(4,5)) else -1 if (i,j) in ((2,1),(5,4)) else 0
            if cell!=((F(wanted),F()),(F(),F())):
                raise ArithmeticError('four-port exact native response identity failed')
            row.append(wanted)
        result.append(tuple(row))
    return tuple(result)


def declare(angle,targets=(0,1),cell=0,level=0):
    if type(angle)not in (int,F):raise ValueError('exact dyadic angle required')
    angle=F(angle)
    if (abs(angle)>1 or angle.denominator & (angle.denominator-1)
            or F(float(angle))!=angle or F(float(angle/2))!=angle/2):
        raise ValueError('bounded binary64 angle and exact half-angle controls required')
    if (type(targets)is not tuple or len(targets)!=2 or len(set(targets))!=2
            or any(type(v)is not int or not 0<=v<12 for v in targets)
            or type(cell)is not int or type(level)is not int or not 0<=level<=2
            or not 0<=cell<(20 << (2*level))):
        raise ValueError('two addressed local occupation wires and exact cell required')
    controls=(F(),)*4+(-angle/2,-angle/2,angle/2,angle/2)+(F(),)*4
    value={'schema':SCHEMA,'cell':cell,'federation_level':level,'local_targets':list(targets),'angle':str(angle),
            'integrated_controls':[str(x)for x in controls],
            'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                             for p in (compiler.RESPONSE,compiler.CONTROL_RECEIPT)},
            'resources':{'native_response_pulses':1,'proper_code_load_events':1,
                         'proper_code_export_events':1,'six_level_blank_preparations':1,
                         'proper_code_dimension':4,'all_export_complement_rows':2,
                         'local_source_work_bits':12,'native_repair_attempts':0,
                         'angle_sign_bits':1,'angle_numerator_bits':max(1,abs(angle.numerator).bit_length()),
                         'angle_denominator_bits':angle.denominator.bit_length(),
                         'integrated_control_payload_bits':12*64,
                         'cell_address_bits':max(1,cell.bit_length()),
                         'local_target_address_bits':8,
                         'level_address_bits':max(1,level.bit_length()),
                         'canonical_declaration_bytes':0},
            'coherent_control_access':'existing M1 integrated port controls and proper4 load/export',
            'frozen_native_program_schema_changed':False}
    # Charge the complete canonical declaration, including its own decimal
    # byte-count field. This small monotone fixed point is unambiguous.
    for _ in range(8):
        size=len(native.canonical(value).encode())
        if size==value['resources']['canonical_declaration_bytes']:return value
        value['resources']['canonical_declaration_bytes']=size
    raise ArithmeticError('canonical description size did not stabilize')


def authenticate(value):
    if type(value)is not dict:raise ValueError('versioned service declaration required')
    try:
        expected=declare(F(value['angle']),tuple(value['local_targets']),value['cell'],value['federation_level'])
    except (KeyError,TypeError,ZeroDivisionError) as error:
        raise ValueError('complete exact service declaration required')from error
    # JSON scalar types remain distinct, including bool/int and float/int.
    encode=lambda x:json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)
    if encode(value)!=encode(expected):raise ValueError('service differs from exact declaration and source bindings')
    return expected


def execute(value):
    value=authenticate(value);theta=F(value['angle']);identity=exact_identity()
    # Exact generator assembly precedes numerical evaluation. Every entry
    # is an exactly representable dyadic, so no target matrix is substituted.
    generator=np.array([[float(theta*x)for x in row]for row in identity],dtype=float)
    full=expm(generator)
    if not np.isfinite(full).all():raise ArithmeticError('finite complete processor response required')
    active={(0,0),(3,3),(1,1),(1,2),(2,1),(2,2),(4,4),(4,5),(5,4),(5,5)}
    if any(full[i,j]!=0 for i in range(6)for j in range(6)if (i,j)not in active):
        raise ArithmeticError('stored response lacks the exact structural return required for native handoff')
    # Compare the complete M6 response to its independently enclosed source
    # exponential, including the rotation on complement coordinates4,5.
    sine,cosine=iv.sincos(theta);zero=iv.exact(0);one=iv.exact(1)
    error=F()
    for i in range(6):
        row=F()
        for j in range(6):
            target=(one if i==j and i in (0,3) else cosine if i==j else
                    sine if (i,j)in ((1,2),(4,5)) else iv.neg(sine)
                    if (i,j)in ((2,1),(5,4)) else zero)
            row+=iv.scalar_error(full[i,j],target)
        error=max(error,row)
    if error>GATE:raise ArithmeticError('literal numerical response exceeds exact full-source error gate')
    # Row and column bounds agree because all two-dimensional blocks have
    # the same symmetric absolute-error shape. Check rather than assume it.
    column_error=F()
    for j in range(6):
        total=F()
        for i in range(6):
            target=(one if i==j and i in (0,3) else cosine if i==j else
                    sine if (i,j)in ((1,2),(4,5)) else iv.neg(sine)
                    if (i,j)in ((2,1),(5,4)) else zero)
            total+=iv.scalar_error(full[i,j],target)
        column_error=max(column_error,total)
    error=max(error,column_error)
    if error>GATE:raise ArithmeticError('full-source row/column gate failed')
    success=tuple(tuple(F(float(full[i,j]))for j in range(4))for i in range(4))
    complements=tuple(tuple(F(float(full[i,j]))for j in range(4))for i in (4,5))
    if any(x for row in complements for x in row):
        raise ArithmeticError('proper-code complement leakage cannot be returned as a native state')
    return {'declaration':value,'full_processor':tuple(tuple(F(float(x))for x in row)for row in full),
            'success':success,'complement_rows':complements,
            'full_complement_rows':tuple(tuple(F(float(full[i,j]))for j in range(6))for i in (4,5)),
            'operator_error_upper':GATE,'complete_halfdiamond_upper':GATE*(2+GATE)/2,
            'complete_map_gain_upper':(1+GATE)**2,
            'exact_stored_number_preservation':True,
            'literal_exponential_error_diagnostic_upper':error}


def success_interval(theta=F(25,32)):
    # Complementary preparations give equal even amplitudes cs exactly,
    # without a balanced single-qubit angle or a pi payload.
    sine,_=iv.sincos(2*theta)
    square=iv.divide(iv.mul(sine,sine),2)
    return F(square[0],iv.Q),F(square[1],iv.Q)


def lift(value,density,carrier):
    """Execute addressed proper4 load/response/export on sparse native masks.

    The source-to-code map is |n_a,n_b> -> |2 n_a+n_b>; all spectators
    and reference indices are unchanged. Complete complement rows are
    returned as separate zero branches, not projected out of a nonzero map.
    """
    result=execute(value);c=native.carrier(carrier)
    if c['kind']!='federation' or c['level']!=result['declaration']['federation_level']:
        raise ValueError('service addresses require their authenticated federation level')
    rho=native._validate_density(density,c)
    a,b=(12*result['declaration']['cell']+v for v in result['declaration']['local_targets'])
    clear=~((1<<a)|(1<<b));operator=result['success'];output=defaultdict(F)
    def column(mask):
        label=2*((mask>>a)&1)+((mask>>b)&1)
        for i in range(4):
            if operator[i][label]:
                target=(mask&clear)|((i>>1)<<a)|((i&1)<<b)
                if target.bit_count()!=mask.bit_count():
                    raise ArithmeticError('addressed source export changed total occupation')
                yield target,operator[i][label]
    for (x,r,y,s),v in rho.items():
        for xx,aa in column(x):
            for yy,bb in column(y):output[xx,r,yy,s]+=v*aa*bb
    success=native._validate_density({k:v for k,v in output.items()if v},c)
    return {'success':success,'complements':({'label':4,'density':{}},{'label':5,'density':{}}),
            'absolute_targets':(a,b),'declaration':result['declaration'],
            'complete_halfdiamond_upper':result['complete_halfdiamond_upper'],
            'complete_map_gain_upper':result['complete_map_gain_upper']}


def coupled_bell(level=0):
    _,edge,spare=bell.layout(level);masks=bell.code_masks(edge,spare)
    first=declare(F(-25,32),(edge[0]%12,spare[0]%12),edge[0]//12,level)
    second=declare(F(25,32),(edge[1]%12,spare[1]%12),edge[1]//12,level)
    command={'op':'repair','selector':{'kind':'fixed','seam':list(edge)},'law':'E','repeat':1}
    program=native.program(federation={'level':level,'total':2},biased=False,commands=[command],capacity=1,overflow='export_oldest')
    # Literal known dual-rail input A0 B1; this is a preparation interface,
    # not a claim that every input native branch prepares it.
    density={(masks[1],0,masks[1],0):F(1)}
    first_run=lift(first,density,program['carrier'])
    second_run=lift(second,first_run['success'],program['carrier'])
    density=second_run['success']
    history=native.execute(program,density)
    rows=history.branches
    if (history.status!='COMPLETE' or history.attempted_steps!=1
            or rows!=native._outcomes(program,command,density)
            or history.resources['maximum_FIFO_occupied']!=1):
        raise ArithmeticError('public native admission/history handoff did not complete')
    total=native.trace(density)
    if sum(native.trace(row.density)for row in rows)!=total:
        raise ArithmeticError('complete native handoff did not preserve its unnormalized input mass')
    accepted=[row for row in rows if row.record.law=='E' and row.record.raw==0]
    success=sum(native.trace(row.density)for row in accepted)
    per_pulse=GATE
    product_error=2*per_pulse+per_pulse**2
    halfdiamond=product_error*(2+product_error)/2
    interval=success_interval()
    if not interval[0]-2*halfdiamond<=success<=interval[1]+2*halfdiamond:
        raise ArithmeticError('complete Bell score exceeds its actual source implementation bound')
    if any(set((x,y)for x,_,y,_ in row.density)-{(masks[i],masks[j])for i in (0,3)for j in (0,3)}for row in accepted):
        raise ArithmeticError('Bell output has an unaccounted cell-charge sector')
    exact_Bell=all(len(set(row.density.values()))<=1 for row in accepted)
    return {'schema':'oph.derived_boundary_bell.v1','level':level,'edge':edge,'spares':spare,
            'source_service_declarations':[first,second],
            'native_resources':program['declared_resources'],
            'public_native_history_resources':history.resources,
            'public_native_status':history.status,
            'ideal_success_interval':interval,'stored_success_mass':success,
            'stored_complete_output_mass':total,'complete_halfdiamond_error_upper':halfdiamond,
            'number_leakage_mass':F(0),'proper_code_complement_rows_retained':4,
            'native_raw_slots_retained':len(rows),'stored_conditional_Bell_exact':exact_Bell,
            'native_handoff_executed':True,'normalization_or_number_projection_inserted':False,
            'known_initial_logical_basis':'A0 B1',
            'local_source_service_pulses':2,'frozen_V2_admission_claimed':False}
