"""Exact-rational M6 control stream for massive, phase-dressed LOCAL12.

Bounds retain every coherent export/first-failure branch of each pulse.
Reset, preparation, readout and physical clock costs are separate services.
"""
from collections import Counter
from fractions import Fraction as F
import importlib.util
import sys

import larger_walk as walk

a=walk.adapter


def coin_template(h,mass):
    h,mass=walk.mass_source.parameters(h,mass); mu=h*mass
    u=a.uniform_template(); word=a.negate_template(u)
    for mask in range(16):
        qpi=F(15,16) if mask==0 else F(-1,16)
        qmu=7*mu/8 if mask==0 else -mu/8
        word+=a.pauli(''.join('Z' if mask>>q&1 else 'I' for q in range(4)),a.angle(qpi,0,qmu))
    return word+u


def phase_plan(size,seam_phase=None,port_potential=None,h=F(1)):
    """Sparse exact-rational controls keyed by seam id or physical port.

    Compute a basis predicate into a clean flag, apply exp(i theta |1><1|),
    uncompute. This declares extra clean workspace and includes its reset
    obligation. Nonzero pad controls are never inferred.
    """
    info=a.layout(size);g=walk.geometry(size);d=info['data_bits']
    flag=info['memory_bits'];scratch=tuple(range(flag+1,flag+d-1))
    blocks=[]
    for controls,kind in ((seam_phase or {},'seam'),(port_potential or {},'port')):
        for key,value in sorted(controls.items()):
            if type(key)is not int or type(value) not in (int,F):raise ValueError('integer keys and exact rational controls required')
            value=F(value)
            if not value:continue
            if kind=='seam':
                if not 0<=key<len(g.graph.seams):raise ValueError('seam id outside graph')
                targets=((int(g.left[key]),value),(int(g.right[key]),-value))
                masks=[({j:(address>>j)&1 for j in range(d)},theta) for address,theta in targets]
            else:
                if not 0<=key<g.graph.ports:raise ValueError('physical port outside graph')
                address=int(g.addresses[key]); masks=[({j+4:(address>>j)&1 for j in range(d-4)},-F(h)*value)]
            for values,theta in masks:
                gates=a.frozen.multi_control(values,flag,scratch)
                blocks.append({'predicate':values,'compute':gates,'theta':theta,'kind':kind,
                               'phase':[('I',(flag,),a.angle(0,0,theta/2)),
                                        ('Z',(flag,),a.angle(0,0,-theta/2))]})
    return {'size':size,'base_memory_bits':info['memory_bits'],
            'memory_bits':flag+d-1 if blocks else flag,'flag':flag,'scratch':scratch,
            'blocks':blocks,'clean_extra_workspace_bits':d-1 if blocks else 0}


def gates_template(gates):
    primitives=a.primitive_templates()
    for name,targets in gates:
        for kind,local,angle in primitives[name]:
            yield kind,tuple(targets[j] for j in local),angle


def phase_template(plan,kind):
    for block in plan['blocks']:
        if block['kind']!=kind:continue
        yield from gates_template(block['compute'])
        yield from block['phase']
        yield from gates_template(reversed(block['compute']))


def template_stream(size,h=F(1),mass=F(0),seam_phase=None,port_potential=None):
    plan=phase_plan(size,seam_phase,port_potential,h)
    # P commutes with C; seam phase is immediately before the source route.
    yield from phase_template(plan,'port')
    yield from coin_template(h,mass)
    yield from phase_template(plan,'seam')
    yield from gates_template(a.reversible_program(size)['gates'])


def pulse_stream(size,h=F(1),mass=F(0),precision_bits=64,seam_phase=None,port_potential=None):
    for kind,targets,angle in template_stream(size,h,mass,seam_phase,port_potential):
        yield a.certified_pulse(kind,targets,angle,precision_bits)


def certificate(size,h=F(1),mass=F(0),precision_bits=64,seam_phase=None,port_potential=None):
    """Product telescoping bound on the complete prepared-code isometry.

    The certificate counts repeated templates; it does not execute an ideal
    gate list and label that execution native. Each exported pulse carries
    all twelve actual rational integrated native controls.
    """
    walk.verify_sources()
    plan=phase_plan(size,seam_phase,port_potential,h)
    counts=Counter(a.reversible_program(size)['gates'][j][0] for j in range(len(a.reversible_program(size)['gates'])))
    primitives=a.primitive_templates()
    bound=a.template_bound(coin_template(h,mass),precision_bits)
    pulses=len(coin_template(h,mass))
    for name,count in counts.items():
        bound+=count*a.template_bound(primitives[name],precision_bits);pulses+=count*len(primitives[name])
    for kind in ('port','seam'):
        word=list(phase_template(plan,kind));bound+=a.template_bound(word,precision_bits);pulses+=len(word)
    return {'source_id':walk.SOURCE_ID,'size':size,'spacing_exact':str(h),'mass_exact':str(mass),
            'precision_bits':precision_bits,'native_pulses':pulses,'load_export_services':2*pulses,
            'whole_word_full_response_operator_error_exact':str(bound),
            'complete_channel_TV_and_half_diamond_upper_exact':str(min(F(1),bound)),
            'memory_bits':plan['memory_bits'],'clean_extra_workspace_bits':plan['clean_extra_workspace_bits'],
            'reset_preparation_clock_detector_included':False,
            'phase_ancestry':'declared external rational seam phases and port potentials; no autonomous field dynamics',
            'native_word':'pulse_stream: actual twelve rational controls per source pulse',
            'serialization':'exact rational adapter; frozen binary64 API unchanged'}


def strict_development_budget(size,steps,h=None,mass=F(0),precision_bits=64):
    """Free massive replacement of the audited all-selector apparatus budget.

    This wrapper admits the frozen two-mode preparation and terminal
    detector only. Gaussian preparation and phase-control workspace have
    separate reset/preparation obligations and are outside this wrapper.
    """
    if type(steps)is not int or steps<0:raise ValueError('nonnegative integer steps required')
    h=F(1,size) if h is None else h
    if h!=F(1,size):raise ValueError('frozen apparatus clock requires h=1/size')
    path=walk.ROOT/'openwave-tuning/next/CHAIN2/particle_field_limit_audit/strict_certificate.py'
    spec=importlib.util.spec_from_file_location('larger_walk_strict_preparation',path)
    strict=importlib.util.module_from_spec(spec);sys.modules[spec.name]=strict;spec.loader.exec_module(strict)
    row=strict.development_budget(size,steps,precision_bits)
    source=certificate(size,h,mass,precision_bits)
    delta=F(source['whole_word_full_response_operator_error_exact'])
    row['base_massless_delta_W_exact']=row['delta_W_exact']
    row['delta_W_exact']=str(delta);row['n_delta_W_exact']=str(steps*delta)
    row['source_compression_TV_exact']=str(min(F(1),steps*delta))
    total=F(row['delta_preparation_exact'])+steps*delta+F(row['clock_and_detector_exact'])+F(row['reset']['timeout_bound_exact'])
    row['native_vs_ideal_complete_law_TV_exact']=str(min(F(1),total))
    row['massive_word_certificate']=source
    row['comparison_target']='same finite free massive step; no physical continuum or species transfer'
    return row
