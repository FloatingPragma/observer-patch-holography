"""Exact dynamic restriction of admitted independent native-cell products."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import replace
from fractions import Fraction as F

from ..native_program import interpreter as source


def _integer(value, minimum=1):
    if type(value) is not int or value < minimum:
        raise ValueError('positive finite integer required')
    return value


def _carrier(shell, cells):
    return source.carrier({'kind':'independent_native_product','shell':shell,'cells':cells})


def lift_program(value, cells):
    """Append unused native factors, preserving every operation and parameter."""
    result = source.admit(value)
    current = result['carrier']
    if current['kind'] not in ('local','independent_native_product'):
        raise ValueError('the embedding is only for local or independent-product carriers')
    _integer(cells)
    if cells < current.get('cells',1):
        raise ValueError('appending factors cannot remove an addressed native cell')
    result['carrier'] = _carrier(current['shell'],cells)
    if current['kind'] == 'local':
        for command in result['commands']:
            if command['op'] == 'repair':
                command['cell'] = 0
    result['declared_resources'] = source._native_cost(result)
    return source.admit(result)


def restrict_operator(density, shell, cells, retained):
    """Actual factor partial trace, keeping low-address native factors and R."""
    carrier = _carrier(shell,cells); _integer(retained)
    if retained > cells:
        raise ValueError('retained factor count exceeds the source carrier')
    value = source._validate_density(density,carrier)
    dimension = int(shell[1:])**retained
    result = defaultdict(F)
    for (i,a,j,b), coefficient in value.items():
        if i//dimension == j//dimension:
            result[i % dimension,a,j % dimension,b] += coefficient
    return {key:value for key,value in result.items() if value}


def restrict_history(history, shell, cells, retained, *, target_program):
    """Trace appended factors and explicitly relabel the known program digest.

    Complete chronological records, FIFO/export placement, stops and runtime
    counts are preserved. Declared allocation accounts and specification
    queries are separate and need not agree.
    """
    if not isinstance(history,source.History):
        raise ValueError('a complete native interpreter history is required')
    target = source.admit(target_program)
    carrier = target['carrier']
    if (carrier['kind'] not in ('local','independent_native_product')
            or carrier['shell'] != shell or carrier.get('cells',1) != retained):
        raise ValueError('target identifier must name the retained carrier')
    if history.program_digest != source.digest(lift_program(target,cells)):
        raise ValueError('history identifier does not match the exact lifted target program')
    return replace(history,program_digest=source.digest(target),branches=tuple(
        replace(branch,density=restrict_operator(branch.density,shell,cells,retained))
        for branch in history.branches))


def _state(shell, cells, retained=1):
    """Pure state entangling active factors, a highest unused factor and R."""
    if cells <= retained:
        raise ValueError('this correlated fixture needs an appended factor')
    states = source.local.basis(shell == 'H66'); dimension = len(states)
    configurations = ((0,4),(1,4)) if shell == 'H66' else ((0,0),(0,5))
    first,second = (states.index(tuple(pair.count(i) for i in range(12))) for pair in configurations)
    active = (first,second+sum(2*dimension**j for j in range(1,retained)))
    high = dimension**(cells-retained-1)
    vectors = ((active[0],0,0,1),(active[1],0,1,1),
               (active[0],high,1,1j),(active[1],high,0,-1j))
    real,imaginary = {},{}
    scale = dimension**retained
    for a,s,r,phase in vectors:
        for b,t,u,other in vectors:
            key = (a+scale*s,r,b+scale*t,u)
            coefficient = phase*complex(other).conjugate()/4
            if coefficient.real:
                real[key] = F(coefficient.real)
            if coefficient.imag:
                imaginary[key] = F(coefficient.imag)
    return real,imaginary


def _program(shell='H66', overflow='export_oldest'):
    return source.program(shell=shell,overflow=overflow,commands=[
        {'op':'repair','selector':{'kind':'fixed','seam':[0,5]},'law':'mixture','repeat':1},
        {'op':'repair','selector':{'kind':'last_current','positive':[1,7],
                                  'otherwise':[2,3]},'law':'E','repeat':1},
        {'op':'hold','repeat':2}])


def compare(value, cells, density, **limits):
    base = source.admit(value); extended = lift_program(base,cells)
    retained = base['carrier'].get('cells',1); shell = base['carrier']['shell']
    small_input = restrict_operator(density,shell,cells,retained)
    small = source.execute(base,small_input,**limits)
    large = source.execute(extended,density,**limits)
    reduced = restrict_history(large,shell,cells,retained,target_program=base)
    if reduced != small:
        raise ArithmeticError('native product restriction does not intertwine complete execution')
    return small,large


def fair_direct_countercontrol():
    """The same one-attempt coin record prevents a direct fair-source lift."""
    command = {'op':'repair','selector':{'kind':'fixed','seam':[0,5]},'law':'C','repeat':1}
    local = source.program(commands=[command])
    fair = source.program(federation={'level':0,'total':2},commands=[command])
    index = source.local.basis().index(tuple(int(i in (0,4)) for i in range(12)))
    local_history = source.execute(local,{(index,0,index,0):F(1)})
    mask = (1 << 0)+(1 << 4)
    fair_history = source.execute(fair,{(mask,0,mask,0):F(1)})
    def coin_one(history):
        return sum(source.trace(branch.density) for branch in history.branches
                   if branch.record.kind == 'repair' and branch.record.coin == 1)
    probabilities = coin_one(local_history),coin_one(fair_history)
    if probabilities != (F(25,33),F(1,2)):
        raise ArithmeticError('the natural one-attempt source-family distinction changed')
    return {'fixed_seam':[0,5],'weighted_coin_one':str(probabilities[0]),
            'fair_coin_one':str(probabilities[1]),
            'complete_output_TV_lower_bound':str(abs(probabilities[0]-probabilities[1])),
            'domain':'one actual repair, same complete coin record; not arbitrary encoded multi-attempt simulation'}


def certificate():
    rows = []
    for shell in ('H66','H78'):
        for cells in (2,3,64):
            value = _program(shell)
            for part,density in zip(('real','imaginary'),_state(shell,cells),strict=True):
                small,large = compare(value,cells,density)
                rows.append({'shell':shell,'cells':cells,'operator_part':part,
                    'complete_blocks':len(small.branches),'attempts':small.attempted_steps,
                    'status':small.status,'known_identifier_changes':small.program_digest != large.program_digest})
    stop_rows = []
    value = _program(); real,_ = _state('H66',2)
    for label,base,limits in (
        ('intrinsic_FIFO_timeout',_program(overflow='halt'),{}),
        ('attempt_limit',value,{'max_steps':1}),
        ('branch_limit',value,{'max_branches':1})):
        small,_ = compare(base,2,real,**limits)
        stop_rows.append({'case':label,'status':small.status,'attempts':small.attempted_steps,
                          'complete_blocks':len(small.branches)})
    huge,_ = _state('H66',65)
    extended = source.execute(lift_program(value,65),huge)
    small = source.execute(value,restrict_operator(huge,'H66',65,1))
    if (extended.status != 'RESOURCE_STOP_UNEVALUATED_CARRIER' or small.status != 'COMPLETE'):
        raise ArithmeticError('the product evaluator boundary changed')
    return {'schema':'oph.native_capacity_restriction.product.v1',
        'complete_correlated_operator_comparisons':rows,
        'matched_stop_comparisons':stop_rows,
        'runtime_domain_countercontrol':{'appended_cells':65,'large_status':extended.status,
                                        'local_status':small.status},
        'direct_fair_federation_countercontrol':fair_direct_countercontrol(),
        'theorem':{'full_quantum_and_arbitrary_reference_restriction':'exact',
                   'record_adaptive_selectors_included':True,
                   'all_raw_slots_FIFO_exports_and_stops_retained':True,
                   'independent_spectator_preparation_required':False,
                   'same_declared_allocation_or_specification_reads':False,
                   'compiled_coherent_source_word_services_in_embedding':False,
                   'entire_record_clock_member_lifted':False,
                   'connected_fair_federation_embedding_established':False,
                   'physical_capacity_space_clock_or_mass_selected':False}}
