"""Actual admitted programs with equal capacity and different coherent action."""
from collections import defaultdict
from fractions import Fraction as F

from model.native_program import interpreter as native
from model.observer import capacity as record_capacity


def selected_program(c_weight):
    c_weight = F(c_weight)
    if not 0 < c_weight < 1:
        raise ValueError('both source laws must have positive original weight')
    program = native.program(biased=False, capacity=1, commands=[{
        'op':'repair', 'selector':{'kind':'fixed','seam':[0,1]},
        'law':'mixture', 'repeat':1}])
    program['parameters']['c_weight'] = str(c_weight)
    program['declared_resources'] = native._native_cost(program)
    return native.admit(program)


def summed(history):
    out = defaultdict(F)
    for branch in history.branches:
        for key, value in branch.density.items():
            out[key] += value
    return {key:value for key,value in out.items() if value}


def certificate():
    basis = native.local.basis(True)
    x = basis.index(tuple(int(k in (0,4)) for k in range(12)))
    y = basis.index(tuple(int(k in (2,4)) for k in range(12)))
    input_state = {(a,r,b,s):F(1,2) for a,r in ((x,0),(y,1))
                   for b,s in ((x,0),(y,1))}
    _, capacity = record_capacity.dependencies()
    programs = [selected_program(c) for c in (F(1,4),F(3,4))]
    results = []
    for program in programs:
        run = native.execute(program,input_state)
        if run.status != 'COMPLETE':
            raise ArithmeticError('the admitted literal program did not finish')
        output = summed(run)
        expected = F(program['parameters']['c_weight'])/4
        if output.get((x,0,y,1),F()) != expected:
            raise ArithmeticError('actual coherent source response changed')
        if native.reference_marginal([output]) != {(0,0):F(1,2),(1,1):F(1,2)}:
            raise ArithmeticError('the complete program changed the untouched reference')
        uniform = {(i,0,i,0):F(1,len(basis)) for i in range(len(basis))}
        census = native.execute(program,uniform)
        records = [str(branch.record) for branch in census.branches
                   if native.trace(branch.density)>0]
        if len(records) != len(set(records)):
            raise ArithmeticError('record atoms are not distinct')
        identity = {'rows':{r:{r:1} for r in records}}
        graph = capacity.compound_confusability_graph(records,[identity])
        code = capacity.maximum_independent_set(graph)
        if len(code) != 10:
            raise ArithmeticError('the actual reachable full-record capacity changed')
        results.append({'c_weight':program['parameters']['c_weight'],
            'program_digest':native.digest(program),
            'actual_allocated_complete_output_blocks':len(census.branches),
            'reachable_full_records':len(records),'correctable_record_count':len(code),
            'terminal_capacity_nats':'ln(10)',
            'source_reference_off_diagonal':str(expected),
            'complete_source_reference_trace':str(native.trace(output)),
            'resources':program['declared_resources']})
    if programs[0]['declared_resources'] != programs[1]['declared_resources']:
        raise ArithmeticError('the source control changed the resource budget')
    defaults = native.program()
    return {'same_carrier':'H66','same_memory_slots':1,
        'same_source_preparation_and_untouched_reference':True,
        'same_complete_capacity_and_service_counts':True,
        'programs':results,
        'literal_default_c_weight':defaults['parameters']['c_weight'],
        'literal_default_repair_probability':defaults['parameters']['repair_probability'],
        'capacity_alone_selects_program':False,
        'fixed_program_predictions_are_undetermined':False,
        'spin_Hamiltonian_or_capacity_trajectory_recovered':False}
