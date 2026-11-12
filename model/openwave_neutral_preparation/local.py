"""Explicit native local eigenlines and finite state-preparation services."""
from functools import lru_cache
import sympy as s
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source,target
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import evaluate,divide,outward
from model.openwave_neutral_stationary.isolation import data


@lru_cache(maxsize=None)
def line(kind):
    if kind not in ('neutral','ground_up','ground_down'):raise ValueError('specified native local eigenline required')
    block=target(2,0) if kind=='neutral' else source(1 if kind=='ground_up' else -1)
    h=block['H'];g=block['gram'];x=s.Symbol('x');e=s.Symbol('e')
    roots=data()['en'] if kind=='neutral' else data()['e0']
    factor=s.Poly((data()['neutral_factor'] if kind=='neutral' else data()['ground_factor']).as_expr().subs(x,e),e)
    q=s.Poly(h.charpoly(x).as_expr(),x).sqf_part()
    coefficients=[s.Poly(sum(q.nth(j)*e**(j-i-1) for j in range(i+1,q.degree()+1)),e).rem(factor)
                  for i in range(q.degree())]
    denominator=s.Poly(q.diff().as_expr().subs(x,e),e).rem(factor)
    powers=[s.eye(h.rows)]
    for _ in range(q.degree()-1):powers.append(h*powers[-1])
    numerator=sum((c.as_expr()*p for c,p in zip(coefficients,powers)),s.zeros(h.rows))
    # Choose a concrete normalized native basis vector with nonzero ground
    # overlap. This choice depends only on the actual local source matrices.
    weights=[]
    for j in range(h.rows):
        weight=s.Poly((g*numerator)[j,j]/g[j,j],e).rem(factor)
        interval=divide(evaluate(weight,roots),evaluate(denominator,roots))
        weights.append((interval,j,weight))
    interval,j,weight=max(weights,key=lambda z:z[0][0]);assert 0<interval[0]<=interval[1]<1
    # P e_j/sqrt(g_jj p_j) is a fully specified algebraic local eigenvector.
    vector=[s.Poly(numerator[i,j],e).rem(factor) for i in range(h.rows)]
    residual=(h-e*s.eye(h.rows))*s.Matrix([v.as_expr() for v in vector])
    assert all(s.Poly(z,e).rem(factor).is_zero for z in residual)
    return {'kind':kind,'native_block_dimension':h.rows,'native_full_local_dimension':5*2**30,
        'selected_energy_factor':str(factor.as_expr()),'selected_root':'smallest real root',
        'seed_basis_column':j,'seed_norm_squared':str(g[j,j]),
        'seed_selected_line_overlap_interval':outward(interval,160),
        'line_numerator_coefficients_in_e':[str(v.as_expr()) for v in vector],
        'projector_denominator':str(denominator.as_expr()),
        'line_normalization':'P(e) e_j / sqrt(Gram_jj * p_j), p_j=<seed|P(e)|seed>',
        'unitary_service':'Let |b> be that normalized seed and |v>=P(e)|b>/sqrt(p_j). U=I-2|b-v><b-v|/[2(1-sqrt(p_j))]; extend by identity on the complete orthogonal complement.',
        'same_block_unitary_preserves_Gauss_raw_L_B3_and_fermion_parity':True,
        'local_line_not_an_unknown_full_interacting_ground':True,
        'primitive_repair_gate_compilation_derived':False}


def certificate():
    rows=[line(kind) for kind in ('neutral','ground_up','ground_down')]
    return {'local_lines':rows,
        'initial_product_service':'Choose one neutral cell uniformly from eight; independently choose each of the seven actual ground spin lines with probability1/2; prepare the explicit local lines and all24 blank links.',
        'source_input_density':'rho0=P/1024=(1/8) sum_y P_y/128',
        'uniform_classical_preparation_bits':10,
        'local_state_services_per_preparation':8,'blank_link_registers':24,
        'classical_recipe_records_retained':1024,
        'local_input_resources':'Eight specified gauge-singlet native seed states in their stated raw-number/parity sectors and blank links are supplied; the explicit local unitaries do not prepare these reservoirs from an arbitrary unknown state.',
        'parity_reference':'The odd source is accompanied by a retained gauge-neutral occupied fermionic reference, giving even total parity. It is idle during the source clock.',
        'external_charge_reference_required':False,
        'initial_source_density_stationary_under_full_H':False,
        'all_preparation_memory_and_external_references_retained':True,
        'reduced_source_target_does_not_assert_decoupling_from_preparation_archive':True}
