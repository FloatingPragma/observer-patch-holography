"""Exact complete local source preparation and finite-polynomial countercontrol."""
from fractions import Fraction as F
import numpy as np
from model.protected_defect_catalog import charge
from model.centered_particle_spectrum.certificate import decimal_bounds
from model.spin_operational_position.bounds import generic_instrument_plan


def build():
    op=charge.operators();g,error,certificate=charge.ground_input()
    identity=np.eye(256,dtype=np.int64)
    matrices={q:np.asarray(a,dtype=np.int64)for q,a in op['twice_creation'].items()}
    # f(Z/2)=(2/5)Z; its exact source charge projections are A_q*2/5.
    # Since matrices[q]=2A_q, the creator is matrices[q]/5.
    if not np.array_equal(sum(a.T@a for a in matrices.values()),4*identity):
        raise ArithmeticError('complete local charged effects differ from4I')
    rows=[];norm=int(g@g)
    for q,a in sorted(matrices.items()):
        if np.any(a@a):raise ArithmeticError('single-site finite creator is not nilpotent')
        if not np.array_equal(op['twice_charge']@a-a@op['twice_charge'],2*q*a):
            raise ArithmeticError('source local charge differs')
        image=a.astype(object)@g
        probability=F(int(image@image),25*norm)
        uncertainty=F(8,25)*error
        rows.append({'internal_charge':q,
                     'original_probability_interval':decimal_bounds(max(F(),probability-uncertainty),min(F(4,25),probability+uncertainty),12),
                     'creator_square_exactly_zero':True,
                     'operator_norm_upper':'2/5'})
    return {'source_ground_certificate':certificate,
            'source_Hamiltonian_and_clock':{'integer_scale':op['scale'],'source_matrix_dimension':256},
            'local_test':'x=Z_0/2; f(x)=x/(1+x^2)',
            'exact_creator':'B_q=(2/5) A_q=twice_creation[q]/5',
            'records':rows,'complete_complement_effect':'21/25 I',
            'source_charge_and_local_support_checked':True,
            'all_source_spin_rows_retained':256,
            'complete_effect_sum_identity':True,
            'finite_source_claimed_to_have_continuum_polynomial_Gram_independence':False,
            'actual_native_instrument_plan':generic_instrument_plan(8,3,F(1,10**10)),
            'full_eight_site_native_tape_executed':False,
            'finite_ground_substituted_for_ordered_continuum':False}
