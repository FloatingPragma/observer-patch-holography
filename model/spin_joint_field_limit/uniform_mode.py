"""Actual uniform Ward concentration and finite-box order certificates."""
from fractions import Fraction as F
from model.spin_fixed_field_locality import bounds as locality
from model.spin_charge_dispersion import order
from . import bounds


def concentration(n, order_floor=F(1,20), sites=None):
    row=bounds.family(n)
    m=bounds.rational(order_floor,'order floor',True)
    if m>F(1,2):raise ValueError('spin-half order exceeds its source norm')
    ratio=6/(m*n)+F(1,n**6)
    out={'n':n,'p':row['p'],'field':row['field'],'order_floor':m,
         'scaled_inverse_weighted_second_moment_upper':ratio,
         'scaled_third_inverse_weighted_energy_square_mean_upper':ratio,
         'scaled_third_inverse_weighted_energy_square_variance_upper':ratio**2/4,
         'selected_scaled_energy_window_upper_squared':2*ratio,
         'weighted_centre_converges_to_zero':True,
         'nonzero_momentum_or_particle_shell_proved':False,
         'original_continuum_residue_lower_proved':False}
    if sites is not None:
        if type(sites)is not int or sites<1:raise ValueError('positive source site count required')
        out['source_sites']=sites
        out['original_preparation_and_window_probability_lower']=m/(8*sites)
    return out


def order_box(n):
    row=bounds.family(n);d=1
    if order.CERTIFIED_M0<F(1,10):raise ArithmeticError('owning source order floor weakened')
    while True:
        radius=384*d;side=2*radius+1
        marginal=locality.marginal(1,row['field'],d,1,6*side**2)
        if marginal['half_trace_upper']<=F(1,40):break
        d*=2
    omitted=1-F((side-2*d)**3,side**3)
    if omitted>F(1,64):raise ArithmeticError('finite order boundary fraction too large')
    floor=(1-omitted)*F(3,40)-omitted/2
    if floor<=F(1,20):raise ArithmeticError('source order floor was lost')
    return {'n':n,'source_radius':radius,'source_sites':side**3,
            'boundary_distance':d,'bulk_one_site_marginal':marginal,
            'excluded_boundary_fraction':omitted,'finite_order_floor':floor,
            'concentration':concentration(n,F(1,20),side**3),
            'ground_matrix_allocated':False,'large_native_tape_executed':False,
            'positive_field_bulk_order_parent':'spin_charge_dispersion/order.py; spin_joint_field_limit/PROOF.md'}


def fixed_box_gap(n, zero_field_gap):
    row=bounds.family(n);g=bounds.rational(zero_field_gap,'zero field gap',True)
    return max(F(),n*n*g-F(8,n**3))


def source_identity(n):
    """Exact all-row Ward and first-moment identities on the actual cube."""
    import numpy as np
    from model.native_ordered_ground import source
    row=bounds.family(n);h,scale=source.matrix(field=row['field'])
    size=len(h);indices=np.arange(size);identity=np.eye(size,dtype=np.int64)
    signs=tuple(-1 if x in source.cube()[1]else 1 for x in range(8))
    m=np.array([sum(signs[x]*(1-2*((word>>x)&1))for x in range(8))for word in indices])
    t=np.zeros_like(h);o=np.zeros_like(h);xy=np.zeros_like(h)
    for x,sign in enumerate(signs):
        t[indices^(1<<x),indices]=-sign*np.array([1-2*((word>>x)&1)for word in indices])
        o[indices^(1<<x),indices]+=1
    for x,y in source.cube()[0]:
        for word in indices:
            if ((word>>x)^(word>>y))&1:xy[word^(1<<x)^(1<<y),word]+=2
    # A=iT/2, M=diag(m/2), O=o/2 in the actual Marshall frame.
    if not np.array_equal(h@t-t@h,int(scale*row['field'])*np.diag(m)):
        raise ArithmeticError('source uniform Ward sign or coefficient changed')
    doubled=-(m[:,None]-m[None,:])**2*h
    expected=4*scale*xy+int(2*scale*row['field'])*o
    if not np.array_equal(doubled,expected):
        raise ArithmeticError('complete source first-moment identity changed')
    return {'n':n,'all_source_rows':size,'uniform_Ward_operator_identity_exact':True,
            'double_commutator_first_moment_identity_exact':True,
            'source_coupling':'1','source_field':str(row['field'])}


def certificate():
    source_order=order.certificate()
    if F(source_order['selected_magnetization_strict_lower'])<F(1,10):
        raise ArithmeticError('owning source order certificate is insufficient')
    return {'owning_order_certificate_replayed':source_order,
            'actual_source_operator_identities':[source_identity(n)for n in(2,3,4)],
            'finite_source_order_and_concentration':[order_box(n)for n in(2,4,16)],
            'uniform_mode_centre_and_weighted_variance_vanish_on_joint_clock':True,
            'fixed_compact_test_field_mode_identified':False,
            'uniform_original_probability_floor_proved':False,
            'nonzero_momentum_or_physical_mass_shell_proved':False}
