"""Exact local selection rules for the source spin-creation mechanism."""
import sympy as s

from model.spin_canonical_field.algebra import spin, comm


def report():
    x,y,z=spin(); identity=s.eye(2)
    sites=[s.kronecker_product(2*x,identity),s.kronecker_product(identity,2*x)]
    total=(sites[0]-sites[1])/2
    local={}
    for site,epsilon in ((0,1),(1,-1)):
        for q in (-1,1):
            a=z-q*s.I*epsilon*y
            if a*a != s.zeros(2):
                raise ArithmeticError('local spin creator lost its occupancy restriction')
            local[site,q]=(s.kronecker_product(a,identity) if site==0
                           else s.kronecker_product(identity,a))
            if comm(total,local[site,q]) != q*local[site,q]:
                raise ArithmeticError('local source charge selection rule failed')
    for q in (-1,1):
        for r in (-1,1):
            a,b=local[0,q],local[1,r]
            if comm(a,b) != s.zeros(4) or comm(total,a*b) != (q+r)*a*b:
                raise ArithmeticError('disjoint creation or additive charge failed')
            if a*b == s.zeros(4):
                raise ArithmeticError('disjoint creation was removed by a false exclusion rule')
    return {'single_site_repeated_creation_is_zero':True,
            'disjoint_site_creators_commute':True,
            'two_site_created_charges':[-2,0,2],
            'composite_charge_rule':'[Q,A_a A_b]=(a+b) A_a A_b',
            'neutral_dynamics_selection_rule':'P_qprime U(t) P_q=0 for qprime != q',
            'nonzero_charge_can_turn_into_only_neutral_products':False,
            'charge_conservation_proves_one_particle_stability':False,
            'internal_spin_half_implies_Lorentz_spin_half':False,
            'bare_repair_records_distinguish_protected_spin_states':False,
            'charge_readout_requires_an_admitted_effect_sensitive_to_charge':True}
