"""Exact local contact inequality behind the uniform physical-smearing bound."""
import sympy as s


def _annihilator(mode):
    result=s.zeros(4)
    for mask in range(4):
        if mask>>mode&1:
            result[mask^(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
    return result


def contact_certificate():
    a,b,c,d=s.symbols('a b c d',real=True)
    u,v=a+s.I*b,c+s.I*d
    fields=[_annihilator(j) for j in range(2)]
    number=[x.adjoint()*x for x in fields]
    contact=number[0]*number[1]
    composite=u*number[1]*fields[0]+v*number[0]*fields[1]
    field=u*fields[0]+v*fields[1]
    if contact*field-field*contact!=-composite:
        raise ArithmeticError('contact commutator changed')
    q=s.simplify(composite*composite.adjoint()+composite.adjoint()*composite)
    norm=a*a+b*b+c*c+d*d
    if s.simplify(q*q-norm*q)!=s.zeros(4):
        raise ArithmeticError('scaled contact projection identity failed')
    if s.simplify(s.trace(q)-2*norm)!=0:
        raise ArithmeticError('contact projection rank changed')
    return {'Fock_dimension':4,'complex_spinor_parameters':4,
            'exact_identity':'Q={B(u),B(u)^*}>=0, Q^2=||u||^2 Q, tr Q=2||u||^2',
            'onsite_operator_bound':'0<=Q<=||u||^2 I',
            'global_operator_bound':'||sum_x B(u_x)||<=||u||_ell2 from disjoint odd anticommutation',
            'contact_commutator_exact':True,
            'state_or_particle_number_projection_used':False}
