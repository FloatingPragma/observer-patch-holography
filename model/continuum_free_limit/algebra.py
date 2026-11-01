"""Exact finite controls for the cubic channel and the hierarchy argument."""
import sympy as s

from ..centered_source_ground import dimer


def cubic_control():
    value=dimer.data(t=s.Integer(1),U=s.Integer(0))
    c=dimer.annihilator(0);other=dimer.annihilator(1)
    d=(2*other.T*other-s.eye(16))*c
    G=value['ground'];L=value['K']-value['ground_energy']*s.eye(16)
    norm=value['ground_norm_squared']
    for field,energy in ((c,1),(d,3)):
        for op in (field,field.T):
            vector=op*G
            if s.simplify(L*vector-energy*vector)!=s.zeros(16,1):
                raise ArithmeticError('free cubic and linear energies were conflated')
            if s.simplify((vector.T*vector)[0]/norm)!=s.Rational(1,2):
                raise ArithmeticError('a charge branch was discarded')
    if s.simplify((G.T*c.T*d*G)[0])!=0:
        raise ArithmeticError('the normal-ordered cubic retains a linear component')
    return {'full_Fock_dimension':16,'reference':'free centered t=1 dimer full ground',
            'linear_addition_and_removal_energy':'1',
            'cubic_addition_and_removal_energy':'3',
            'each_charge_weight':'1/2','linear_cubic_overlap':'0',
            'cubic_channel_discarded':False}


def vacuum_only_countercontrol():
    # Empty finite Fock vacuum: an interaction is invisible in all one-field
    # vacuum propagators but changes the doubly occupied level.
    c=dimer.annihilator(0)[:4,:4];other=dimer.annihilator(1)[:4,:4]
    n0=c.T*c;n1=other.T*other;H=n0+n1+n0*n1
    vacuum=s.Matrix([1,0,0,0]);pair=s.Matrix([0,0,0,1])
    defect=H*c-c*H+c
    if defect*vacuum!=s.zeros(4,1) or defect.T*vacuum!=s.zeros(4,1):
        raise ArithmeticError('vacuum-only negative control failed')
    if defect*pair==s.zeros(4,1) or H[3,3]!=3:
        raise ArithmeticError('nonlinear two-particle interaction disappeared')
    return {'empty_vacuum_one_field_motion':'free with energy 1',
            'two_particle_energy':'3 rather than 2',
            'vacuum_motion_alone_proves_full_free_dynamics':False,
            'outside_cyclicity_required_in_hierarchy_argument':True}


def certificate():
    return {'cubic_control':cubic_control(),
            'vacuum_only_countercontrol':vacuum_only_countercontrol()}
