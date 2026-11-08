"""Actual centered-ground nodal response and dissipative-memory lower bounds."""
from fractions import Fraction as F
from math import isqrt
import sympy as s
from .local_matter import annihilator
from . import momentum_response as momentum
from . import uv_tail


def rational(value,positive=False):
    if type(value)not in(int,F):raise ValueError('exact rational source input required')
    value=F(value)
    if value<0 or(positive and not value):raise ValueError('positive source input required')
    return value


def sqrt_lower(value,bits=40):
    value=rational(value)
    if type(bits)is not int or not 1<=bits<=128:raise ValueError('integer square-root grid required')
    scale=1<<bits
    return F(isqrt(value.numerator*scale*scale//value.denominator),scale)


def star_algebra():
    """Full four-mode source, including both singly occupied and eta blocks."""
    a=tuple(annihilator(2,j)for j in range(4));identity=s.eye(16)
    def generators(site,epsilon):
        u,d=a[2*site:2*site+2];nu=u.H*u;nd=d.H*d
        return ((u.H*d,(nu-nd)/2),(epsilon*u.H*d.H,(nu+nd-identity)/2))
    left,right=generators(0,1),generators(1,-1)
    q=sum((z*w+(p*v.H+p.H*v)/2 for(p,z),(v,w)in zip(left,right,strict=True)),s.zeros(16))
    t=sum((a[j].H*a[j+2]+a[j+2].H*a[j]for j in(0,1)),s.zeros(16))
    if q*t!=t*q:raise ArithmeticError('actual star invariant does not commute with hopping')
    expected={s.Integer(0):10,s.Integer(2):5,s.Integer(6):1}
    for sign in(-1,1):
        if (identity-4*q+sign*t).eigenvals()!=expected:
            raise ArithmeticError('full star source operator inequality failed')
    # On the one-fermion sector Q=0 and |T|=1: coefficient four is sharp.
    vector=s.zeros(16,1);vector[1]=1;vector[4]=1
    if q*vector!=s.zeros(16,1)or t*vector!=vector:
        raise ArithmeticError('source sharpness witness changed')
    y=momentum.bond_derivative(0,1)
    b=(2*a[1].H*a[1]-identity)*a[0]
    force=b*t-t*b
    if (force.H*force).eigenvals()!={s.Integer(0):8,s.Integer(1):4,s.Integer(5):4}:
        raise ArithmeticError('source local derivative norm changed')
    return {'full_Fock_dimension':16,'both_signs_PSD_spectra':{'0':10,'2':5,'6':1},
            'star_comparison_coefficient':F(4),'coefficient_sharp':True,
            'D_operator_norm_squared':F(5),'all_occupation_sectors_retained':True,
            'source_D_polynomial_terms':len(y)}


def dimer_energy(c,u):
    c,u=rational(c,True),rational(u,True)
    delta=(sqrt_lower(u*u+4*c*c)-u)/4
    if delta<=0:raise ValueError('square-root resolution does not certify positive dimer gain')
    return delta


def source_budget(c,u,*,s_max=F(0),eta=None):
    """Positive source residual bound on the specified invariant ground average."""
    c,u=rational(c,True),rational(u,True);ss=rational(s_max)
    if ss>3:raise ValueError('kinetic S lies in [0,3]')
    delta=dimer_energy(c,u)
    node=2*sqrt_lower(F(6))*c*delta
    lower=node-12*c*c*ss
    current_bound=322*c*(u/2+12*c)+23*c*c
    row={'c':c,'U':u,'dimer_gain_lower':delta,'S_upper':ss,
         'nodal_residual_Gram_lower':node,'neighbourhood_residual_Gram_lower':max(F(0),lower),
         'Lperp_R_norm_upper':current_bound,'source_state':'proper_cubic_average_of_invariant_centered_full_ground',
         'source_clock_unchanged':True,'zero_frequency_atom_proved':False,
         'physical_mass_assigned':None}
    if eta is not None:
        eta=rational(eta,True)
        if lower<=0:raise ValueError('positive residual required for the frequency lower bound')
        row['imaginary_frequency']=eta
        row['negative_imaginary_memory_lower']=eta*lower*lower/(eta*eta*lower+current_bound*current_bound)
    return row


def default_budget(c=F(1,8),eta=F(1,8)):
    c=rational(c,True);eta=rational(eta,True)
    node=F(36,125)*c*c;lower=F(18,125)*c*c;ss=F(3,250)
    audit=source_budget(c,8*c,s_max=ss,eta=eta)
    if audit['nodal_residual_Gram_lower']<=node or audit['neighbourhood_residual_Gram_lower']<=lower:
        raise ArithmeticError('actual dimer energy does not imply the fixed nodal enclosure')
    current=5175*c*c
    if audit['Lperp_R_norm_upper']!=current:raise ArithmeticError('default source current scale changed')
    return {'c':c,'U':8*c,'S_upper':ss,'node_Gram_strict_lower':node,
            'neighbourhood_Gram_lower':lower,'original_fourth_moment_defect_lower':16*c*c*lower,
            'Lperp_R_norm_upper':current,'imaginary_frequency':eta,
            'negative_imaginary_memory_lower':eta*lower*lower/(eta*eta*lower+current*current),
            'every_original_coefficient_packet_in_region_covered':True,
            'zero_frequency_residue_inferred':False,'continuum_particle_mass':None}


def certificate():
    if uv_tail.ball(4)!=129:raise ArithmeticError('source locality ball differs')
    star=star_algebra()
    if sqrt_lower(F(17))<=F(103,25)or sqrt_lower(F(6))<=F(12,5):
        raise ArithmeticError('default rational square-root admission failed')
    observables=momentum.invariant_observables()
    return {'same_centered_full_CAR_source':True,'full_star_operator_audit':star,
            'local_qstar_score_plan':momentum.score_plan(observables['qstar']),
            'default_source_result':default_budget(),
            'local_derivative_support_radius':2,'Gram_range_ball_sites':129,
            'exact_positive_memory_for_nonzero_imaginary_frequency':True,
            'ground_marginal_or_native_services_implicitly_supplied':False,
            'source_ground_replaced_by_dimer_trial':False,
            'physical_species_identified':[]}
