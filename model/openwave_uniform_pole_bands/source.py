"""A three-coordinate finite Schur shortcut for every real momentum."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as sy
from model.openwave_baryon_poles.source import exact_symbol
from model.openwave_baryon_dispersion.source import repair_matrix
from model.openwave_hadron_refinement_flow.analytic import constants,full_bounds
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_hadron_refinement_flow.certified_rest import eigen_intervals,sqrt_upper
from model.openwave_capacity_stress.source import enclosure
LOW=(0,2,4);HIGH=(1,3);CELLS=32


def norm_frobenius(A,left_gram,right_gram):
    square=sum(F(sy.expand(sy.conjugate(A[i,j])*A[i,j]))/F(left_gram[i]*right_gram[j])for i in range(A.rows)for j in range(A.cols))
    return sqrt_upper(square)


@lru_cache(None)
def member(proton=False):
    A,G=exact_symbol(F(1),F(0),proton);J,_,_=repair_matrix(proton)
    grams=list(G.diagonal());g=[grams[i]for i in LOW]
    j=F(J[0,0])/6;small=J.extract(LOW,LOW)-sy.Rational(j)*G.extract(LOW,LOW)
    small_norm=norm_frobenius(small,g,g)
    _,_,_,p,_,unit=constants(*degree_coefficients());bounds=full_bounds(*degree_coefficients())
    C=bounds['Q_floor']-(p['delta']if proton else 0)
    # The two high compact coordinates have common diagonal c and positive
    # correlated addition. Their mutual native entry is -t/2.
    high_compact=F(A[1,1])-p['t']/2
    b=sqrt_upper(5)*p['t'];D=min(C,high_compact)-b
    native_cross=norm_frobenius(A.extract(LOW,HIGH),g,[grams[i]for i in HIGH])
    kinetic_cross=norm_frobenius(J.extract(LOW,HIGH),g,[grams[i]for i in HIGH])
    cross2=(native_cross+4*kinetic_cross)**2+5*p['t']**2
    # After subtracting f(k)*j I, only one light-frame hopping pair depends
    # on k. Its derivative has norm sqrt(3/2)t; f'(k) has magnitude<=2.
    derivative=sqrt_upper(F(3,2))*p['t']+2*small_norm
    return {'A':A,'G':G,'J':J,'scalar':j,'Q_floor':D,'coupling_square':cross2,
            'reduced_derivative':derivative,'unit':unit,'high_compact_floor':high_compact,
            'native_Q_floor':C,'original_high_Q_cross_norm_upper':b,
            'small_inverse_rest_norm_upper':small_norm}


def halfangle(t,sign=1):
    t=F(t)
    if sign not in(-1,1)or not 0<=t<=1:raise ValueError('two half-circle charts require t in[0,1] and sign+-1')
    return sign*(1-t*t)/(1+t*t),sign*2*t/(1+t*t)


def levels3(A,G):
    if any(sy.im(x)!=0 for x in A):
        R,I=A.applyfunc(sy.re),A.applyfunc(sy.im);R=R.row_join(-I).col_join(I.row_join(R));G=sy.diag(G,G)
        levels=eigen_intervals(R,G);assert levels[::2]==levels[1::2];return levels[::2]
    return eigen_intervals(A,G)


def cell(index,sign,proton=False):
    d=member(proton);left,right=F(index,CELLS),F(index+1,CELLS);mid=(left+right)/2
    c,s=halfangle(mid,sign);A,G=exact_symbol(c,s,proton)
    small=A.extract(LOW,LOW)-sy.Rational(2*(1-c)*d['scalar'])*G.extract(LOW,LOW)
    levels=levels3(small,G.extract(LOW,LOW))
    # k(t)=2atan(t) (plus pi) has derivative<=2.
    radius=2*d['reduced_derivative']*(right-left)/2
    ca,_=halfangle(left,sign);cb,_=halfangle(right,sign)
    scalar_lo,scalar_hi=sorted((2*(1-ca)*d['scalar'],2*(1-cb)*d['scalar']))
    compact=[(lo-radius+scalar_lo,hi+radius+scalar_hi)for lo,hi in levels]
    top=compact[-1][1];assert top<d['Q_floor']
    eps=d['coupling_square']/(d['Q_floor']-top)
    # The common scalar shift cancels in the local eigenvalue gap.
    gap=levels[1][0]-levels[0][1]-2*radius-eps
    assert gap>0,(index,sign,proton,float(gap))
    q=d['coupling_square']/(d['Q_floor']-compact[0][1])**2
    assert 0<=q<1
    energies=[((lo-eps)*d['unit'],hi*d['unit'])for lo,hi in compact]
    return {'half_angle_chart_sign':sign,'half_angle_tangent_interval':[str(left),str(right)],
        'rational_center_phase':[str(c),str(s)],'reduced_symbol_radius_source':str(radius),
        'scalar_shift_source_interval':[str(scalar_lo),str(scalar_hi)],
        'all_three_compact_source_intervals':[list(map(str,x))for x in compact],
        'all_three_full_source_GeV_intervals':[list(map(str,x))for x in energies],
        'complete_ground_isolation_gap_lower_source':str(gap),'uniform_Schur_error_source':str(eps),
        'full_high_complement_probability_upper':str(q),
        'three_field_ground_pole_residue_trace_interval':[str(1-q),'1']}


def run():
    rows=[]
    for proton in(False,True):
        records=[cell(i,sign,proton)for sign in(1,-1)for i in range(CELLS)]
        d=member(proton);gap=min(F(r['complete_ground_isolation_gap_lower_source'])for r in records)
        q=max(F(r['full_high_complement_probability_upper'])for r in records)
        energy=(min(F(r['all_three_full_source_GeV_intervals'][0][0])for r in records),max(F(r['all_three_full_source_GeV_intervals'][0][1])for r in records))
        rows.append({'assignment':'proton'if proton else'neutron','all_momentum_GeV_interval':list(map(str,energy)),
            'all_momentum_GeV_display':enclosure(energy),'all_momentum_gap_lower_source':str(gap),
            'all_momentum_high_probability_upper':str(q),'all_momentum_three_field_residue_trace_lower':str(1-q),
            'light_coordinates':list(LOW),'high_compact_coordinates':list(HIGH),
            'complete_high_complement_floor_source':str(d['Q_floor']),
            'complete_light_to_high_coupling_square_upper':str(d['coupling_square']),
            'removed_common_inverse_rest_scalar':str(d['scalar']),
            'remaining_inverse_rest_norm_upper':str(d['small_inverse_rest_norm_upper']),
            'records':records})
    return {'schema':'oph.uniform_native_pole_bands.v1','verdict':'PASS_COMPLETE_MOMENTUM_CIRCLE_ISOLATED_NATIVE_POLES',
        'members':rows,'circle_cover':{'charts':'z=(1+it)/(1-it), -z, and their conjugates, t in[0,1]',
            'interval_cells_per_chart':CELLS,'every_real_momentum_covered':True,'complex_conjugation_preserves_spectrum':True},
        'finite_theorem_shortcut':'Separate the three light compact coordinates from both high compact coordinates and the complete infinite relative tail; remove the shared scalar inverse-rest propagation before interval transport.',
        'full_source_spectrum_replaced_by_three_coordinates':False,
        'mass_or_current_parameter_refitted':False,
        'ancestry':'Inherited electron GeV anchor, strong and axial mass calibration, supplied speed and strong length. The native source is unchanged.',
        'physical_four_dimensional_particle_dictionary_proved':False,
        'scope':'Every fixed-spacing momentum fibre of the complete infinite one-dimensional interacting baryon source has an isolated simple lowest pole and a positive native three-field residue. Normalizable packets require the separate direct-integral construction.'}
