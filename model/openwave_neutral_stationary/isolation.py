"""Full-source neutral band isolation from complete local and flux bounds."""
from functools import lru_cache
import sympy as s
from model.coherent_defect_wave.chiral_joint.local_reference.spectrum import certificate as local_spectrum
from model.coherent_defect_wave.chiral_joint.atomic_reference_control.lattice import hopping_norm,local_norm
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import outward

LAMBDA=s.Rational(1,2**34)
KAPPA=s.Integer(2)
DIMENSION=(5*2**30)**8*119**24


@lru_cache(maxsize=4)
def data(kappa=KAPPA):
    kappa=s.Rational(kappa)
    if kappa<=0:raise ValueError('positive declared link penalty required')
    x=s.Symbol('x');u=x*(x-1)
    ground=s.Poly(u**3-313*u**2+31708*u-1041468,x)
    neutral=s.Poly(u**3-313*u**2+31738*u-1041088,x)
    roots=lambda p:tuple(map(s.Rational,outward(tuple(map(s.Rational,
        p.intervals(eps=s.Rational(1,10**110))[0][0])),384)))
    e0=roots(ground);en=roots(neutral);local=local_spectrum()
    for exact,name in ((e0,'global_physical_ground_energy_interval'),(en,'first_excited_energy_interval')):
        lo,hi=map(s.Rational,local[name]);assert lo<=exact[0]<exact[1]<=hi
    assert local['global_physical_ground_rank']==local['first_excited_rank']==2
    assert {z['raw_L'] for z in local['complete_maximal_module_energy_records'][1]['sectors']}=={0,2}
    low=s.Rational(-1171,100);loss=e0[1]-low
    gap=(en[0]-e0[1],en[1]-e0[0])
    next_local=min(s.Rational(local['complete_maximal_module_energy_records'][2]['energy_interval'][0]),
                   s.Rational(local['excluded_physical_sectors_strict_energy_lower_bound']))
    next_gap=next_local-e0[1]
    flux_gap=kappa-2*loss
    separation=min(next_gap,2*gap[0],flux_gap)-gap[1]
    if separation<=0:raise ValueError('the complete support bound does not isolate the neutral band at this penalty')
    tau=hopping_norm();assert tau==360
    v=LAMBDA*tau;assert v<separation/5
    width=v*v/(separation-v)
    ground_shift=v*v/(gap[0]-v)
    eta=v/separation
    energy=(7*e0[0]+en[0],7*e0[1]+en[1])
    band=(energy[0]-width,energy[1])
    outside=energy[0]+separation-v
    rest=(gap[0]-width,gap[1]+ground_shift)
    assert rest[0]>0 and outside>band[1]
    norm=8*local_norm()+24*kappa+v
    assert norm==456+24*kappa+v
    return {'ground_factor':ground,'neutral_factor':neutral,'e0':e0,'en':en,'gap':gap,
        'local_other_lower':next_local,'unrestricted_local_lower':low,'flux_gap':flux_gap,
        'old_flux_gap':1-2*loss,'separation':separation,'v':v,'tau':tau,'width':width,
        'ground_shift':ground_shift,'eta':eta,'atomic_band_energy':energy,'band':band,
        'outside_lower':outside,'rest':rest,'norm':norm}


def pair(values,bits=160):return outward(tuple(values),bits)


def certificate():
    d=data()
    return {'source':'N2 full PW119, all local Yukawa/scalar coefficients one, periodic fermions, link penalty two, lambda=2^-34',
        'original_link_penalty_one_source_identified':False,
        'original_penalty_one_bound_overlap':{'one_link_lower_gap':str(d['old_flux_gap']),
            'two_link_lower_gap':str(2*d['old_flux_gap']),
            'neutral_atomic_gap_interval':pair(d['gap']),
            'these_coarse_bounds_prove_isolation_at_penalty_one':False,
            'actual_penalty_one_spectral_collision_proved':False},
        'complete_local_ground_factor':str(d['ground_factor'].as_expr()),
        'complete_local_neutral_factor':str(d['neutral_factor'].as_expr()),
        'local_ground_energy_interval':pair(d['e0']),
        'local_neutral_energy_interval':pair(d['en']),
        'other_local_physical_energy_lower':str(d['local_other_lower']),
        'all_nonblank_support_gap_lower':str(d['flux_gap']),
        'complete_atomic_band_separation_lower':str(d['separation']),
        'atomic_band':'one raw-L0 spin-zero local eigenstate, seven actual local ground doublets, every link blank',
        'raw_L_B3':[7,0],'band_rank':1024,'all_eight_momentum_sector_ranks':[128]*8,
        'atomic_band_full_hopping_compression':'P_nu T P_nu=0 by every actual blank-link matrix element',
        'full_hopping_norm_upper':str(d['tau']),'hopping_perturbation_norm_upper':str(d['v']),
        'full_source_band_energy_interval':pair(d['band']),
        'every_other_same_sector_energy_lower':str(d['outside_lower']),
        'band_width_upper':str(d['width']),
        'equal_rank_projector_distance_upper':str(d['eta']),
        'actual_ground_energy_interval':pair((8*d['e0'][0]-d['ground_shift'],8*d['e0'][1])),
        'actual_ground_relative_neutral_band_interval':pair(d['rest']),
        'ground_relative_decimal_interval':[float(z) for z in d['rest']],
        'band_width_decimal_upper':float(d['width']),
        'original_unique_ground_or_preparation_transferred':False,
        'lowest_band_projector_is_explicit_known_interval_not_unknown_minimum':True,
        'full_physical_and_nonblank_complements_retained':True,
        'full_source_dimension':str(DIMENSION),'full_source_operator_norm_upper':str(d['norm'])}
