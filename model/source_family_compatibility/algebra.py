"""Exact finite CAR source, grading and conjugate energy control."""
from functools import lru_cache
import copy
import sympy as s

from model.strong_coupling_refinement.algebra import operators
from model.standard_model_source_map.algebra import source_components,gauge_and_higgs
from model.gauge_defect_transport.transport import neutralized_control


def creation(mode):
    if type(mode)is not int or mode not in range(4):raise ValueError('one of four actual dimer modes required')
    out=s.zeros(16)
    for mask in range(16):
        if not mask>>mode&1:
            sign=(-1)**(mask&((1<<mode)-1)).bit_count()
            out[mask|(1<<mode),mask]=sign
    return out


@lru_cache(maxsize=1)
def _dimer():
    c=s.Rational(1,8);u=s.Integer(1);h=s.zeros(16)
    for number in range(5):
        masks,t2,d=operators(2,((0,1,1),),number)
        block=c*s.Matrix(t2.toarray())/2+u*s.diag(*map(int,d))-u*number*s.eye(len(masks))/2
        for a,x in enumerate(masks):
            for b,y in enumerate(masks):h[x,y]=block[a,b]
    n=s.diag(*[mask.bit_count()for mask in range(16)])
    q=n-2*s.eye(16);parity=s.diag(*[(-1)**mask.bit_count()for mask in range(16)])
    protected=s.diag(*[int((mask&3).bit_count()==1 and (mask>>2).bit_count()==1)for mask in range(16)])
    comm=lambda a,b:a*b-b*a
    assert comm(h,q)==s.zeros(16)and comm(h,parity)==s.zeros(16)
    checks=0
    for j in range(4):
        a=creation(j)
        assert comm(q,a)==a and parity*a+a*parity==s.zeros(16)
        assert protected*a*protected==s.zeros(16)
        for k in range(4):
            b=creation(k)
            assert a*b+b*a==s.zeros(16)
            assert a.H*b+b*a.H==(s.eye(16)if j==k else s.zeros(16))
        checks+=1
    e0=-(u+s.sqrt(u*u+4*c*c))/2
    eplus=(-u+s.sqrt(u*u+4*c*c))/2
    expected_spectrum={s.Integer(0):3,-u:3,-u/2-c/2:4,-u/2+c/2:4,e0:1,eplus:1}
    assert h.eigenvals()==expected_spectrum
    assert all(s.simplify(energy-e0).is_nonnegative for energy in expected_spectrum)
    ground=(h-e0*s.eye(16)).nullspace()
    assert len(ground)==1
    v=ground[0];norm=s.simplify((v.H*v)[0])
    assert comm(h,q)==s.zeros(16)and q*v==s.zeros(16,1)
    weights=[]
    for j in range(4):
        a=creation(j)
        plus=s.simplify((v.H*a.H*a*v)[0]/norm)
        minus=s.simplify((v.H*a*a.H*v)[0]/norm)
        assert plus==minus==s.Rational(1,2)
        weights.append({'mode':j,'addition_weight':str(plus),'removal_weight':str(minus)})
    expected_gap=(s.sqrt(u*u+4*c*c)-c)/2
    triplet_gap=(s.sqrt(u*u+4*c*c)-u)/2
    sectors=[]
    for number in (1,3):
        indices=[i for i in range(16)if i.bit_count()==number]
        values=h.extract(indices,indices).eigenvals()
        minimum=-u/2-c/2
        assert values=={minimum:2,-u/2+c/2:2}
        assert s.simplify(minimum-e0-expected_gap)==0
        sectors.append({'occupation':number,'lowest_excitation':str(expected_gap)})
    assert len((h+u*s.eye(16)).nullspace())==3
    return {'full_source_dimension':16,'source_c':str(c),'source_U':str(u),
            'source_hopping_is_c_over_two':True,'exact_CAR_modes_checked':checks,
            'source_odd_creator_compression_to_single_occupancy_zero':True,
            'unique_centered_ground_occupation':2,'ground_energy':str(e0),
            'complete_source_eigenvalues':{str(k):int(v)for k,v in expected_spectrum.items()},
            'charge_sectors':sectors,'spin_triplet_excitation':str(triplet_gap),
            'complete_addition_removal_weights':weights,
            'finite_odd_and_spin_excitations_both_nonzero':True,
            'finite_dimer_identified_as_relativistic_vacuum':False}


def dimer():return copy.deepcopy(_dimer())


def internal_and_transport():
    components=source_components()
    assert len(components)==15 and all(row['occupation_mask'].bit_count()in (2,4)for row in components)
    return {'internal_source_components':list(components),
            'all_fifteen_internal_occupation_parities_even':True,
            'actual_weak_color_Higgs_current_identities':gauge_and_higgs(),
            'same_singlet_gauge_transport_countercontrol':neutralized_control(-4),
            'finite_internal_parity_is_physical_fermion_statistics':False,
            'representation_covariance_selects_transport_parameter':False}
