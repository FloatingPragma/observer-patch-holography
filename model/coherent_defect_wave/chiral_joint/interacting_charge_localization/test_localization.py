"""Independent native hopping, cut geometry and full spectral projection tests."""
from itertools import combinations
import pytest
import sympy as s
from . import geometry,localization,screening
from ..periodic_family.model import Configuration,kinetic_terms
from ..peter_weyl_links.model import hopping_column
from ..interacting_charged_response.model import FAMILIES


def test_all_minimal_native_flux_supports_and_disconnected_reference_cut():
    result=geometry.certificate();edges,distance=geometry.graph()
    assert result['distance_shell_multiplicities']=={'1':3,'2':3,'3':1}
    assert [r['minimal_support_count'] for r in result['all_nonorigin_cells']]==[2,2,8,2,8,8,48]
    # All fewer-than-distance supports fail connectivity, independently of
    # the minimal-path enumeration used by the receipt.
    for y,ell in distance.items():
        if y==0:continue
        for size in range(ell):
            assert all(y not in geometry.component(chosen) for chosen in combinations(edges,size))
    for mask in range(1,256,2):
        region={i for i in range(8) if mask>>i&1}
        for y in range(1,8):
            if y in region:continue
            # A component without outgoing flux contains nonzero net source
            # plus reference charge under the single-unscreened-charge rule.
            assert sum((1 if x==0 else -1 if x==y else 0) for x in region)==1


def test_independent_full_car_compiler_remote_pair_norm():
    term=next(t for t in kinetic_terms(Configuration()) if t.path==(9,))
    for family in ('ec','L'):
        modes=[15*spin+a for spin in range(2) for a in FAMILIES[family]]
        count=len(modes);total=s.Integer(0)
        for occupancy in range(1<<(2*count)):
            mask=sum(1<<(30*site+mode) for j,(site,mode) in enumerate(
                [(site,mode) for site in (3,7) for mode in modes]) if occupancy>>j&1)
            source_count=((mask>>(30*3))&((1<<30)-1)).bit_count()
            image=hopping_column((mask,(0,)*24,(0,)*8),(term,))
            for (out,links,scalars),value in image.items():
                if ((out>>(30*3))&((1<<30)-1)).bit_count()!=source_count-1:continue
                assert [i for i,x in enumerate(links) if x]==[9] and scalars==(0,)*8
                total+=s.conjugate(value)*value
        # Product one-body contraction equals that of the actual spin-mixed
        # local ground, whose within-family density is I/2.
        assert s.simplify(total/2**(2*count))==s.Rational(len(FAMILIES[family]),32)
    record=screening.certificate()
    assert all(r['excited_link_count']<r['remote_hole_distance_from_origin'] for r in record['records'])
    assert all(not r['flux_path_from_origin_to_remote_hole_exists'] for r in record['records'])


def test_ordered_projection_bound_with_a_dark_state_inside_the_full_bin():
    cosine,sine=s.Rational(99,101),s.Rational(20,101)
    rotation=s.eye(4)
    rotation[0,0]=rotation[2,2]=cosine
    rotation[0,2]=-sine;rotation[2,0]=sine
    atomic=s.diag(0,0,2,5);full=rotation*atomic*rotation.T
    perturbation=full-atomic;v=2*sine
    # Exact norm certificate without an eigensolver.
    assert perturbation**2==s.diag(v*v,0,v*v,0)
    pi=s.diag(0,0,1,1);energy=rotation*s.diag(1,1,0,0)*rotation.T
    assert energy**2==energy and full*energy==s.zeros(4)
    a,b=s.Integer(2),s.Rational(1,100)
    bound=localization.ordered_projection_bound(a,b,v)
    assert sine*sine<=bound<1
    # Both a visible low vector and a dark degenerate low vector are retained.
    assert s.trace(energy)==2
    assert pi*energy*pi==s.diag(0,0,sine*sine,0)
    assert atomic*pi*energy-pi*energy*full==-pi*perturbation*energy
    for args in ((0,0,1),(0,1,1),(2,1,-1),(2,1,0.1)):
        with pytest.raises(ValueError):localization.ordered_projection_bound(*args)


def test_actual_energy_bounds_and_complete_sequential_law():
    for family,ceiling in (('ec',s.Rational(866,10**17)),('L',s.Rational(170,10**17))):
        d=localization.data(family)
        assert 0<d['gamma'] and d['delta']+d['charged_gap']<1
        assert 0<d['remote_upper']<ceiling
        for row in d['rows']:
            ell=row['distance']
            # The nonminimal support bound dominates the constrained-path bound.
            assert ell+1-d['delta']>ell+d['charged_gap']
            assert row['spectral_separation_lower']>0
        result=localization.certificate(family)
        records=result['sequential_complete_original_preparation_records']
        assert len(records['records'])==records['coarse_record_count']==13
        assert s.Rational(records['unconditional_energy_bin_success_lower'])>0
        assert not result['energy_filtered_localization']['origin_probability_lower_inferred_by_discarding_screened_complement']
    q,a,b=s.symbols('q a b');w=s.symbols('w:9')
    total=(1-q)+q*(1-a)+q*a/2+q*(a/2-b)+q*b*sum(w)
    assert s.expand(total.subs(w[-1],1-sum(w[:-1])))==1


def test_missing_localization_dependency_is_rejected(monkeypatch):
    from . import certificate as custody
    monkeypatch.setattr(custody,'DEPENDENCIES',(*custody.DEPENDENCIES,'missing-localization-source.json'))
    with pytest.raises(FileNotFoundError,match='missing declared charge localization'):
        custody.input_hashes()
