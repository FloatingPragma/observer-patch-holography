"""Certified whole-domain perturbation, source inverse and scope checks."""
from copy import deepcopy
from fractions import Fraction as F
import json
import mpmath as mp
import pytest

if __package__:
    from . import capacity_nonselection as c
else:
    import capacity_nonselection as c


@pytest.fixture(scope='module')
def parent():
    return json.loads(c.owner.OUTPUT.read_text())


@pytest.fixture(scope='module')
def certificate():
    return c.certificate()


def test_actual_strict_images_survive_and_shapes_are_distinct(certificate):
    assert len(certificate['rows'])==2
    for row in certificate['rows']:
        domain=list(map(F,row['original_alpha_domain']))
        image=list(map(F,row['perturbed_alpha_image']))
        assert domain[0]<image[0]<image[1]<domain[1]
        assert F(row['closure_uniform_perturbation_upper'])==F(1,25*10**14)
        old,new=map(lambda key:list(map(F,row[key])),('original_chi_interval','perturbed_chi_interval'))
        assert new[0]>old[1]
        q0,qs=map(lambda key:list(map(F,row[key])),('original_Q_interval','perturbed_Q_interval'))
        assert qs[1]<q0[0]<q0[1]<F(2,3)


def test_complete_native_inverse_recomputed_for_perturbed_boxes(certificate,parent):
    for row,old in zip(certificate['rows'],parent['closures']):
        inverse=row['perturbed_complete_memory_inverse']
        assert inverse['cubic_root_count']==3 and inverse['quadratic_root_count']==2
        assert inverse['all_complete_record_probabilities_retained'] is True
        for branch in inverse['inverse_branches']:
            p,q=map(lambda key:list(map(F,branch[key])),('p','q'))
            assert 0<p[0]<p[1]<F(1,12) and 0<q[0]<q[1]<F(1,6)
        for after,before in zip(row['perturbed_masses_over_v'],old['implicit_and_mass_enclosures']['masses_over_v']):
            lo,hi=map(F,after)
            assert lo<F(before['lo'])<F(before['hi'])<hi


@pytest.mark.parametrize('shift',[True,1e-14,0,-1,F(1,1000),F(1,10**30),F(1,10**6)])
def test_shift_must_meet_all_stated_admission_conditions(parent,shift):
    with pytest.raises(ValueError):c.audit_row(parent['closures'][0],shift)


@pytest.mark.parametrize('mutation',['image','lipschitz','response','implicit','unique'])
def test_missing_owning_premise_is_rejected(parent,mutation):
    row=deepcopy(parent['closures'][0])
    if mutation=='image':row['fixed_point_alpha']['lo']=row['uniqueness_domain_alpha']['lo']
    elif mutation=='lipschitz':row['lipschitz_upper']='1'
    elif mutation=='response':row['implicit_and_mass_enclosures']['alpha_inv']['lo']='100'
    elif mutation=='implicit':row['implicit_and_mass_enclosures']['alpha_u_verification']['R_u_sign_definite']=False
    else:row['conditional_unique_root_in_domain']=False
    with pytest.raises(ValueError):c.audit_row(row)


@pytest.mark.parametrize('design',['mssm_beta4_proxy','sm_color_beta7'])
def test_exact_kernel_bound_at_actual_source_coordinates(parent,design):
    """Independent high-precision quadrature, without a perturbed root solve."""
    row=next(r for r in parent['closures'] if r['design']==design)
    with mp.workdps(70):
        p,u,mz,v=(mp.mpf(row['display'][k]) for k in ('p','alpha_u','mz','v'))
        phi=(1+mp.sqrt(5))/2
        chi=u*u/(512*(1-u/77))
        delta=mp.mpf(2)/9+(p-phi)**2/(21*mp.pi*(1-u/27-u*u/135))
        kappa=(u/50)/(1+u/31+u*u/310)
        def masses(shift):
            roots=[1+mp.sqrt(2)*mp.exp(-chi-shift)*mp.cos(delta+2*mp.pi*j/3) for j in range(3)]
            squares=[r*r for r in roots]
            scale=mp.root(mp.exp(3*kappa)/(2*6**14*mp.fprod(squares)),3)
            return [v*scale*s for s in squares],roots
        def response(mass):
            z=(mz/mass)**2
            return 2/mp.pi*mp.quad(lambda x:x*(1-x)*mp.log1p(z*x*(1-x)),[0,mp.mpf('.5'),1])
        old,roots=masses(0)
        assert min(roots)>mp.mpf('.01')
        for shift in (mp.mpf('1e-14'),mp.mpf('1e-4')):
            new,changed_roots=masses(shift)
            assert min(changed_roots)>mp.mpf('.01')
            assert all(abs(mp.log(b/a))<600*shift for a,b in zip(old,new))
            assert abs(mp.fsum(response(b)-response(a) for a,b in zip(old,new)))<400*shift
            assert abs(mp.fprod(new)/mp.fprod(old)-1)<mp.mpf('1e-65')


def test_selection_and_uniqueness_claims_are_separate(certificate):
    json.dumps(certificate,allow_nan=False)
    assert certificate['pixel_self_consistency_selects_attenuation_law'] is False
    assert certificate['old_fixed_law_unique_root_disputed'] is False
    assert certificate['original_source_law_or_receipts_modified'] is False
    assert certificate['parent_full_root_producer_cold_replayed_here'] is False
    assert certificate['physical_particle_assignments']==certificate['experimental_mass_inputs']==[]
    for row in certificate['rows']:
        assert row['perturbed_fixed_point_existence_certified'] is True
        assert row['source_program_is_identical_to_original'] is False
        assert row['perturbed_numeric_root_computed'] is False
        assert row['explicit_shift_root_uniqueness_certified'] is False
