"""Independent full-CAR, tile geometry, spectral and source controls."""
from copy import deepcopy
from fractions import Fraction as F
from math import comb

import numpy as np
import pytest
from scipy import sparse

from . import model as m


def _car(modes,i,j):
    # Independent Jordan-Wigner tensor construction.
    one=sparse.eye(2,format='csr');z=sparse.diags([1,-1],format='csr',dtype=float)
    annihilator=sparse.csr_matrix([[0,1],[0,0]])
    def c(k):
        result=sparse.csr_matrix([[1.]])
        for site in reversed(range(modes)):
            result=sparse.kron(result,z if site<k else annihilator if site==k else one,format='csr')
        return result
    return c(i).T@c(j)


def test_exact_integer_blocks_match_independent_full_car_operators():
    identity=sparse.eye(256,format='csr')
    D=4*identity
    for site in (0,4):
        for j in (0,1):D+=_car(8,site+j,site+j+2)+_car(8,site+j+2,site+j)
    Q=sum((_car(8,j,j) for j in (0,1,4,5)),sparse.csr_matrix((256,256)))
    for a in (0,1):
        for b in (0,1):Q-=2*_car(8,a,b)@_car(8,4+b,4+a)
    represented=np.zeros((256,256));representedD=np.zeros_like(represented)
    blocks=m.q6_integer_blocks()
    assert len(blocks)==25 and sum(len(row[1]) for row in blocks)==256
    assert max(len(row[1]) for row in blocks)==36
    for sector,states,d,q in blocks:
        represented[np.ix_(states,states)]=q;representedD[np.ix_(states,states)]=d
        assert all(((state&15).bit_count(),(state>>4).bit_count())==sector for state in states)
    np.testing.assert_array_equal(represented,Q.toarray())
    np.testing.assert_array_equal(representedD,D.toarray())
    assert sparse.linalg.norm(D@Q-Q@D)>0


@pytest.mark.parametrize('kind',['all','q6'])
def test_exact_local_pivots_enclose_the_literal_all_occupation_minimum(kind):
    value=m.local_allocation(kind);row=value['native_blocks'][-1]
    onsite=row['onsite_per_bond_endpoint'];w=row['inverse_electric_margin']
    lower,upper=row['dyadic_bracket'];minimum=0.
    for _,_,D,Q in m.q6_integer_blocks():
        matrix=float(onsite)*np.array(D)-float(w)*np.array(Q)
        minimum=min(minimum,float(np.linalg.eigvalsh(matrix)[0]))
    assert float(lower)-1e-12<=-minimum<float(upper)+1e-12
    assert m.q6_positive_shift(onsite,w,upper)
    assert not m.q6_positive_shift(onsite,w,lower)
    assert upper<row['coarse_occupation_loss_upper']
    assert value['nonconstant_link_margin_lower']>0
    assert row['full_Fock_dimension']==256


@pytest.mark.parametrize('side',[2,4,8])
def test_tiled_geometry_has_a_private_constant_link_on_every_cross_plaquette(side):
    def tile(x):return tuple(k//2 for k in x)
    def step(x,j):return tuple(k+(axis==j) for axis,k in enumerate(x))
    bonds=[];plaquettes=[]
    for x in np.ndindex((side,side,side)):
        for j in range(3):
            if x[j]+1<side:bonds.append((x,step(x,j)))
        for i in range(3):
            for j in range(i+1,3):
                if x[i]+1>=side or x[j]+1>=side:continue
                xi=step(x,i);xj=step(x,j);xij=step(xi,j)
                plaquettes.append(((x,xi),(xi,xij),(xj,xij),(x,xj)))
    inner_bonds=sum(tile(x)==tile(y) for x,y in bonds)
    inner_plaquettes=0
    for plaquette in plaquettes:
        crosses=[tile(x)!=tile(y) for x,y in plaquette]
        if not any(crosses):inner_plaquettes+=1
        else:
            # Every selected edge occurs once in this elementary loop.
            chosen=plaquette[crosses.index(True)]
            assert sum(edge==chosen for edge in plaquette)==1
            # Every native matrix entry has nonzero integer central charge.
            assert all(charge!=0 for _,_,charge in m.NATIVE)
    tiles=side**3//8
    assert inner_bonds==12*tiles and inner_plaquettes==6*tiles
    assert len(bonds)==3*side*side*(side-1)
    assert len(plaquettes)==3*side*(side-1)**2


def test_complete_collective_spin_decomposition_retains_all_twelve_pair_fock_states():
    pairs=12;dimension=0;maximum=F()
    # Each pair is one spin-doublet plus two spin-singlets. The binomial
    # difference gives all irreducible spins in the k active doublets.
    for active in range(pairs+1):
        inactive_multiplicity=comb(pairs,active)*2**(pairs-active)
        for q in range(active//2+1):
            twice_j=active-2*q
            multiplicity=comb(active,q)-(comb(active,q-1) if q else 0)
            dimension+=inactive_multiplicity*multiplicity*(twice_j+1)
            j=F(twice_j,2)
            for twice_m in range(-twice_j,twice_j+1,2):
                magnetic=F(twice_m,2)
                maximum=max(maximum,j*(j+1)-magnetic*(magnetic+1))
    assert dimension==4**12 and maximum==42


@pytest.mark.parametrize('native_dimension',[1,2])
def test_collective_commutator_and_norm_on_the_literal_full_fock_space(native_dimension):
    pairs=2*native_dimension;modes=2*pairs;dimension=1<<modes
    raising=sum((_car(modes,j,j+pairs) for j in range(pairs)),sparse.csr_matrix((dimension,dimension)))
    beta=sum((_car(modes,j,j)-_car(modes,j+pairs,j+pairs) for j in range(pairs)),sparse.csr_matrix((dimension,dimension)))
    D=beta+pairs*sparse.eye(dimension)
    np.testing.assert_allclose((raising.T@raising-raising@raising.T).toarray(),
                               (pairs*sparse.eye(dimension)-D).toarray(),atol=1e-14)
    largest=np.linalg.eigvalsh((raising.T@raising).toarray())[-1]
    assert abs(largest-native_dimension*(native_dimension+1))<1e-12


def test_complete_collective_instrument_keeps_reference_and_aggregate_failure():
    C=(_car(4,0,2)+_car(4,1,3)).toarray()/np.sqrt(2)
    Cs=(C,C.T);effect=sum(c.T@c for c in Cs)/2
    eig,vec=np.linalg.eigh(np.eye(16)-effect)
    assert min(eig)>-1e-13
    failure=(vec*np.sqrt(np.maximum(0,eig)))@vec.T
    kraus=[c/np.sqrt(2) for c in Cs]+[failure]
    np.testing.assert_allclose(sum(k.T@k for k in kraus),np.eye(16),atol=1e-13)
    psi=np.zeros(32);psi[8]=np.sqrt(.4);psi[21]=np.sqrt(.6)
    outputs=[np.kron(k,np.eye(2))@psi for k in kraus]
    assert abs(sum(v@v for v in outputs)-1)<1e-13
    success=sum(v@v for v in outputs[:-1])
    assert abs(success-psi@np.kron(effect,np.eye(2))@psi)<1e-13
    assert 0<success<1


def test_polynomial_passivity_and_centered_moment_on_a_mixed_source():
    energies=[F(0),F(2),F(7)];beta=F(4);scale=16;order=3;U=F(8)
    def taylor(x,k):
        result=term=F(1)
        for j in range(1,k+1):term*=x/j;result+=term
        return result
    weights=[taylor(beta*(U-e)/(2*scale),order)**(2*scale) for e in energies]
    probabilities=[x/sum(weights) for x in weights]
    B=[[F(),F(1,5),F(1,7)],[F(3,5),F(),F()],[F(1,4),F(),F()]]
    records=[(energies[j]-energies[i],probabilities[i]*B[j][i]**2)
             for i in range(3) for j in range(3)]
    chi=sum((probabilities[i]-probabilities[j])*B[j][i]**2 for i in range(3) for j in range(3))
    assert chi>0 and any(w>0 for e,w in records if e<0)
    total=sum(w for _,w in records);moment=sum((e-2)**2*w for e,w in records)
    low=F(1,100);high=F(20)
    result=m.passive_window(chi,beta,2,moment,low,high,mass_upper=1)
    original=sum(w for e,w in records if low<=e<=high)
    assert original>=result['original_transfer_weight_lower']>0
    assert total<1
    for energy in energies:
        x=beta*(U-energy)/(2*scale)
        slope=beta*taylor(x,order-1)/taylor(x,order)
        assert 0<slope<=beta


def test_positive_polynomial_cooling_and_populated_moving_band_do_not_imply_fixed_energy_survival():
    # Independent two-level countercontrol, not an actual primary output.
    # Q(E)=1+4(U-E) is a positive Taylor amplitude with beta=8.
    for n in (1,3,7,12):
        U=F(2**n);Z=(1+4*U)**2+1
        excited=1/Z;ground=1-excited
        assert excited*U<F(1,16)/U
        assert ground>F(9,10)
        fixed_band_weight=ground if F(1,4)<=U<=4 else F()
        if n>=3:assert fixed_band_weight==0


@pytest.fixture(scope='module')
def members():return tuple(m.contract.member(n) for n in (0,1,2))


def test_actual_uniform_density_and_level_one_q6_response(members):
    assert F('9.657')<m.uniform_density_bound()<F('9.658')
    for member in members:
        value=m.bind_member(member)['result'];level=member['source']['level']
        assert F(value['all_block_defect_density_upper'])<F(97,10)
        assert F(value['collective_commutator_lower'])>=F(1,20)
        assert F(value['collective_original_transfer_weight_lower'])==F(3,100)
        assert F(value['final_excitation_weight_lower'])==F(1,50)
        assert F(value['finite_record_transfer_weight_lower'])>=F(1,50)
        assert value['complete_probe_outcomes']==member['source']['side']**3+1
        assert not value['fixed_finite_energy_window_uniform_weight_proved']
        assert not value['clock_histogram_executed'] and value['physical_mass_table_GeV']==[]
        assert not value['tiled_ground_upper']['trial_replaces_selected_state']
        if level==1:
            assert F(value['q6_only_defects']['defect_density_upper'])<F('1.82')
            assert F(value['q6_averaged_matched_pair_weight_lower'])>F(9,100)
            assert F(value['q6_scalar_original_transfer_window']['original_transfer_weight_lower'])>F(7,100)
        if level==2:
            assert F(value['q6_averaged_matched_pair_weight_lower'])==0
            assert value['q6_scalar_original_transfer_window'] is None


def test_tiled_upper_bound_is_extensive_and_stays_in_the_actual_source(members):
    seed=m.tile_seed()
    assert F('-723.767820')<seed['dimensionless_tile_energy_upper']<F('-723.767819')
    for member in members:
        tile=m.tiled_upper(member);V=member['source']['side']**3
        assert tile['tile_count']==V//8
        assert tile['dimensionless_ground_upper']/V==seed['dimensionless_tile_energy_upper']/8
        assert tile['physical_ground_upper']*F(member['source']['spacing'])==tile['dimensionless_ground_upper']


def test_actual_bindings_reject_changed_source_state_clock_and_result(members):
    member=members[1];receipt=m.bind_member(member);m.verify_bound(member,receipt)
    wrong=deepcopy(receipt);wrong['result']['collective_original_transfer_weight_lower']='1'
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
    wrong=deepcopy(member);wrong['state']['beta']='1'
    with pytest.raises(ValueError):m.bind_member(wrong)
    wrong=deepcopy(member);wrong['clock']['GeV_or_proper_time_conversion']='fit'
    with pytest.raises(ValueError):m.bind_member(wrong)
    with pytest.raises(ValueError):m.verify_bound(members[0],receipt)
    with pytest.raises(ValueError):m.local_allocation('unproved')


def test_all_scientific_evidence_paths_exist():
    assert all((m.contract.ROOT/name).is_file() for name in m.SOURCE_INPUTS)
    assert m.native_control()==m.NATIVE
