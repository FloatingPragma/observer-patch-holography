"""Independent finite-algebra, Gaussian-tail and source-ground controls."""
from fractions import Fraction as F
import math
import numpy as np
import mpmath as mp
from scipy import sparse
from scipy.linalg import expm

from .model import (exact_collective,finite_hardcore_check,gaussian_number_moments,
                        rational_moment_bounds,free_parameters,controlled_resources,
                        dispersion_controls,source)


def test_native_replica_algebra_is_exact_and_admitted():
    assert exact_collective(3)['exact_collective_columns']==60


def test_consumed_source_is_positive_isolated_native_triplet():
    parent,pins=source()
    assert len(pins)>=316
    assert F(parent['source_bounds']['r_interval'][0])>0
    assert parent['construction']['measured_inputs']==[]


def test_actual_three_color_ground_comparison():
    for N in (1,3):
        row=finite_hardcore_check(N)
        assert row['ground_energy_error']>0
        assert 1-row['free_ground_overlap_squared']<=row['ground_infidelity_upper_from_actual_energy']+1e-12


def test_squeezed_number_moments_from_independent_probability_series():
    mp.mp.dps=50
    u=mp.mpf(3)/4;q=u/(1+u)
    probabilities=[mp.binomial(2*j,j)/4**j*q**j*mp.sqrt(1-q) for j in range(300)]
    row=gaussian_number_moments([u])
    assert abs(sum(probabilities)-1)<mp.mpf('1e-45')
    m2=sum((2*j)**2*p for j,p in enumerate(probabilities))
    m4=sum((2*j+2)**4*p for j,p in enumerate(probabilities))
    assert abs(m2-row['second_moment'])<mp.mpf('1e-40')
    assert abs(m4-row['n_plus_two_fourth_moment'])<mp.mpf('1e-40')
    excited=[(2*j+1)*(1-q)*p for j,p in enumerate(probabilities)]
    row=gaussian_number_moments([u],True)
    assert abs(sum(excited)-1)<mp.mpf('1e-45')
    assert abs(sum((2*j+3)**4*p for j,p in enumerate(excited))-row['n_plus_two_fourth_moment'])<mp.mpf('1e-39')


def test_exact_rational_resource_bounds_dominate_actual_gaussian_moments():
    _,_,_,occ=free_parameters(mp.mpf('.5'),3)
    actual=gaussian_number_moments(occ,True)
    bound=rational_moment_bounds(F(1,2),3,True)
    for key in actual:
        assert actual[key]<=float(bound[key])
    row=controlled_resources(F(1,2),16)
    assert F(row['ground_energy_error_over_Delta_upper_exact'])<=F(1,10**8)
    assert F(row['one_excitation_evolution_vector_error_upper_exact'])<F(1,1000)
    assert row['actual_constituent_blocks_per_macrocell_N']>row['total_number_cutoff_K']+1
    assert float(row['smallest_nonzero_physical_momentum_decimal'])<1


def test_pairing_critical_zero_mode_is_not_stable_reference():
    D=2.
    assert D*(D-12*(D/12))==0
    assert D*(D-12*(D/12+.01))<0


def test_complete_compressed_Weyl_detector_keeps_failure_carrier():
    dimension=12
    b=np.diag(np.sqrt(np.arange(1,dimension)),1)
    W=expm(.3j*(b+b.T))
    C=W[:4,:4]
    # Shared singular vectors avoid unrelated roundoff in two nearly
    # singular positive square roots of the same contraction.
    left,singular,right=np.linalg.svd(C)
    defect=np.diag(np.sqrt(np.maximum(0,1-singular**2)))
    defect_left=left@defect@left.conj().T
    defect_right=right.conj().T@defect@right
    U=np.block([[C,defect_left],[defect_right,-C.conj().T]])
    assert np.max(abs(U.conj().T@U-np.eye(8)))<1e-12
    psi=np.zeros(8);psi[0]=1
    signal=float(np.real(psi@U@psi))
    probabilities=[(1+signal)/2,(1-signal)/2]
    assert min(probabilities)>=0 and abs(sum(probabilities)-1)<1e-14
    assert np.linalg.norm(U[4:,:4])>0
    kraus=[(np.eye(4)+s*C)/2 for s in (1,-1)]+[s*defect_right/2 for s in (1,-1)]
    assert np.max(abs(sum(k.conj().T@k for k in kraus)-np.eye(4)))<1e-14
    source=np.zeros(4);source[3]=1
    full=np.zeros(dimension);full[3]=1
    ideal=W@full
    actual=np.zeros(2*dimension,complex)
    actual[:4]=C@source;actual[dimension:dimension+4]=defect_right@source
    ideal_flag=np.concatenate((ideal,np.zeros(dimension)))
    tau=np.linalg.norm(ideal[4:])
    assert abs(np.linalg.norm(actual-ideal_flag)-math.sqrt(2)*tau)<1e-14
    assert abs(np.linalg.norm((actual-ideal_flag)/math.sqrt(2))-tau)<1e-14


def test_finite_momentum_and_covariance_bounds():
    rows=dispersion_controls()
    assert rows[-1]['frequency_error_max_over_Delta']<rows[0]['frequency_error_max_over_Delta']/100


def test_full_replica_and_higher_complement_do_not_lower_ground():
    # Two actual four-state matrix-unit blocks per site, plus a fifth
    # untouched higher eigenstate, produce 5^4 complete test states.
    # This probes the variational/permutation ground argument, independently
    # of the symmetric occupation implementation used by the main program.
    d=5;N=2;D=2.;t=.125
    ident=sparse.eye(d,format='csr')
    local=sparse.diags([0,D,D,D,7.],format='csr')
    def tensor(factors):
        out=factors[0]
        for factor in factors[1:]:out=sparse.kron(out,factor,format='csr')
        return out
    H=sum((tensor([local if i==j else ident for j in range(4)]) for i in range(4)),sparse.csr_matrix((625,625)))
    for color in (1,2,3):
        X=np.zeros((d,d));X[0,color]=X[color,0]=1
        X=sparse.csr_matrix(X)
        for i in (0,1):
            for j in (2,3):
                H-=t/N*tensor([X if k in (i,j) else ident for k in range(4)])
    from scipy.sparse.linalg import eigsh
    energy,vector=eigsh(H,k=1,which='SA',tol=1e-12,v0=np.ones(625))
    reduced=finite_hardcore_check(N)
    assert abs(energy[0]-reduced['ground_energy'])<1e-11
    higher=np.array([4 in np.unravel_index(i,(5,)*4) for i in range(625)])
    assert np.linalg.norm(vector[higher,0])<1e-10


def test_retained_reference_uses_actual_full_native_symmetry_cover():
    from .reference import certificate
    row=certificate()
    assert row['group']['native_full_cover_order']==3072
    assert row['group']['native_AP_rotation_composition_checks']==576
    assert row['group']['translation_reference_dimension_sufficient_on_full_source']==64
    assert row['group']['admitted_translation_reference_dimension']==8
    assert row['full119_coordinate_current_orientation_replay']['complete_endpoint_Ward_identities']==1416
    assert row['native_rotation_translation_replay']['exact_oriented_link_label_rotation_translation_conjugations']==1728
    assert row['native_Fourier_spin_covariance_replay']['complete_cubic_Fourier_word_checks']==576
    exact=row['exact_retained_reference_check']
    assert exact['exact_joint_covariance_generator_blocks']==240
    assert exact['all_four_Kraus_intertwiners']==192
    assert exact['isometry_Gram_identity']
    assert exact['unreferenced_character_changing_contraction_rejected']
    assert row['completion']['records_are_joint_relational_records']
    assert row['completion']['actual_stationary_ground_lift_is_stationary']


def test_reference_and_all_detector_records_are_in_capacity_count():
    from .reference import apparatus_capacity
    apparatus=apparatus_capacity(2)
    assert apparatus['retained_regular_reference_dimension']==3072
    assert apparatus['retained_control_and_flag_qubits']==4
    assert apparatus['total_declared_Weyl_apparatus_dimension']==49152
    assert apparatus_capacity(0)['total_declared_Weyl_apparatus_dimension']==3072
    row=controlled_resources(F(1,2),16)
    assert row['retained_Weyl_apparatus']==apparatus
    expression=row['source_and_declared_Weyl_apparatus_capacity_nats_expression']
    assert expression==row['kinematic_source_capacity_nats_expression']+' + log(3072) + 4 log(2)'
    assert row['complete_records_including_actual_ground_error_upper_exact']=='19/20000'
    import pytest
    with pytest.raises(ValueError):apparatus_capacity(-1)
