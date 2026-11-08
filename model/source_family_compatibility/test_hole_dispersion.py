"""Independent source, sign, complete-record and continuum controls."""
from fractions import Fraction as F
from itertools import combinations, product

import numpy as np
import pytest
from scipy.linalg import expm

from . import hole_dispersion as h


@pytest.mark.parametrize('sides',[(2,1,1),(1,3,2),(3,4,5),(2,2,2)])
def test_actual_pi_flux_square_and_exact_open_box_norm(sides):
    audit=h.box_identity(sides)
    a=sum(h.directional_hopping(sides)).toarray()/2
    exact=np.sqrt(sum(np.cos(np.pi/(n+1))**2 for n in sides))
    values=np.linalg.eigvalsh(a)
    assert abs(values[-1]-exact)<3e-14
    assert abs(values[0]+exact)<3e-14
    assert audit['h_squared_separates_exactly']


def test_positive_hopping_substitution_is_detected(monkeypatch):
    monkeypatch.setattr(h.symmetry,'link_sign',lambda x,j:1)
    with pytest.raises(ArithmeticError,match='anticommute'):
        h.box_identity((2,2,2))


@pytest.mark.parametrize('sides',[(2,1,1),(3,1,1),(2,2,1)])
def test_every_hole_number_reduces_full_source_without_contact(sides):
    for q in range(np.prod(sides)+1):
        audit=h.hole_lift(sides,q)
        assert audit['full_source_invariant']
        assert audit['contact_restriction_zero']


def test_full_cube_one_hole_source_reduction():
    audit=h.hole_lift((2,2,2),1)
    assert audit['whole_number_dimension']==11440
    assert audit['aligned_dimension']==8


def test_hole_exchange_has_actual_car_minus_sign_and_determinant_propagation():
    sites=4;sides=(2,2,1);q=2
    _,edges=h.algebra.geometry(sides)
    masks,two_hopping,_=h.algebra.operators(sites,edges,sites-q)
    labels=tuple(combinations(range(sites),q));lookup={mask:i for i,mask in enumerate(masks)}
    inclusion=np.zeros((len(masks),len(labels)))
    for j,label in enumerate(labels):
        mask,sign=h._hole_column(sites,label);inclusion[lookup[mask],j]=sign
        reversed_mask,reversed_sign=h._hole_column(sites,tuple(reversed(label)))
        assert reversed_mask==mask and reversed_sign==-sign
    time=.37;speed=.125
    full=expm(-1j*time*speed*two_hopping.toarray()/2)
    one=-speed*sum(h.directional_hopping(sides)).toarray()/2
    u=expm(-1j*time*one)
    wedge=np.array([[np.linalg.det(u[np.ix_(y,x)])for x in labels]for y in labels])
    np.testing.assert_allclose(full@inclusion,inclusion@wedge,atol=4e-15,rtol=0)
    wrong=np.array([[np.prod(np.diag(u[np.ix_(y,x)]))+np.prod(np.diag(np.fliplr(u[np.ix_(y,x)])))
                     for x in labels]for y in labels])
    assert np.linalg.norm(wrong-wedge)>1e-4


def test_original_weyl_frame_and_bare_repair_non_covariance():
    row=h.frame_audit()
    assert row['actual_original_Weyl_frame_identity']
    assert (row['original_input_mask'],row['bare_original_repair_output_mask'],row['raw'])==(8,2,6)
    assert row['hidden_aligned_output_weight']==0


def test_site_position_and_actual_magnetic_translations_share_covariance():
    audit=h.translation_identity(7)
    assert audit['interior_columns']==64
    assert audit['one_step_magnetic_translations_anticommute']
    assert audit['two_step_translations_commute']
    assert not audit['finite_boundary_discarded']


@pytest.mark.parametrize('sites', [2,3,4])
def test_complete_tree_records_resolve_all_q_configurations_with_reference(sites):
    edges=tuple((x,x+1)for x in range(sites-1))
    for q in range(sites+1):
        labels=tuple(combinations(range(sites),q))
        masks=[h._hole_column(sites,label)[0]for label in labels]
        rows=h.configuration_scan(sites,q,edges)
        counts={mask:0 for mask in masks}
        # Reference entanglement: each source column has a distinct reference.
        # Every complete branch retains exactly its matching reference column.
        for word,row in rows.items():
            assert len(row)==1
            initial,(terminal,sign)=next(iter(row.items()))
            assert h.decode_configuration(sites,q,edges,word)==initial
            assert terminal.bit_count()==sites-q and sign in(-1,1)
            counts[initial]+=1
        assert set(counts.values())=={2**len(edges)}
        # Nonzero rows plus the zero-default padded alphabet form the full PVM.
        assert len(rows)<=32**len(edges)


@pytest.mark.parametrize('theta',[F(0),F(1,3),F(1)])
def test_many_hole_complete_dilation_gram_all_sectors(theta):
    edges=((0,1),(0,2),(1,3),(2,3))
    for q in range(5):
        row=h.many_hole_dilation(4,q,edges,theta)
        assert row['maximum_squared_step_difference']<=(1+theta)*F(2*q,4)
        if q in(0,4):assert row['maximum_squared_step_difference']==0


def test_native_energy_analyzer_retains_full_inputs_and_failure_records():
    row=h.native_dimer_analyzer()
    assert row['native_pulses']==78
    assert row['full_source_input_columns']==16
    assert row['all_first_failure_slots']==288
    assert row['terminal_E_raw_slots']==32
    assert row['numerical_admission_passed']
    assert F(row['full_dilation_operator_error_interval'][1])<F(1,10**8)
    assert row['lower_and_upper_one_hole_mode_position_probabilities']==[[1,0],[0,1]]


def test_deleted_native_export_is_rejected(monkeypatch):
    from model.field_thermal import thermal
    original=thermal.execute_with_exports
    def corrupt(*args,**kwargs):
        output,events=original(*args,**kwargs)
        return output,events[:-1]
    monkeypatch.setattr(thermal,'execute_with_exports',corrupt)
    with pytest.raises(ArithmeticError,match='omitted'):
        h.native_dimer_analyzer()


def test_bloch_matrix_is_actual_parity_cell_action():
    k=np.array([.21,-.33,.17]);gamma=h.gamma_matrices()
    matrix=h.bloch_h(k)
    np.testing.assert_allclose(matrix@matrix,sum(np.cos(k)**2)*np.eye(8),atol=2e-15,rtol=0)
    values=np.linalg.eigvalsh(matrix);r=np.linalg.norm(np.cos(k))
    np.testing.assert_allclose(values,[-r]*4+[r]*4,atol=3e-15,rtol=0)
    # Independently evaluate forward/backward signed plane waves at every cell.
    parities=tuple(product(range(2),repeat=3));lookup={x:i for i,x in enumerate(parities)}
    for column in range(8):
        for row,x in enumerate(parities):
            value=0j
            for j in range(3):
                for step in(-1,1):
                    y=list(x);y[j]+=step
                    residue=tuple(a%2 for a in y)
                    if lookup[residue]==column:
                        value+=h.symmetry.link_sign(x,j)*np.exp(1j*step*k[j])/2
            assert abs(value-matrix[row,column])<1e-15
    assert gamma[0].shape==(8,8)


def test_canonical_band_rotor_and_original_position_derivative_bound():
    rng=np.random.default_rng(101)
    for _ in range(12):
        k=rng.normal(size=3);k*=.45/np.linalg.norm(k)
        w,p,p0=h.band_identification(k)
        np.testing.assert_allclose(w.conj().T@w,np.eye(8),atol=8e-16,rtol=0)
        np.testing.assert_allclose(w@p0@w.conj().T,p,atol=8e-16,rtol=0)
        for j in range(3):
            delta=np.zeros(3);delta[j]=1e-6
            derivative=(h.band_identification(k+delta)[0]-h.band_identification(k-delta)[0])/(2e-6)
            assert np.linalg.norm(derivative,2)<2*np.linalg.norm(k)


def test_controlled_parabolic_dispersion_uses_gap_subtraction_not_relativistic_mass():
    rng=np.random.default_rng(901);c=1/8
    for spacing in(.125,.0625,.03125):
        for _ in range(20):
            xi=rng.normal(size=3);xi*=2/np.linalg.norm(xi)
            k=spacing*xi
            scaled=c*(np.sqrt(3)-np.linalg.norm(np.cos(k)))/spacing**2
            expected=c*np.dot(xi,xi)/(2*np.sqrt(3))
            bound=c*spacing**2*np.linalg.norm(xi)**4/3
            assert abs(scaled-expected)<bound
    row=h.continuum_budget(F(1,16),2,3,holes=2)
    assert row['q_hole_dynamics_operator_error']==2*row['one_hole_dynamics_operator_error']
    assert not row['relativistic_or_interacting_limit']


def test_reference_gap_rational_enclosure_and_finite_size_bound():
    for sides in((2,2,2),(7,8,9),(17,17,17)):
        row=h.finite_gap_interval(sides)
        bulk=.5-np.sqrt(3)/8
        finite=.5-np.sqrt(sum(np.cos(np.pi/(n+1))**2 for n in sides))/8
        assert float(row['bulk_gap_lower'])<=bulk<=float(row['bulk_gap_upper'])
        assert 0<=finite-bulk<=float(row['finite_box_gap_minus_bulk_upper'])
        assert row['physical_mass']is None


def test_default_source_mutation_is_not_silently_relabelled(monkeypatch):
    original=h.FieldGeometry
    monkeypatch.setattr(h,'FieldGeometry',lambda:original(speed=.25))
    with pytest.raises(ArithmeticError,match='defaults'):
        h.certificate()
