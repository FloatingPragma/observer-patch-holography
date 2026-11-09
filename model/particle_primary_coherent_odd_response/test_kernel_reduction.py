"""Independent source rotation, Fourier normalization and cancellation tests."""
from itertools import product
import pytest
import sympy as s

from . import kernel_reduction as k


def test_actual_source_ordered_pair_orbits_and_reversals():
    x=k.orbit_control()
    assert x['proper_spatial_rotations']==24
    assert x['ordered_pair_orbit_sizes']==[8,24,24,8]
    assert x['every_ordered_pair_has_a_proper_rotation_swap']


def test_exact_fourier_shell_transform_inverse_and_both_menu_normalizations():
    # Separate direct sums over all ordered site pairs check every Walsh row.
    shells=(s.Rational(3),s.Rational(1,5),s.Rational(-1,11),s.Rational(1,17))
    eigen=k.walsh_eigenresponses(shells)
    assert k.invert_eigenresponses(eigen)==shells
    for q in k.CELLS:
        value=sum((-1)**sum(a*b for a,b in zip(q,x))
                  *(-1)**sum(a*b for a,b in zip(q,y))*shells[k.distance(x,y)]
                  for x,y in product(k.CELLS,repeat=2))/8
        assert value==eigen[sum(q)]
        assert value/48==(value/6)/8
    assert sum(n*v/48 for n,v in zip(k.MULTIPLICITIES,eigen))==shells[0]/6
    assert min(eigen)>0
    with pytest.raises(ValueError):k.radial_matrix((1,2,3))
    with pytest.raises(ValueError):k.invert_eigenresponses((1,2))


def test_positive_populated_rows_do_not_force_contrast_or_energy_shape():
    x=k.control()
    for label,row in x['diagonal_no_contrast_countercontrol'].items():
        values=tuple(map(s.Rational,row['four_eigenresponses']))
        assert len(set(values))==1 and values[0]>0
        floor=s.Rational(1,128) if label=='G1' else s.Rational(9,512)
        assert s.Rational(row['each_complete_menu_row_probability'])>floor
    contrast=x['intensity_only_contrast_countercontrol']
    assert len(set(contrast['four_eigenresponse_intensities']))==4
    assert not x['actual_primary_off_diagonal_shell_values_computed']
    assert not x['nonidentical_normalized_energy_shapes_proved']
    assert not x['full_component_or_cross_irrep_kernel_reduced_to_four_scalars']


def test_prepared_cubic_transport_can_cancel_in_faithful_original_response():
    x=k.prepared_transport_countercontrol()
    assert x['normalizer']==2240 and x['polynomial_amplitude_lower']==1
    assert x['prepared_inter_site_moments_orders_zero_to_four'][:3]==[0,0,0]
    assert abs(x['prepared_inter_site_moments_orders_zero_to_four'][3])==6
    assert set(x['original_off_diagonal_transfer_moments_orders_zero_to_four'])=={'0'}
    assert x['invariant_unitary_forces_all_off_diagonal_spectral_events_zero']
    assert not x['comparison_state_or_H_substituted_for_primary_source']
