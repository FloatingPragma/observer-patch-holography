"""Guard exact native-weight multiplicities and complete sector partitions."""
from .sector_counts import sector_counts


def test_exact_uniform_two_cell_sector_counts():
    result = sector_counts()
    neutral, compensated = result['sectors']
    assert result['native_internal_indices'] == [9, 10, 11]
    assert neutral['fixed_right_boundary_hypercharge'] == 0
    assert compensated['fixed_right_boundary_hypercharge'] == -1
    for sector in (neutral, compensated):
        assert sector['active_source_dimension_before_Gauss_constraints'] == 614400
        assert sum(row['dimension'] for row in sector['blocks']) == sector['diagonal_Gauss_support_dimension']
        assert len(sector['blocks']) == sector['conserved_block_count']
        assert all(row['dimension'] > 0 for row in sector['blocks'])
        # Spin-singlet interactions and the complete two spin components
        # retain reflected spin multiplicities, including charged sectors.
        table = {(row['L'], row['two_Sz']): row['dimension'] for row in sector['blocks']}
        assert all(table[ell, -spin] == count for (ell, spin), count in table.items())
    assert (neutral['diagonal_Gauss_support_dimension'], neutral['conserved_block_count'],
            neutral['largest_diagonal_block_dimension']) == (4140, 33, 642)
    assert (compensated['diagonal_Gauss_support_dimension'], compensated['conserved_block_count'],
            compensated['largest_diagonal_block_dimension']) == (2496, 25, 354)
    assert neutral['largest_block_labels'] == [{'L': 2, 'two_Sz': 0}]
    assert compensated['largest_block_labels'] == [{'L': 1, 'two_Sz': -1}, {'L': 1, 'two_Sz': 1}]
    assert not result['full_weak_Gauss_projection_computed']
    assert not result['Hamiltonian_or_ground_computed']
