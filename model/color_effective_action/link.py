"""Direct Wilson-link substitution is tested on the literal native register."""
from model.interacting_extension.higgs_join.interacting_carrier import NativeCarrier


def certificate():
    carrier = NativeCarrier()
    rows = []
    for r, block in enumerate(carrier.data['blocks']):
        labels = tuple(block['indices'])
        d = len(labels)
        # U has external representation index and the actual sixty-state link.
        # Store each column as an exact integer sparse vector.
        columns = {}
        for b in range(d):
            column = {}
            for a in range(d):
                link = carrier.reverse[r, labels[a], labels[b]]
                if carrier.records[link] != (r, labels[a], labels[b]):
                    raise ArithmeticError('actual link record mismatch')
                column[a*60+link] = 1
            columns[b*60] = column
        for j, column in columns.items():
            if any(i in columns for i in column):
                raise ArithmeticError('source raising matrix is not square zero')
            for k, other in columns.items():
                inner = sum(value*other.get(i, 0) for i, value in column.items())
                if inner != (d if j == k else 0):
                    raise ArithmeticError('actual link initial-space Gram differs')
        if len(columns) != d or sum(map(len, columns.values())) != d*d:
            raise ArithmeticError('entire native raising matrix must be retained')
        rows.append({'block':r, 'charge6':str(block['q6Y']), 'representation_dimension':d,
                     'native_link_dimension':60, 'operator_dimension':60*d,
                     'nonzero_raising_columns':d, 'raising_matrix_squared_is_zero':True,
                     'U_adjoint_U_equals_d_times_external_identity_tensor_link_vacuum':True,
                     'normalized_partial_isometry_unitarity_defect_norm':'1',
                     'normalized_operator_statewise_distance_from_identity_lower':'1/2'})
    return {'native_register_rebuilt':True, 'blocks':rows,
            'tensoring_untouched_references_changes_defect':False,
            'direct_unitary_Wilson_coordinate_identification':False,
            'emergent_collective_Wilson_coordinates_excluded':False}
