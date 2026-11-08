"""Exact rational reflector certificate for the selected source cube."""
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np

from model.spin_covariant_localization import localization
from model.centered_particle_spectrum.certificate import sqrt_bounds

HERE = Path(__file__).resolve().parent


def reflect(vector, normal):
    if len(vector) != len(normal) or not normal or not any(normal):
        raise ValueError('a nonzero complete reflection normal is required')
    denominator = sum(n*n for n in normal)
    factor = 2*sum(x*n for x, n in zip(vector, normal, strict=True))/F(denominator)
    return tuple(x-factor*n for x, n in zip(vector, normal, strict=True))


def proposal():
    """Untrusted numerical proposal; certify() uses exact rational residuals."""
    data = localization.source_dictionary()
    columns = np.column_stack([np.asarray(v, float)/np.sqrt(float(n))
                               for v, n in zip(data['columns'], data['norms'], strict=True)])
    normals = []
    for k in range(8):
        v = columns[:, k].copy()
        v[k] -= 1
        v[:k] = 0
        normal = tuple(int(round(float(x)*(1 << 50))) for x in v)
        if not any(normal):
            raise ArithmeticError('zero proposal normal')
        normals.append(normal)
        u = np.asarray(normal, float)
        columns -= 2*np.outer(u, u@columns)/(u@u)
    return {'normal_integer_columns': normals}


def witness():
    value = json.loads((HERE/'reflectors.json').read_text())
    if (set(value) != {'normal_integer_columns'}
            or type(value['normal_integer_columns']) is not list
            or len(value['normal_integer_columns']) != 8):
        raise ValueError('eight complete reflector normals are required')
    rows = value['normal_integer_columns']
    if any(type(v) is not list or len(v) != 256
           or any(type(n) is not int or abs(n) > 1 << 53 for n in v)
           or not any(v) for v in rows):
        raise ValueError('bounded complete integer reflector data are required')
    return tuple(tuple(v) for v in rows)


def certify(normals=None):
    data = localization.source_dictionary()
    normals = witness() if normals is None else normals
    if len(normals) != 8:
        raise ValueError('eight reflector normals required')
    residuals = []
    for k, (column, norm) in enumerate(zip(data['columns'], data['norms'], strict=True)):
        vector = tuple(F(int(x)) for x in column)
        for normal in normals:
            vector = reflect(vector, normal)
        # Orthogonal rational reflectors preserve the squared norm exactly.
        if sum(x*x for x in vector) != norm:
            raise ArithmeticError('rational reflector lost the complete norm')
        lower, upper = sqrt_bounds(norm, 180)
        if vector[k] <= 0:
            raise ArithmeticError('reflected source mode has the wrong target phase')
        residuals.append(max(F(), 2-2*vector[k]/upper))
    isometry_error = sqrt_bounds(sum(residuals, F()), 180)[1]
    if isometry_error >= F(1, 10**10):
        raise ArithmeticError('source character-to-register map exceeds its certified tolerance')
    # PVM Stinespring difference: rank-one rows <=2d in square sum;
    # complementary projection <=2d; sqrt(8)d<3d.
    return {'reflectors': 8, 'dimension': 256,
            'rational_orthogonal_transform_exact': True,
            'mode_residual_squared_bounds': [str(x) for x in residuals],
            'mode_isometry_error': str(isometry_error),
            'complete_detector_error_from_reflectors': str(3*isometry_error),
            'complete_detector_error_from_true_ground': str(data['instrument_error']),
            'actual_ground_preparation_vector_error': str(data['state_error'])}


def matrix(normals=None):
    """Numerical diagnostic of the certified rational operator, not its proof."""
    value = np.eye(256)
    for normal in witness() if normals is None else normals:
        v = np.asarray(normal, float)
        value -= 2*np.outer(v, v@value)/(v@v)
    walsh = np.array([[localization.character(i, j) for j in range(8)]
                     for i in range(8)])/np.sqrt(8)
    return np.kron(np.eye(32), walsh)@value
