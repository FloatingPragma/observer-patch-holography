"""Complete charge/spin spectrum of the unchanged six-site joined field.

The exact SU(2) weight count fixes every multiplicity.  Numerical charts and
energies are admitted by residuals, not promoted to interval or continuum
certificates.  Every Fock sector and every highest-weight eigenvalue remains.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy import linalg, sparse

from ..full_fock_field.field import FullFock
from . import spin


SITES = 6
DIMENSION = 4**SITES
ADMISSION = 2e-10
ADDRESSES = np.arange(DIMENSION)
NUMBERS = np.array([int(i).bit_count() for i in ADDRESSES])
TWICE_MAGNETIZATIONS = np.array([sum(((int(i) >> (2*j)) & 1)-
                                    ((int(i) >> (2*j+1)) & 1) for j in range(SITES))
                                 for i in ADDRESSES])


def weight_indices(number, twice_m):
    if (type(number) is not int or not 0 <= number <= 12
            or type(twice_m) is not int or abs(twice_m) > 6 or (number-twice_m) % 2):
        raise ValueError('admitted integer number and compatible doubled spin weight required')
    return np.flatnonzero((NUMBERS == number) & (TWICE_MAGNETIZATIONS == twice_m))


def multiplicity(number, twice_spin):
    if (type(number) is not int or not 0 <= number <= 12 or type(twice_spin) is not int
            or not 0 <= twice_spin <= min(number, 12-number) or (number-twice_spin) % 2):
        raise ValueError('admitted number and total doubled spin required')
    up, down = (number+twice_spin)//2, (number-twice_spin)//2
    weight = math.comb(SITES, up)*math.comb(SITES, down)
    above = math.comb(SITES, up+1)*math.comb(SITES, down-1) if up < SITES and down else 0
    return weight-above


@lru_cache(maxsize=1)
def raising():
    rows, columns = [], []
    for mask in range(DIMENSION):
        for site in range(SITES):
            if mask >> (2*site) & 3 == 2:
                rows.append(mask ^ (3 << (2*site))); columns.append(mask)
    return sparse.coo_matrix((np.ones(len(rows)), (rows, columns)), shape=(DIMENSION,)*2).tocsr()


@lru_cache(maxsize=1)
def hubbard_generator():
    ops = spin.operators()
    return (ops.frame@ops.K@ops.frame.adjoint()).matrix(ops.scale).real.tocsr()


@dataclass(frozen=True)
class Chart:
    number: int
    twice_spin: int
    indices: tuple[np.ndarray, ...]
    columns: tuple[np.ndarray, ...]
    energies: np.ndarray
    energy_vectors: np.ndarray
    generator: np.ndarray
    largest_residual: float

    @property
    def dimension(self):
        return self.columns[0].shape[1]

    def complete_columns(self, weight=0, *, original=False, eigenbasis=False):
        if type(weight) is not int or not 0 <= weight <= self.twice_spin:
            raise ValueError('weight slot from highest to lowest required')
        result = np.zeros((DIMENSION, self.dimension), complex)
        values = self.columns[weight]
        result[self.indices[weight]] = values@self.energy_vectors if eigenbasis else values
        return spin.operators().frame.adjoint().matrix()@result if original else result


@lru_cache(maxsize=32)
def chart(number, twice_spin):
    expected = multiplicity(number, twice_spin)
    highest_indices = weight_indices(number, twice_spin)
    above = np.flatnonzero((NUMBERS == number) & (TWICE_MAGNETIZATIONS == twice_spin+2))
    matrix = raising()[above][:, highest_indices].toarray()
    highest = linalg.null_space(matrix, rcond=1e-13) if len(above) else np.eye(len(highest_indices))
    if highest.shape[1] != expected:
        raise ArithmeticError('numerical highest-weight rank differs from the exact count')
    columns, indices = [highest], [highest_indices]
    residuals = [np.linalg.norm(matrix@highest), np.linalg.norm(highest.T@highest-np.eye(expected))]
    for slot in range(twice_spin):
        next_indices = weight_indices(number, twice_spin-2*(slot+1))
        coefficient = math.sqrt((twice_spin-slot)*(slot+1))
        next_columns = raising().T[next_indices][:, indices[-1]]@columns[-1]/coefficient
        residuals.append(np.linalg.norm(next_columns.T@next_columns-np.eye(expected)))
        columns.append(next_columns); indices.append(next_indices)
    block = hubbard_generator()[highest_indices][:, highest_indices]
    module = highest.T@(block@highest)
    residuals.append(np.linalg.norm(module-module.T))
    energies, vectors = np.linalg.eigh(module)
    residuals.append(np.linalg.norm(module@vectors-vectors*energies))
    for where, value in zip(indices, columns, strict=True):
        block = hubbard_generator()[where][:, where]
        residuals.append(np.linalg.norm(block@value-value@module))
    largest = float(max(residuals))
    if not np.isfinite(largest) or largest > ADMISSION:
        raise ArithmeticError('same multiplicity generator fails its complete spin ladder')
    return Chart(number, twice_spin, tuple(indices), tuple(columns), energies, vectors, module, largest)


def census():
    rows = []
    for number in range(13):
        for twice_spin in range(number % 2, min(number, 12-number)+1, 2):
            value = chart(number, twice_spin)
            rows.append({'number': number, 'charge_from_half_filling': number-6,
                         'twice_spin': twice_spin, 'multiplicity': value.dimension,
                         'irrep_dimension': twice_spin+1,
                         'total_sector_dimension': value.dimension*(twice_spin+1),
                         'energies': value.energies.tolist(),
                         'largest_chart_Frobenius_residual': value.largest_residual})
        if sum(row['total_sector_dimension'] for row in rows if row['number'] == number) != math.comb(12, number):
            raise ArithmeticError('number-sector dimension does not include all states')
    if sum(row['total_sector_dimension'] for row in rows) != DIMENSION:
        raise ArithmeticError('highest-weight census is not the full Fock carrier')
    return rows


def ground():
    selected = chart(6, 0)
    state = selected.complete_columns(eigenbasis=True)[:, 0]
    return float(selected.energies[0]), state


def probe_spectrum(kind):
    """All spectral weights of a declared local spinor or triplet-pair probe.

    Operators and states are both in the exact Hubbard frame.  No spectral
    threshold drops small/dark lines.  These are numerical finite response
    weights about a supplied pure ground, not physical detector execution.
    """
    if type(kind) is not str or kind not in ('hole', 'particle', 'pair_hole', 'pair_particle'):
        raise ValueError('one declared charge-one or triplet-pair probe required')
    energy, vacuum = ground()
    field = FullFock(12)
    first, second = field.annihilator(0), field.annihilator(2)
    if kind == 'hole':
        operator, number, total_spin, weight = first, 5, 1, 1
    elif kind == 'particle':
        operator, number, total_spin, weight = first.conj().T, 7, 1, 0
    elif kind == 'pair_hole':
        operator, number, total_spin, weight = second@first, 4, 2, 2
    else:
        operator, number, total_spin, weight = first.conj().T@second.conj().T, 8, 2, 0
    response = operator@vacuum
    retained = chart(number, total_spin)
    basis = retained.complete_columns(weight, eigenbasis=True)
    coefficients = basis.conj().T@response
    weights = abs(coefficients)**2
    remainder = float(np.linalg.norm(response-basis@coefficients))
    norm = float(np.vdot(response, response).real)
    if remainder > ADMISSION or abs(float(sum(weights))-norm) > ADMISSION:
        raise ArithmeticError('declared probe response left its complete spin/charge sector')
    return {'kind': kind, 'number': number, 'twice_spin': total_spin,
            'hamiltonian_model_gaps': (retained.energies-energy).tolist(),
            'all_highest_weight_residues': weights.tolist(),
            'lowest_multiplet_gap': float(retained.energies[0]-energy),
            'lowest_multiplet_residue': float(weights[0]),
            'probe_norm_squared': norm,
            'outside_spin_charge_sector_norm': remainder,
            'physical_detector_executed': False,
            'repair_access': ('outside outcome only' if number in (5, 7)
                              else 'inside and outside outcomes both retained')}


def report():
    rows = census()
    vacuum_energy, _ = ground()
    minima = sorted((row['energies'][0], row['number'], row['twice_spin']) for row in rows)
    if (minima[0][1:] != (6, 0) or minima[1][0]-vacuum_energy <= ADMISSION
            or chart(6, 0).energies[1]-vacuum_energy <= ADMISSION):
        raise ArithmeticError('the declared half-filled unique ground changed')
    trace = sum((row['twice_spin']+1)*sum(row['energies']) for row in rows)
    square = sum((row['twice_spin']+1)*sum(x*x for x in row['energies']) for row in rows)
    operator = hubbard_generator()
    trace_error = abs(trace-float(operator.diagonal().sum()))
    square_error = abs(square-float((operator.multiply(operator)).sum()))
    if max(trace_error, square_error) > 2e-8:
        raise ArithmeticError('complete spectral moments disagree with the full action')
    probes = {kind: probe_spectrum(kind) for kind in ('hole', 'particle', 'pair_hole', 'pair_particle')}
    for kind in ('hole', 'pair_hole'):
        partner = kind.replace('hole', 'particle')
        if not np.allclose(probes[kind]['hamiltonian_model_gaps'],
                           probes[partner]['hamiltonian_model_gaps'], atol=ADMISSION, rtol=0):
            raise ArithmeticError('same centered particle/hole spectrum failed')
    return {'reference': 'unique finite six-site half-filled singlet ground of K=H-N/2',
            'reference_is_faithful': False,
            'ground_number': 6, 'ground_twice_spin': 0, 'ground_energy': vacuum_energy,
            'complete_Fock_dimension': DIMENSION,
            'full_charge_spin_spectrum': rows,
            'trace_residual': trace_error, 'second_trace_moment_residual': square_error,
            'numerical_chart_admission': ADMISSION,
            'complete_local_probe_responses': probes,
            'lowest_neutral_triplet_gap': float(chart(6, 2).energies[0]-vacuum_energy),
            'full_Hamiltonian_gap': float(min([chart(6, 0).energies[1]]+
                                             [row['energies'][0] for row in rows
                                              if (row['number'], row['twice_spin']) != (6, 0)])-vacuum_energy),
            'absolute_physical_masses_predicted': 0,
            'finite_gap_is_a_continuum_mass_claim': False,
            'interacting_spectrum_is_stationary_under_repair': False}
