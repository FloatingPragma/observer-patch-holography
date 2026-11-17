"""Finite loop-flux derivatives, complete spectral sums and gauge covariance.

Flux is a dimensionless supplied Peierls phase. Derivatives use Hamiltonian
model time/energy, not a measured electromagnetic flux or physical current.
The exact identities stated in README hold algebraically; their evaluations
here are floating diagnostics on every state of the 256-dimensional Fock space.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json

from ..full_fock_field.field import FullFock
from ..reference_observation.observation import frozen


DIMENSION = 256
TOLERANCE = 5e-11
# Only the replay projection uses clustering. All dynamics and Kubo sums use
# the unclustered complete eigensystem. The retained fixtures have distinct
# gaps at least 2.7e-6 apart, well beyond this numerical resolution.
SPECTRAL_CLUSTER_WIDTH = 2e-10
SERIALIZATION_UNIT = 1e-11


def _real(value, name, *, nonnegative=False):
    if (type(value) not in (int, float) or not math.isfinite(value)
            or (nonnegative and value < 0)):
        raise ValueError(f'finite {"nonnegative " if nonnegative else ""}{name} required')
    return float(value)


class FluxRing:
    """Same four-site z-Weyl/contact action with one supplied loop flux.

    Bond x points from x to (x+1)%4 and has phase weights[x]*flux.
    The weights sum to one; changing their distribution changes the gauge.
    Zero hopping is admitted for the atomic control, not for a unique-ground
    response. No particle-number or low-energy projection is applied.
    """
    def __init__(self, flux=0., *, speed=.125, interaction=1., weights=(.25,)*4):
        self.flux = _real(flux, 'dimensionless flux')
        self.speed = _real(speed, 'coordinate speed', nonnegative=True)
        self.interaction = _real(interaction, 'onsite interaction', nonnegative=True)
        if type(weights) not in (tuple, list) or len(weights) != 4:
            raise ValueError('four finite flux-distribution weights required')
        self.weights = tuple(_real(w, 'flux-distribution weight') for w in weights)
        if not math.isclose(sum(self.weights), 1., rel_tol=0., abs_tol=1e-14):
            raise ValueError('flux-distribution weights must sum to one')
        if max(abs(self.flux), self.speed, self.interaction, *map(abs, self.weights)) > 100:
            raise ValueError('finite numerical domain bounds every input magnitude by 100')
        space = FullFock(8)
        a = tuple(space.annihilator(j) for j in range(8))
        self.occupations = tuple(frozen(np.diag(
            sum(((space.addresses >> (2*x+s)) & 1) for s in range(2)))) for x in range(4))
        self.number = frozen(sum(self.occupations))
        contact = sum(((space.addresses >> (2*x)) & 1)*
                      ((space.addresses >> (2*x+1)) & 1) for x in range(4))
        h = np.diag(self.interaction*contact).astype(complex)
        first = np.zeros_like(h); second = np.zeros_like(h); currents = []
        for x, weight in enumerate(self.weights):
            directed = sum((1j*self.speed*sign/2)*
                (a[2*x+s].conj().T@a[2*((x+1) % 4)+s])
                for s, sign in enumerate((1, -1))).toarray()
            directed *= np.exp(1j*weight*self.flux)
            bond = directed+directed.conj().T
            current = 1j*(directed-directed.conj().T)
            h += bond; first += weight*current; second -= weight**2*bond
            currents.append(frozen(current))
        self.hamiltonian = frozen(h)
        self.generator = frozen(h-self.interaction*self.number/2)
        self.current = frozen(first)
        self.diamagnetic = frozen(second)
        self.bond_currents = tuple(currents)
        values, vectors = np.linalg.eigh(self.generator)
        self.energies, self.vectors = frozen(values), frozen(vectors)
        self.scale = max(1., float(np.max(abs(values))))
        self.ground_degeneracy = int(np.count_nonzero(values-values[0] <= TOLERANCE*self.scale))
        self.residual = max(float(np.linalg.norm(self.generator@vectors-vectors*values, 2)),
                            float(np.linalg.norm(vectors.conj().T@vectors-np.eye(DIMENSION), 2)))
        if not math.isfinite(self.residual) or self.residual > TOLERANCE*self.scale:
            raise ArithmeticError('complete flux eigensystem failed numerical admission')

    def gauge_to(self, other):
        """G with H_other=G H_self G*: phase redistribution, not flux removal."""
        if (not isinstance(other, FluxRing) or self.flux != other.flux
                or self.speed != other.speed or self.interaction != other.interaction):
            raise ValueError('the same loop flux and local action are required')
        chi = [0.]
        for x in range(3):
            chi.append(chi[-1]+(self.weights[x]-other.weights[x])*self.flux)
        phases = sum(chi[x]*np.diag(self.occupations[x]) for x in range(4))
        return frozen(np.diag(np.exp(1j*phases)))

    def ground_response(self):
        """Complete nondegenerate perturbation response; never divides by a cutoff gap."""
        if self.ground_degeneracy != 1:
            raise ValueError('unique numerically resolved ground required; degeneracy is retained')
        ground = self.vectors[:, 0]
        gaps = self.energies-self.energies[0]
        amplitudes = self.vectors.conj().T@self.current@ground
        weights = abs(amplitudes)**2
        first = float(np.vdot(ground, self.current@ground).real)
        dia = float(np.vdot(ground, self.diamagnetic@ground).real)
        para = float(2*np.sum(weights[1:]/gaps[1:]))
        commutator = self.generator@self.current-self.current@self.generator
        double_commutator = self.current@commutator-commutator@self.current
        moment = float(.5*np.vdot(ground, double_commutator@ground).real)
        norm = float(np.vdot(self.current@ground, self.current@ground).real)
        residual = max(abs(float(weights.sum())-norm), abs(float(weights@gaps)-moment),
                       abs(float(weights[0])-first**2))
        if residual > TOLERANCE*self.scale**3:
            raise ArithmeticError('complete current spectral moment failed')
        return dict(ground_energy=float(self.energies[0]),
                    ground_number=float(np.vdot(ground, self.number@ground).real),
                    first_excitation_gap=float(gaps[1]),
                    flux_derivative=first, diamagnetic_expectation=dia,
                    paramagnetic_sum=para, flux_curvature=dia-para,
                    complete_gaps=frozen(gaps), complete_current_weights=frozen(weights),
                    current_zeroth_moment=norm, current_first_moment=moment,
                    stationary_bond_currents=tuple(float(np.vdot(ground, j@ground).real)
                                                  for j in self.bond_currents),
                    numerical_moment_residual=residual)


def _summary(ring):
    response = ring.ground_response()
    return dict(flux=ring.flux, flux_weights=ring.weights,
                **response,
                retained_spectral_entries=len(response['complete_gaps']))


def report():
    rows = [_summary(FluxRing(flux)) for flux in (0., .4, 1.2)]
    uniform = FluxRing(.4); seam = FluxRing(.4, weights=(0., 0., 0., 1.))
    change = uniform.gauge_to(seam)
    a, b = uniform.ground_response(), seam.ground_response()
    gauge = dict(hamiltonian_residual=float(np.linalg.norm(
        seam.generator-change@uniform.generator@change.conj().T, 2)),
        spectrum_residual=float(np.max(abs(uniform.energies-seam.energies))),
        flux_derivative_difference=abs(a['flux_derivative']-b['flux_derivative']),
        flux_curvature_difference=abs(a['flux_curvature']-b['flux_curvature']),
        uniform_diamagnetic=a['diamagnetic_expectation'],
        seam_diamagnetic=b['diamagnetic_expectation'],
        uniform_paramagnetic=a['paramagnetic_sum'], seam_paramagnetic=b['paramagnetic_sum'])
    free = FluxRing(.4, interaction=0.)
    free_response = free.ground_response()
    alpha = free.flux/4
    free_control = dict(flux=free.flux, ground_degeneracy=free.ground_degeneracy,
        computed_energy=free_response['ground_energy'],
        independent_Bloch_energy=-2*free.speed*(math.cos(alpha)+math.sin(alpha)),
        computed_curvature=free_response['flux_curvature'],
        independent_Bloch_curvature=free.speed*(math.cos(alpha)+math.sin(alpha))/8,
        zero_flux_ground_degeneracy=FluxRing(interaction=0.).ground_degeneracy)
    atom = FluxRing(.4, speed=0.)
    root = Path(__file__).resolve().parents[2]
    source_files = ('model/flux_response/response.py', 'model/flux_response/README.md',
                    'model/full_fock_field/field.py', 'model/reference_observation/observation.py',
                    'tools/strict_json.py')
    return dict(schema='oph.complete_finite_flux_response.v1',
        source_sha256={name: hashlib.sha256((root/name).read_bytes()).hexdigest()
                       for name in source_files},
        scope='SUPPLIED_FOUR_SITE_WEYL_CONTACT_LOOP_RESPONSE',
        complete_Fock_dimension=DIMENSION, coordinate_speed=.125, onsite_interaction=1.,
        chemical_potential=.5, flux_is_dimensionless_supplied_Peierls_phase=True,
        interacting_responses=rows, gauge_redistribution=gauge, free_control=free_control,
        atomic_control=dict(ground_degeneracy=atom.ground_degeneracy,
            current_norm=float(np.linalg.norm(atom.current, 2)),
            diamagnetic_norm=float(np.linalg.norm(atom.diamagnetic, 2))),
        reference='unique full-Fock ground of H(phi)-U*N/2 when numerically resolved',
        full_spectrum_and_zero_weight_entries_retained=True,
        exact_identity_numerical_evaluation=True,
        rigorous_numerical_error_bound=False,
        source_pulse_execution_added=False,
        dynamical_photon_or_physical_QED_identified=False,
        physical_clock_or_current_calibrated=False,
        thermodynamic_conductivity_or_Lorentz_limit_proved=False,
        integrated_completion_gate_accepted=False)


def _json_value(value):
    if isinstance(value, dict):
        return {k: _json_value(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError('nonfinite flux response diagnostic')
        return round(value, 11)+0.
    return value



def spectral_subspaces(row):
    """Basis-independent replay projection, retaining every state and weight.

    Equal-eigenvalue eigenvectors have no canonical basis. Their complete
    projector weight <0|J P_E J|0>, multiplicity and energy are invariant.
    Numerically ambiguous near-degeneracies are rejected, not merged into
    a claim of exact spectral degeneracy.
    """
    gaps, weights = row['complete_gaps'], row['complete_current_weights']
    if (type(gaps) is not list or type(weights) is not list or
            len(gaps) != DIMENSION or len(weights) != DIMENSION or
            row['retained_spectral_entries'] != DIMENSION):
        raise ValueError('complete 256-state spectrum and every current weight required')
    for values in (gaps, weights):
        if any(type(x) not in (int, float) or not math.isfinite(x) or x < 0 for x in values):
            raise ValueError('finite nonnegative spectral entries required')
    if gaps[0] != 0 or any(a > b for a, b in zip(gaps, gaps[1:])):
        raise ValueError('ordered ground-relative complete spectrum required')
    groups = []
    start = 0
    while start < DIMENSION:
        end = start+1
        while end < DIMENSION and gaps[end]-gaps[start] <= SPECTRAL_CLUSTER_WIDTH:
            end += 1
        if end < DIMENSION and gaps[end]-gaps[end-1] <= 100*SPECTRAL_CLUSTER_WIDTH:
            raise ValueError('unresolved spectral cluster separation')
        groups.append(dict(minimum_gap=gaps[start], maximum_gap=gaps[end-1],
                           multiplicity=end-start,
                           current_weight=math.fsum(weights[start:end])))
        start = end
    return groups


def _same_json_schema(expected, actual, path='root'):
    """Validate exact JSON types/shape and finiteness before spectral replay."""
    if type(expected) is not type(actual):
        raise ValueError(f'flux receipt type mismatch at {path}')
    if type(expected) is dict:
        if (any(type(key) is not str for key in expected) or
                expected.keys() != actual.keys()):
            raise ValueError(f'flux receipt keys mismatch at {path}')
        for key in expected:
            _same_json_schema(expected[key], actual[key], f'{path}.{key}')
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError(f'flux receipt length mismatch at {path}')
        for index, (a, b) in enumerate(zip(expected, actual, strict=True)):
            _same_json_schema(a, b, f'{path}[{index}]')
    elif type(expected) is float:
        if not math.isfinite(expected) or not math.isfinite(actual):
            raise ValueError(f'nonfinite flux receipt value at {path}')
    elif type(expected) not in (str, bool, int, type(None)):
        raise ValueError(f'unsupported flux receipt type at {path}')


def verify_receipt(expected, actual):
    """Compare physical spectral measures, never a degenerate basis choice.

    All nonspectral fields (including source hashes, budgets, aggregate
    responses and scope flags) still compare exactly. Group weight budgets
    are the sum of both lists' eleven-decimal serialization roundings plus
    two serialization units for binary64 spectral reconstruction; this does
    not enlarge an eigenvector-residue tolerance to hide a basis change.
    """
    _same_json_schema(expected, actual)
    left, right = copy.deepcopy(expected), copy.deepcopy(actual)
    try:
        if len(left['interacting_responses']) != len(right['interacting_responses']):
            raise ValueError('different retained flux domain')
        for a, b in zip(left['interacting_responses'], right['interacting_responses'], strict=True):
            ag, bg = spectral_subspaces(a), spectral_subspaces(b)
            if len(ag) != len(bg):
                raise ValueError('different complete spectral subspace count')
            for x, y in zip(ag, bg, strict=True):
                if x['multiplicity'] != y['multiplicity']:
                    raise ValueError('different spectral multiplicity')
                if any(abs(x[k]-y[k]) > 2*SERIALIZATION_UNIT
                       for k in ('minimum_gap', 'maximum_gap')):
                    raise ValueError('different spectral energy')
                weight_budget = (x['multiplicity']+2)*SERIALIZATION_UNIT
                if abs(x['current_weight']-y['current_weight']) > weight_budget:
                    raise ValueError('different complete spectral projector weight')
            for key in ('complete_gaps', 'complete_current_weights'):
                a.pop(key); b.pop(key)
    except (KeyError, TypeError, OverflowError) as error:
        raise ValueError('malformed complete spectral receipt') from error
    if left != right:
        raise ValueError('flux response receipt differs outside the spectral basis choice')
    return True


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', type=Path)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.write and args.check:
        parser.error('choose write or check')
    result = _json_value(report())
    rendered = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
    if args.check:
        try:
            verify_receipt(strict_json.load(args.check), result)
        except ValueError as error:
            raise SystemExit(str(error)) from error
        print('complete finite flux response replay passed')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as output:
            output.write(rendered)
    else:
        print(rendered, end='')


if __name__ == '__main__':
    main()
