"""One-history native occupation fields and quantitative Green-limit controls.

The observation-count CLT is distinct from spatial refinement. No coherent
Hamiltonian, independent replica or physical time unit enters either proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import math
from pathlib import Path
from types import MappingProxyType

import numpy as np
import sympy as sp

from . import construction as passage
from . import source_green
from . import ancestry_green
from . import uniform_green
from . import causal_placement
from ..observer import native
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
FORMAL_SOURCE_PATHS = (
    'Lean/EventAlgebra/Basic.lean',
    'Lean/EventAlgebra/Lueders.lean',
    'Lean/EventAlgebra/PartitionAverage.lean',
    'Lean/EventAlgebra/PartitionPinching.lean',
    'Lean/EventAlgebra/PublicRecordAlgebra.lean',
    'Lean/EventAlgebra/TwoScalePublicRepair.lean',
    'Lean/Geometry/CanonicalLorentzModule.lean',
    'Lean/Geometry/CausalOrderComposition.lean',
    'Lean/Geometry/CelestialNullCone.lean',
    'Lean/Geometry/EventGermDisplacement.lean',
    'Lean/Geometry/EventPopulationChartInterface.lean',
    'Lean/Geometry/LorentzOverlapCocycle.lean',
    'Lean/Geometry/MassShellKinematics.lean',
    'Lean/Geometry/NativeRepairSemanticPlacement.lean',
    'Lean/Geometry/ObserverFrameHyperboloid.lean',
    'Lean/Geometry/ObserverRestSpace.lean',
    'Lean/Geometry/SourceDerivedSpacetimeCarrier.lean',
    'Lean/ObserverPatchHolography/Locality/NoSignalling.lean',
    'Lean/ObserverPatchHolography/Provenance/CausalInterval.lean',
    'Lean/ObserverPatchHolography/Provenance/HistoryCausalInvariance.lean',
    'Lean/ObserverPatchHolography/Provenance/MismatchProvenance.lean',
    'Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean',
    'Lean/QFT/CoverageReceipt.lean',
    'Lean/QFT/FiniteCausalObserverNet.lean',
    'Lean/QFT/NoncommutativeWitness.lean',
    'Lean/QFT/ObserverAccessCut.lean',
    'Lean/QFT/ObserverNetDescent.lean',
    'Lean/QFT/RichFibreRegionalNet.lean',
    'Lean/QFT/RichFibreWitness.lean',
    'Lean/Screen/A5PortAction.lean',
    'Lean/Screen/A5PortModule.lean',
    'Lean/Screen/A5PrimitivePortPrediction.lean',
    'Lean/Screen/PrimitivePortTranslationBridge.lean',
    'Lean/Time/TimeOrderLedger.lean',
    'Lean/Tower/ConsensusTower.lean',
    'Lean/lake-manifest.json',
    'Lean/lakefile.lean',
    'Lean/lean-toolchain')
SOURCE_PATHS = tuple(sorted(set((*passage.SOURCE_PATHS, *source_green.SOURCE_PATHS,
    *ancestry_green.SOURCE_PATHS, *uniform_green.SOURCE_PATHS,
    *causal_placement.SOURCE_PATHS, *FORMAL_SOURCE_PATHS,
    'model/physical_space_time/occupation_field.py',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py',
    'model/__init__.py', 'model/common_source/__init__.py',
    'model/common_source/experiment.py', 'model/common_source/gauge_matter.py',
    'model/common_source/native_join.py', 'model/common_source/records.py',
    'model/common_source/scalar_coupling.py', 'model/observer/__init__.py',
    'model/observer/p_response/__init__.py', 'model/observer/p_response/experiment.py',
    'model/observer/scalar_response/__init__.py',
    'model/observer/source_scalar/__init__.py', 'model/observer/source_scalar/bounds.py',
    'model/observer/source_scalar/compiler.py', 'model/observer/source_scalar/cutoff.py',
    'model/observer/source_scalar/experiment.py', 'model/observer/source_scalar/field.py',
    'model/observer_reconstruction/__init__.py', 'model/physical_space_time/__init__.py',
    'model/relative_observer_geometry/__init__.py'))))
CHI = sp.Rational(5, 33)
P0 = sp.ImmutableMatrix(sp.eye(12)-sp.ones(12)/12)
MODES = ((1, 0), (0, 1), (1, 1), (2, 0))


def _integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('an original bounded integer is required')
    return value


def _rows(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def _laplacian(size, edges):
    value = sp.zeros(size)
    for a, b in edges:
        value[a, a] += 1; value[b, b] += 1
        value[a, b] -= 1; value[b, a] -= 1
    return sp.ImmutableMatrix(value)


@lru_cache(maxsize=1)
def native_counts():
    """Execute all 3960 native branches and decode every field by currents."""
    states = tuple(native.states(binary=True))
    index = {state: i for i, state in enumerate(states)}
    seams, laplacian = passage.recovered_carrier()
    for bound in (laplacian-2*P0, 10*P0-laplacian):
        _, certificate = bound[:11, :11].LDLdecomposition(hermitian=False)
        if any(certificate[i, i] <= 0 for i in range(11)):
            raise ArithmeticError('the exact source spectral bounds changed')
    a, b = sp.zeros(66), sp.zeros(66)
    checked = 0
    for column, before in enumerate(states):
        for edge in seams:
            for coin in (0, 1):
                after, (_, current) = native.native(before, edge, coin)
                decoded = tuple(before[i]+current*(int(i == edge[1])-int(i == edge[0]))
                                for i in range(12))
                if decoded != after:
                    raise ArithmeticError('native current disclosure lost occupation readback')
                (b if 0 in edge else a)[index[after], column] += 1
                checked += 1
    z = sp.Matrix(12, 66, lambda i, j: sp.Rational(states[j][i])-sp.Rational(1, 6))
    if z*z.T/66 != CHI*P0:
        raise ArithmeticError('canonical fixed-total covariance changed')
    if a != a.T or b != b.T or sp.ones(1, 66)*a != 50*sp.ones(1, 66):
        raise ArithmeticError('remote native counts lost symmetry or exposure')
    if sp.ones(1, 66)*b != 10*sp.ones(1, 66):
        raise ArithmeticError('incident exposure changed')
    if z*(a+b)/60 != (sp.eye(12)-laplacian/60)*z:
        raise ArithmeticError('native density closure failed')
    return states, sp.ImmutableMatrix(a), sp.ImmutableMatrix(b), sp.ImmutableMatrix(z), checked


@lru_cache(maxsize=2, typed=True)
def field(clock='attempt'):
    """Exact 12-dimensional density map and CLT covariance in the named clock."""
    if type(clock) is not str or clock not in ('attempt', 'incident'):
        raise ValueError('attempt or incident clock required')
    seams, laplacian = passage.recovered_carrier()
    green = (laplacian+sp.ones(12)/12).inv()-sp.ones(12)/12
    incident = seams if clock == 'attempt' else tuple(e for e in seams if 0 in e)
    m = len(incident)
    remote = _laplacian(12, tuple(e for e in seams if e not in incident))
    d = (2*m*sp.eye(12)+remote)/60
    b = sp.Rational(m, 30)*sp.eye(12)-(laplacian-remote)/60
    k = b*d.inv()
    covariance = CHI*(4*m*green+remote*green+green*remote-P0)
    inverse = (sp.eye(12)-k+sp.ones(12)/12).inv()-sp.ones(12)/12
    if covariance != CHI*(inverse+inverse.T-P0):
        raise ArithmeticError('stopped density Green sum disagrees with the clock correction')
    if inverse*(CHI*(P0-k*k.T))*inverse.T != covariance:
        raise ArithmeticError('the martingale Poisson covariance differs from the Green sum')
    # Principal minors certify positivity on the constant-free space exactly.
    _, pivots = covariance[:11, :11].LDLdecomposition(hermitian=False)
    if any(pivots[i, i] <= 0 for i in range(11)):
        raise ArithmeticError('limiting field covariance is not strictly positive')
    precision = (covariance+sp.ones(12)/12).inv()-sp.ones(12)/12
    if covariance*precision != P0:
        raise ArithmeticError('the field action is not the inverse covariance')
    return MappingProxyType({'clock': clock, 'm': m, 'laplacian': laplacian,
            'remote_laplacian': remote, 'green': sp.ImmutableMatrix(green),
            'kernel': sp.ImmutableMatrix(k), 'covariance': sp.ImmutableMatrix(covariance),
            'precision': sp.ImmutableMatrix(precision),
            'positive_LDL_pivots': tuple(pivots[i, i] for i in range(11))})


@lru_cache(maxsize=2, typed=True)
def configuration_kernel(clock='attempt'):
    """Independent complete 66-state kernel, retaining all native alternatives."""
    _, a, b, z, _ = native_counts()
    if clock == 'attempt':
        k = (a+b)/60
    elif clock == 'incident':
        k = b*(60*sp.eye(66)-a).inv(method='DM')
    else:
        raise ValueError('attempt or incident clock required')
    if z*k != field(clock)['kernel']*z or sp.ones(1, 66)*k != sp.ones(1, 66):
        raise ArithmeticError('the complete state law and reduced field law differ')
    if k*sp.ones(66, 1) != sp.ones(66, 1) or any(x < 0 for x in k):
        raise ArithmeticError('the uniform stopped state is not stationary')
    return sp.ImmutableMatrix(k)


def covariance_at_count(clock, count):
    """Covariance of count**(-1/2) sum_{j=0}^{count-1} eta_j; no CLT used."""
    _integer(count, 1, 128)
    k = field(clock)['kernel']
    answer, power = count*P0, P0
    for lag in range(1, count):
        power = k*power
        answer += (count-lag)*(power+power.T)
    return sp.ImmutableMatrix(CHI*answer/count)


def covariance_error_bound(clock, count):
    """Rational Euclidean operator bound, valid for every positive count."""
    _integer(count, 1, 10**12)
    field(clock)
    # L>=2 P0 and L<=10 P0; for incident D between (1/6)I and (1/3)I.
    # Similarity condition number <=sqrt(2)<2; full-attempt D=I.
    r, condition = (F(29, 30), F(1)) if clock == 'attempt' else (F(9, 10), F(2))
    return 2*F(CHI)*condition*r/(count*(1-r)**2)


def scalar_moments(clock, count):
    """Exact first four moments of the SAME history's port-zero sum.

    Integer score 6*n_0-1 is converted back to centered occupation. The
    recursion retains every terminal configuration and all counter moments.
    """
    _integer(count, 1, 64)
    states, _, _, _, _ = native_counts()
    kernel = configuration_kernel(clock)
    denominator = int(sp.ilcm(*(x.q for x in kernel)))
    counts = [[int(x*denominator) for x in row] for row in kernel.tolist()]
    destinations = [[(i, counts[i][j]) for i in range(66) if counts[i][j]] for j in range(66)]
    score = [6*state[0]-1 for state in states]
    values = [[s**order for s in score] for order in range(5)]
    for _ in range(1, count):
        propagated = [[0]*66 for _ in range(5)]
        for j in range(66):
            for i, weight in destinations[j]:
                for order in range(5):
                    propagated[order][i] += weight*values[order][j]
        values = [[sum(math.comb(order, old)*score[i]**(order-old)*propagated[old][i]
                       for old in range(order+1)) for i in range(66)] for order in range(5)]
    normalization = 66*denominator**(count-1)
    moments = tuple(F(sum(values[order]), normalization*6**order) for order in range(5))
    if moments[0] != 1 or moments[1] != 0:
        raise ArithmeticError('the full counter law lost normalization or centering')
    if moments[2]/count != covariance_at_count(clock, count)[0, 0]:
        raise ArithmeticError('independent complete-history second moment disagrees')
    return {'clock': clock, 'count': count, 'sum_moments': [str(x) for x in moments],
            'normalized_fourth_cumulant': str((moments[4]-3*moments[2]**2)/count**2),
            'every_terminal_configuration_retained': True}


def incomplete_tape_control():
    """Same initial state/local incident tape, different remote field values."""
    before = tuple(int(i in (1, 2)) for i in range(12))
    seams, _ = passage.recovered_carrier()
    if (1, 5) not in seams or (0, 10) not in seams:
        raise ArithmeticError('the native access counterexample changed its seams')
    outputs, tapes = [], []
    for coin in (0, 1):
        after, _ = native.native(before, (1, 5), coin)
        after, record = native.native(after, (0, 10), 0)
        outputs.append(after); tapes.append(record)
    if tapes[0] != tapes[1] or outputs[0] == outputs[1] or outputs[0][0] != outputs[1][0]:
        raise ArithmeticError('the incomplete-tape control no longer separates access')
    return {'initial_occupied_ports': [1, 2], 'remote_seam': [1, 5],
            'first_incident_seam': [0, 10], 'identical_incident_record': list(tapes[0]),
            'terminal_occupied_ports': [[i for i, x in enumerate(out) if x] for out in outputs],
            'port_only_tape_determines_the_full_field': False}


def _expanded_product(*matrices):
    answer = matrices[0]
    for matrix in matrices[1:]:
        answer = (answer*matrix).applyfunc(sp.expand)
    return sp.ImmutableMatrix(answer)


@lru_cache(maxsize=1)
def carrier_band():
    """Joint spatial/clock field on the actual record-recovered rank-three band.

    No coordinate frame is selected. The radial rank-one projector is
    determined by the observed port zero and the same invariant spatial Gram.
    """
    laplacian = field('attempt')['laplacian']
    checks, encoded = passage.ng._band_certificate(laplacian.tolist())
    if not all(checks.values()):
        raise ArithmeticError('record-recovered source band certificate failed')
    p = sp.ImmutableMatrix([[sp.Rational(a)+sp.sqrt(5)*sp.Rational(b)
                            for a, b in row] for row in encoded])
    radial = _expanded_product(4*p[:, 0], p[0, :])
    tangent = p-radial
    eigenvalue = 5-sp.sqrt(5)
    li = laplacian-field('incident')['remote_laplacian']
    if _expanded_product(p, li, p) != ((3-sp.sqrt(5))*radial/2+tangent/2).applyfunc(sp.expand):
        raise ArithmeticError('source incident-clock axial decomposition failed')
    result = {'projector': p, 'radial': radial, 'tangent': tangent,
              'laplacian_eigenvalue': eigenvalue}
    for clock in ('attempt', 'incident'):
        data = field(clock)
        covariance = _expanded_product(p, data['covariance'], p)
        radial_value = sp.simplify(sp.trace(covariance*radial))
        tangent_value = sp.simplify(sp.trace(covariance*tangent)/2)
        if covariance != (radial_value*radial+tangent_value*tangent).applyfunc(sp.expand):
            raise ArithmeticError('projected field has an unaccounted spatial component')
        precision = (radial/radial_value+tangent/tangent_value).applyfunc(sp.radsimp)
        if _expanded_product(covariance, precision) != p:
            raise ArithmeticError('projected Gaussian action is not its marginal precision')
        k = data['kernel']
        memory = _expanded_product(p, k, sp.eye(12)-p, k, p)
        result[clock] = MappingProxyType({'covariance': covariance,
            'precision': sp.ImmutableMatrix(precision),
            'radial_variance': radial_value, 'tangent_variance': tangent_value,
            'two_tick_closure_defect': memory})
    if result['attempt']['two_tick_closure_defect'] != sp.zeros(12):
        raise ArithmeticError('the native attempt band lost exact temporal closure')
    if result['incident']['two_tick_closure_defect'] == sp.zeros(12):
        raise ArithmeticError('the incident memory control lost discrimination')
    return MappingProxyType(result)


def carrier_band_receipt():
    data = carrier_band()
    return {'spatial_readout': 'record-recovered unit Gram 4P; P has rank three',
        'coordinate_frame_selected': False, 'radial_direction': 'P e_0',
        'laplacian_eigenvalue': str(data['laplacian_eigenvalue']),
        'same_stationary_native_history_and_access_contract_as_full_field': True,
        'clock_fields': [{
            'clock': clock,
            'radial_variance': str(data[clock]['radial_variance']),
            'each_tangential_variance': str(data[clock]['tangent_variance']),
            'marginal_action_radial_precision': str(sp.simplify(1/data[clock]['radial_variance'])),
            'marginal_action_each_tangential_precision': str(sp.simplify(1/data[clock]['tangent_variance'])),
            'two_tick_projected_linear_closure_defect_00': str(data[clock]['two_tick_closure_defect'][0, 0]),
            'projected_linear_transition_closes': clock == 'attempt'}
            for clock in ('attempt', 'incident')],
        'three_dimensional_physical_bulk_or_Lorentzian_metric_established': False}


def torus_edges(size):
    _integer(size, 4, 256)
    return tuple(sorted({tuple(sorted((size*x+y, size*((x+dx) % size)+(y+dy) % size)))
                         for x, y in product(range(size), repeat=2) for dx, dy in ((1, 0), (0, 1))}))


def torus_native_control(size=8):
    """Execute the unchanged native rule; topology is explicitly supplied."""
    _integer(size, 4, 8)
    edges = torus_edges(size)
    ports, m = size**2, len(edges)
    # Single-particle enumeration audits the complete density matrix. The local
    # four-state identity proves its extension to every fixed-number sector.
    matrix = sp.zeros(ports)
    cases = 0
    for before in product((0, 1), repeat=2):
        outputs = [native.native(before, (0, 1), coin)[0] for coin in (0, 1)]
        if sorted(outputs) != sorted((before, before[::-1])):
            raise ArithmeticError('native binary exchange identity failed')
    for source in range(ports):
        before = tuple(int(i == source) for i in range(ports))
        for edge in edges:
            for coin in (0, 1):
                after, _ = native.native(before, edge, coin)
                matrix[after.index(1), source] += 1
                cases += 1
    if matrix != 2*m*sp.eye(ports)-_laplacian(ports, edges):
        raise ArithmeticError('declared torus changed the actual native density law')
    return {'supplied_torus_side': size, 'ports': ports, 'seams': m,
            'native_branch_checks': cases, 'arbitrary_binary_sector_identity_cases': 4,
            'this_is_the_existing_spherical_source_screen': False}


def torus_bounds(size):
    """Certified fixed-band Green and incident-clock error, using pi**2<10.

    Each listed integer wavevector supplies real cos and sin coordinates.
    The spatial scale is a=N**2/(4*pi**2). This is a declared torus chart.
    """
    _integer(size, 8, 65536)
    errors, positive = [], []
    for k in MODES:
        energy = sum(x*x for x in k)
        loss = F(10*sum(x**4 for x in k), 3*size**2)
        lower = energy-loss
        if lower <= 0:
            raise ArithmeticError('fixed Fourier band exceeds the certified positive domain')
        errors.append(loss/(energy*lower)); positive.append(lower)
    green = max(errors)
    # ||W^T L_inc W|| <= 16*pi**2 sum_real_modes |k|**2 / N**4.
    # W is an exact orthonormal sampling of the selected real Fourier modes.
    total_energy = 2*sum(sum(x*x for x in k) for k in MODES)
    incident = F(10, 4*size**2)+F(8*10*total_energy, 4*size**4)/min(positive)
    attempt = F(5, size**4)
    return {'size': size, 'real_mode_count': 2*len(MODES),
            'canonical_total': size**2//3,
            'susceptibility': str(F((size**2//3)*(size**2-size**2//3), size**2*(size**2-1))),
            'green_operator_error_bound': str(green),
            'attempt_clock_correction_bound': str(attempt),
            'incident_clock_correction_bound': str(incident),
            'attempt_total_covariance_error_bound': str(green+attempt),
            'incident_total_covariance_error_bound': str(green+incident),
            'clock_to_clock_covariance_error_bound': str(attempt+incident)}


def torus_diagnostic(size):
    """Independent floating band-matrix check; rational bounds certify errors."""
    bounds = torus_bounds(size)
    points = np.asarray(tuple(product(range(size), repeat=2)), dtype=float)
    columns, energies = [], []
    for k in MODES:
        phase = 2*np.pi*(points@np.asarray(k))/size
        columns.extend((np.sqrt(2)/size*np.cos(phase), np.sqrt(2)/size*np.sin(phase)))
        energy = (size/np.pi)**2*sum(math.sin(math.pi*x/size)**2 for x in k)
        energies.extend((energy, energy))
    w = np.asarray(columns).T
    if np.linalg.norm(w.T@w-np.eye(len(columns))) > 1e-12 or np.linalg.norm(w.sum(axis=0)) > 1e-11:
        raise ArithmeticError('the declared Fourier observation map lost orthonormality')
    target = np.diag([1/sum(x*x for x in k) for k in MODES for _ in range(2)])
    green = np.diag(1/np.asarray(energies))
    local = np.zeros_like(green)
    for a, b in torus_edges(size):
        if 0 in (a, b):
            delta = w[a]-w[b]
            local += np.outer(delta, delta)
    a_scale = size**2/(4*np.pi**2)
    attempt = green-np.eye(len(columns))/(4*(2*size**2)*a_scale)
    incident = green+np.eye(len(columns))/(16*a_scale)-(local@green+green@local)/16
    errors = {'attempt': np.linalg.norm(attempt-target, 2),
              'incident': np.linalg.norm(incident-target, 2)}
    for clock, error in errors.items():
        if error > float(F(bounds[clock+'_total_covariance_error_bound']))+1e-12:
            raise ArithmeticError('Fourier covariance diagnostic exceeded its analytic allowance')
    return {**bounds, 'attempt_error_diagnostic': float(errors['attempt']),
            'incident_error_diagnostic': float(errors['incident']),
            'incident_cross_mode_diagnostic': float(incident[1, 5])}


def subordination_budget(rate, jump_cutoff, count_cutoff):
    """Positive native-lag mixture for exp(-rate*sqrt(I-K)), with finite caps.

    ``rate`` is a declared mixture parameter. It is neither the native attempt
    count nor a distance. Omitted histories retain their missing mass.
    """
    if type(rate) not in (int, F) or rate <= 0:
        raise ValueError('an original positive rational mixture rate is required')
    _integer(jump_cutoff, 1, 4096)
    _integer(count_cutoff, 0, 64)
    rate = F(rate)
    try:
        numeric_rate = float(rate)
    except OverflowError as error:
        raise ValueError('the mixture rate must admit a finite float diagnostic') from error
    if not math.isfinite(numeric_rate) or numeric_rate == 0:
        raise ValueError('the mixture rate must admit a positive finite float diagnostic')
    weights = tuple(F(math.comb(2*n, n), 4**n*(2*n-1))
                    for n in range(1, jump_cutoff+1))
    tail = F(math.comb(2*jump_cutoff, jump_cutoff), 4**jump_cutoff)
    if sum(weights, F(0))+tail != 1:
        raise ArithmeticError('the exact native-lag probability tail changed')
    # P(Poisson(rate)>J) <= E[binom(Poisson(rate),J+1)].
    poisson_bound = rate**(count_cutoff+1)/math.factorial(count_cutoff+1)
    error_bound = min(F(1), rate*tail+poisson_bound)
    log_retained_rate = math.log(numeric_rate)+math.log1p(-float(tail))
    mass = sum(math.exp(-numeric_rate+j*log_retained_rate-math.lgamma(j+1))
               for j in range(count_cutoff+1))
    if not 0 <= 1-mass <= float(error_bound)+2e-15:
        raise ArithmeticError('finite subordination omitted mass exceeded its budget')
    return {
        'rate': str(rate), 'jump_cutoff': jump_cutoff,
        'count_cutoff': count_cutoff,
        'maximum_native_attempt_lag': jump_cutoff*count_cutoff,
        'exact_single_jump_tail': str(tail),
        'poisson_count_tail_upper': str(poisson_bound),
        'operator_error_upper': str(error_bound),
        'retained_probability_mass_diagnostic': mass,
    }


def _subordinated_polynomial(kernel, rate, jump_cutoff, count_cutoff):
    """Evaluate the capped positive polynomial without conditioning on success."""
    budget = subordination_budget(rate, jump_cutoff, count_cutoff)
    identity = np.eye(len(kernel))
    power, jumps = identity.copy(), np.zeros_like(kernel)
    for n in range(1, jump_cutoff+1):
        power = power@kernel
        jumps += float(F(math.comb(2*n, n), 4**n*(2*n-1)))*power
    term, result = identity.copy(), identity.copy()
    for j in range(1, count_cutoff+1):
        term = (float(rate)/j)*(term@jumps)
        result += term
    return math.exp(-float(rate))*result, budget


def subordinated_native_readout():
    """Audit the Poisson extension as positive processing of native lag records."""
    _, a, b, exact_z, checked = native_counts()
    complete = np.asarray((a+b)/60, dtype=float)
    reduced = np.asarray(field('attempt')['kernel'], dtype=float)
    z = np.asarray(exact_z, dtype=float)
    projector = np.asarray(P0, dtype=float)
    reduced_generator = np.eye(12)-reduced
    values, vectors = np.linalg.eigh(reduced_generator)
    if abs(values[0]) > 1e-14 or np.min(values[1:]) <= 0:
        raise ArithmeticError('native density subordination lost its stationary mode')
    values[0] = 0.0
    rate, jump_cutoff, count_cutoff = F(1, 2), 64, 8
    full, budget = _subordinated_polynomial(
        complete, rate, jump_cutoff, count_cutoff)
    reduced_poly, _ = _subordinated_polynomial(
        reduced, rate, jump_cutoff, count_cutoff)
    target = (vectors*np.exp(-float(rate)*np.sqrt(values)))@vectors.T
    error = float(np.linalg.norm(reduced_poly-target, 2))
    covariance = z@full@z.T/66
    covariance_error = float(np.linalg.norm(
        covariance-float(CHI)*reduced_poly@projector, 2))
    mass = budget['retained_probability_mass_diagnostic']
    column_error = float(np.max(abs(full.sum(axis=0)-mass)))
    # Unit jumps give heat flow, not the square-root Poisson extension.
    unit_jump_target = (vectors*np.exp(-float(rate)*values))@vectors.T
    wrong_jump_error = float(np.linalg.norm(unit_jump_target-target, 2))
    if (float(np.min(full)) < -1e-15 or column_error > 2e-14
            or covariance_error > 2e-14
            or error > float(F(budget['operator_error_upper']))+2e-14
            or wrong_jump_error < 0.05):
        raise ArithmeticError('native positive-lag subordination control failed')
    return {
        **budget, 'native_branches_checked': checked,
        'complete_configurations': len(complete),
        'complete_history_covariance_intertwiner_error': covariance_error,
        'column_mass_error': column_error,
        'density_Poisson_extension_error_diagnostic': error,
        'unit_jump_heat_substitution_error_diagnostic': wrong_jump_error,
        'missing_mass_is_not_renormalized': True,
        'finite_native_record_polynomial_executed': True,
        'mixture_rate_identified_as_physical_radius_or_clock': False,
        'finite_sampling_error_included_in_truncation_bound': False,
    }


def run():
    _, _, _, z, checks = native_counts()
    rows = []
    for clock in ('attempt', 'incident'):
        data = field(clock)
        complete = configuration_kernel(clock)
        # Independent full-state lag laws, including noncommuting incident map.
        power = sp.eye(66)
        for lag in range(4):
            if z*power*z.T/66 != data['kernel']**lag*(CHI*P0):
                raise ArithmeticError('complete native lag covariance disagrees with the field')
            power = complete*power
        rows.append({'clock': clock, 'incident_seams': data['m'],
                     'covariance': _rows(data['covariance']),
                     'action_precision': _rows(data['precision']),
                     'positive_LDL_pivots': [str(x) for x in data['positive_LDL_pivots']],
                     'finite_count_covariance_error_bounds': {
                         str(n): str(covariance_error_bound(clock, n)) for n in (128, 4096, 65536)}})
    attempt, incident = field('attempt'), field('incident')
    correction = incident['covariance']/5-attempt['covariance']/30
    if correction == sp.zeros(12):
        raise ArithmeticError('the negative control lost its random-clock sampling defect')
    return {'schema': 'oph.physical_space_time.occupation_field.v1',
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
            'contract': {'stationary_binary_total': 2, 'ports': 12,
                         'stationary_configurations': 66, 'native_branch_checks': checks,
                         'clock_names': ['native attempts', 'port-zero incident attempts'],
                         'initial_occupation_disclosure_and_all_seam_currents_supplied': True,
                         'same_history_no_independent_replicas': True,
                         'physical_seconds_or_Hamiltonian_supplied': False},
            'field_laws': rows,
            'clock_negative_control': {'mean_attempts_per_incident_tick': 6,
                'equal_mean_duration_rescaling_makes_fields_equal': False,
                'scaled_covariance_difference_00': str(correction[0, 0]),
                'scaled_covariance_difference_11': str(correction[1, 1]),
                'scaled_covariance_difference_01': str(correction[0, 1])},
            'observer_access_negative_control': incomplete_tape_control(),
            'actual_carrier_spatial_clock_field': carrier_band_receipt(),
            'actual_source_screen_Green': source_green.run(),
            'actual_source_ancestry_covariances': ancestry_green.run(),
            'actual_source_uniform_Green_and_subsequential_field': uniform_green.run(),
            'actual_carrier_causal_placement': causal_placement.run(),
            'native_positive_lag_Poisson_extension': subordinated_native_readout(),
            'finite_history_controls': [scalar_moments(clock, n) for clock, ns in
                (('attempt', (1, 4, 16, 64)), ('incident', (1, 4, 16))) for n in ns],
            'finite_cap_control': {'requested_incident_ticks': 64,
                'attempt_cap_per_interval': 256,
                'complete_output_timeout_bound': str(64*F(5, 6)**256),
                'retaining_timeouts_without_success_renormalization_required': True,
                'full_cap_256_history_distribution_executed': False},
            'declared_torus_native_control': torus_native_control(),
            'declared_torus_green_convergence': [torus_diagnostic(n) for n in (8, 16, 32, 64)],
            'acceptance': {'one_history_Gaussian_field_limit_and_action_proved': True,
                'actual_incident_clock_correction_derived': True,
                'actual_rank_three_carrier_marginal_clock_actions_proved': True,
                'actual_source_screen_Green_elimination_and_uniform_approximation_proved': True,
                'actual_source_fixed_ancestry_covariances_certified_at_levels_zero_through_three': True,
                'actual_source_all_level_Green_bound_and_vanishing_clock_correction_proved': True,
                'actual_source_nonzero_compatible_subsequential_Gaussian_law_proved': True,
                'native_positive_lag_Poisson_extension_with_finite_cap_bound': True,
                'conditional_graph_Green_convergence_criterion_proved': True,
                'supplied_torus_fixed_band_iterated_limit_proved': True,
                'uniform_finite_count_Gaussian_error_bound_proved': False,
                'existing_spherical_source_continuum_criterion_verified': False,
                'instantaneous_native_occupations_are_Gaussian': False,
                'interacting_CFT_or_Lorentzian_spacetime_derived': False,
                'observed_Universe_identified': False},
            'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--output', type=Path)
    modes.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native occupation field, clock correction and Green-limit controls replayed.')
    elif args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
