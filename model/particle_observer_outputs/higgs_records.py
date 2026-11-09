"""Complete observer outputs on the actual finite quantum-Higgs sector.

The preparation pulse, source coefficient, number energy and program clock
are declared architecture data. No architectural-uniqueness premise or
measured particle input is needed for these finite output laws.
"""
from __future__ import annotations

from functools import lru_cache
import json

import sympy as s

from model.particle_calibration_selection import phase_attenuation as source
from model.particle_candidate_atlas import attenuation


@lru_cache(maxsize=32768)
def _clean(value):
    return s.simplify(s.expand(value))


def _matrix(value):
    value = s.ImmutableSparseMatrix(value)
    return s.ImmutableSparseMatrix(value.rows, value.cols,
        {key: item for key, entry in value.todok().items()
         if (item := _clean(entry)) != 0})


def _zero(value):
    return not _matrix(value).todok()


def _exact(value):
    value = s.sympify(value)
    if value.free_symbols or value.atoms(s.Float) or value.is_finite is not True:
        raise ValueError('a finite exact coefficient is required')
    return value


@lru_cache(maxsize=12)
def _sector_data(y, d):
    h, basis = source.zero_charge_block(y, pair_cost=d)
    n = len(basis)
    vacuum = s.ImmutableSparseMatrix(n, 1, {(0, 0): 1})
    p0 = vacuum*vacuum.H
    v = h*vacuum
    g2 = _clean(4*s.trace(y.H*y))
    if g2.is_positive is not True:
        raise ValueError('a nonzero family coupling is required for a bright line')
    if _clean((v.H*v)[0]-g2) != 0:
        raise ArithmeticError('the actual source bright norm differs')
    gap = s.sqrt(d*d+4*g2)
    bright = _matrix(v*v.H/g2)
    dark = _matrix(s.eye(n)-p0-bright)
    vertex = _matrix(v*vacuum.H+vacuum*v.H)
    ground = _matrix((1+d/gap)*p0/2+(1-d/gap)*bright/2-vertex/gap)
    excited = _matrix(p0+bright-ground)
    ground_energy, excited_energy = (d-gap)/2, (d+gap)/2
    if not _zero(h-(vertex+d*(bright+dark))):
        raise ArithmeticError('the full source block has an omitted term')
    if not _zero(h-ground_energy*ground-excited_energy*excited-d*dark):
        raise ArithmeticError('the spectral projectors do not resolve the source action')
    return {'Y': y, 'pair_cost': d, 'basis': basis, 'H': h,
            'vacuum': vacuum, 'vacuum_projector': p0, 'bright_vector': v,
            'bright_projector': bright, 'dark_projector': dark,
            'vertex': vertex, 'occupation': _matrix(s.eye(n)-p0),
            'ground_projector': ground, 'excited_projector': excited,
            'g_squared': g2, 'gap': gap,
            'ground_energy': ground_energy, 'excited_energy': excited_energy,
            'dark_energy': d}


def sector_data(y, pair_cost=1):
    """Exact full307-state Q*=N_F-2N_H=0 sector for a chosen family matrix."""
    y = s.ImmutableMatrix(y)
    if y.shape != (3, 3):
        raise ValueError('a three-by-three family matrix is required')
    y = s.ImmutableMatrix(y.applyfunc(_exact))
    d = _exact(pair_cost)
    if d.is_real is not True:
        raise ValueError('a real pair cost is required')
    return dict(_sector_data(y, d))


def coherent_sector(pass_amplitude='1', phase=s.S.One, pair_cost=1):
    """Use the complete coherent analyzer's configured family coefficient."""
    y = attenuation.family(pass_amplitude, phase)['Y']
    return sector_data(y, pair_cost)


def spectral_projectors(data):
    return {'ground': data['ground_projector'],
            'excited': data['excited_projector'], 'dark': data['dark_projector']}


def measurement_projectors(data, refined=False):
    """Complete, gauge-invariant occupation pointer, optionally refined."""
    if refined:
        return {'empty': data['vacuum_projector'],
                'bright': data['bright_projector'], 'dark': data['dark_projector']}
    return {'empty': data['vacuum_projector'], 'occupied': data['occupation']}


def preparation_angle(data):
    return s.acos(s.sqrt((1+data['pair_cost']/data['gap'])/2))


def preparation_unitary(data):
    """exp(+i eta V(iY)/g) on the full retained sector, fixing every dark input.

    The number term is paused or refocused during this declared preparation
    pulse. Free observation evolution uses H afterwards on the same carrier.
    """
    alpha = (1+data['pair_cost']/data['gap'])/2
    beta = (1-data['pair_cost']/data['gap'])/2
    v, vacuum = data['bright_vector'], data['vacuum']
    return _matrix(data['dark_projector']+s.sqrt(alpha)*(data['vacuum_projector']+data['bright_projector'])
                   +s.sqrt(beta/data['g_squared'])*(vacuum*v.H-v*vacuum.H))


def gauge_invariant_dark_probe(data):
    """A retained neutral pair state obtained from an orthogonal source vertex.

    Every W(X) is gauge invariant. Hilbert-Schmidt orthogonality to Y makes
    W(X)|Omega> orthogonal to the bright state, without assigning a charged
    occupation basis vector the status of a physical Gauss-kernel state.
    """
    y = data['Y']
    norm2 = s.trace(y.H*y)
    x = None
    for row in range(3):
        for col in range(3):
            seed = s.MutableDenseMatrix.zeros(3)
            seed[row, col] = 1
            candidate = _matrix(seed-s.trace(y.H*seed)*y/norm2)
            if _clean(s.trace(candidate.H*candidate)).is_positive is True:
                x = candidate
                break
        if x is not None:
            break
    if x is None:
        raise ArithmeticError('the family matrix complement disappeared')
    local, basis = source.zero_charge_block(x, pair_cost=0)
    if basis != data['basis']:
        raise ArithmeticError('the dark probe changed the actual source carrier')
    vector = local*data['vacuum']
    vector = _matrix(vector/s.sqrt((vector.H*vector)[0]))
    if not _zero(data['dark_projector']*vector-vector):
        raise ArithmeticError('the orthogonal gauge-invariant probe is not dark')
    return {'family_probe': x, 'vector': vector, 'density': _matrix(vector*vector.H)}


def propagator(data, time):
    """The same full Hamiltonian clock, including the relative dark phase."""
    time = _exact(time)
    if time.is_real is not True:
        raise ValueError('real common evolution time required')
    # The spectral representation preserves the entire dark block and avoids
    # an unproved two-level replacement of the retained carrier.
    return _matrix(s.exp(-s.I*data['ground_energy']*time)*data['ground_projector']
                   +s.exp(-s.I*data['excited_energy']*time)*data['excited_projector']
                   +s.exp(-s.I*data['dark_energy']*time)*data['dark_projector'])


def sequential_kraus(data, time, refined=False):
    """K_(a,b)=P_b exp(-itH) P_a; every four or nine pointer outcome is kept."""
    pointers = measurement_projectors(data, refined)
    u = propagator(data, time)
    return {(a, b): _matrix(pb*u*pa)
            for a, pa in pointers.items() for b, pb in pointers.items()}


def transition_probability(data, time):
    time = _exact(time)
    if time.is_real is not True:
        raise ValueError('real common evolution time required')
    return _clean(4*data['g_squared']/data['gap']**2*s.sin(data['gap']*time/2)**2)


def occupation_joint_law(data, time, *, empty, bright, dark):
    """Unconditional joint probabilities for arbitrary admitted sector states.

    The three inputs are Tr(P_empty rho), Tr(P_bright rho), Tr(P_dark rho).
    Coherences are retained in the quantum instrument but cancel from these
    four pointer probabilities. No event probability is renormalized.
    """
    alpha, beta, delta = map(_exact, (empty, bright, dark))
    if any(v.is_nonnegative is not True for v in (alpha, beta, delta)) or _clean(alpha+beta+delta-1) != 0:
        raise ValueError('nonnegative original sector weights summing to one required')
    transition = transition_probability(data, time)
    return {('empty', 'empty'): _clean(alpha*(1-transition)),
            ('empty', 'occupied'): _clean(alpha*transition),
            ('occupied', 'empty'): _clean(beta*transition),
            ('occupied', 'occupied'): _clean(beta*(1-transition)+delta)}


def ground_joint_law(data, time):
    return occupation_joint_law(data, time,
        empty=(1+data['pair_cost']/data['gap'])/2,
        bright=(1-data['pair_cost']/data['gap'])/2, dark=0)


def ground_correlation(data, time):
    """Ordered quantum NH correlator, distinct from projective joint records."""
    time = _exact(time)
    if time.is_real is not True:
        raise ValueError('real common evolution time required')
    mean = _clean((1-data['pair_cost']/data['gap'])/2)
    residue = _clean(data['g_squared']/data['gap']**2)
    connected = residue*s.exp(-s.I*data['gap']*time)
    return {'mean': mean, 'residue': residue, 'connected': connected,
            'ordered': mean*mean+connected,
            'zero_frequency_weight': mean*mean,
            'positive_frequency': data['gap'], 'dark_weight': s.S.Zero}


def gap_interval_from_records(probability_t, probability_2t, time):
    """Exact gap interval from complete J_(empty,occupied) count intervals.

    The experiment uses identical preparation and pointer services at t and
    2t. Its declared unaliased window is 0<Delta*t<pi. Bounds on probabilities
    include sampling and complete implementation errors supplied by the
    apparatus contract. Conditioning on the first empty result is absent.
    """
    if len(probability_t) != 2 or len(probability_2t) != 2:
        raise ValueError('lower and upper probability endpoints required')
    p_lo, p_hi, q_lo, q_hi = map(_exact, (*probability_t, *probability_2t))
    time = _exact(time)
    if time.is_positive is not True or not 0 < p_lo <= p_hi <= 1 or not 0 <= q_lo <= q_hi <= 1:
        raise ValueError('positive time and a strictly positive denominator interval required')
    r_lo, r_hi = _clean(q_lo/(4*p_hi)), _clean(q_hi/(4*p_lo))
    if not 0 < r_lo <= r_hi < 1:
        raise ValueError('probability intervals do not isolate the admitted unaliased phase')
    lo = _clean(2*s.acos(s.sqrt(r_hi))/time)
    hi = _clean(2*s.acos(s.sqrt(r_lo))/time)
    return {'ratio_interval': (r_lo, r_hi), 'gap_interval': (lo, hi),
            'assumed_alias_window': '0<Delta*time<pi',
            'conditioning_on_first_record': False,
            'common_time_parameter': time}


def _string_records(records):
    return {a+','+b: str(value) for (a, b), value in records.items()}


def certificate():
    data = sector_data(s.diag(1, 2, 3), 1)
    projections = spectral_projectors(data)
    for key, projection in projections.items():
        if projection != projection.H or not _zero(projection*projection-projection):
            raise ArithmeticError('a source spectral pointer is not a projector')
        for other_key, other in projections.items():
            if key != other_key and not _zero(projection*other):
                raise ArithmeticError('spectral outcomes overlap')
    if not _zero(sum(projections.values(), s.zeros(307))-s.eye(307)):
        raise ArithmeticError('the complete retained sector was truncated')
    prep = preparation_unitary(data)
    if not _zero(prep.H*prep-s.eye(307)) or not _zero(prep*data['vacuum_projector']*prep.H-data['ground_projector']):
        raise ArithmeticError('the declared complete preparation failed')
    if not _zero(prep*data['dark_projector']-data['dark_projector']):
        raise ArithmeticError('preparation changed the retained dark inputs')
    nh = data['occupation']
    direct_residue = _clean(s.trace(data['ground_projector']*nh*data['excited_projector']*nh))
    if direct_residue != data['g_squared']/data['gap']**2:
        raise ArithmeticError('the full-state ground response residue differs')
    if not _zero(data['dark_projector']*nh*data['ground_projector']):
        raise ArithmeticError('the ground occupation pointer couples to omitted dark lines')
    time = s.pi/45
    first, second = ground_joint_law(data, time), ground_joint_law(data, 2*time)
    p, q = first['empty', 'occupied'], second['empty', 'occupied']
    readout = gap_interval_from_records((p, p), (q, q), time)
    if readout['gap_interval'] != (15, 15):
        raise ArithmeticError('the unconditional count inverse changed the gap')
    coherent = coherent_sector()
    coherent_response = ground_correlation(coherent, 0)
    dark_probe = gauge_invariant_dark_probe(data)
    if not _zero(data['H']*dark_probe['vector']-data['dark_energy']*dark_probe['vector']):
        raise ArithmeticError('the retained neutral dark probe changed its energy')
    return {'declared_fixture': {'Y': ['1', '2', '3'], 'pair_cost': '1',
                                'measured_inputs': [], 'dimension': 307},
            'actual_source': 'finite scalar vacuum plus conjugate doublet; charged full-CAR vertex',
            'complete_sector': {'ground_rank': int(s.trace(data['ground_projector'])),
                                'excited_rank': int(s.trace(data['excited_projector'])),
                                'dark_rank': int(s.trace(data['dark_projector'])),
                                'energies': ['-7', '8', '1']},
            'ground_preparation': {'unitary': 'exp(+i eta V(iY)/g)',
                'fixes_all_dark_inputs': True, 'number_term_paused_or_refocused': True,
                'autonomous_ground_relaxation_claimed': False},
            'occupation_spectral_response': {'gap': '15', 'mean': '7/15',
                'connected_residue': str(direct_residue), 'dark_line_weight_in_ground': '0',
                'ordered_connected_correlator': '56 exp(-15 i t)/225'},
            'two_time_complete_records': {'time': str(time), 'first': _string_records(first),
                'double_time': _string_records(second), 'J_01_2t_divided_by_4J_01_t': '3/4',
                'recovered_gap': '15', 'all_pointer_outcomes_retained': True,
                'fixture_alias_window_from_declared_coefficients': 'Delta*time=pi/3',
                'postselection': False},
            'dark_negative_control': {'input': 'normalized W(X)Omega, Tr(Y^dagger X)=0',
                'gauge_invariant': True, 'occupied,occupied': '1',
                'trace_after_invalid_bright_only_projection': '0'},
            'coherent_source_control': {'phase': '1', 'pass_amplitude': '1', 'pair_cost': '1',
                'g_squared': str(coherent['g_squared']), 'gap': str(coherent['gap']),
                'ground_energy': str(coherent['ground_energy']),
                'mean_occupation': str(coherent_response['mean']),
                'connected_residue': str(coherent_response['residue']),
                'fixed_macro_time': '1/8',
                'joint_records_at_macro_time': _string_records(ground_joint_law(coherent, s.Rational(1, 8))),
                'reference_mass_inputs': []},
            'clock': 'the same declared U(t) and repeated macro evolution in both protocols',
            'scope': {'finite_interacting_gauge_invariant_line_derived': True,
                'complete_Q_zero_sector_retained': True,
                'all_other_charge_sector_spectra_computed': False,
                'native_pulse_execution_claimed_here': False,
                'relativistic_particle_identification': False,
                'physical_mass_or_GeV_conversion': False}}


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2, sort_keys=True))
