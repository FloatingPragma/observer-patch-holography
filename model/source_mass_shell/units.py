"""Exact conditional SI conversion with an explicit Hamiltonian clock.

Energies refer to exp(-i H_model tau). Matching light speed fixes the length
to time ratio. Provenance fields declare premises, without validating physics.
"""
from dataclasses import dataclass
from fractions import Fraction as F

PLANCK_H_J_S = F('6.62607015e-34')
ELEMENTARY_CHARGE_C = F('1.602176634e-19')
GEV_IN_J = 10**9*ELEMENTARY_CHARGE_C
CESIUM_FREQUENCY_HZ = F(9192631770)
LIGHT_SPEED_M_S = F(299792458)
CESIUM_TRANSITION_ENERGY_GEV = PLANCK_H_J_S*CESIUM_FREQUENCY_HZ/GEV_IN_J
SI_SOURCE = 'https://www.nist.gov/pml/special-publication-330/sp-330-section-2'
CESIUM_PROCESS = 'unperturbed cesium-133 ground-state hyperfine transition'


def exact(value):
    """Reject implicit binary floating input in an exact interval claim."""
    if type(value) not in (int, str, F):
        raise ValueError('an exact integer, rational string or Fraction is required')
    try:
        return F(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError('invalid exact rational input') from exc


def interval(value, *, positive=False, nonnegative=False):
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise ValueError('an ordered pair of exact interval endpoints is required')
    lo, hi = map(exact, value)
    if lo > hi or (positive and lo <= 0) or (nonnegative and lo < 0):
        raise ValueError('interval endpoints violate the declared domain')
    return lo, hi


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(name+' must be an explicit nonempty string')
    return value


def _provenance(value, name):
    if not isinstance(value, tuple) or not value:
        raise ValueError(name+' must be a nonempty tuple of references or premises')
    return tuple(_text(x, name) for x in value)


def _encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: _encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_encode(item) for item in value]
    return value


def pi_interval(terms=24):
    """Machin identity and exact alternating-series outward bounds."""
    if type(terms) is not int or not 2 <= terms <= 128:
        raise ValueError('integer term count from 2 to 128 required')

    def atan_bounds(x):
        total = sum(((-1)**k*x**(2*k+1)/(2*k+1)
                     for k in range(terms)), F())
        adjacent = total+(-1)**terms*x**(2*terms+1)/(2*terms+1)
        return min(total, adjacent), max(total, adjacent)

    a, b = atan_bounds(F(1, 5)), atan_bounds(F(1, 239))
    return 16*a[0]-4*b[1], 16*a[1]-4*b[0]


def hbar_gev_seconds_interval():
    lo, hi = pi_interval()
    coefficient = PLANCK_H_J_S/(2*GEV_IN_J)
    return coefficient/hi, coefficient/lo


@dataclass(frozen=True)
class CsClockInput:
    """Conditional same-generator clock premises with no inferred identity.

    angular_gap_interval is a gap of H in exp(-i H tau). The three provenance
    groups identify that gap, the physical process and the common time.
    Nonempty references document inputs; this class does not prove them.
    """
    generator_id: str
    angular_gap_interval: tuple
    source_gap_provenance: tuple
    physical_process_provenance: tuple
    same_time_provenance: tuple
    physical_process: str

    def __post_init__(self):
        _text(self.generator_id, 'generator identity')
        if self.physical_process != CESIUM_PROCESS:
            raise ValueError('the SI frequency refers to the specified cesium process')
        object.__setattr__(self, 'angular_gap_interval',
                           interval(self.angular_gap_interval, positive=True))
        for name in ('source_gap_provenance', 'physical_process_provenance',
                     'same_time_provenance'):
            object.__setattr__(self, name, _provenance(getattr(self, name), name))


def c_normalized_units(seconds_per_source_time, model_signal_speed=1):
    """Conditional light-speed matching; no default absolute time unit."""
    tlo, thi = interval(seconds_per_source_time, positive=True)
    velocity = exact(model_signal_speed)
    if velocity <= 0:
        raise ValueError('positive derived model signal speed required')
    hlo, hhi = hbar_gev_seconds_interval()
    energy = hlo/thi, hhi/tlo
    return _encode({
        'status': 'CONDITIONAL_ON_SUPPLIED_SECONDS_PER_SOURCE_TIME',
        'model_evolution': 'exp(-i H_model tau)',
        'seconds_per_source_time_interval': (tlo, thi),
        'model_signal_speed': velocity,
        'matched_physical_signal_speed_m_s': LIGHT_SPEED_M_S,
        'length_unit_m_interval': (LIGHT_SPEED_M_S*tlo/velocity,
                                   LIGHT_SPEED_M_S*thi/velocity),
        'energy_unit_GeV_interval': energy,
        'energy_equivalent_kg_per_model_energy_interval': tuple(
            x*GEV_IN_J/LIGHT_SPEED_M_S**2 for x in energy),
        'length_relation': 'length_unit_m=c_SI*seconds_per_source_time/v_model',
        'energy_relation': 'energy_unit_GeV=hbar_GeV_s/seconds_per_source_time',
        'light_speed_matching_fixes_absolute_scale': False,
        'rest_mass_identification_from_units_alone': False,
        'physical_seconds_per_source_time_derived': False,
        'SI_source': SI_SOURCE,
    })


def energy_conversion(model_gap_interval, *, generator_id, clock=None,
                      model_signal_speed=1):
    """Nonnegative energy conversion, retaining a free scale absent a clock.

    A Cs input cancels 2*pi: E_GeV=h*nu_Cs*gap/(GeV_J*epsilon).
    """
    glo, ghi = interval(model_gap_interval, nonnegative=True)
    identity = _text(generator_id, 'target generator identity')
    velocity = exact(model_signal_speed)
    if velocity <= 0:
        raise ValueError('positive derived model signal speed required')
    hlo, hhi = hbar_gev_seconds_interval()
    row = {
        'generator_id': identity,
        'model_gap_interval': (glo, ghi),
        'model_evolution': 'exp(-i H_model tau)',
        'model_signal_speed': velocity,
        'matched_physical_signal_speed_m_s': LIGHT_SPEED_M_S,
        'free_scale_parameter': 'T=seconds_per_source_time>0',
        'length_relation': 'length_unit_m=c_SI*T/v_model',
        'energy_GeV_relation': 'E_GeV=(hbar_GeV_s*gap_model)/T',
        'energy_GeV_times_seconds_per_source_time_interval': (hlo*glo, hhi*ghi),
        'Cs_energy_GeV_relation': 'E_GeV=(h*nu_Cs/GeV_J)*gap_model/epsilon_Cs_model',
        'Cs_energy_numerator_GeV_interval': (
            CESIUM_TRANSITION_ENERGY_GEV*glo, CESIUM_TRANSITION_ENERGY_GEV*ghi),
        'physical_rest_mass_established': False,
        'absolute_scale_derived_by_this_conversion': False,
        'SI_source': SI_SOURCE,
    }
    if clock is None:
        row.update({
            'status': 'UNRESOLVED_ABSOLUTE_SCALE',
            'clock_input': None,
            'energy_GeV_interval': None,
            'seconds_per_source_time_interval': None,
            'energy_equivalent_kg_interval': None,
            'unique_physical_mass_GeV': None,
        })
        return _encode(row)
    if not isinstance(clock, CsClockInput):
        raise ValueError('a typed same-generator CsClockInput is required')
    if clock.generator_id != identity:
        raise ValueError('clock and excitation have different generator normalizations')
    elo, ehi = clock.angular_gap_interval
    plo, phi = pi_interval()
    time = elo/(2*phi*CESIUM_FREQUENCY_HZ), ehi/(2*plo*CESIUM_FREQUENCY_HZ)
    energy = (CESIUM_TRANSITION_ENERGY_GEV*glo/ehi,
              CESIUM_TRANSITION_ENERGY_GEV*ghi/elo)
    row.update({
        'status': 'CONDITIONAL_ON_SUPPLIED_SAME_GENERATOR_CS_IDENTIFICATION',
        'clock_input': {
            'generator_id': clock.generator_id,
            'angular_gap_interval': clock.angular_gap_interval,
            'physical_process': clock.physical_process,
            'source_gap_provenance': clock.source_gap_provenance,
            'physical_process_provenance': clock.physical_process_provenance,
            'same_time_provenance': clock.same_time_provenance,
            'provenance_truth_verified_by_arithmetic': False,
        },
        'energy_GeV_interval': energy,
        'seconds_per_source_time_interval': time,
        'length_unit_m_interval': tuple(LIGHT_SPEED_M_S*x/velocity for x in time),
        'energy_equivalent_kg_interval': tuple(
            x*GEV_IN_J/LIGHT_SPEED_M_S**2 for x in energy),
        'physical_clock_identity_derived_by_this_conversion': False,
        'unique_physical_mass_GeV': None,
    })
    return _encode(row)


def certificate():
    from model.radial_interacting_core.core import core_bounds
    radial = core_bounds()
    return _encode({
        'exact_SI_definitions': {
            'Planck_h_J_s': PLANCK_H_J_S,
            'elementary_charge_C': ELEMENTARY_CHARGE_C,
            'GeV_in_J': GEV_IN_J,
            'cesium_frequency_Hz': CESIUM_FREQUENCY_HZ,
            'light_speed_m_s': LIGHT_SPEED_M_S,
            'cesium_transition_energy_GeV': CESIUM_TRANSITION_ENERGY_GEV,
            'hbar_GeV_s_interval': hbar_gev_seconds_interval(),
            'pi_interval': pi_interval(),
            'source': SI_SOURCE,
        },
        'radial_charged_gap_scale_family': energy_conversion(
            (radial['actual_charged_gap_lower'], radial['actual_charged_gap_upper']),
            generator_id='bounded-radial-H:a=1,v=1,m=64,U=1/65536,radius=1',
            model_signal_speed=1),
        'radial_exterior_one_body_rest_energy_scale_family': energy_conversion(
            (64, 64), generator_id='bounded-radial-one-body-H:v=1,m=64,continuum-exterior-limit',
            model_signal_speed=1),
        'radial_exterior_rest_energy_domain': 'one-body exterior packets in the radial continuum scaling; no interacting charged-pole dispersion',
        'localized_charged_gap_is_relativistic_rest_energy': False,
        'c_matching': 'v_model=1 is identified with c_SI; length_unit_m=c_SI*T',
        'scale_nonidentifiability': {
            'domain': 'lambda>0, T>0',
            'transformation': 'T->lambda*T; length_unit->lambda*length_unit; energy_unit->energy_unit/lambda',
            'dimensionless_Hamiltonian_and_observer_histories_unchanged': True,
            'matched_physical_signal_speed_unchanged': True,
            'absolute_GeV_energies_unchanged': False,
        },
        'accepted_same_generator_physical_clock_inputs': [],
        'host_runtime_identified_with_source_time': False,
        'G_or_observed_particle_mass_used_to_choose_a_clock': False,
        'actual_absolute_particle_mass_prediction': False,
    })
