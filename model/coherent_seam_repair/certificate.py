"""Regenerate the finite observer, controller and apparatus certificate."""
from functools import lru_cache
import hashlib
from pathlib import Path
import sympy as s
from model.common_source import slow_geometry
from model.particle_common_contract.provenance import import_closure
from . import apparatus, coherence, committed, gauge_walk, holonomy, walk
from .common import require
from .program import Observer, PROGRAM_ID

ROOT = Path(__file__).resolve().parents[2]


def input_hashes():
    paths = set(Path(__file__).parent.glob('*.py'))
    paths = {p for p in paths if not p.name.startswith('test_')}
    paths.update(ROOT/p for p in (
        'model/observer/native.py', 'model/observer/coherent_repair.py',
        'model/common_source/slow_geometry.py',
        'model/particle_common_contract/provenance.py', 'model/requirements.txt',
        'Lean/Screen/PrimitivePortFrameQuotient.lean', 'Lean/Screen/PortFrameGram.lean',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/__init__.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/__init__.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/__init__.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    paths = {ROOT/name for name in import_closure({str(path.relative_to(ROOT)) for path in paths})}
    # On a case-insensitive filesystem an imported class named Session can
    # spuriously resolve as a module Session.py beside the real session.py.
    # Keep the actual directory-entry spelling so Linux replay has the same pins.
    paths = {path for path in paths if path.name in {p.name for p in path.parent.iterdir()}}
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}


def execution_check():
    observer = Observer()
    initial = observer.embedded_port(0)
    once = observer.step(initial)
    twice = observer.step(once)
    require(initial.norm_squared() == once.norm_squared() == twice.norm_squared() == 1,
            'public observer preserves normalization')
    require(observer.inverse_step(once) == initial, 'actual inverse walk')
    require(observer.source_populations(twice)[0] == s.Rational(643, 900),
            'public step matches exact coherent receipt')
    phased, record = observer.quarter_turn(initial, observer.index(0, 0, 0))
    require(phased.norm_squared() == 1 and record['kind'] == 'quarter_turn', 'optional phase')
    branches = observer.register_measurement(phased, observer.index(0, 0, 0))
    require(sum(st.norm_squared() for _, st in branches) == 1, 'complete optional detector')
    filters = observer.preparation_filter(initial)
    require(sum(st.norm_squared() for _, st in filters) == 1, 'complete preparation branches')
    central = tuple(slow_geometry.source().seams()[0])
    outputs, probabilities = [], []
    for edge in (None, central):
        gauge = Observer(fibre=5, central_edge=edge)
        state = gauge.embedded_port(central[0], internal=3)
        for _ in range(3):
            state = gauge.step(state)
        outputs.append(gauge.embedded_return_amplitude(state, central[0], internal=3))
        require(state.norm_squared() == 1, 'actual coupled walk normalization')
        detector = gauge.embedded_return_measurement(state, central[0], internal=3)
        probabilities.append(detector[0][1].norm_squared())
        require(sum(st.norm_squared() for _, st in detector) == 1, 'complete source-return instrument')
        require(probabilities[-1] == outputs[-1]**2, 'amplitude attached to complete detector')
    require(outputs[1]-outputs[0] == -s.Rational(1, 6750), 'executed central-link signal')
    return {'exact_callable_observer_checked': True,
            'idle_inverse_checked': True, 'complete_optional_probe_and_filter_checked': True,
            'scalar_two_step_return_probability': '643/900',
            'coupled_three_step_return_amplitudes': list(map(str, outputs)),
            'coupled_amplitude_difference': '-1/6750',
            'coupled_complete_return_probabilities': list(map(str, probabilities)),
            'coupled_probability_difference': str(probabilities[1]-probabilities[0]),
            'flux_experiment_preparation': 'A|central-edge endpoint, parity0, internal3>; separate from the stationary reference preparation',
            'all_return_detector_outcomes_retained': True}


@lru_cache(maxsize=1)
def certificate():
    scalar = walk.certificate()
    return {
        'schema': 'oph.coherent_seam_repair.v1',
        'program_id': PROGRAM_ID,
        'ancestry': {
            'native': ['pinned twelve-port carrier', 'reference endpoint population kernel'],
            'declared': ['complex amplitude and instrument framework', 'owned coherent repair record',
                         'uniform coherent seam-address reflection', 'heralded reference preparation',
                         'optional basis quarter-turn apparatus', 'five-level encoded fibre',
                         'static encoded link background and typed overlap controller'],
            'calibrated_from_particle_data': [],
            'committed_native_rule_replaced_in_place': False,
        },
        'input_sha256': input_hashes(),
        'committed_pump': {'finite_state_and_effects': committed.certificate(),
                           'coherence_obstruction': coherence.certificate()},
        'walk': scalar,
        'reference': dict(scalar['stationary_filter'],
                          success_probability=scalar['stationary_filter']['success_probability_from_maximally_mixed'],
                          algebra='reduced walk-register algebra; conditioning on preparation labels changes the state'),
        'apparatus': apparatus.certificate(),
        'holonomy': holonomy.certificate(),
        'gauged_walk': gauge_walk.certificate(),
        'program_execution': execution_check(),
        'scope': {
            'finite_reversible_repair': True,
            'same_walk_stationary_reference': True,
            'optional_complex_tomography_1440': True,
            'kinematic_gauge_coupling': True,
            'full_7200_dimensional_tomography': False,
            'global_gauge_premises_for_complete_apparatus': False,
            'physical_particle_identification': False,
            'dynamical_gauge_field': False,
            'continuum_limit': False,
            'GeV_mass_prediction': False,
        },
    }
