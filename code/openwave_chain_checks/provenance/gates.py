"""Typed source, preparation, clock, observable and error-domain joins."""
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import marshal
import math
from types import FunctionType
import numpy as np
from strict_facade import density


@dataclass(frozen=True)
class Source:
    identity: str
    carrier: str
    dimension: int
    regulator: str
    generator: str


@dataclass(frozen=True)
class Clock:
    law: str
    duration: Fraction
    rate: Fraction
    idle_rate: Fraction
    units: str
    source_paused_during_apparatus: bool


@dataclass(frozen=True)
class Endpoint:
    source: Source
    preparation: str
    clock: Clock | None
    observable: str
    complete_counts: tuple
    complete_failures: bool


@dataclass(frozen=True)
class ActualQuantum:
    source_h: str
    preparation: str
    dimension: int
    clock: Clock
    evolved_state: str


def array_hash(array, dtype):
    return hashlib.sha256(np.asarray(array, dtype=dtype).tobytes()).hexdigest()


def actual_quantum(raw, preparation, count_states=None, *, duration=Fraction(37,100),
                   rate=Fraction(4,5), uniform_rate=Fraction(8,5)):
    density(preparation, raw.n, normalized=True, sectors=raw.sectors)
    actual, _ = raw.preparation()
    if not np.array_equal(preparation, actual):
        raise ValueError('actual complete reservoir preparation required')
    if raw.law != 'MIX' or (duration,rate,uniform_rate)!=(Fraction(37,100),Fraction(4,5),Fraction(8,5)):
        raise ValueError('actual declared complete MIX common-clock parameters required')
    if count_states is None:
        count_states=raw.count_split(preparation,float(duration),float(rate))
    if len(count_states)!=3:
        raise ValueError('all three actual count states required')
    for state in count_states:density(state,raw.n,sectors=raw.sectors)
    total=sum(count_states)
    expected=raw.semigroup(preparation,float(duration),float(rate))
    weights=[float(np.trace(state).real) for state in count_states]
    mean=float(duration*rate)
    expected_weights=[math.exp(-mean),mean*math.exp(-mean),1-(1+mean)*math.exp(-mean)]
    if np.max(abs(total-expected))>2e-12 or np.max(abs(np.asarray(weights)-expected_weights))>2e-12:
        raise ValueError('actual complete common-clock output required')
    clock=Clock(raw.law,duration,rate,uniform_rate-rate,'supplied_inverse_C10_H',True)
    return ActualQuantum(array_hash(raw.h, '<f8'), array_hash(preparation, '<c16'), raw.n,
                         clock,array_hash(total,'<c16'))


def quantum_endpoint(contract, clock_contract):
    if contract.get('clock_id') != clock_contract.get('clock_id'):
        raise ValueError('particle and actual clock identifiers differ')
    if 'particle_instrument_id' in contract:
        if contract.get('count_classes') != ['N=0','N=1','N>=2 complete infinite-history complement']:
            raise ValueError('particle contract must retain the complete count complement')
        if set(contract.get('outcome_partition',{})) != {'first_detector_stops','second_detector_stops','conclusive_primary','fine_flags'}:
            raise ValueError('particle contract must retain all stop and fine-flag families')
    numerical = clock_contract['numerical_member']
    return Endpoint(
        Source(contract['source_id'], 'complex_density', contract['source_dimension'],
               'Z6_periodic', contract['source_H_array_sha256_le_float64']),
        contract['preparation_array_sha256_le_complex128'],
        Clock(numerical['law'], Fraction(str(numerical['duration'])),
              Fraction(str(numerical['native_attempt_rate'])),
              Fraction(str(numerical['uniformization_rate']))-Fraction(str(numerical['native_attempt_rate'])),
              'supplied_inverse_C10_H', True),
        'complete_quantum_terminal_state', ('N=0','N=1','N>=2'), True)


def join_quantum(left, right, witness):
    if not isinstance(witness, ActualQuantum):
        raise ValueError('labels alone do not supply an executed source/preparation witness')
    if left != right:
        raise ValueError('source/preparation/clock/observable/complement join mismatch')
    if left.clock != witness.clock:
        raise ValueError('regenerated clock mismatch')
    if left.source.carrier != 'complex_density' or left.source.regulator != 'Z6_periodic':
        raise ValueError('quantum384 Z6 source required')
    if (left.source.generator, left.preparation, left.source.dimension) != (witness.source_h,witness.preparation,witness.dimension):
        raise ValueError('regenerated source or preparation mismatch')
    if left.complete_counts != ('N=0','N=1','N>=2') or not left.complete_failures:
        raise ValueError('all count and apparatus complements required')
    return True


def field_bound(dynamic, modulus, tolerance, *, particle_modulus=None, claimed=None):
    if type(modulus) is not int or modulus < 6:
        raise ValueError('declared periodic integer modulus >= 6 required')
    if not isinstance(tolerance, Fraction) or not 0 < tolerance < 1:
        raise ValueError('nonvacuous rational error tolerance required')
    result = dynamic.bound(modulus, 5)
    if claimed is not None and claimed != result:
        raise ValueError('claimed bound differs from complete Poisson/Dyson tail')
    if result > tolerance:
        raise ValueError('complete-history error bound exceeds observable tolerance')
    if particle_modulus is not None and modulus != particle_modulus:
        raise ValueError('changed regulator has no transferred particle instrument')
    return result


def require_same_source(left, right):
    if left != right:
        raise ValueError('different native carriers/generators/regulators require an explicit transfer')
    return True


def require_trajectory_clock(endpoint):
    if endpoint.clock is None:
        raise ValueError('static field refinement supplies no trajectory clock')
    return True


def classical_witness(fixed, edges, native):
    """Bind domain-checked calls to a regenerated pinned native source."""
    actual_edges, _ = fixed.graph()
    canonical = fixed.literal_native()
    if list(edges) != list(actual_edges):
        raise ValueError('actual frozen source seam order required')
    if (not isinstance(native, FunctionType) or native.__closure__ is not None
            or marshal.dumps(native.__code__) != marshal.dumps(canonical.__code__)
            or set(native.__globals__) != set(canonical.__globals__)
            or native.__globals__.get('__builtins__') is not canonical.__globals__.get('__builtins__')
            or native.__globals__.get('native') is not native
            or native.__defaults__ != canonical.__defaults__
            or native.__kwdefaults__ != canonical.__kwdefaults__):
        raise ValueError('actual frozen literal native endpoint function required')
    return {'seams':len(edges), 'ports':96, 'joint_positions':9216,
            'edge_order_sha256':hashlib.sha256(np.asarray(edges,dtype='<i8').tobytes()).hexdigest(),
            'native_function':'Exact compiled frozen native AST, with the canonical globals namespace.',
            'source_regenerated':True}
