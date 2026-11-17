#!/usr/bin/env python3
"""One full native CCF-2 polynomial, spatial carrier and field descriptor.

This builds the actual joint source. No full 2**(96 V) state vector is
allocated. Matrix witnesses and operator-polynomial checks have different
scopes. All original calibration ancestry is retained.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
import sympy as sy
from . import core
from . import spatial

PARENT_DIGEST = '591949feb5120a050c49a8b4edf52db9a920bfaf0db920d8c7937089b4613013'
CANDIDATE = 'CCF-LT-2'


@lru_cache(None)
def parent_source():
    result = core.compiler.native_export()
    raw = json.dumps(result, sort_keys=True, separators=(',', ':')) + '\n'
    digest = hashlib.sha256(raw.encode()).hexdigest()
    assert digest == PARENT_DIGEST, 'parent snapshot differs from CCF-1 export'
    assert result['source_terms'] == 114 and result['mode_count'] == 192
    gauged = {(tuple((int(i), bool(c)) for i, c in row['word']), tuple(row['link_shift'])):
              sy.sympify(row['coefficient_exact_as_supplied']) for row in result['terms']}
    return core.forget_links(gauged), gauged, result['mode_charge_native_thirds']


def parent_parts(parent):
    """Split actual native-vacuum scalar equally, preserving energy zero.

    Pure partial vacuum restrictions each retain the entire scalar, so
    subtracting/replicating those naively changes the scalar per cell. Here
    the convention is explicit and recovers the original two-cell operator.
    """
    scalar = parent.get((), sy.Integer(0))
    nonconstant = {w: c for w, c in parent.items() if w}
    left = core.add(core.vacuum_restrict(nonconstant, range(96)), core.mono(coefficient=scalar/2))
    right = core.add(core.vacuum_restrict(nonconstant, range(96, 192)), core.mono(coefficient=scalar/2))
    translated = core.reindex(right, {i: i-96 for i in range(96, 192)})
    assert left == translated, 'parent local blocks differ; no homogeneous replication permitted'
    bond = core.add(parent, core.scale(left, -1), core.scale(right, -1))
    assert not bond.get((), 0)
    assert core.add(left, right, bond) == parent
    return left, bond


def replicated_parent(graph, parent):
    local, bond = parent_parts(parent)
    output = {}
    for v in range(graph.nvertices):
        p = core.reindex(local, {i: i+96*v for i in range(96)})
        output = core.gadd(output, core.lift(p, len(graph.edges)))
    coordinate_charges = core.native_charge_coordinates()
    for edge_index, (u, v) in enumerate(graph.edges):
        mapping = {i: i+96*u for i in range(96)}
        mapping.update({i: i-96+96*v for i in range(96, 192)})
        p = core.reindex(bond, mapping)
        dressed = {}
        for word, coefficient in p.items():
            charge_change_u = sum(coordinate_charges[(i % 48)//3]*(1 if creation else -1)
                                  for i, creation in word if i//96 == u)
            shift = [0]*len(graph.edges)
            shift[edge_index] = charge_change_u
            dressed[(word, tuple(shift))] = coefficient
        output = core.gadd(output, dressed)
    return output


@lru_cache(None)
def assemble(graph, *, include_dirac=True, include_interactions=True,
             speed=1, wilson_r=1, electric=1, magnetic=1):
    """Return all actual components of ONE finite source, with no exclusions."""
    from .dirac import dirac_terms
    for label, value in [('electric', electric), ('magnetic', magnetic)]:
        if sy.sympify(value).is_positive is not True:
            raise ValueError(label + ' field coefficient must be positive')
    parent, _, _ = parent_source()
    components = {'parent': replicated_parent(graph, parent)}
    if include_dirac:
        directions, spacings = graph.edge_geometry()
        components['dirac'] = dirac_terms(graph.nvertices, graph.edges, directions,
                                          spacings=spacings, c=speed, wilson_r=wilson_r)
    if include_interactions:
        additions = replicated_interactions(graph.nvertices)
        components.update({name: core.lift(p, len(graph.edges)) for name, p in additions.items()})
    combined = core.gadd(*components.values())
    B, C = graph.incidence()
    return {'candidate': CANDIDATE, 'graph': graph, 'components': components,
            'matter': combined,
            'field': {'electric': sy.sympify(electric), 'magnetic': sy.sympify(magnetic),
                      'B': B, 'C': C,
                      'Hamiltonian': 'electric/2 sum E_e^2 + magnetic sum (I-(W_f+W_f†)/2)'},
            'scope': 'Complete operator-polynomial source; no many-body spectrum or continuum claim.'}


@lru_cache(None)
def local_interactions():
    from .interactions import interaction_components
    return interaction_components(1)


def replicated_interactions(nvertices):
    return {name: core.add(*(core.reindex(poly, {i: i+96*v for i in range(96)})
                              for v in range(nvertices)))
            for name, poly in local_interactions().items()}


def describe_joint(source):
    """The complete local-capacity law executed by apply_total_column."""
    return {'matter_digest': core.polynomial_digest(source['matter']),
            'field_coefficients': {key: str(source['field'][key]) for key in ('electric', 'magnetic')},
            'local_energy_operator': 'k_v = g_v + N_v (h_v-e_v I) + e_v I',
            'geometry_operator': 'sum_v [A_v + M_v (k_v-F I) + F I]',
            'capacity_auxiliaries_per_vertex': 1,
            'default_capacity_dimension_per_vertex': 16,
            'stress_includes': list(source['components']) + ['electric', 'magnetic'],
            'source_energy_reference': 'e_v=0 and common F=0; scalar source terms retained',
            'no_separate_detector_energy_or_field_source': True,
            'full_source_state_vector_executed': False}


def apply_total_column(source, key, cap=6, capacity_config=None, local_config=None):
    """Complete H column, key=(capacity_tuple,(Phi_index,(occupation,flux)))."""
    from .geometry import GeometryConfig
    from .local_geometry import LocalGeometryConfig
    from .local_capacity import apply_native_column
    return apply_native_column(source, key, cap=cap,
        capacity_config=capacity_config or GeometryConfig(),
        local_config=local_config or LocalGeometryConfig())


def apply_total_phase_column(source, key, edge, other_edge=None, cap=6,
                             capacity_config=None, local_config=None, *, matter_only=True):
    """Exact current/contact sum_v M_v N_v times the source-shard derivative.

    Matter phase derivatives define the conserved charge current. Set
    matter_only=False to include magnetic loop response in the derivative
    of the entire Hamiltonian. Auxiliary energies and references are fixed.
    """
    from .geometry import GeometryConfig
    from .local_geometry import LocalGeometryConfig, _exact
    from .local_capacity import phi_labels, geometry_shard_column
    from .joint_geometry import auxiliary_columns
    capacity_config = capacity_config or GeometryConfig()
    local_config = local_config or LocalGeometryConfig()
    capacities, (local_index, (mask, flux)) = key
    n = source['graph'].nvertices
    labels = phi_labels(local_index, n, local_config)
    geometry_shard_column(local_index, 0, source['graph'], local_config)
    capacity_columns = auxiliary_columns(capacity_config)[1]
    if len(capacities) != n or any(not isinstance(i,int) or not 0 <= i < len(capacity_columns) for i in capacities):
        raise ValueError('one valid capacity index per vertex required')
    if not isinstance(cap, int) or cap < 0 or any(not isinstance(e,int) or abs(e)>cap for e in flux):
        raise ValueError('integer flux outside carrier')
    if not isinstance(mask,int) or mask<0 or mask.bit_length()>96*n or len(flux)!=len(source['graph'].edges):
        raise ValueError('native basis label outside carrier')
    output = defaultdict(lambda: sy.Integer(0))
    for vertex in range(n):
        lapse = 1+_exact(local_config.lapse_coupling)*_exact(local_config.phi_values[labels[vertex]])
        derivative = phase_derivative(source, edge, other_edge, vertex=vertex, matter_only=matter_only)
        for (word, shift), coefficient in derivative.items():
            acted = core.apply_word(word, mask)
            if acted is None:
                continue
            target_mask, sign = acted
            target_flux = tuple(x+y for x,y in zip(flux,shift))
            if any(abs(x)>cap for x in target_flux):
                continue
            for row, value in capacity_columns[capacities[vertex]]:
                changed = list(capacities); changed[vertex] = row
                target = (tuple(changed),(local_index,(target_mask,target_flux)))
                output[target] += coefficient*sign*lapse*value
    return core.clean(output)


def term_support(word, shift, edges):
    sites = {i//96 for i, _ in word}
    for amount, edge in zip(shift, edges):
        if amount:
            sites.update(edge)
    return sites


def energy_shards(source):
    """Hermitian gauge-invariant h_v summing H_S, an explicit stress choice."""
    graph = source['graph']
    shards = [defaultdict(lambda: sy.Integer(0)) for _ in range(graph.nvertices)]
    for (word, shift), coefficient in source['matter'].items():
        sites = term_support(word, shift, graph.edges) or set(range(graph.nvertices))
        for v in sites:
            shards[v][(word, shift)] += coefficient/sy.Integer(len(sites))
    return [core.clean(p) for p in shards]


def phase_derivative(source, edge, other_edge=None, *, matter_only=False, vertex=None):
    """Exact first/second external link-phase derivative at zero phase.

    U_e -> exp(i A_e) U_e is applied to every link-dependent term of H_S.
    Set matter_only=True for the matter current separately from loop response.
    Local allocation uses the same supports as the source energy shards.
    Electric energy and every source scalar have zero phase derivative.
    """
    graph = source['graph']
    selected = (edge,) if other_edge is None else (edge, other_edge)
    if any(not isinstance(e, int) or not 0 <= e < len(graph.edges) for e in selected):
        raise ValueError('phase derivative edge outside carrier')
    if vertex is not None and (not isinstance(vertex, int) or not 0 <= vertex < graph.nvertices):
        raise ValueError('energy shard vertex outside carrier')
    polynomial = dict(source['matter'])
    if not matter_only:
        for row in source['field']['C']:
            for sign in (-1, 1):
                shift = tuple(sign*int(s) for s in row)
                key = ((), shift)
                polynomial[key] = polynomial.get(key, 0)-source['field']['magnetic']/2
    output = {}
    for (word, shift), coefficient in polynomial.items():
        factor = sy.I*shift[edge] if other_edge is None else -shift[edge]*shift[other_edge]
        if not factor:
            continue
        if vertex is not None:
            sites = term_support(word, shift, graph.edges)
            if vertex not in sites:
                continue
            factor /= sy.Integer(len(sites))
        output[(word, shift)] = coefficient*factor
    return core.clean(output)


def phase_response_audit(source):
    """Check complete first/second derivatives, continuity and local sums."""
    graph = source['graph']
    charges, _ = core.mode_metadata(graph.nvertices)
    currents = [phase_derivative(source, e) for e in range(len(graph.edges))]
    result = {'first_derivative_terms': [], 'second_derivative_terms': []}
    for e, current in enumerate(currents):
        core.check_gauged(current, graph.nvertices, graph.edges)
        result['first_derivative_terms'].append(len(current))
        assert core.add(*(phase_derivative(source, e, vertex=v)
                          for v in range(graph.nvertices))) == current
        for f in range(e, len(graph.edges)):
            contact = phase_derivative(source, e, f)
            assert contact == phase_derivative(source, f, e)
            core.check_gauged(contact, graph.nvertices, graph.edges)
            result['second_derivative_terms'].append({'edges': [e, f], 'terms': len(contact)})
    for v in range(graph.nvertices):
        divergence = core.add(*(core.scale(j, int(source['field']['B'][v,e]))
                                for e, j in enumerate(currents)))
        time_charge = {}
        for (word, shift), coefficient in source['matter'].items():
            change = sum(charges[i]*(1 if creation else -1)
                         for i, creation in word if i//96 == v)
            if change:
                time_charge[(word,shift)] = -sy.I*change*coefficient
        assert core.add(time_charge, divergence) == {}
    result.update({'all_passed': True, 'exact_total_phase_current': True,
                   'exact_second_phase_contact': True,
                   'local_current_shards_sum': True,
                   'exact_native_charge_continuity': True,
                   'matter_current_and_loop_response_separately_available': True})
    return result


def apply_source_column(source, occupation, flux, cap=6, vertex=None):
    """Execute a literal full-source basis column without Fock truncation.

    All 96*V native modes are represented by the integer bit mask. The
    source is the explicitly hard-wall finite member. No reachable matter
    occupation is removed, and failed boundary shifts act as zero as the
    specified hard-wall law requires.
    """
    graph = source['graph']
    if not isinstance(cap, int) or cap < 0:
        raise ValueError('flux cap must be a nonnegative integer')
    if vertex is not None and (not isinstance(vertex, int) or not 0 <= vertex < graph.nvertices):
        raise ValueError('energy shard vertex outside carrier')
    if not isinstance(occupation, int) or occupation < 0 or occupation.bit_length() > 96*graph.nvertices:
        raise ValueError('occupation outside native Fock carrier')
    if len(flux) != len(graph.edges) or any(not isinstance(e, int) or abs(e) > cap for e in flux):
        raise ValueError('flux outside declared finite carrier')
    out = defaultdict(lambda: sy.Integer(0))
    for (word, shift), coefficient in source['matter'].items():
        if vertex is not None:
            sites = term_support(word, shift, graph.edges) or set(range(graph.nvertices))
            if vertex not in sites:
                continue
            coefficient /= sy.Integer(len(sites))
        applied = core.apply_word(word, occupation)
        if applied is None:
            continue
        destination_flux = tuple(e+s for e, s in zip(flux, shift))
        if any(abs(e) > cap for e in destination_flux):
            continue
        destination, sign = applied
        out[(destination, destination_flux)] += sign*coefficient
    electric = source['field']['electric']
    magnetic = source['field']['magnetic']
    electric_sum = sum(e*e*(sy.Rational(1, 2) if vertex is not None else 1)
                       for e, edge in zip(flux, graph.edges) if vertex is None or vertex in edge)
    out[(occupation, tuple(flux))] += electric*electric_sum/2
    for face, row in zip(graph.faces, source['field']['C']):
        if vertex is not None and vertex not in face:
            continue
        weight = sy.Rational(1, len(face)) if vertex is not None else sy.Integer(1)
        out[(occupation, tuple(flux))] += weight*magnetic
        for sign in (-1, 1):
            destination_flux = tuple(e+sign*int(s) for e, s in zip(flux, row))
            if all(abs(e) <= cap for e in destination_flux):
                out[(occupation, destination_flux)] -= weight*magnetic/2
    return core.clean(out)


def column_audit(source):
    graph = source['graph']
    B = source['field']['B']
    charges, _ = core.mode_metadata(graph.nvertices)
    columns = []
    for name, initial in [('native_empty', 0), ('native_neutral_excitation', 1 << core.index(0, 0, 10, 0))]:
        flux = (0,)*len(graph.edges)
        action = apply_source_column(source, initial, flux)
        for (mask, field), coefficient in action.items():
            local = [0]*graph.nvertices
            for i in range(mask.bit_length()):
                if mask >> i & 1:
                    local[i//96] += charges[i]
            assert [sum(int(B[v, e])*field[e] for e in range(len(field)))
                    for v in range(graph.nvertices)] == local
        mean = sy.expand(action.get((initial, flux), 0))
        variance = sy.expand(sum(sy.conjugate(a)*a for a in action.values()) - mean*mean)
        assert sy.im(mean) == 0 and float(sy.re(variance)) >= -1e-12
        columns.append({'preparation': name, 'nonzero_complete_source_outputs': len(action),
                        'mean_source_energy_exact': str(mean), 'energy_variance_exact': str(variance),
                        'every_output_gauss_consistent': True,
                        'no_matter_occupation_projection': True})
    assert columns[0]['mean_source_energy_exact'] != columns[1]['mean_source_energy_exact']
    empty_action = apply_source_column(source, 0, (0,)*len(graph.edges))
    shard_action = core.add(*(apply_source_column(source, 0, (0,)*len(graph.edges), vertex=v)
                              for v in range(graph.nvertices)))
    assert shard_action == empty_action
    # A complete finite-source short-time theorem, not a sector truncation.
    cap = 6
    norm_bound = sum(int(sy.ceiling(sy.Abs(sy.re(c))+sy.Abs(sy.im(c))))
                     for c in source['matter'].values())
    norm_bound += int(sy.ceiling(source['field']['electric']*len(graph.edges)*cap**2/2
                                + 2*source['field']['magnetic']*len(graph.faces)))
    candidates = [(key, value) for key, value in empty_action.items()
                  if any(key[1]) and key[0].bit_count() == 2
                  and all((i % 48)//3 in (9, 11) for i in range(key[0].bit_length()) if key[0] >> i & 1)]
    assert candidates, 'no actual joint-source charged transition'
    (mask, field), amplitude = candidates[0]
    reduced = dict(source)
    reduced['matter'] = core.gadd(source['matter'], core.scale(source['components']['dirac'], -1))
    assert apply_source_column(reduced, 0, (0,)*len(graph.edges)).get((mask, field), 0) == 0
    # |a| has positive exact modulus; tau=|a|/(2 B^2) gives lower |A| >= 3 tau |a|/4.
    amplitude_squared = sy.simplify(sy.conjugate(amplitude)*amplitude)
    assert amplitude_squared.is_positive
    time_squared = amplitude_squared/(4*norm_bound**4)
    probability_lower = sy.Rational(9, 16)*time_squared*amplitude_squared
    return {'columns': columns, 'flux_cap': cap,
            'local_energy_shard_columns_sum_to_full_source': True,
            'joint_source_short_time_witness': {
                'native_final_occupation': str(mask), 'final_flux': list(field),
                'matrix_element': str(amplitude), 'source_operator_norm_upper_integer': norm_bound,
                'time_squared_exact': str(time_squared),
                'unconditional_probability_lower_exact': str(probability_lower),
                'dirac_deleted_same_matrix_element': '0',
                'meaning': 'Strict positive finite-source transition at a very small specified time; not a photon/decay-rate identification.'},
            'execution': 'literal full Hamiltonian columns, not time propagation',
            'basis_states_are_native_occupations_not_claimed_physical_particles': True,
            'neutral_energy_enters_same_source_operator_used_by_geometry': True}


def audit(graph):
    started = time.perf_counter()
    parent, original_gauged, charges = parent_source()
    original_rebuilt = replicated_parent(spatial.edge(), parent)
    assert original_rebuilt == original_gauged
    local, bond = parent_parts(parent)
    from .mass_check import parent_mass_compatibility
    mass_join = parent_mass_compatibility(local)
    result = assemble(graph)
    components = {name: core.check_gauged(p, graph.nvertices, graph.edges)
                  for name, p in result['components'].items()}
    joint = core.check_gauged(result['matter'], graph.nvertices, graph.edges)
    shards = energy_shards(result)
    shard_checks = [core.check_gauged(p, graph.nvertices, graph.edges) for p in shards]
    assert core.gadd(*shards) == result['matter']
    assert core.gadd(result['matter'], *[core.scale(p, -1) for name, p in result['components'].items()
                                      if name != 'parent']) == result['components']['parent']
    # Native charge/flux labels can represent the failure. Do not project it out.
    charged_key = next(k for k in result['matter'] if any(k[1]))
    mutation = dict(result['matter'])
    word, shift = charged_key
    value = mutation.pop(charged_key)
    bad_shift = list(shift)
    bad_shift[0] += 1
    mutant = {(word, tuple(bad_shift)): value}
    original_pair = core.gadd({charged_key: value}, core.gadjoint({charged_key: value}))
    mutation = core.gadd(result['matter'], core.scale(original_pair, -1), mutant, core.gadjoint(mutant))
    assert mutation == core.gadjoint(mutation)
    try:
        core.check_gauged(mutation, graph.nvertices, graph.edges)
    except AssertionError:
        wrong_link_rejected = True
    else:
        raise AssertionError('wrong-link mutation not rejected')
    B, C = graph.incidence()
    return {'all_passed': True, 'candidate': CANDIDATE,
            'graph': graph.name, 'vertices': graph.nvertices,
            'edges': len(graph.edges), 'faces': len(graph.faces),
            'native_CAR_modes': 96*graph.nvertices,
            'parent_export_digest': PARENT_DIGEST, 'parent_factors': 114,
            'parent_monomials': len(parent), 'parent_local_monomials': len(local),
            'parent_bond_monomials': len(bond),
            'exact_two_cell_parent_recovery': True,
            'remove_all_additions_recovers_replicated_parent': True,
            'boundary_of_boundary_zero': not bool((B @ C.T).any()),
            'components': components, 'joint': joint,
            'local_matter_energy_shards': shard_checks,
            'local_matter_shards_sum_to_full_polynomial': True,
            'actual_parent_mass_kinetic_join': mass_join,
            'complete_source_column_execution': column_audit(result),
            'phase_current_and_contact': phase_response_audit(result),
            'wrong_link_mutation_rejected': wrong_link_rejected,
            'combined_source': describe_joint(result),
            'elapsed_seconds': time.perf_counter()-started,
            'physical_rows_passed_by_this_check': [],
            'calibrations_and_free_parameters_retained': True}

