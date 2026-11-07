"""Exact tree reduction of the configured native quantum-link Hamiltonian.

The full carrier retains M15 matter, every sixty-state link and one declared
static dual reference. Gauss invariance selects its complete one-fermion
sector. No relativistic, chiral-continuum or physical-mass claim is made.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import math

import sympy as s

from ..native_quantum_gauss import algebra as native


@dataclass(frozen=True)
class Tree:
    """Root zero; each other vertex points to an earlier parent.

    A finite tree admits this numbering. The explicit order fixes the link
    tensor factors and orients every literal raising link child to parent.
    """
    parents: tuple[int, ...] = (-1, 0, 0, 1)

    def __post_init__(self):
        if (type(self.parents) is not tuple or not self.parents
                or self.parents[0] != -1 or type(self.parents[0]) is not int
                or any(type(p) is not int or not 0 <= p < v
                       for v, p in enumerate(self.parents) if v)):
            raise ValueError('a nonempty rooted tree with parent earlier than child is required')

    @property
    def size(self):
        return len(self.parents)

    @property
    def edges(self):
        return tuple((child, parent) for child, parent in enumerate(self.parents) if child)

    def path(self, vertex):
        if type(vertex) is not int or not 0 <= vertex < self.size:
            raise ValueError('an actual tree vertex is required')
        path = [vertex]
        while path[-1]:
            path.append(self.parents[path[-1]])
        return tuple(reversed(path))

    def distance(self, vertex):
        return len(self.path(vertex))-1


@lru_cache(maxsize=1)
def link_dictionary():
    """The full literal link basis, including all five native blocks."""
    records = [None]
    for charge, indices in native.native_blocks():
        records.extend((charge, a, b) for a in indices for b in indices)
    if len(records) != 60:
        raise ArithmeticError('the complete native link dimension changed')
    return tuple(records), {record: i for i, record in enumerate(records)}


def _block(charge):
    blocks = dict(native.native_blocks())
    if type(charge) is not str or charge not in blocks:
        raise ValueError('an exact native central-eigenvalue label is required')
    return blocks[charge]


@lru_cache(maxsize=1)
def representation_certificate():
    """Actual-current irreducibility and distinctness for the leaf proof."""
    generators = native.native_currents()
    blocks = native.native_blocks()
    rows = []
    for charge, indices in blocks:
        q, d = s.Rational(charge), len(indices)
        if q == 0:
            raise ArithmeticError('an uncharged nonvacuum link invalidates leaf pruning')
        restricted = tuple(g.extract(indices, indices) for g in generators)
        if restricted[-1] != s.I*q*s.eye(d):
            raise ArithmeticError('the actual central generator differs on a native block')
        equations = s.Matrix.vstack(*(s.kronecker_product(g, s.eye(d))
                                     -s.kronecker_product(s.eye(d), g.T) for g in restricted))
        rank = equations.rank()
        if rank != d*d-1:
            raise ArithmeticError('actual native block irreducibility failed')
        rows.append({'central_eigenvalue': charge, 'dimension': d,
                     'commutator_constraint_rank': rank, 'endomorphism_dimension': 1,
                     'nonvacuum_single_endpoint_invariant_dimension': 0})
    if len({charge for charge, _ in blocks}) != len(blocks):
        raise ArithmeticError('distinct native central charges are required')
    return {'actual_native_currents': len(generators), 'matter_dimension': 15,
            'link_dimension': len(link_dictionary()[0]), 'blocks': rows,
            'Hom_dimensions': [[int(a == b) for b, _ in blocks] for a, _ in blocks],
            'off_diagonal_Hom_proof': 'nonzero difference of the actual central eigenvalues',
            'diagonal_Hom_proof': 'exact rank of all twelve actual commutator equations'}


def string_state(tree, charge, vertex):
    """Unnormalized A_v; its squared norm is d**(distance(v)+1).

    Keys retain the matter vertex and full M15 label, the full label of the
    static dual reference, and one full sixty-state address on every edge.
    The root position is included. No external spin register is supplied
    by the owning native_quantum_gauss model.
    """
    indices, path = _block(charge), tree.path(vertex)
    _, address = link_dictionary()
    result = {}
    for labels in product(indices, repeat=len(path)):
        links = [0]*len(tree.edges)
        for position, child in enumerate(path[1:], 1):
            links[child-1] = address[charge, labels[position], labels[position-1]]
        result[vertex, labels[-1], labels[0], tuple(links)] = s.Integer(1)
    return result


def gauss_action(tree, state, generator, vertex, *, include_reference=True, link_sign=1):
    """Full native endpoint generators acting on sparse ambient states."""
    if generator.shape != (15, 15):
        raise ValueError('an actual fifteen-state generator is required')
    tree.path(vertex)
    records, address = link_dictionary()
    entries = tuple((i, j, x) for (i, j), x in generator.todok().items() if x)
    output = {}
    for (site, matter, reference, links), amplitude in state.items():
        for i, j, value in entries:
            if site == vertex and j == matter:
                native.add(output, (site, i, reference, links), amplitude*value)
            if vertex == 0 and include_reference and i == reference:
                native.add(output, (site, matter, j, links), -amplitude*value)
            for edge, (child, parent) in enumerate(tree.edges):
                record = records[links[edge]]
                if record is None:
                    continue
                charge, a, b = record
                if vertex == child and i == a:
                    changed = list(links)
                    changed[edge] = address[charge, j, b]
                    native.add(output, (site, matter, reference, tuple(changed)), -link_sign*amplitude*value)
                if vertex == parent and j == b:
                    changed = list(links)
                    changed[edge] = address[charge, a, i]
                    native.add(output, (site, matter, reference, tuple(changed)), link_sign*amplitude*value)
    return native.normalized(output)


def hamiltonian_action(tree, state, *, hopping=1, link_energy=1):
    """Same literal full-link action as native_quantum_gauss on every edge.

    Electric energy is epsilon times each nonvacuum link projector. This
    is the owning model's electric law; no Casimir weight is substituted.
    """
    hopping = native.exact_real(hopping)
    link_energy = native.exact_real(link_energy, nonnegative=True)
    records, address = link_dictionary()
    partners = {i: (charge, indices) for charge, indices in native.native_blocks() for i in indices}
    output = {}
    for (site, matter, reference, links), amplitude in state.items():
        key = site, matter, reference, links
        native.add(output, key, link_energy*sum(value != 0 for value in links)*amplitude)
        for edge, (child, parent) in enumerate(tree.edges):
            record = records[links[edge]]
            if site == parent and record is None:
                charge, indices = partners[matter]
                for new in indices:
                    changed = list(links)
                    changed[edge] = address[charge, new, matter]
                    native.add(output, (child, new, reference, tuple(changed)), hopping*amplitude)
            elif site == child and record is not None and matter == record[1]:
                changed = list(links)
                changed[edge] = 0
                native.add(output, (parent, record[2], reference, tuple(changed)), hopping*amplitude)
    return native.normalized(output)


def reduced_model(tree, charge, *, hopping=1, link_energy=1):
    """Exact coordinates in A_v/sqrt(d), including the nontrivial Gram metric."""
    d = len(_block(charge))
    hopping = native.exact_real(hopping)
    link_energy = native.exact_real(link_energy, nonnegative=True)
    depths = tuple(tree.distance(v) for v in range(tree.size))
    gram = s.diag(*(d**depth for depth in depths))
    h = s.diag(*(link_energy*depth for depth in depths))
    normalized = h.copy()
    for child, parent in tree.edges:
        h[child, parent], h[parent, child] = hopping, d*hopping
        normalized[child, parent] = normalized[parent, child] = hopping*s.sqrt(d)
    if gram*h != h.H*gram:
        raise ArithmeticError('the exact reduced Hamiltonian is not Gram self-adjoint')
    return h, gram, normalized


def central_sector_census(tree=Tree()):
    """Exhaust every matter/link block assignment on the checked finite tree.

    The central equations alone force the reference block along the unique
    path. Exact non-Abelian Hom dimensions then determine all tensor states
    within those surviving assignments.
    """
    if len(tree.edges) > 5:
        raise ValueError('the explicit census is limited to five links; the leaf proof is general')
    blocks = native.native_blocks()
    # Twice the actual central eigenvalue is integral on all five blocks.
    charges = {name: int(2*s.Rational(name)) for name, _ in blocks}
    if any(s.Integer(value) != 2*s.Rational(name) for name, value in charges.items()):
        raise ArithmeticError('the exact integer central census convention changed')
    choices = (None, *charges)
    checked, admitted = 0, []
    for reference in charges:
        for site in range(tree.size):
            legal = []
            for matter in charges:
                for links in product(choices, repeat=len(tree.edges)):
                    weights = [0]*tree.size
                    weights[site] += charges[matter]
                    weights[0] -= charges[reference]
                    for (child, parent), link in zip(tree.edges, links, strict=True):
                        value = 0 if link is None else charges[link]
                        weights[child] -= value
                        weights[parent] += value
                    checked += 1
                    if not any(weights):
                        legal.append((matter, links))
            path = set(tree.path(site)[1:])
            expected = (reference, tuple(reference if child in path else None for child, _ in tree.edges))
            if legal != [expected]:
                raise ArithmeticError('the exhaustive central census found an omitted block sector')
            admitted.append({'reference': reference, 'matter_position': site,
                             'matter_block': reference, 'link_blocks': list(expected[1])})
    return {'complete_block_assignments_checked': checked,
            'admitted_assignments': admitted,
            'one_assignment_per_reference_block_and_position': True,
            'no_extra_off_path_flux_or_other_matter_blocks': True}


def one_edge_identity():
    """Compare every matter/link column with the owning source implementation.

    Hamiltonians act as identity on the reference, so its fixed label here
    represents every reference column and arbitrary untouched extensions.
    """
    tree = Tree((-1, 0))
    records, address = link_dictionary()
    blocks = tuple(indices for _, indices in native.native_blocks())

    def owning_key(key):
        site, matter, reference, links = key
        record = records[links[0]]
        a, b = (-1, -1) if record is None else record[1:]
        return 1-site, matter, reference, a, b

    columns = 0
    for site, matter, link in product(range(2), range(15), range(60)):
        key = site, matter, 0, (link,)
        derived = hamiltonian_action(tree, {key: s.Integer(1)})
        actual = native.hamiltonian_action({owning_key(key): s.Integer(1)}, 15, blocks=blocks)
        if {owning_key(k): value for k, value in derived.items()} != actual:
            raise ArithmeticError('the one-edge Hamiltonian differs from its owning source')
        columns += 1
    return {'all_matter_and_link_columns_compared': columns,
            'reference_factor_untouched_on_both_sides': True,
            'owning_action': 'model.native_quantum_gauss.algebra.hamiltonian_action',
            'electric_projector_and_literal_raising_amplitudes_identical': True}


def exact_tree_certificate(tree=Tree()):
    representation_certificate()
    rows, checks = [], 0
    for charge, indices in native.native_blocks():
        d = len(indices)
        states = [string_state(tree, charge, vertex) for vertex in range(tree.size)]
        h, gram, normalized = reduced_model(tree, charge)
        for vertex, state in enumerate(states):
            norm = sum(abs(value)**2 for value in state.values())
            if norm != d**(tree.distance(vertex)+1):
                raise ArithmeticError('the complete contracted-string norm failed')
            for generator in native.native_currents():
                for endpoint in range(tree.size):
                    if gauss_action(tree, state, generator, endpoint):
                        raise ArithmeticError('a literal native Gauss equation failed')
                    checks += 1
            expected = native.sum_states(*(native.scale(other, h[index, vertex])
                                           for index, other in enumerate(states)))
            if hamiltonian_action(tree, state) != expected:
                raise ArithmeticError('the same full Hamiltonian leaks from the complete string sector')
        # Different matter positions are orthogonal on the full ambient
        # carrier. Kernel completeness itself uses the exact leaf/Schur proof.
        if any(set(states[i]) & set(states[j]) for i in range(tree.size) for j in range(i)):
            raise ArithmeticError('different complete positions have overlapping ambient supports')
        rows.append({'central_eigenvalue': charge, 'dimension': d,
                     'complete_Gauss_kernel_dimension': tree.size,
                     'ambient_one_fermion_reference_link_dimension': tree.size*15*d*60**len(tree.edges),
                     'string_support_sizes': list(map(len, states)),
                     'relative_Gram_diagonal': list(map(int, gram.diagonal())),
                     'normalized_H': [[str(value) for value in normalized.row(i)] for i in range(tree.size)]})
    return {'parents': list(tree.parents), 'all_root_and_nonroot_positions_retained': True,
            'all_native_blocks_retained_before_Gauss_kernel': True,
            'all_sixty_link_states_retained_before_Gauss_kernel': True,
            'exact_Gauss_equations_checked': checks,
            'complete_kernel_proof': 'orthogonal block sectors, uncharged-leaf pruning, exact native Hom dimensions',
            'full_H_reduces_without_complement_leakage': True,
            'electric_law': 'epsilon*sum_edges(I-|0><0|); native_quantum_gauss projector law',
            'normalized_reduction': 'J*sqrt(d_r)*adjacency_tree + epsilon*distance_from_root',
            'blocks': rows}


def complete_history(tree, charge, *, hopping=1, link_energy=1, duration=s.Rational(1, 8), rounds=3):
    """All position records of the configured Cayley dynamics, including holds."""
    duration = native.exact_real(duration, nonnegative=True)
    if type(rounds) is not int or not 1 <= rounds <= 6 or tree.size > 8:
        raise ValueError('one through six rounds on at most eight vertices are required')
    h, gram, _ = reduced_model(tree, charge, hopping=hopping, link_energy=link_energy)
    eye = s.eye(tree.size)
    unitary = (eye+s.I*duration*h/2).inv()*(eye-s.I*duration*h/2)
    unitary = unitary.applyfunc(s.cancel)
    if (unitary.H*gram*unitary-gram).applyfunc(s.simplify) != s.zeros(tree.size):
        raise ArithmeticError('the configured Cayley step lost exact unitarity')
    if (unitary.H*gram*h*unitary-gram*h).applyfunc(s.simplify) != s.zeros(tree.size):
        raise ArithmeticError('the same-Hamiltonian Cayley step changed its energy')
    transition = s.Matrix(tree.size, tree.size,
                          lambda i, j: s.cancel(s.conjugate(unitary[i, j])*unitary[i, j]*gram[i, i]/gram[j, j]))
    if any(sum(transition[:, j]) != 1 for j in range(tree.size)):
        raise ArithmeticError('the complete position transition lost probability')
    outcomes = []
    for word in product(range(tree.size), repeat=rounds):
        before, probability = 0, s.Integer(1)
        for vertex in word:
            probability *= transition[vertex, before]
            before = vertex
        outcomes.append({'word': list(word), 'probability': str(s.cancel(probability))})
    if sum(s.Rational(row['probability']) for row in outcomes) != 1:
        raise ArithmeticError('the complete recorded law is not normalized')
    distribution = s.eye(tree.size)[:, 0]
    energies = s.Matrix([h[i, i] for i in range(tree.size)])
    ledger, total_work = [], s.Integer(0)
    for attempt in range(rounds):
        before = (energies.T*distribution)[0]
        distribution = transition*distribution
        after = (energies.T*distribution)[0]
        work = s.cancel(after-before)
        total_work += work
        ledger.append({'attempt': attempt+1, 'energy_before_readout': str(s.cancel(before)),
                       'energy_after_readout': str(s.cancel(after)), 'detector_work': str(work)})
    if s.cancel(total_work-(energies.T*distribution)[0]) != 0:
        raise ArithmeticError('complete detector work did not telescope')
    return {'duration': str(duration), 'rounds': rounds, 'initial_position': 0,
            'complete_position_alphabet': list(range(tree.size)), 'outcomes': outcomes,
            'same_position_holding_records_retained': True,
            'zero_probability_records_retained': True,
            'Cayley_step_exactly_preserves_same_H_energy': True,
            'energy_ledger': ledger, 'total_detector_work': str(s.cancel(total_work)),
            'detector_reservoir_or_reset_cost_derived': False}


def tail_bound(depth, initial_depth, degree, dimension, time, *, hopping=s.Integer(1)):
    """Rational majorant of the continuous-Hamiltonian finite-volume error.

    The reference root ball and position records are explicit inputs. The
    native declared Cayley evolution has a separate finite-step error.
    """
    if any(type(x) is not int for x in (depth, initial_depth, degree, dimension)):
        raise ValueError('exact integer graph and representation bounds are required')
    if not 0 <= initial_depth <= depth or degree < 1 or dimension < 1:
        raise ValueError('valid root-ball, degree and representation bounds are required')
    time = native.exact_real(time, nonnegative=True)
    hopping = native.exact_real(hopping)
    if not time.is_Rational or not hopping.is_Rational:
        raise ValueError('rational time and hopping are required for this rational majorant')
    # ceil(sqrt(d)) is an exact upper bound, avoiding a rounded threshold.
    root = math.isqrt(dimension)
    root += root*root < dimension
    z, first = degree*abs(hopping)*root*time, depth-initial_depth+1
    if z >= first+1:
        return {'first_omitted_order': first, 'z_upper': str(z),
                'complete_record_total_variation_upper': '1'}
    bound = 2*z**first/s.factorial(first)/(1-z/(first+1))
    return {'first_omitted_order': first, 'z_upper': str(z),
            'complete_record_total_variation_upper': str(min(s.Integer(1), bound))}
