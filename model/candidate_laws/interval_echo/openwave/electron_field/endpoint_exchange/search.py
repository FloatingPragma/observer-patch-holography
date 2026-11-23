"""Exchange a pair of signed charges through a junction with resolved field endpoints."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
import json
from math import prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[6]))
from model.candidate_laws.interval_echo.openwave.electron_field.flux_string.search import Patch, phase
from model.executable_model import federation


def occupied(a, b):
    return tuple(sorted((a, b)))


class Exchange:
    def __init__(self, patch, sites, hops, places, contacts, sign):
        self.patch = patch
        self.sites = sites
        self.hops = hops
        self.places = places
        self.sign = sign
        self.paths = {(f, g): patch.path(f, g) for f, g in combinations(range(20), 2)}
        # Every pair in this assay has the same external reservoirs and fixed divergence.
        sector = defaultdict(set)
        sector[self.places[0]].add(patch.seed)
        todo = deque([(self.places[0], patch.seed)])

        def visit(place, word):
            if word not in sector[place]:
                sector[place].add(word)
                todo.append((place, word))

        while todo:
            place, word = todo.popleft()
            for _, other in patch.ring_moves(word):
                visit(place, other)
            for edge in contacts:
                owners = [x for x in edge if x in place]
                if len(owners) != 1:
                    continue
                source, = owners
                target, after = patch.move(source, word, edge, sign)
                if target != source:
                    spectator, = (x for x in place if x != source)
                    visit(occupied(target, spectator), after)
        assert set(self.places) <= set(sector)
        self.words = [sector[place] for place in self.places]
        self.counts = [len(words) for words in self.words]
        self.active_faces = [{f for word in words for f, _ in patch.ring_moves(word)}
                             for words in self.words]
        self.enabled = []
        self.toggled = 0
        for k, (source, target) in enumerate(self.hops):
            edge = occupied(source, target)
            self.toggled ^= patch.bit[edge]
            enabled = []
            for word in self.words[k]:
                pos, after = patch.move(source, word, edge, sign)
                if pos == target:
                    assert after in self.words[(k+1) % 6] and after == word ^ patch.bit[edge]
                    assert patch.move(target, after, edge, sign) == (source, word)
                    enabled.append(word)
            self.enabled.append(enabled)

    def mask(self, faces):
        return self.paths[occupied(*faces)]

    def resolved(self, side, old_mask, new_mask):
        """Integer numerator between two normalized parity-string states."""
        edge = self.patch.bit[occupied(*self.hops[side])]
        return phase(edge, new_mask) * self.character_sum(side, old_mask ^ new_mask)

    @lru_cache(maxsize=4096)
    def character_sum(self, side, mask):
        return sum(phase(word, mask) for word in self.enabled[side])

    def position_only(self):
        accepted = 0
        for initial in self.words[0]:
            word = initial
            for source, target in self.hops:
                pos, word = self.patch.move(source, word, occupied(source, target), self.sign)
                if pos != target:
                    break
            else:
                assert word == initial ^ self.toggled and word in self.words[0]
                accepted += 1
        return accepted

    def verify_projection(self, assignments, numerator):
        """Contract each accepted transfer directly in the words, including its adjoint."""
        at = dict(zip(self.sites, assignments))
        masks = [self.mask((at[x], at[y])) for x, y in self.places]
        values = []
        for k, (source, target) in enumerate(self.hops):
            edge = occupied(source, target)
            following = (k + 1) % 6
            forward = 0
            for word in self.words[k]:
                new_pos, moved = self.patch.move(source, word, edge, self.sign)
                if new_pos == target:
                    forward += phase(word, masks[k]) * phase(moved, masks[following])
            backward = 0
            for word in self.words[following]:
                new_pos, moved = self.patch.move(target, word, edge, self.sign)
                if new_pos == source:
                    backward += phase(word, masks[following]) * phase(moved, masks[k])
            assert forward == backward
            values.append(forward)
        assert prod(values) == numerator

    def scan(self):
        denom = prod(self.counts)
        data = {group: {'counts': Counter(), 'moving_vs_parked': Counter(), 'max': {}}
                for group in ('all', 'disjoint_at_each_read')}
        fixed_inputs = {group: defaultdict(dict) for group in data}
        qualified = 0
        zero_energy = 0
        distinct = 0
        incident = [[i for i, face in enumerate(self.patch.faces) if p in face] for p in self.sites]
        bare = prod(len(w) for w in self.enabled)
        for assignments in product(*incident):
            if len(set(assignments)) != len(self.sites):
                continue
            distinct += 1
            at = dict(zip(self.sites, assignments))
            endpoints = [(at[x], at[y]) for x, y in self.places]
            if any(f not in self.active_faces[k] for k, pair in enumerate(endpoints) for f in pair):
                zero_energy += 1
                continue
            qualified += 1
            masks = [self.mask(pair) for pair in endpoints]
            numerators = [self.resolved(k, masks[k], masks[(k+1) % 6]) for k in range(6)]
            value = prod(numerators)
            stationary_phase = phase(self.toggled, masks[0])
            assert prod(self.resolved(k, masks[0], masks[0]) for k in range(6)) == bare * stationary_phase
            outcome = 'negative' if value < 0 else 'positive' if value > 0 else 'zero'
            disjoint = all(set(self.patch.faces[f]).isdisjoint(self.patch.faces[g])
                           for f, g in endpoints)
            for group in ('all', 'disjoint_at_each_read') if disjoint else ('all',):
                result = data[group]
                result['counts'][outcome] += 1
                stationary_outcome = 'negative' if stationary_phase < 0 else 'positive'
                result['moving_vs_parked'][f'{outcome}_vs_{stationary_outcome}'] += 1
                initial_faces = tuple(sorted(endpoints[0]))
                records = fixed_inputs[group][initial_faces]
                if outcome not in records or abs(value) > abs(records[outcome][0]):
                    records[outcome] = (value, assignments)
                best = result['max'].get(outcome)
                if best is None or abs(value) > abs(best['signed_numerator']):
                    result['max'][outcome] = {
                        'face_ids_in_site_order': assignments, 'per_hop_numerators': numerators,
                        'signed_numerator': value, 'amplitude': float(Fraction(value, denom)),
                        'six_read_success': float(Fraction(value, denom)**2),
                        'parked_field_phase_with_position_checks': stationary_phase,
                    }
        assert distinct == zero_energy + qualified
        for group, result in data.items():
            for row in result['max'].values():
                self.verify_projection(row['face_ids_in_site_order'], row['signed_numerator'])
            paired = [(key, outcomes) for key, outcomes in fixed_inputs[group].items()
                      if 'negative' in outcomes and 'positive' in outcomes]
            result['initial_face_pairs_with_both_signs'] = len(paired)
            if paired:
                key, outcomes = max(paired, key=lambda entry: min(abs(entry[1][s][0])
                                                                   for s in ('negative', 'positive')))
                examples = {'initial_faces': key}
                for sign in ('negative', 'positive'):
                    value, assignments = outcomes[sign]
                    self.verify_projection(assignments, value)
                    amplitude = Fraction(value, denom)
                    examples[sign] = {'faces': assignments, 'amplitude': float(amplitude),
                                      'six_read_success': float(amplitude**2)}
                result['same_initial_state_examples'] = examples
        self.character_sum.cache_clear()
        accepted = self.position_only()
        return {'sign': self.sign, 'sites_in_assignment_order': self.sites, 'hops': self.hops,
                'starting_ports': self.places[0], 'words_per_position': self.counts,
                'normalization_denominator': denom, 'assignments_considered': prod(map(len, incident)),
                'distinct_face_assignments': distinct, 'distinct_but_missing_endpoint_energy': zero_energy,
                'resolved': data, 'uniform_or_stationary_endpoint_success': float(Fraction(bare, denom)**2),
                'position_only_accepted_words': accepted, 'position_only_exchange_probability': accepted / self.counts[0],
                'position_only_xor_mask_with_unit_transition_amplitude': self.toggled}


def run():
    patch = Patch()
    gateways = {p: q for p, q in federation(0).external if p < 12}
    results = []
    for centre, neighbours in enumerate(patch.cell.neighbours):
        for a, b in combinations(sorted(set(neighbours) & gateways.keys()), 2):
            for spare in sorted(set(neighbours) - {a, b}):
                for sign in (-1, 1):
                    sites = (a, b, centre, spare)
                    hops = ((a, centre), (centre, spare), (b, centre),
                            (centre, a), (spare, centre), (centre, b))
                    places = (occupied(a, b), occupied(centre, b), occupied(spare, b),
                              occupied(spare, centre), occupied(spare, a), occupied(centre, a))
                    contacts = [occupied(centre, arm) for arm in (a, b, spare)]
                    results.append(Exchange(patch, sites, hops, places, contacts, sign).scan())
    return {'initial_head_word': patch.seed, 'face_ports': patch.faces,
            'external_compensator_by_gateway': gateways, 'fixtures': results}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
