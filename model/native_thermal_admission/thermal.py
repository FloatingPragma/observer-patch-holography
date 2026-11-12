"""Exact recurrent classes and all-positive-temperature stationarity obstruction."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from model.native_abundance_epochs import model as native
from model.particle_common_contract.provenance import import_closure

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def encode(value):
    if isinstance(value, F): return str(value)
    raise TypeError(type(value).__name__)


def digest(value):
    return hashlib.sha256(json.dumps(value, default=encode, sort_keys=True,
                                    separators=(',', ':')).encode()).hexdigest()


def source_bindings():
    paths = {str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py', '.md')}
    paths.update({'model/native_abundance_epochs/model.py', 'model/native_abundance_epochs/PROOF.md',
        'model/native_abundance_epochs/results.json', 'model/defect_integer_lift/model.py',
        'model/defect_integer_lift/federation.json', 'Lean/Dynamics/ProtectedCharge.lean',
        'Lean/Dynamics/PublicMarkov.lean', 'Lean/Dynamics/ChargeConjugation.lean',
        'Lean/Dynamics/StationaryLyapunov.lean'})
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def strongly_connected_components(matrix):
    """Deterministic iterative Kosaraju on the exact positive-edge graph."""
    states = tuple(sorted(matrix))
    reverse = {s: [] for s in states}
    for s, row in matrix.items():
        for t, weight in row.items():
            if weight <= 0: raise ValueError('only strictly positive edges belong in the kernel')
            reverse[t].append(s)
    visited, order = set(), []
    for start in states:
        if start in visited: continue
        visited.add(start); stack = [(start, iter(sorted(matrix[start])))]
        while stack:
            node, adjacent = stack[-1]
            target = next(adjacent, None)
            if target is None:
                stack.pop(); order.append(node)
            elif target not in visited:
                visited.add(target); stack.append((target, iter(sorted(matrix[target]))))
    visited, components = set(), []
    for start in reversed(order):
        if start in visited: continue
        visited.add(start); stack = [start]; part = []
        while stack:
            node = stack.pop(); part.append(node)
            for target in reverse[node]:
                if target not in visited: visited.add(target); stack.append(target)
        components.append(tuple(sorted(part)))
    return tuple(sorted(components, key=lambda part: part[0]))


@lru_cache(maxsize=1)
def classification():
    matrix = native.kernel()
    components = strongly_connected_components(matrix)
    closed = tuple(part for part in components if all(t in part for s in part for t in matrix[s]))
    return dict(matrix=matrix, components=components, closed=closed,
                class_of={s: i for i, part in enumerate(closed) for s in part})


def stationary_push(distribution):
    result = defaultdict(F)
    for s, mass in distribution.items():
        for t, weight in native.kernel()[s].items(): result[t] += mass*weight
    return {s: v for s, v in result.items() if v}


def drift_obstruction(matrix, cost):
    if any(sum(row.values()) != 1 or any(p <= 0 for p in row.values()) for row in matrix.values()):
        raise ValueError('a normalized positive-edge kernel is required')
    drifts = {s: sum(p*(cost[t]-cost[s]) for t,p in row.items()) for s,row in matrix.items()}
    return all(v <= 0 for v in drifts.values()) and any(v < 0 for v in drifts.values())


def certificate():
    pins = source_bindings()
    source = classification(); matrix = source['matrix']; closed = source['closed']
    records = []
    for state in native.STATES:
        for seam in range(3):
            for coin in range(2):
                target, mark = native.event(state, seam, coin)
                records.append(dict(source=list(state), target=list(target), probability=F(1, 6),
                                    energy_change=native.energy(target)-native.energy(state), **mark))
    drift = {s: sum(p*(native.energy(t)-native.energy(s)) for t, p in row.items())
             for s, row in matrix.items()}
    recurrent = set(source['class_of'])
    classes = []
    for index, part in enumerate(closed):
        equilibrium = {s: F(1, len(part)) for s in part}
        if stationary_push(equilibrium) != equilibrium: raise ArithmeticError('class measure is not stationary')
        if any(equilibrium[s]*matrix[s].get(t, 0) != equilibrium[t]*matrix[t].get(s, 0)
               for s in part for t in part): raise ArithmeticError('class detailed balance fails')
        energies = {native.energy(s) for s in part}
        if len(energies) != 1: raise ArithmeticError('closed class is not cost-flat')
        if len(part) == 3 and any(matrix[s].get(t, 0) != (F(2,3) if s==t else F(1,6))
                                  for s in part for t in part):
            raise ArithmeticError('three-state recurrent mixing law differs')
        classes.append(dict(id=index, load=sum(part[0][:3]), states=[list(s) for s in part],
            state_probability=F(1, len(part)), energy=energies.pop(), face_class=part[0][3],
            transition_matrix=[[matrix[s].get(t, F(0)) for t in part] for s in part],
            exact_stationarity=True, detailed_balance=True,
            nonconstant_eigenvalues=[F(1, 2)]*(len(part)-1)))
    sectors = []
    for load in range(16):
        states = [s for s in native.STATES if sum(s[:3]) == load]
        strict = [(s, t) for s in states for t in matrix[s] if native.energy(t) < native.energy(s)]
        witness = None
        if strict:
            s, t = strict[0]
            witness = dict(source=list(s), target=list(t), source_energy=native.energy(s),
                target_energy=native.energy(t), forward_probability=matrix[s][t],
                reverse_probability=matrix[t].get(s, F(0)), source_drift=drift[s])
            if witness['reverse_probability'] != 0: raise ArithmeticError('strict descent has reverse edge')
        weights = {s: F(1, 2**native.energy(s)) for s in states}
        normalization = sum(weights.values()); gibbs = {s: v/normalization for s, v in weights.items()}
        pushed = stationary_push(gibbs)
        gibbs_drift = sum(gibbs[s]*drift[s] for s in states)
        residual = sum(abs(pushed.get(s, 0)-gibbs[s]) for s in states)/2
        admitted = not strict
        if drift_obstruction({s: matrix[s] for s in states}, {s: native.energy(s) for s in states}) != (not admitted):
            raise ArithmeticError('pointwise drift obstruction classification differs')
        if (gibbs_drift == 0) != admitted or (residual == 0) != admitted:
            raise ArithmeticError('exact Gibbs diagnostic contradicts drift theorem')
        selected = [row for row in classes if row['load'] == load]
        sectors.append(dict(total_load=load, state_count=len(states),
            recurrent_state_count=sum(len(c['states']) for c in selected),
            recurrent_class_ids=[c['id'] for c in selected], stationary_simplex_dimension=len(selected)-1,
            strict_descent_edge_count=len(strict), strict_descent_witness=witness,
            drift_min=min(drift[s] for s in states), drift_max=max(drift[s] for s in states),
            full_sector_all_positive_temperature_Gibbs_stationary=admitted,
            full_sector_all_positive_temperature_detailed_balance=admitted,
            temperature_selected_by_native_law=False,
            exact_beta_log2_control=dict(mean_cost_drift=gibbs_drift, half_l1_stationarity_defect=residual),
            recurrent_only_Gibbs_stationary=True,
            recurrent_only_class_weights='size(C)*exp(-energy(C)/T) divided by the recurrent partition sum; optional supplied choice',
            limiting_native_preparation_selects_these_weights=False))
    checks = dict(all_1296_rows_normalized=all(sum(row.values()) == 1 for row in matrix.values()),
        all_7776_original_records_retained=len(records) == 7776,
        all_edges_nonpositive_cost_change=all(r['energy_change'] <= 0 for r in records),
        pointwise_cost_drift_nonpositive=all(v <= 0 for v in drift.values()),
        closed_class_union_equals_native_terminal_set=recurrent == {s for s in native.STATES if native.extinct_or_frozen(s)},
        recurrent_class_count_46=len(closed) == 46,
        recurrent_state_count_66=len(recurrent) == 66,
        absorbing_singletons_36=sum(len(c) == 1 for c in closed) == 36,
        permutation_classes_10=sum(len(c) == 3 for c in closed) == 10,
        all_nontrivial_load_sectors_obstructed=all(not row['full_sector_all_positive_temperature_Gibbs_stationary']
            for row in sectors if 0 < row['total_load'] < 15),
        endpoint_load_sectors_are_identity=all(matrix[s] == {s: F(1)} for s in native.STATES if sum(s[:3]) in (0, 15)))
    if not all(checks.values()): raise ArithmeticError('native thermal admission check failed')
    if pins != source_bindings(): raise ArithmeticError('source changed during thermal admission')
    return dict(schema='oph.native_thermal_admission.v1', all_checks_pass=True,
        source_bindings=pins, source_state_reference_clock=dict(
            source='unchanged closed three-port/two-face native repair/register kernel',
            preparation='arbitrary distributions on the1296 states; full-support sector Gibbs is a tested comparison',
            clock='complete attempted native events including rejected/no-operation records',
            physical_temperature=None, physical_time=None),
        complete_attempt_records=records, complete_attempt_record_sha256=digest(records),
        recurrent_classes=classes, sectors=sectors, checks=checks,
        stationary_simplex=dict(extreme_measures=46, affine_dimension=45,
            statement='all stationary laws are arbitrary nonnegative normalized mixtures of the46 class-uniform measures',
            weights_selected_by_kernel=False),
        theorem=dict(domain='every strictly positive distribution on each conserved-load sector',
            Gibbs_family='exp(-E/T)/Z for every finite T>0, E=5*sum(load^2)+2*(C(k)+C(-k))',
            obstruction='stationarity forces zero expected cost drift; full support weights a strictly negative drift in sectors1..14',
            minimum_missing_reverse_moves='at least one positive reverse edge for every strict native cost-descent edge, or removal of full support; either changes the proposed thermal model',
            missing_physical_services=['bath energy exchange and detailed-balance ratios', 'temperature readout', 'physical energy and clock conversion']),
        acceptance_flags=dict(exact_finite_native_equilibrium_question_resolved=True,
            substituted_thermal_sampler_used=False, physical_thermal_abundances=False,
            expansion_or_cosmological_freezeout_derived=False))


def main():
    parser=argparse.ArgumentParser();m=parser.add_mutually_exclusive_group(required=True)
    m.add_argument('--write',action='store_true');m.add_argument('--check',action='store_true');args=parser.parse_args()
    result=json.loads(json.dumps(certificate(),default=encode,allow_nan=False))
    path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    elif json.loads(path.read_text()) != result:raise ArithmeticError('thermal admission receipt differs')
    print('PASS:46 stationary classes; all finite positive temperatures obstructed in load sectors1..14')


if __name__=='__main__':main()
