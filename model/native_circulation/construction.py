"""Exact current counters and a single-history Gaussian limit.

Native attempts index the history. The uniform scheduler, binary total-two
sector, stationary preparation and chosen cycle readout are declared.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp

from ..observer import native
from ..primitive_observer import runtime, transfer
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
VALUES = (-1, 0, 1)
SOURCE_PATHS = tuple(dict.fromkeys((
    'model/native_circulation/__init__.py', 'model/native_circulation/__main__.py',
    'model/native_circulation/construction.py', 'model/common_source/protocol.py',
    *transfer.SOURCE_PATHS,
    'Lean/Thermodynamics/GreenKubo.lean',
    'Lean/Thermodynamics/NativeRepairFourCycle.lean',
    'Lean/Screen/PrimitivePortFrameQuotient.lean',
    'Lean/Screen/PortFrameGram.lean',
)))


def bounded(value, name):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 2e-11:
        raise ArithmeticError(f'{name} failed numerical comparison')
    return value


@lru_cache(maxsize=1)
def source():
    kernel = transfer.native_counts()
    transfer.exact_audit(kernel)
    index = {state: i for i, state in enumerate(kernel.states)}
    edges = tuple(kernel.seams)
    # Apparatus designation uses only the source graph, never measurement data.
    triangle = next(t for t in combinations(range(12), 3)
                    if all(tuple(sorted(e)) in edges for e in combinations(t, 2)))
    oriented = tuple(zip(triangle, triangle[1:]+triangle[:1]))
    weights = tuple(sum(int((a,b)==pair)-int((b,a)==pair) for pair in oriented)
                    for a,b in edges)
    branches = []
    for e,seam in enumerate(edges):
        for coin in (0,1):
            for i,state in enumerate(kernel.states):
                target,record = native.native(state,seam,coin)
                if record != (coin,state[seam[0]]-target[seam[0]]):
                    raise ArithmeticError('native record grammar changed')
                current = record[1]
                branches.append((i,index[target],e,coin,current,weights[e]*current))
    return kernel, triangle, weights, tuple(branches)


def exact_geometry():
    kernel,triangle,weights,branches = source()
    b=sp.zeros(12,30)
    for e,(a,c) in enumerate(kernel.seams):
        b[a,e]=-1; b[c,e]=1
    lap=b*b.T; j=sp.ones(12)/12
    inverse=(lap+j).inv()-j
    projection=sp.eye(30)-b.T*inverse*b
    if b*projection != sp.zeros(12,30) or projection*projection != projection or projection.T != projection:
        raise ArithmeticError('cycle/gradient decomposition failed')
    if sp.trace(projection)!=19 or b*sp.Matrix(weights)!=sp.zeros(12,1):
        raise ArithmeticError('the source cycle sector changed')
    # A current changes loads by B j. Its conditional mean is -B^T x/60.
    conditional=[[F(0) for _ in range(30)] for _ in kernel.states]
    second=[[F(0) for _ in range(30)] for _ in range(30)]
    for i,out,e,coin,current,g in branches:
        conditional[i][e]+=F(current,60)
        second[e][e]+=F(current*current,3960)
        if tuple(kernel.states[out][v]-kernel.states[i][v] for v in range(12)) != tuple(int(b[v,e])*current for v in range(12)):
            raise ArithmeticError('current does not reproduce the actual load change')
    for i,state in enumerate(kernel.states):
        expected=-(b.T*sp.Matrix(state))/60
        if sp.Matrix(conditional[i])!=expected or projection*expected!=sp.zeros(30,1):
            raise ArithmeticError('projected currents are not martingale differences')
    if sp.Matrix(second)!=sp.eye(30)/198:
        raise ArithmeticError('stationary current covariance changed')
    sigma=projection/198
    return {'chosen_triangle':list(triangle), 'oriented_edge_weights':list(weights),
            'incidence_rank':int(b.rank()), 'cycle_rank':19,
            'cycle_projection':[[str(x) for x in row] for row in projection.tolist()],
            'asymptotic_covariance':[[str(x) for x in row] for row in sigma.tolist()],
            'all_native_branch_continuity_checks':len(branches),
            'projected_conditional_mean_exactly_zero_on_all_states':True,
            'stationary_single_edge_second_moment':'1/198',
            'native_count_parameter':'integer uniform native attempts',
            'spatial_dimension_inferred_from_cycle_rank':False}


def axial_weights():
    """The recovered rank-three frame's axial seam readout, exactly over Q(sqrt5)."""
    kernel, _, _, _ = source()
    frame = [sp.Matrix([sp.Rational(a)+sp.Rational(b)*sp.sqrt(5) for a,b in row])
             for row in transfer.slow_geometry.raw_frame()]
    radius_squared = sp.simplify(frame[0].dot(frame[0]))
    retained_gram = transfer.slow_geometry.source().intrinsic_gram_exact()
    for i in range(12):
        for j in range(12):
            a,b = retained_gram[i][j]
            if sp.simplify(frame[i].dot(frame[j])/radius_squared-a-b*sp.sqrt(5)) != 0:
                raise ArithmeticError('axial frame differs from the recovered rank-three Gram')
    weights = sp.zeros(3, 30)
    for e,(a,b) in enumerate(kernel.seams):
        weights[:,e] = frame[a].cross(frame[b])/radius_squared
    return sp.ImmutableMatrix(weights.applyfunc(sp.simplify))


def axial_geometry():
    """Audit every native increment and a closed-state, nonzero-circulation cycle."""
    kernel, triangle, _, branches = source()
    weights = axial_weights()
    incidence = sp.zeros(12, 30)
    conditional = sp.zeros(66, 30)
    second = [0]*30
    for e,(a,b) in enumerate(kernel.seams):
        incidence[a,e] = -1
        incidence[b,e] = 1
    for i,_,e,_,current,_ in branches:
        conditional[i,e] += current
        second[e] += current*current
    if ((weights*incidence.T).applyfunc(sp.simplify) != sp.zeros(3,12)
            or (weights*conditional.T).applyfunc(sp.simplify) != sp.zeros(3,66)):
        raise ArithmeticError('axial current is not a native martingale difference')
    if (weights*weights.T).applyfunc(sp.simplify) != 8*sp.eye(3):
        raise ArithmeticError('axial frame covariance is not isotropic')
    covariance = (weights*sp.diag(*second)*weights.T/sp.Integer(len(branches))).applyfunc(sp.simplify)
    if covariance != sp.Rational(4,99)*sp.eye(3):
        raise ArithmeticError('native axial current covariance changed')
    # One token goes around an actual face while the other stays outside it.
    # Returning the carrier state does not reset the accessible counter.
    other = next(p for p in range(12) if p not in triangle)
    initial = tuple(int(p in (triangle[0],other)) for p in range(12))
    state = initial
    total = sp.zeros(3,1)
    events = []
    for a,b in zip(triangle, triangle[1:]+triangle[:1]):
        seam = tuple(sorted((a,b)))
        e = kernel.seams.index(seam)
        desired = tuple(state[p]-int(p==a)+int(p==b) for p in range(12))
        candidates = [(coin,*native.native(state,seam,coin)) for coin in (0,1)]
        coin,target,record = next(row for row in candidates if row[1] == desired)
        total += weights[:,e]*record[1]
        events.append({'seam':list(seam),'coin':coin,'current':record[1]})
        state = target
    total = total.applyfunc(sp.simplify)
    if state != initial or total == sp.zeros(3,1):
        raise ArithmeticError('closed native cycle lost its nonzero axial readout')
    return {'rank':3, 'native_branches_checked':len(branches),
            'frame_gram_entries_checked':144,
            'weights_by_oriented_seam':[[str(x) for x in weights[:,e]] for e in range(30)],
            'cycle_orthogonality_exact':True,'all_conditional_drifts_zero':True,
            'weight_gram_diagonal':'8','covariance_per_attempt_diagonal':'4/99',
            'off_diagonal_covariances':'0',
            'closed_cycle':{'initial_and_final_loads':list(initial),'events':events,
                            'axial_counter_increment':[str(x) for x in total],
                            'increment_norm_squared':str(sp.simplify(total.dot(total))),
                            'repeatable_without_resetting_carrier':True},
            'interpretation':'axial circulation in the recovered frame; no spatial translation or bulk inferred'}


def axial_radial_limit_control():
    """Finite hypotheses for the same counter's all-initial-law radial limit.

    Initial laws are classical distributions on the 66 binary configurations.
    Native attempt count supplies the scaling parameter. No assertion about
    the convergence of quantum terminal/reference states is made here.
    """
    kernel, _, _, branches = source()
    weights = axial_weights()
    size = len(kernel.states)
    counts = [[0]*size for _ in range(size)]
    means = sp.zeros(size, 30)
    squares = [[0]*30 for _ in range(size)]
    for i, j, edge, _, current, _ in branches:
        counts[i][j] += 1
        means[i, edge] += current
        squares[i][edge] += current*current
    if (tuple(tuple(row) for row in counts) != kernel.counts
            or any(sum(row) != 60 for row in counts)
            or any(sum(counts[i][j] for i in range(size)) != 60 for j in range(size))):
        raise ArithmeticError('the complete native kernel lost its uniform invariant law')
    if (weights*means.T).applyfunc(sp.simplify) != sp.zeros(3, size):
        raise ArithmeticError('the axial counter has a conditional drift')
    diameters = []
    for start in range(size):
        distance, frontier = {start: 0}, {start}
        depth = 0
        while frontier:
            next_frontier = {j for i in frontier for j, count in enumerate(counts[i])
                             if count and j not in distance}
            depth += 1
            for j in next_frontier:
                distance[j] = depth
            frontier = next_frontier
        if len(distance) != size:
            raise ArithmeticError('the native binary configuration chain is reducible')
        diameters.append(max(distance.values()))
    if any(counts[i][i] <= 0 for i in range(size)):
        raise ArithmeticError('a native self-loop is required for this aperiodicity certificate')
    jump_norms = tuple(sp.simplify(weights[:, e].dot(weights[:, e])) for e in range(30))
    if set(jump_norms) != {sp.Rational(4, 5)}:
        raise ArithmeticError('the native axial jump bound changed')
    stationary_covariance = (weights*sp.diag(*(sum(row[e] for row in squares)
        for e in range(30)))*weights.T/sp.Integer(60*size)).applyfunc(sp.simplify)
    if stationary_covariance != sp.Rational(4, 99)*sp.eye(3):
        raise ArithmeticError('the predictable covariance average changed')
    radial = {'adjacent': set(), 'nonadjacent': set()}
    for state, row in zip(kernel.states, squares, strict=True):
        occupied = tuple(i for i, value in enumerate(state) if value)
        kind = 'adjacent' if occupied in kernel.seams else 'nonadjacent'
        radial[kind].add(sum(jump_norms[e]*row[e] for e in range(30))/60)
    if radial != {'adjacent': {sp.Rational(8, 75)}, 'nonadjacent': {sp.Rational(2, 15)}}:
        raise ArithmeticError('the same-radius hidden-state control changed')
    # Exact radial generator identities on polynomial probes. The angular
    # term is tested on f(z)=z_3: dropping r^-2 leaves a nonzero drift.
    fourth = sp.simplify(sum(jump_norms[e]**2*current**4
                            for _, _, e, _, current, _ in branches)/(60*size))
    gaussian_fourth = sp.trace(stationary_covariance)**2+2*sp.trace(stationary_covariance**2)
    if fourth != sp.Rational(16, 165) or fourth == gaussian_fourth:
        raise ArithmeticError('the finite radial law falsely became Gaussian')
    diffusion = stationary_covariance[0, 0]/2
    radius = sp.Symbol('r', positive=True)
    radial_generator = lambda f: sp.simplify(diffusion*(sp.diff(f, radius, 2)
                                                       +2*sp.diff(f, radius)/radius))
    if (radial_generator(radius**2) != sp.Rational(4, 33)
            or radial_generator(radius**4) != sp.Rational(40, 99)*radius**2
            or sp.simplify(radial_generator(radius)-2*diffusion/radius) != 0):
        raise ArithmeticError('the three-dimensional radial generator changed')
    theta = sp.Symbol('theta', real=True)
    coordinate = radius*sp.cos(theta)
    angular_part = sp.simplify(sp.diff(coordinate, theta, 2)
                              +sp.cot(theta)*sp.diff(coordinate, theta))
    if sp.simplify(radial_generator(coordinate)+diffusion*angular_part/radius**2) != 0:
        raise ArithmeticError('the polar form fails the Cartesian linear coordinate test')
    unwarped = sp.simplify(radial_generator(coordinate)+diffusion*angular_part)
    if sp.simplify(unwarped-2*diffusion*(1/radius-radius)*sp.cos(theta)) != 0 or unwarped == 0:
        raise ArithmeticError('a constant-rate angular product passed as the native polar law')
    return {'initial_domain': 'every classical law on the 66 binary total-two configurations; counter starts at zero',
        'complete_native_branches_checked': len(branches),
        'irreducible_configuration_graph': True,
        'configuration_graph_diameter': max(diameters),
        'Doeblin_attempt_block': max(diameters),
        'Doeblin_uniform_reset_probability_lower': str(F(size, 60**max(diameters))),
        'minimum_self_loop_probability': str(F(min(counts[i][i] for i in range(size)), 60)),
        'uniform_stationary_law': '1/66',
        'conditional_drift_zero_at_every_configuration': True,
        'maximum_increment_squared_norm': '4/5',
        'ergodic_predictable_covariance_per_attempt': '(4/99) I3',
        'stationary_initial_preparation_required_for_functional_limit': False,
        'exact_finite_covariance_formula_requires_stationary_preparation': True,
        'same_zero_counter_one_attempt_radial_second_moment':
            {kind: str(next(iter(values))) for kind, values in radial.items()},
        'stationary_one_attempt_radial_second_moment': '4/33',
        'stationary_one_attempt_radial_fourth_moment': str(fourth),
        'Gaussian_one_attempt_radial_fourth_moment_control': str(gaussian_fourth),
        'scaling': 'Z_floor(N t)/sqrt(N); t is scaled native attempt count',
        'limit_heat_generator': '(2/99) Delta_R3',
        'radial_generator_for_r_positive': '(2/99) (d_r^2+2 d_r/r)',
        'polar_generator_for_r_positive': '(2/99) (d_r^2+2 d_r/r+Delta_S2/r^2)',
        'radial_volume_density_up_to_constant': 'r^2',
        'radial_squared_generator': '4/33',
        'radial_fourth_generator': '(40/99) r^2',
        'dropped_radial_Jacobian_squared_generator_control': '4/99',
        'unwarped_product_z3_drift_control': '(4/99) (1/r-r) cos(theta)',
        'radius_entrance_value': '0; angular direction is undefined there',
        'finite_radius_is_autonomous_for_every_initial_configuration': False,
        'physical_spatial_transport_or_independent_product_radius_identified': False,
        'new_Lean_limit_or_uniform_finite_count_error_bound': False}


def marked():
    kernel,_,_,branches=source()
    matrices=[[[0]*66 for _ in range(66)] for _ in VALUES]
    for i,j,_,_,_,g in branches:
        matrices[g+1][i][j]+=1
    if tuple(tuple(sum(matrices[k][i][j] for k in range(3)) for j in range(66)) for i in range(66))!=kernel.counts:
        raise ArithmeticError('marked current maps changed the native transition')
    for i in range(66):
        if sum(g*sum(matrices[g+1][i]) for g in VALUES)!=0:
            raise ArithmeticError('cycle current has a conditional drift')
    return tuple(tuple(tuple(row) for row in m) for m in matrices)


def instrument():
    _,_,weights,branches=source()
    groups={}
    for i,j,e,coin,current,g in branches:
        key=(e,coin,current)
        groups.setdefault(key,np.zeros((66,66),complex))[j,i]=1/math.sqrt(60)
    return runtime.Instrument(kraus=tuple(tuple(groups[key] for key in sorted(groups)
        if weights[key[0]]*key[2]==g) for g in VALUES),visible=True,
        resources=(('native_attempts',1),('circulation_reads',1)),
        binding={'candidate':'C','observable':'oriented triangular circulation from seam and native current',
                 'values':list(VALUES),'all_hidden_labels':'probability summed',
                 'native_channel_semantics':True,'native_pulse_words_emitted':False})


def horizon(value):
    if type(value) is not int or not 0<=value<=512:
        raise ValueError('integer attempt count from zero through512 required')
    return value


def count_laws(horizons=(1,4,16,64,256)):
    """Exact retained counter and terminal-state law; no trajectory cutoff."""
    horizons=tuple(horizon(n) for n in horizons)
    if not horizons or len(set(horizons))!=len(horizons):
        raise ValueError('distinct nonempty attempt counts required')
    m=marked()
    sparse=[[(i,j,count) for i,row in enumerate(matrix) for j,count in enumerate(row) if count]
            for matrix in m]
    rows=[[1]*66] # all66 stationary inputs, integer denominator66
    denominator=66
    saved={}
    for n in range(max(horizons)+1):
        if n in horizons:
            saved[n]=(tuple(tuple(row) for row in rows),denominator)
        if n==max(horizons): break
        updated=[[0]*66 for _ in range(len(rows)+2)]
        for g,transitions in enumerate(sparse):
            for s,row in enumerate(rows):
                out=updated[s+g]
                for i,j,c in transitions:
                    if row[i]: out[j]+=row[i]*c
        rows=updated; denominator*=60
        if sum(map(sum,rows))!=denominator:
            raise ArithmeticError('counter coarsening lost a history or terminal state')
    return saved


def runtime_check():
    operation=instrument()
    program=runtime.Program('one_observer_triangle_current',66,np.eye(66)/66,
        tuple(runtime.Event(f'cycle_read_{n}',instruments=(operation,)) for n in range(3)),
        binding={'stationary_preparation':'I66/66','parameter':'native attempts'})
    output=runtime.execute(program)
    matrices=marked(); residual=0.; records=[]
    aggregate=defaultdict(lambda:np.zeros((66,66),complex))
    for branch in output.branches:
        counts=[1]*66
        word=tuple(VALUES[k] for k in branch.record)
        for k in branch.record:
            counts=[sum(counts[i]*matrices[k][i][j] for i in range(66)) for j in range(66)]
        exact=[F(n,66*60**3) for n in counts]
        expected=np.diag(list(map(float,exact)))
        residual=max(residual,float(np.linalg.norm(branch.density-expected,'fro')))
        aggregate[sum(word)]+=branch.density
        if dict(branch.resources)!={'native_attempts':3,'circulation_reads':3}:
            raise ArithmeticError('recorded circulation resource accounting changed')
        records.append({'word':list(word),'probability':str(sum(exact))})
    law,denominator=count_laws((3,))[3]
    for s,row in enumerate(law):
        expected=np.diag([float(F(value,denominator)) for value in row])
        residual=max(residual,float(np.linalg.norm(aggregate[s-3]-expected,'fro')))
    prefix=runtime.execute(program,stop_after=1)
    restored=runtime.checkpoint_from_json(program,runtime.checkpoint_to_json(prefix))
    continued=runtime.execute(program,checkpoint=restored)
    return {'complete_records':records,'full_terminal_density_dimension':66,
            'maximum_complete_output_residual':bounded(residual,'native current output'),
            'checkpoint_complete_distance':bounded(runtime.complete_distance(output,continued),'checkpoint'),
            'exact_counter_coarsening_checked_against_all_runtime_records':True}


def run():
    geometry=exact_geometry(); laws=count_laws(); summaries=[]
    for n,(rows,denominator) in sorted(laws.items()):
        probabilities=[F(sum(row),denominator) for row in rows]
        moments=[sum((F(s-n)**k*p for s,p in enumerate(probabilities)),F(0)) for k in range(1,5)]
        if moments[0]!=0 or moments[1]!=F(n,66) or moments[2]!=0:
            raise ArithmeticError('native exact counter moments disagree with the martingale identity')
        variance=moments[1]
        empirical=sum(float(p)*math.cos((s-n)/math.sqrt(float(variance))) for s,p in enumerate(probabilities))
        summaries.append({'attempts':n,'complete_counter_law':[str(p) for p in probabilities],
            'score_min':-n,'score_max':n,'mean':str(moments[0]),'variance':str(variance),
            'fourth_moment':str(moments[3]),'excess_kurtosis':str(moments[3]/variance**2-3),
            'standardized_characteristic_at_one_diagnostic':empirical,
            'Gaussian_characteristic_at_one':math.exp(-.5),
            'terminal_counter_joint_integer_counts_sha256':hashlib.sha256(canonical(rows).encode()).hexdigest(),
            'all_scores_and_terminal_states_retained':True})
    m=marked()
    square_correlation=F(0)
    for i in range(66):
        for j in range(66):
            square_correlation+=F((m[0][i][j]+m[2][i][j])*(sum(m[0][j])+sum(m[2][j])),66*60**2)
    if square_correlation!=F(1,1980):
        raise ArithmeticError('the nonindependent squared-current control changed')
    return {'schema':'oph.native_circulation.v1','all_pass':True,
        'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'source_geometry':geometry,'axial_circulation':axial_geometry(),
        'axial_radial_limit':axial_radial_limit_control(),
        'complete_native_execution':runtime_check(),
        'one_history_counter_predictions':summaries,
        'temporal_dependence_control':{'successive_signed_covariance':'0',
            'successive_squared_product':str(square_correlation),'independent_squared_product':str(F(1,66)**2),
            'increments_independent':False},
        'limit':{'triangle_variance_per_native_attempt':'1/66',
            'vector_covariance':'(I-B^T L^+ B)/198','stationary_finite_ergodic_martingale_limit':True,
            'axial_three_vector_covariance':'(4/99) I3',
            'axial_functional_limit_from_every_binary_initial_law':True,
            'independent_replicas_required':False,'new_Lean_CLT':False,
            'uniform_Berry_Esseen_or_trajectory_error_certificate':False},
        'scope':{'stationary_binary_preparation_uniform_scheduler_and_cycle_readout_supplied':True,
            'physical_clock_spatial_dimension_or_interacting_QFT_derived':False,
            'native_quantum_C_E_selected':False,'new_empirical_agreement':False}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--write',type=Path);mode.add_argument('--check',type=Path)
    args=parser.parse_args(); result=run()
    if args.check:
        compare_receipt(strict_json.load(args.check),result);print('Native circulation receipt reproduced')
    elif args.write:
        with args.write.open('x',encoding='utf-8') as stream:stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else: print(canonical(result),end='')


if __name__=='__main__':main()
