"""Complete local capacity response with one auxiliary per vertex."""
from __future__ import annotations
import math
import sympy as sy

from .geometry import GeometryConfig
from .joint_geometry import auxiliary_columns
from .local_geometry import LocalGeometryConfig, _graph, _references, _exact

CANDIDATE = 'A-ccf2-local-capacity-v2'


def phi_labels(index, n, config):
    count = len(config.phi_values)
    if not isinstance(index, int) or not 0 <= index < count**n:
        raise ValueError('local geometry index outside carrier')
    labels = [0]*n
    for v in reversed(range(n)):
        index, labels[v] = divmod(index, count)
    return tuple(labels)


def phi_index(labels, count):
    answer = 0
    for value in labels:
        answer = count*answer+value
    return answer


def geometry_shard_column(index, vertex, graph, config=LocalGeometryConfig()):
    """Complete g_v|index>, with sum_v g_v equal to local H_geometry.

    Uses local formulas rather than allocating the full exponential matrix.
    """
    n, edges = _graph(graph)
    labels = phi_labels(index, n, config)
    if not 0 <= vertex < n:
        raise ValueError('invalid shard vertex')
    parameters = (*config.phi_values, config.hopping, config.regulator, config.lapse_coupling, config.stiffness, *config.boundary_weights)
    if not all(math.isfinite(x) for x in parameters):
        raise ValueError('finite local geometry labels and coefficients required')
    values = tuple(_exact(x) for x in config.phi_values)
    count = len(values)
    if count < 2 or any(values[i+1] <= values[i] for i in range(count-1)):
        raise ValueError('ordered local geometry labels required')
    weights = config.boundary_weights or (0.,)*n
    if len(weights) != n or any(w < 0 for w in weights):
        raise ValueError('one nonnegative boundary weight per vertex')
    if min(config.hopping, config.regulator, config.lapse_coupling) < 0 or config.stiffness <= 0:
        raise ValueError('invalid local geometry coefficients')
    if any(1+_exact(config.lapse_coupling)*x <= 0 for x in values):
        raise ValueError('lapse must be positive')
    output = {}
    j = labels[vertex]
    hopping, stiffness = _exact(config.hopping), _exact(config.stiffness)
    diagonal = sy.Integer(0)
    for neighbour in (j-1, j+1):
        if 0 <= neighbour < count:
            target = list(labels); target[vertex] = neighbour
            output[phi_index(target, count)] = -hopping
            diagonal += hopping
    diagonal += stiffness*(_exact(weights[vertex])+_exact(config.regulator))*values[j]**2/2
    for u, v, weight in edges:
        if vertex in (u, v):
            diagonal += stiffness*_exact(weight)*(values[labels[u]]-values[labels[v]])**2/4
    output[index] = sy.expand(diagonal)
    return {key: value for key, value in output.items() if value != 0}


def local_energy_column(shard_action, key, vertex, graph, config=LocalGeometryConfig()):
    """k_v on (Phi_index, native_source_key), including all source outputs."""
    n, _ = _graph(graph)
    index, native = key
    labels = phi_labels(index, n, config)
    reference = _exact(_references(n, config)[vertex])
    lapse = 1+_exact(config.lapse_coupling)*_exact(config.phi_values[labels[vertex]])
    output = {(row, native): c for row, c in geometry_shard_column(index, vertex, graph, config).items()}
    for target, coefficient in shard_action(native).items():
        k = (index, target)
        output[k] = output.get(k, 0)+lapse*sy.sympify(coefficient)
    output[key] = output.get(key, 0)+(1-lapse)*reference
    return {k: sy.expand(c) for k, c in output.items() if sy.expand(c) != 0}


def apply_local_capacity_column(shard_actions, key, graph,
                                local_config=LocalGeometryConfig(),
                                capacity_config=GeometryConfig()):
    """Full H_v2 column in order (capacity_indices,(Phi_index,source_key))."""
    n, _ = _graph(graph)
    if len(shard_actions) != n:
        raise ValueError('every source shard must be supplied')
    capacities, local_key = key
    if len(capacities) != n:
        raise ValueError('one complete capacity auxiliary per vertex required')
    acols, ncols = auxiliary_columns(capacity_config)
    if any(not isinstance(i, int) or not 0 <= i < len(acols) for i in capacities):
        raise ValueError('invalid capacity index')
    output = {}
    def add(k, c):
        output[k] = output.get(k, sy.Integer(0))+c
    reference = _exact(capacity_config.source_energy_reference)
    for v in range(n):
        for target, coefficient in acols[capacities[v]]:
            changed = list(capacities); changed[v] = target
            add((tuple(changed), local_key), coefficient)
        k_action = local_energy_column(shard_actions[v], local_key, v, graph, local_config)
        for target, coefficient in ncols[capacities[v]]:
            changed = list(capacities); changed[v] = target
            for local_target, c in k_action.items():
                add((tuple(changed), local_target), coefficient*c)
            add((tuple(changed), local_key), -coefficient*reference)
        add(key, reference)
    return {k: sy.expand(c) for k, c in output.items() if sy.expand(c) != 0}


def apply_native_column(source, key, cap=6, **kwargs):
    from . import source as master
    graph = source['graph']
    actions = [lambda s, v=v: master.apply_source_column(source, s[0], s[1], cap=cap, vertex=v)
               for v in range(graph.nvertices)]
    return apply_local_capacity_column(actions, key, graph, **kwargs)


def carrier_dimension(vertices, edges, cap=6, local_config=LocalGeometryConfig(),
                      capacity_config=GeometryConfig()):
    """Intrinsic kinematic dimension for a fixed reference-charge sector."""
    if any(not isinstance(x,int) for x in (vertices,edges,cap)) or vertices<1 or edges<0 or cap<0:
        raise ValueError('positive vertex count and nonnegative edge count/cutoff required')
    capacity=len(auxiliary_columns(capacity_config)[0])
    return (1<<(96*vertices))*(2*cap+1)**edges*len(local_config.phi_values)**vertices*capacity**vertices
