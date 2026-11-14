"""Independent finite gravity-source checks and their observation domains.

Uniform native seam attempts give E[n'] = (I-L/(2M))n. A tagged
hard-core load has hop rate (1-K/11)/12 on the twelve-port carrier,
while its attempted-contact rate is 1/6. For ds²=-N²dt²+dx²,
G_tt=0 and G_ij=(delta_ij Delta N-d_i d_j N)/N.
"""
from __future__ import annotations

import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations
import json
import math
from pathlib import Path
import sys

import numpy as np
import sympy as sy


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(workspace):
    d3 = workspace / 'openwave-push/issues/D3_gravity'
    kernel = workspace / 'openwave-push/kernel/ophk.py'
    native_file = workspace / 'reverse-engineering-reality/model/observer/native.py'
    names = [kernel, native_file, d3/'parts.py', d3/'d3lib.py',
             d3/'results_native.json', d3/'results_clock.json',
             d3/'results_lookup.json', d3/'results_bridge.json',
             workspace/'reverse-engineering-reality/evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
             workspace/'reverse-engineering-reality/evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/federation.py']
    before = {str(p.relative_to(workspace)): digest(p) for p in names}
    spec = importlib.util.spec_from_file_location('gravity_audit_geometry', kernel)
    geometry = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = geometry
    spec.loader.exec_module(geometry)
    tree = ast.parse(native_file.read_text())
    ns = {}
    body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'native']
    assert len(body) == 1
    exec(compile(ast.Module(body=body, type_ignores=[]), str(native_file), 'exec'), ns)
    native = ns['native']
    car = geometry.carrier()
    edges = car.seams
    assert len(edges) == 30 and all(sum(v in e for e in edges) == 5 for v in range(12))
    local_checks = 0
    for a in range(-8, 9):
        for b in range(-8, 9):
            out = [native((a, b), (0, 1), c) for c in (0, 1)]
            assert sum(t[1][1] for t in out) == a-b
            assert all(sum(t[0]) == a+b for t in out)
            assert all(t == geometry.native((a, b), (0, 1), c) for c, t in enumerate(out))
            local_checks += 2

    # Only existing development totals one, two, three; total four is reserved.
    nr = json.loads((d3/'results_native.json').read_text())
    existing = {r['K']: r for r in nr['carrier_exact']}
    rates = []
    for k in (0, 1, 2):
        moved = touched = count = 0
        for tag in range(12):
            for occupied in combinations([v for v in range(12) if v != tag], k):
                n = tuple(int(v == tag or v in occupied) for v in range(12))
                for e in edges:
                    for coin in (0, 1):
                        nn, _ = native(n, e, coin)
                        moved += nn[tag] == 0
                        touched += tag in e
                        count += 1
        hop, touch = F(moved, count), F(touched, count)
        assert hop == F(11-k, 132) == F(existing[k]['hop_rate'])
        assert touch == F(1, 6) == F(existing[k]['touch_rate'])
        rates.append(dict(other_loads=k, branches=count, hop=str(hop), contact=str(touch)))

    # Rebuild the existing pin preparation; use our own adjacency/BFS, no T eigensolver.
    fed = geometry.flat(3, 3)
    pos = np.array(fed.cell_points)
    delta = pos-pos[0]
    v = delta[:, 1]/1.5
    u = (delta[:, 0]-v*math.sqrt(3)/2)/math.sqrt(3)
    u -= 3*np.round(u/3)
    v -= 3*np.round(v/3)
    radii = np.hypot(u*math.sqrt(3)+v*math.sqrt(3)/2, 1.5*v)
    carriers = np.flatnonzero(radii < 1.01)
    rng = np.random.default_rng(7)
    pinned = set()
    for c in carriers:
        for p in range(12):
            if rng.random() < .5:
                pinned.add(int(12*c+p))
    adj = {v: set() for v in range(fed.ports) if v not in pinned}
    for a, b in fed.seams:
        if a in adj and b in adj:
            adj[a].add(b)
            adj[b].add(a)
    pending = set(adj)
    components = []
    while pending:
        todo = [pending.pop()]
        component = set(todo)
        while todo:
            x = todo.pop()
            for y in adj[x] & pending:
                pending.remove(y)
                component.add(y)
                todo.append(y)
        components.append(sorted(component))
    sizes = sorted(map(len, components))
    assert sizes == [1, 5, 7, 180]
    saved = nr['pinned_patch_exact']
    assert saved['components'] == len(components)
    assert saved['isolated_available_ports'] == sum(not a for a in adj.values())
    assert saved['component_sizes'] == sizes
    assert saved['stationary_probability_per_main_component_port'] == str(F(1, max(sizes)))
    patch = dict(component_sizes=sizes, isolated_vertices=sum(not a for a in adj.values()),
                 outside_largest_component=len(adj)-max(sizes), available_vertices=len(adj),
                 largest_component_stationary_probability=str(F(1, max(sizes))),
                 saved_isolated_available_ports=saved['isolated_available_ports'])

    # Independent linearized curvature, derived from h rather than a nonlinear metric solver.
    x, y, z, t, mu, gamma = sy.symbols('x y z t mu gamma', real=True)
    coords = [t, x, y, z]
    eta = [-1, 1, 1, 1]
    potential = mu/sy.sqrt(x*x+y*y+z*z)
    h = sy.diag(2*potential, 2*gamma*potential, 2*gamma*potential, 2*gamma*potential)
    trace = sum(eta[a]*h[a, a] for a in range(4))
    ric = sy.zeros(4)
    for a in range(4):
        for b in range(4):
            ric[a, b] = sy.simplify((sum(eta[c]*(sy.diff(h[c, b], coords[c], coords[a])
                + sy.diff(h[c, a], coords[c], coords[b])-sy.diff(h[a, b], coords[c], 2))
                for c in range(4))-sy.diff(trace, coords[a], coords[b]))/2)
    scalar = sum(eta[a]*ric[a, a] for a in range(4))
    einstein = ric-sy.diag(*eta)*scalar/2
    residuals = {}
    br = json.loads((d3/'results_bridge.json').read_text())
    for g in (0, 1):
        axis = sy.simplify(einstein[1, 1].subs({gamma:g, x:3, y:0, z:0})/mu)
        zero = sy.simplify(einstein[0, 0]/mu)
        assert axis == sy.Rational(2*(1-g), 27) and zero == 0
        assert str(axis) == br['T10_exterior_residual'][f'gamma_{g}']['first_order_G11_at_3_0_0_per_mu']
        residuals[str(g)] = str(axis)

    lookup = json.loads((d3/'results_lookup.json').read_text())['strength']
    kappa = F(4835703278458516698824704, 642754164150366239279661)
    energy = F(51099895069, 10**14)
    length = F(236221477831373991, 312500000000000000)
    hc = F('0.1973269804')
    gobs = F('6.70883e-39')
    coupling = 1/(4*math.pi*float(kappa))
    speed = float(energy*length/hc)
    grav = coupling*float(length/(energy*hc))
    invariant = coupling/float(energy**2)
    assert math.isclose(grav/speed, invariant, rel_tol=2e-15)
    assert math.isclose(grav, lookup['G_model_GeV_m2'], rel_tol=2e-15)
    assert math.isclose(invariant/float(gobs), lookup['ratio_with_light_matched'], rel_tol=2e-15)

    clocks = json.loads((d3/'results_clock.json').read_text())
    m, lam = clocks['m'], 6.0
    e = math.sqrt(m*m+lam)
    assert math.isclose(e, clocks['E'], rel_tol=2e-15)
    clock_outputs = []
    for row in clocks['clocks']:
        n = row['N']
        assert abs(row['fitted_omega_over_E']-n) < 2e-14
        assert abs(row['fitted_depth']-lam/e**2) < 2e-14
        clock_outputs.append(dict(carrier=row['carrier'], sin_argument_rate=n*e,
                                  population_angular_frequency=2*n*e, lapse_ratio=n))
    assert before == {str(p.relative_to(workspace)): digest(p) for p in names}
    return dict(schema='oph.openwave.gravity_source_audit.v1',
        checks=dict(native_endpoint_mean_exact=True, literal_tagged_rates_match=True,
                    independently_rebuilt_disconnected_domain=True,
                    linearized_exterior_metric_match=True, common_unit_invariant_match=True,
                   population_frequency_factor_two_explicit=True, source_pins_unchanged=True),
        local_endpoint_branches=local_checks, tagged_rates=rates, pinned_domain=patch,
        exterior_Gxx_over_mu_at_x3=residuals,
        common_units=dict(G=grav, speed_over_c=speed, G_over_speed=invariant,
                          G_over_observed_with_light_matched=invariant/float(gobs)),
        population_clocks=clock_outputs, source_sha256=before,
        auditor_sha256=digest(Path(__file__)), reserved_acceptance_executed=False,
        scope='Finite native C events, supplied static pins and supplied lapse comparison. No derived gravitational field, physical Newton constant, universal event-clock lapse or self-starting clock.',
        official_rows_changed=[])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.workspace.resolve())
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['checks','pinned_domain','common_units']}))
