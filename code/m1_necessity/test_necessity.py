"""Semantic, analytic-control and actual CLI rejection tests."""

from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import pytest

from . import model, check, verify, audit, balanced, balanced_check

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def packet():
    return check.strict_load(Path(__file__).with_name("receipt.json"))


def test_complete_reference_execution_and_counts(packet):
    actual = verify.verify(packet)
    assert {k:actual[k] for k in ('routes','route_steps','executed_events','executed_reads','consumer_claims')} == {
        "routes":298,"route_steps":2215,"executed_events":149760,"executed_reads":13479912,"consumer_claims":8}
    assert actual['balanced_routes']==161 and actual['balanced_steps']>0
    assert actual['downstream_claims']==133


def test_operational_clock_descendants_do_not_inherit_dense_m1_premises():
    downstream = audit.downstream_audit()
    inputs = {
        'OPH-M1-MASSIVE-FLIGHT-SPECTRUM': {
            'declared_charged_coherent_flight_grammar', 'declared_massive_tetrahedral_process'},
        'OPH-M1-RESOLVED-LOCAL-READ-CONE': {
            'declared_massive_tetrahedral_process', 'declared_smooth_local_car_operations'},
        'OPH-M1-OBSERVABLE-MASS-CLOCK': {
            'declared_massive_tetrahedral_process', 'declared_two_mass_clock_preparation'}}
    for name, assumptions in inputs.items():
        row = downstream['claims'][name]
        assert set(row['assumptions']) == assumptions
        assert row['disposition'] == 'existing_contract_retained_without_automatic_transfer'
        assert not (set(row['assumptions']) & audit.KEYS)
    # The massive result is attached through the actual tetrahedral parent;
    # reachability records that use, without transferring a dense read law.
    assert any(edge['from'] == 'OPH-M1-MINIMAL-TETRAHEDRAL-TRANSPORT' and
               edge['to'] == 'OPH-M1-MASSIVE-FLIGHT-SPECTRUM'
               for edge in downstream['links'])


@pytest.mark.parametrize("bad", [True,False,0,-1,25,1.0,"2",None])
def test_invalid_levels_fail(bad):
    with pytest.raises(ValueError):
        model.sparse(bad)
    with pytest.raises(ValueError):
        check.expected_sparse(bad)


@pytest.mark.parametrize("mutation", ["empty","drop","duplicate","too_fast","fixed_axes","boolean","asymmetric"])
def test_substituted_stencil_fails(mutation):
    vectors = model.sparse(2)
    if mutation == "empty": vectors = []
    elif mutation == "drop": vectors.pop()
    elif mutation == "duplicate": vectors.append(vectors[0])
    elif mutation == "too_fast": vectors[0] = (-5,0,0)
    elif mutation == "fixed_axes": vectors = [v for v in vectors if sum(x!=0 for x in v) <= 1]
    elif mutation == "boolean": vectors[0] = (True,0,0)
    else: vectors = [v for v in vectors if v[0] >= 0]
    with pytest.raises(ValueError):
        check.check_sparse(2,vectors)


@pytest.mark.parametrize("mutation", ["empty","missing","duplicate","target","repeat_bool","zero_repeat",
                                      "negative_repeat","huge_repeat","illegal_step","dirty_endpoint","count","digest"])
def test_routes_cannot_be_truncated_relabelled_or_resealed(packet,mutation):
    rows = deepcopy(packet["routes"])
    if mutation == "empty": rows = []
    elif mutation == "missing": rows.pop()
    elif mutation == "duplicate": rows[-1] = rows[0]
    elif mutation == "target": rows[1]["displacement"] = [0,0,0]
    elif mutation == "count": rows[1]["steps"] -= 1
    elif mutation == "digest": rows[1]["points_sha256"] = "0"*64
    elif mutation in ("repeat_bool","zero_repeat","negative_repeat","huge_repeat"):
        rows[1]["instructions"][0][0] = {"repeat_bool":True,"zero_repeat":0,
                                         "negative_repeat":-1,"huge_repeat":10001}[mutation]
    elif mutation == "illegal_step": rows[1]["instructions"][0][1:] = [3,0,0]
    else: rows[1]["instructions"].pop()
    with pytest.raises(ValueError):
        verify.verify_routes(rows)


@pytest.mark.parametrize("mutation", ["schema","missing_level","fake_stencil","cut","bool_cut",
                                      "moment","zero_action","fake_stability","missing_consumer","scope",
                                      "graph","missing_intervention","source_pin"])
def test_false_or_incomplete_receipts_fail(packet,mutation):
    bad = deepcopy(packet)
    row = bad["cuts"]["2"]
    if mutation == "schema": bad["trusted"] = True
    elif mutation == "missing_level": del bad["cuts"]["20"]
    elif mutation == "fake_stencil": row["stencil_sha256"] = "0"*64
    elif mutation == "cut": row["cuts"]["axial_diagonal"] += 1
    elif mutation == "bool_cut": row["cuts"]["sparse_degree"] = True
    elif mutation == "moment": row["actions"]["sparse"]["second_moments"][0][1] = 1
    elif mutation == "zero_action": row["actions"]["sparse"]["alpha_times_h_squared"] = "0"
    elif mutation == "fake_stability": row["actions"]["sparse"]["tick_squared_eigenvalues"] = ["1"]*3
    elif mutation == "missing_consumer": bad["consumers"].popitem()
    elif mutation == "scope": bad["consumers"]["OPH-BH-CROSSING-READ-AREA-LAW"]["disposition"] = "preserved"
    elif mutation == "graph": bad["graphs"]["2:sparse"]["baseline"]["reads"] -= 1
    elif mutation == "missing_intervention": del bad["graphs"]["2:sparse"]["intermediate"]
    else: bad["sources"].popitem()
    with pytest.raises(ValueError):
        verify.verify(bad)


def test_short_legal_route_with_wrong_boundary_behavior_is_rejected():
    # All steps are legal, the endpoint is exact and the time bound passes.
    # Its off-axis detour still leaves the proved tube; resealing a points
    # digest cannot turn this path into a valid boundary certificate.
    steps = [(192,168,0)]*160+[(192,-168,0)]*160
    with pytest.raises(ValueError,match="boundary tube"):
        check.check_route(8,(61440,0,0),steps)


def test_forged_intermediate_states_fail_even_with_correct_event_counts(packet):
    bad = deepcopy(packet)
    bad["graphs"]["2:sparse"]["intermediate"]["states"][2] = check.digest([0]*4096)
    with pytest.raises(ValueError,match="false graph execution"):
        verify.verify(bad)


def test_zero_boundary_energy_includes_missing_neighbors():
    sites,outgoing,_ = check.graph_replay(1,"sparse")
    degree = len(check.expected_sparse(1))-1
    values = [((13*i)%11)-5 for i in range(len(sites))]
    action = [degree*values[i]-sum(values[j] for j in targets if j != i)
              for i,targets in enumerate(outgoing)]
    twice_energy = sum(a*b for a,b in zip(values,action))
    internal = sum((values[i]-values[j])**2 for i,targets in enumerate(outgoing) for j in targets)
    boundary = sum((degree-(len(targets)-1))*values[i]**2 for i,targets in enumerate(outgoing))
    assert twice_energy == F(internal,2)+boundary > 0
    assert boundary > 0


def test_consumer_inventory_detects_an_added_dependency():
    registry = check.strict_load(ROOT/"claims/claim_registry.yaml")
    extra = deepcopy(registry["claims"][0])
    extra["claim_id"] = "NEW-M1-CONSUMER"
    extra["assumptions"] = [next(iter(audit.KEYS))]
    registry["claims"].append(extra)
    with pytest.raises(ValueError,match="coverage"):
        audit.consumer_audit(registry)


def test_downstream_closure_rejects_omission_and_new_dependency(packet):
    expected=check.strict_load(ROOT/'code/m1_necessity/downstream.json')
    with pytest.raises(ValueError,match='coverage'):
        audit.downstream_audit(expected=expected[:-1])
    registry=check.strict_load(ROOT/'claims/claim_registry.yaml')
    graph=check.strict_load(ROOT/'claims/dependency_graph.json')
    other=next(row['claim_id'] for row in registry['claims'] if row['claim_id'] not in set(expected)|audit.OWN)
    graph['edges'].append({'from':next(iter(audit.EXPECTED)),'to':other,'role':'new explicit dependency'})
    with pytest.raises(ValueError,match='coverage'):
        audit.downstream_audit(registry,graph)
    bad=deepcopy(packet);bad['downstream']['claims'].popitem()
    with pytest.raises(ValueError,match='downstream'):
        verify.verify(bad)


def test_all_small_cut_rectangles_by_direct_lattice_enumeration():
    for q in (4,7):
        for x,y in product(range(-3,4),repeat=2):
            if x+y <= 0:
                continue
            actual = sum(1 for i,j in product(range(q),repeat=2)
                         if 0<=i+x<q and 0<=j+y<q and i+j<q<=i+j+x+y)*q
            assert check.cut(q,[(x,y,0)],True) == actual
            assert model.cut_count(q,[(x,y,0)],True) == actual


def test_exact_axial_sums_and_actual_pair_census():
    for r in range(1,65):
        assert sum(j*j for j in range(1,r+1)) == r*(r+1)*(2*r+1)//6
        assert sum(j**4 for j in range(1,r+1)) == r*(r+1)*(2*r+1)*(3*r*r+3*r-1)//30
    for family in ("sparse","axial","dense"):
        for diagonal in (False,True):
            expected = model.cuts(1)[family+("_diagonal" if diagonal else "_coordinate")]
            assert check.direct_cut_pairs(1,family,diagonal) == expected


def test_adjacent_layer_deletion_and_lossless_read_counterexample():
    sparse = set(model.stencil(2,"sparse"))
    dense = set(model.stencil(2,"dense"))
    missing = sorted(dense-sparse)
    assert missing
    # One-layer time forbids a substitute two-step path, even if spatial sums agree.
    removed = missing[0]
    point = (8,8,8)
    parent = tuple(x-v for x,v in zip(point,removed))
    assert all(0<=i<16 for i in parent)
    a = {v:0 for v in dense}
    b = dict(a); b[removed] = 1
    assert [a[v] for v in sorted(sparse)] == [b[v] for v in sorted(sparse)]
    assert [a[v] for v in sorted(dense)] != [b[v] for v in sorted(dense)]


def test_exact_operator_on_all_monomials_through_degree_three():
    point = (1,-2,3)
    for family in ("sparse","axial","dense"):
        vs = model.stencil(2,family)
        alpha = F(model.action(2,family)["alpha_times_h_squared"])
        for powers in product(range(4),repeat=3):
            if sum(powers)>3:
                continue
            def value(x):
                result = 1
                for coordinate,power in zip(x,powers): result *= coordinate**power
                return result
            laplacian = 0
            for axis,power in enumerate(powers):
                if power>=2:
                    term = power*(power-1)
                    for i in range(3): term *= point[i]**(powers[i]-(2 if i==axis else 0))
                    laplacian += term
            actual = alpha*sum(value(point)-value(tuple(x+y for x,y in zip(point,v))) for v in vs)
            assert actual == -laplacian


def test_parity_is_an_actual_unstable_periodic_mode():
    r,q = model.parameters(1)
    vectors = [v for v in model.sparse(1) if any(v)]
    alpha = F(model.action(1,"sparse")["alpha_times_h_squared"])
    mode = lambda x: -1 if x[0]%2 else 1
    eigen = F(model.action(1,"sparse")["tick_squared_eigenvalues"][0])
    assert eigen > 4
    for x in product(range(q),repeat=3):
        actual = r*r*alpha*sum(mode(x)-mode(tuple((a+b)%q for a,b in zip(x,v))) for v in vectors)
        assert actual == eigen*mode(x)


def test_verifier_runs_without_importing_producer():
    script = '''
import sys
sys.modules['m1_necessity.model'] = None
sys.modules['m1_necessity.build'] = None
sys.modules['m1_necessity.balanced'] = None
from m1_necessity import verify, check
print(verify.verify(check.strict_load(verify.HERE/'receipt.json')))
'''
    run = subprocess.run([sys.executable,"-c",script],cwd=ROOT,
                         env=dict(os.environ,PYTHONPATH=str(ROOT/"code")),capture_output=True,text=True)
    assert run.returncode == 0,run.stdout+run.stderr


@pytest.mark.parametrize("text", ['{"schema":1,"schema":2}', '{"x":NaN}', '{"x":Infinity}'])
def test_noncanonical_json_rejected(tmp_path,text):
    path = tmp_path/"bad.json";path.write_text(text,encoding="utf-8")
    with pytest.raises(ValueError): check.strict_load(path)


def test_actual_optimized_cli_rejects_junk(tmp_path):
    path = tmp_path/"junk.json";path.write_text('{}',encoding="utf-8")
    run = subprocess.run([sys.executable,"-O","-m","m1_necessity.verify","--receipt",str(path)],
                         cwd=ROOT,env=dict(os.environ,PYTHONPATH=str(ROOT/"code")),capture_output=True,text=True)
    assert run.returncode != 0
    assert "exact receipt schema" in run.stderr


@pytest.mark.parametrize('mutation',['empty','missing_level','radius','spacing','degree','cut','moment','drop_route','relabel','zero_cost','bad_edge','no_stabilizer','hide_aliases'])
def test_balanced_certificate_cannot_hide_geometry_or_work(packet,mutation):
    candidate=deepcopy(packet['balanced'])
    if mutation=='empty': candidate={}
    elif mutation=='missing_level': del candidate['cuts']['3']
    elif mutation in ('radius','spacing','degree'): candidate['cuts']['1'][mutation]+=1
    elif mutation=='cut': candidate['cuts']['1']['coordinate']+=1
    elif mutation=='moment': candidate['cuts']['1']['second_norm_moment']=0
    elif mutation=='drop_route': candidate['routes'].pop()
    elif mutation=='relabel': candidate['routes'][1]['displacement']=[0,0,0]
    elif mutation=='zero_cost': candidate['routes'][1]['instructions'][0][0]=0
    elif mutation=='no_stabilizer': candidate['cuts']['1']['stabilized']['nearest_extra_weight_times_h_squared']='0'
    elif mutation=='hide_aliases': candidate['cuts']['1']['stabilized']['naive_alias_modes']=1
    else: candidate['routes'][1]['instructions'][0][1:]=[129,0,0]
    with pytest.raises(ValueError): balanced_check.verify(candidate)


@pytest.mark.parametrize('mutation,message', [
    ('endpoint', 'balanced endpoint'),
    ('charged_time', 'balanced charged time bound'),
    ('tube', 'balanced boundary tube'),
    ('digest', 'balanced path commitment'),
])
def test_balanced_resealed_routes_fail_semantic_checks(packet, mutation, message):
    candidate = deepcopy(packet['balanced'])
    row = candidate['routes'][1]
    assert row['t'] == 1 and row['displacement'] == [255, 0, 0]
    steps = [tuple(item[1:]) for item in row['instructions'] for _ in range(item[0])]
    if mutation == 'endpoint':
        steps.append((1, 0, 0))  # Legal edge, wrong endpoint.
    elif mutation == 'charged_time':
        steps = [(0, 0, 0)]*20 + steps  # Legal waits, all inside the tube.
    elif mutation == 'tube':
        # Legal radius flights, unchanged endpoint and still within the time
        # budget; the off-axis excursion alone violates the 2m tube.
        steps = [(0, 128, 0)] + steps + [(0, -128, 0)]
        assert len(steps) <= 16
    assert all(balanced_check.admissible(1, step) for step in steps)
    point, points = [0, 0, 0], [[0, 0, 0]]
    for step in steps:
        point = [x+y for x, y in zip(point, step)]
        points.append(point)
    row.update(instructions=[[1, *step] for step in steps], steps=len(steps),
               points_sha256=check.digest(points))
    if mutation == 'digest':
        row['points_sha256'] = '0'*64
    with pytest.raises(ValueError, match=message):
        balanced_check.verify(candidate)


def test_balanced_routes_at_half_ties_corners_and_all_supported_levels():
    from random import Random
    rng = Random(981)
    for t in range(1, 9):
        q, m = 1 << (8*t), 1 << (4*t)
        targets = list(product((-m//2, 0, m//2), repeat=3))
        targets += list(product((-q+1, q-1), repeat=3))
        targets += [tuple(rng.randrange(-q+1, q) for _ in range(3)) for _ in range(16)]
        for target in targets:
            steps = balanced.route(t, target)
            result = balanced_check.check_route(t, target, steps)
            assert result['steps'] == len(steps)


def test_balanced_small_complete_pair_census_and_moments():
    # A miniature member of the same generic spaced-ball-plus-dyadic rule.
    q,m,k=16,2,2
    coarse={tuple(m*x for x in v) for v in product(range(-k,k+1),repeat=3) if sum(x*x for x in v)<=k*k}
    vectors=coarse|{tuple(sign if i==axis else 0 for i in range(3)) for axis in range(3) for sign in (-1,1)}
    for diagonal in (False,True):
        projection=lambda v: 2*(v[0]+v[1] if diagonal else v[0])
        threshold=2*q-1 if diagonal else q-1
        count=0
        for x in product(range(q),repeat=3):
            if projection(x)>=threshold: continue
            for v in vectors:
                y=tuple(a+b for a,b in zip(x,v))
                if all(0<=a<q for a in y) and projection(y)>threshold: count+=1
        assert count==model.cut_count(q,vectors,diagonal)==check.cut(q,vectors,diagonal)
    matrix=[[sum(v[i]*v[j] for v in vectors) for j in range(3)] for i in range(3)]
    assert all(matrix[i][j]==(matrix[0][0] if i==j else 0) for i,j in product(range(3),repeat=2))


def test_balanced_t1_column_sums_against_explicit_stencil():
    q,m,k=balanced.parameters(1)
    vectors={tuple(m*x for x in v) for v in product(range(-k,k+1),repeat=3) if sum(x*x for x in v)<=k*k}
    vectors.update(tuple(sign*2**b if i==axis else 0 for i in range(3))
                   for axis,sign,b in product(range(3),(-1,1),range(4)))
    row=balanced.certificate(1)
    assert row==balanced_check.certificate(1)
    assert len(vectors)==row['degree']
    assert sum(sum(x*x for x in v) for v in vectors)==row['second_norm_moment']
    assert sum(sum(x*x for x in v)**2 for v in vectors)==row['fourth_norm_moment']
    assert check.cut(q,vectors)==row['coordinate']
    assert check.cut(q,vectors,True)==row['diagonal']


def test_fixed_radius_deletions_pay_the_proved_cut_deficit():
    q,r=16,4
    dense=set(model.stencil(2,'dense'))
    for cutoff in (0,1,2,3,4):
        removed={v for v in dense if sum(x*x for x in v)<=cutoff*cutoff}
        deficit=sum(check.cut(q,[tuple(v[(i+j)%3] for j in range(3)) for v in removed]) for i in range(3))
        weight=sum(sum(abs(x) for x in v) for v in removed)
        assert 2*deficit>=(q-r)**2*weight


def test_all_periodic_t1_alias_modes_and_their_correction():
    import numpy as np
    q,m,_=balanced.parameters(1)
    row=balanced.certificate(1)
    ell=np.array(list(product(range(-m//2,m//2),repeat=3)),dtype=float)
    phases=2*np.pi*ell/m
    raw=np.zeros(len(ell))
    alpha=6*q*q/row['second_norm_moment']
    for bit in range(4):
        raw+=2*alpha*np.sum(1-np.cos(phases*2**bit),axis=1)
    nearest=2*q*q*np.sum(1-np.cos(phases),axis=1)
    corrected=(raw+nearest)/2
    wave_squared=np.sum((2*np.pi*(q/m)*ell)**2,axis=1)
    nonzero=wave_squared>0
    assert len(ell)==row['stabilized']['naive_alias_modes']==4096
    assert max(raw)<=float(F(row['stabilized']['naive_alias_eigenvalue_times_L_squared_upper']))+1e-10
    assert max(raw)<min((2/np.pi**2)*wave_squared[nonzero])
    assert np.all(corrected>=(2/np.pi**2)*wave_squared-1e-7)
    assert min(corrected[nonzero])>=row['stabilized']['stabilized_nonzero_alias_eigenvalue_times_L_squared_lower']-1e-7


def test_entire_small_clipped_spectrum_has_the_proved_lower_bound():
    import numpy as np
    q,m,k=8,2,2
    vectors={tuple(m*x for x in v) for v in product(range(-k,k+1),repeat=3) if 0<sum(x*x for x in v)<=k*k}
    nearest={tuple(sign if i==axis else 0 for i in range(3)) for axis in range(3) for sign in (-1,1)}
    vectors|=nearest
    moment=sum(sum(x*x for x in v) for v in vectors)
    weights={v:3*q*q/moment+(q*q/2 if v in nearest else 0) for v in vectors}
    sites=list(product(range(q),repeat=3));index={x:i for i,x in enumerate(sites)}
    operator=np.eye(len(sites))*sum(weights.values())
    for i,x in enumerate(sites):
        for v,w in weights.items():
            y=tuple(a+b for a,b in zip(x,v))
            if y in index: operator[i,index[y]]-=w
    assert np.array_equal(operator,operator.T)
    eigenvalues=np.linalg.eigvalsh(operator)
    nn=sorted(4*q*q*sum(np.sin(np.pi*j/(2*(q+1)))**2 for j in ks)
              for ks in product(range(1,q+1),repeat=3))
    assert np.all(eigenvalues>=np.asarray(nn)/2-1e-8)
    assert np.all(eigenvalues>=np.arange(1,len(sites)+1)**(2/3)/2-1e-8)
    assert eigenvalues[0]>0


def test_all_public_lean_results_have_transitive_audits_and_ci():
    source = (ROOT/"Lean/Geometry/M1Necessity.lean").read_text(encoding="utf-8")
    gate = (ROOT/"Lean/Geometry/M1NecessityAxiomAudit.lean").read_text(encoding="utf-8")
    names = re.findall(r'^theorem (\w+)',source,re.M)
    assert len(names) == 17
    assert {"OPH.M1Necessity."+n for n in names} == set(re.findall(r'^audit_m1_necessity (OPH\.\S+)',gate,re.M))
    assert '#guard_msgs in\naudit_m1_necessity sorryAx' in gate
    assert '#guard_msgs in\naudit_m1_necessity Lean.ofReduceBool' in gate
    ci = (ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8")
    assert re.search(r'^\s*"Geometry.M1NecessityAxiomAudit"\s*$',ci,re.M)
