"""Data-free checks of new direct twelve-port falsification proposals.

The response-rank witness is a new elementary consequence of A1, not an
existing source-selected laboratory identification. The SU(4) example below
is an algebraic negative control, not a measurement or an OPH counterexample.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
from model.observer import native
from model.record_laws import family
from audit_dispersion import geometry


def population_check():
    vertices, _, _ = geometry()
    adjacency = sp.Matrix(12,12,lambda i,j: int(i != j and abs(float(vertices[i]@vertices[j])-1/math.sqrt(5)) < 1e-12))
    laplacian = 5*sp.eye(12)-adjacency
    seams = [(i,j) for i,j in itertools.combinations(range(12),2) if adjacency[i,j]]
    states = native.states(binary=True)
    expected = sp.eye(12)-laplacian/60
    branch_checks = 0
    for state in states:
        summed = sp.zeros(12,1)
        for seam in seams:
            a,b = seam
            for coin in (0,1):
                # Independently specified fair endpoint-balancing oracle.
                oracle = list(state)
                total = state[a]+state[b]
                low,high = total//2,(total+1)//2
                oracle[a],oracle[b] = (low,high) if coin==0 else (high,low)
                actual,_ = native.native(state,seam,coin)
                assert tuple(oracle)==actual
                summed += sp.Matrix(actual)/60
                branch_checks += 1
        assert summed == expected*sp.Matrix(state)
    # Twelve preparation ensembles within the actual binary total-two sector.
    preparations = sp.ones(12,12)/11+sp.eye(12)*sp.Rational(10,11)
    assert preparations.det() != 0
    outputs = expected*preparations
    assert outputs*preparations.inv()==expected
    eta = sp.symbols("eta", positive=True)
    held = sp.eye(12)-eta*laplacian/60
    low,mid,high=(5-sp.sqrt(5))*eta/60,eta/10,(5+sp.sqrt(5))*eta/60
    assert sp.simplify(high/low-(3+sp.sqrt(5))/2)==0
    assert sp.simplify(mid/low-3*(5+sp.sqrt(5))/10)==0
    return {"native_branch_checks":branch_checks,"binary_states":len(states),
            "preparation_matrix_rank":int(preparations.rank()),
            "one_attempt_population_matrix":"I-eta*L/60",
            "deficit_spectrum":"0 [1], eta*(5-sqrt(5))/60 [3], eta/10 [5], eta*(5+sqrt(5))/60 [3]",
            "high_low_deficit_ratio":str(sp.simplify(high/low)),
            "middle_low_deficit_ratio":str(sp.simplify(mid/low))}


def interference_check():
    real,imag = sp.symbols("real imag", real=True)
    p = [sp.Rational(3,8)+real/2,sp.Rational(3,8)+imag/2,
         sp.Rational(3,8)-real/2,sp.Rational(3,8)-imag/2]
    assert sp.simplify(p[0]+p[2])==sp.Rational(3,4)
    assert sp.simplify(p[1]+p[3])==sp.Rational(3,4)
    assert sp.expand((p[0]-p[2])**2+(p[1]-p[3])**2)==real**2+imag**2
    # The repository's exact coherent preparations and probability readouts.
    # Convert sympy rationals to the source's exact number domain explicitly.
    from fractions import Fraction as F
    c = family.table(1,F(1,2),F(1,2),0,F(1,2))
    names = {"C":c,"E":family.table(1,0,0,0,0),"T":family.table(0,0,0,0,0)}
    names.update({n:t for n,t in family.generic_tables().items() if n!="outside"})
    values={}
    for name,table in names.items():
        assert family.positive_semidefinite(table)
        probs={setting:family.plus_probability(family.apply_table(table,family.prepared(setting,family.oi.Z)),family.oi.Z)
               for setting in ("phase_0","phase_pi","phase_pi_over_2","diagonal_control")}
        assert probs["phase_0"]+probs["phase_pi"]==family.S(F(3,4))
        assert probs["diagonal_control"]==family.S(F(3,8))
        values[name]={setting:str(value) for setting,value in probs.items()}
    # A prospective, separately labelled illustrative Hoeffding design.
    # Two independent settings, 100,000 attempts each, per-setting b=.0091.
    n,b=100000,sp.Rational(91,10000)
    alpha=4*math.exp(-2*n*float(b)**2)
    assert alpha < 2.87e-7
    # A .01 independently certified physical error allowance per setting.
    radius=2*(b+sp.Rational(1,100))
    assert radius==sp.Rational(191,5000)
    return {"sum_rule":"p(0)+p(pi)=3/4; p(pi/2)+p(3pi/2)=3/4",
            "disk_rule":"[p(0)-p(pi)]^2+[p(pi/2)-p(3pi/2)]^2 <= 1/2",
            "source_exact_examples":values,
            "illustrative_design":{"attempts_per_setting":n,"settings":2,
                "per_setting_Hoeffding_radius":str(b),"sampling_error_upper_bound":alpha,
                "per_setting_physical_probability_allowance":"1/100",
                "sum_rule_total_radius":str(radius),
                "calibration_failure_probability_must_be_added":True,
                "actual_laboratory_data":False}}


def rank_witness_check():
    I=sp.eye(2)
    X=sp.Matrix([[0,1],[1,0]])
    Y=sp.Matrix([[0,-sp.I],[sp.I,0]])
    Z=sp.diag(1,-1)
    paulis={a+b:sp.kronecker_product(A,B) for (a,A),(b,B) in itertools.product(zip("IXYZ",[I,X,Y,Z]),repeat=2) if a+b!="II"}
    included=[key for key in paulis if not key.startswith("Z") or key=="ZI"]
    assert len(included)==12
    # The bracket of two of the twelve putative port generators is a thirteenth.
    ax=-sp.I*paulis["XI"]/2
    ayx=-sp.I*paulis["YX"]/2
    bracket=ax*ayx-ayx*ax
    assert bracket==-sp.I*paulis["ZX"]/2
    columns=included+["ZX"]
    rows=[]
    settings=[]
    for target in columns:
        P=paulis[target]
        found=None
        for aname,A in paulis.items():
            for bname,B in paulis.items():
                if B*A-A*B==2*sp.I*P:
                    found=(aname,A,bname,B)
                    break
            if found:
                break
        assert found is not None
        aname,A,bname,B=found
        rho=(sp.eye(4)+B)/4
        effect=(sp.eye(4)+A)/2
        assert sp.trace(rho)==1
        assert all(t>=0 for t in rho.eigenvals())
        row=[]
        for key in columns:
            G=-sp.I*paulis[key]/2
            row.append(sp.simplify(sp.trace(effect*(G*rho-rho*G))))
        rows.append(row)
        settings.append({"target":target,"preparation_Pauli":bname,"readout_Pauli":aname})
    J=sp.Matrix(rows)
    assert J==sp.eye(13)/2
    assert J.rank()==13
    # Any exact twelve-dimensional generator factorization obeys rank <= 12.
    rng=np.random.default_rng(20261102)
    left=rng.normal(size=(19,12)); right=rng.normal(size=(12,15))
    positive=left@right
    sv=np.linalg.svd(positive,compute_uv=False)
    assert sv[12]<1e-12
    # Finite response rank is not a Lie-dimension witness: one U(1) generator
    # can have many Fourier harmonics. Four generic sine/cosine columns suffice.
    angles=np.array([.13,.37,.79,1.21])
    finite=np.column_stack([np.cos(angles),np.sin(angles),np.cos(3*angles),np.sin(3*angles)])
    assert np.linalg.matrix_rank(finite)==4
    return {"A1_consequence":"infinitesimal port-generated response J factors through R^12, so rank(J)<=12",
            "rejection_rule":"sigma_13(J_observed)>epsilon when ||J_observed-J_true||_2<=epsilon",
            "exact_negative_control":"12 selected SU(4) Pauli generators plus their XI/YX bracket ZX",
            "negative_control_response_matrix":"I13/2", "exact_negative_control_rank":int(J.rank()),
            "negative_control_minimum_singular_value":"1/2",
            "tomography_settings":settings,
            "synthetic_rank12_factorization_thirteenth_singular_value":float(sv[12]),
            "finite_angle_rank_countercontrol":4,
            "laboratory_single_carrier_identification_established":False,
            "proof_status":"elementary factorization argument with exact finite controls; no Lean formalization"}


def main():
    paths=["docs/AXIOM_REFERENCE.md","model/record_laws/README.md","model/record_laws/family.py",
           "model/record_laws/INSTRUMENT_CONTRACT.md","model/observer/native.py",
           "model/observer/source_discrimination/README.md","model/observer/source_discrimination/DERIVATION.md",
           "model/observer_time_fields/REFINEMENT.md","model/physical_space_time/OCCUPATION_FIELD.md",
           "model/relative_observer_geometry/README.md","model/relativistic_scalar/README.md",
           "model/continuum_free_limit/README.md","docs/OPH_EXECUTABLE_HYPOTHESIS.md",
           "model/observer_time_fields/test_continuum_dynamics.py","model/physical_space_time/test_observation_process.py",
           "model/moving_defect/README.md","model/moving_defect/transport.py",
           "model/defect_field_propagation/README.md",
           "model/relative_observer_geometry/test_construction.py",
           "model/relative_observer_geometry/test_clock_reconstruction.py",
           "model/relative_observer_geometry/test_temporal_records.py",
           "model/moving_defect/test_transport.py",
           "Lean/Geometry/CanonicalLorentzModule.lean","Lean/Geometry/EventGermDisplacement.lean",
           "Lean/Time/WorldlineRealization.lean","Lean/Time/ProperTimeCalibration.lean",
           "Lean/Time/ClockComparison.lean"]
    result={"schema":"oph.ipi.direct_tests_audit.v1","date":"2026-11-02",
            "population":population_check(),"interference":interference_check(),
            "response_rank":rank_witness_check(),
            "input_sha256":{p:hashlib.sha256((REPO/p).read_bytes()).hexdigest() for p in paths},
            "new_script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "status":"PASS_DATA_FREE_CONTROLS_ONLY"}
    target=HERE/"direct_tests_results.json"
    target.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],"output":str(target),
                      "native_branches":result["population"]["native_branch_checks"],
                      "rank_witness":result["response_rank"]["exact_negative_control_rank"]}))


if __name__=="__main__":
    main()
