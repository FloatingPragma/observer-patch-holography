"""Full accepted native population response, projection memory and clock.

Repeating the existing completion operation until its first terminal timeout
is an explicitly declared continuation experiment. It is not a derived
autonomous matter law. No response rate is identified as a particle mass.
The complete source here is the66-state POPULATION instrument; coherent
matrix-unit response has its separate construction.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as s

from .. import slow_geometry
from ..hierarchy_diagnosis import native_family_dynamics as source
from ...native_modular_reference import local_preparation
from ...observer import native
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[3]))
    from tools.pin_scope import bound_pins


ROOT=Path(__file__).resolve().parents[3]
SOURCE_FILES=("model/common_source/full_repair_response/memory.py",
              "model/common_source/hierarchy_diagnosis/native_family_dynamics.py")


def _clean(matrix):
    # All source entries lie in Q(sqrt5). Rationalizing that number field is
    # sufficient; a general symbolic simplifier is substantially slower.
    return s.ImmutableSparseMatrix(s.SparseMatrix(matrix).applyfunc(lambda x:s.expand(s.radsimp(x))))


def _zero(matrix):
    return all(s.expand(x)==0 or s.simplify(x)==0 for x in matrix)


def _serialize(matrix):
    matrix=s.SparseMatrix(matrix)
    return {"shape":list(map(int,matrix.shape)),
            "entries":[[int(i),int(j),str(x)] for (i,j),x in sorted(matrix.todok().items())]}


@lru_cache(maxsize=1)
def _data():
    packet=source.report()
    tau=s.Rational(packet["complete_reduced_instrument"]["terminal_timeout_probability"])
    total=s.ImmutableSparseMatrix([[s.sympify(x) for x in row]
                                  for row in packet["full_preparation_context"]["actual_T12"]])
    accepted=total-tau*s.eye(12)
    state_total=source.proposal_operator()
    state_accepted=state_total-tau*s.eye(66)
    q=1-tau
    c=slow_geometry.source()
    basis=s.ImmutableMatrix([[s.Rational(a)+s.Rational(b)*s.sqrt(5) for a,b in row]
                            for row in slow_geometry.raw_frame()])
    gram=_clean(basis.T*basis)
    left=_clean(gram.inv()*basis.T)
    projector=_clean(basis*left)
    return {"packet":packet,"tau":tau,"q":q,"T":total,"A":s.ImmutableSparseMatrix(accepted),
            "T66":state_total,"A66":s.ImmutableSparseMatrix(state_accepted),
            "C":s.ImmutableSparseMatrix(accepted/q),"C66":s.ImmutableSparseMatrix(state_accepted/q),
            "B":basis,"left":left,"P":projector,"Q":s.ImmutableSparseMatrix(s.eye(12)-projector)}


def source_port_map():
    """One-proposal work marginal T=A+tauI; do not iterate past timeout."""
    return _data()["T"]


def accepted_port_map():
    """Substochastic accepted-proposal first-moment transfer A12."""
    return _data()["A"]


def conditional_port_map():
    """C12=A12/q, conditioned on survival; q is retained separately."""
    return _data()["C"]


def source_state_map():
    return _data()["T66"]


def accepted_state_map():
    """Substochastic66-state population transfer before the first timeout."""
    return _data()["A66"]


def _map(kind):
    names={"source":"T","accepted":"A","conditional":"C"}
    if kind not in names:
        raise ValueError("map_kind must be source, accepted or conditional")
    return _data()[names[kind]]


def proposal_history(proposals,initial=None):
    """Live distribution and every separate terminal timeout block.

    Work vectors are unnormalized. The timeout at proposal j+1 has work
    tau*A66^j*initial. It never receives a later repair. Initial defaults to
    the actual native reference; arbitrary exact population laws are allowed.
    """
    if type(proposals) is not int or not 0<=proposals<=32:
        raise ValueError("a proposal horizon from0 through32 required")
    data=_data()
    if initial is None:
        initial=s.ImmutableMatrix(list(map(s.Rational,local_preparation.reference())))
    else:
        initial=s.ImmutableMatrix(initial)
    if initial.shape!=(66,1) or any(not x.is_Rational or x<0 for x in initial) or sum(initial)!=1:
        raise ValueError("an exact normalized66-state population column required")
    live=initial; failures=[]
    for number in range(proposals):
        failures.append({"accepted_before_timeout":number,"timeout_proposal":number+1,
                         "work":s.ImmutableMatrix(data["tau"]*live)})
        live=data["A66"]*live
    mass=sum(live)+sum(sum(row["work"]) for row in failures)
    if mass!=1:
        raise ArithmeticError("a terminal failure or live branch was lost")
    return {"live":s.ImmutableMatrix(live),"terminal_timeouts":tuple(failures),"total_mass":mass,
            "live_mass":sum(live),"parameter":"completion proposals, stopped on timeout"}


def memory_kernel(orders=4,*,map_kind="accepted"):
    """R and memory coefficients M_l=left*M*Q*(QMQ)^l*Q*M*B.

    For zero initial complement, y_(n+1)=R y_n+sum_(j=0)^(n-1)
    M_(n-1-j)y_j. Initial complement adds its explicit inhomogeneous term.
    Full source and complement are retained in the returned matrices.
    """
    if type(orders) is not int or not 0<=orders<=12:
        raise ValueError("a memory order from0 through12 required")
    d=_data();m=_map(map_kind);q=d["Q"]
    r=_clean(d["left"]*m*d["B"])
    v=_clean(q*m*d["B"]);right=_clean(d["left"]*m*q)
    complement=_clean(q*m*q)
    coefficients=[];current=v
    for _ in range(orders):
        coefficients.append(_clean(right*current))
        current=_clean(complement*current)
    return {"instantaneous":r,"coefficients":tuple(coefficients),"complement_map":complement,
            "input_to_complement":v,"complement_to_output":right,
            "map_kind":map_kind}


def family_resolvent(z=2,*,map_kind="conditional"):
    """Exact rational-frequency family resolvent, retaining the full12-map."""
    z=s.sympify(z)
    if not z.is_Rational:
        raise ValueError("an exact rational resolvent point required")
    return _family_resolvent(z,map_kind)


@lru_cache(maxsize=12)
def _family_resolvent(z,map_kind):
    d=_data();matrix=z*s.eye(12)-_map(map_kind)
    try:
        inverse=matrix.inv()
    except (s.matrices.exceptions.NonInvertibleMatrixError,ZeroDivisionError) as exc:
        raise ValueError("resolvent point lies on a full-source pole") from exc
    return _clean(d["left"]*inverse*d["B"])


def schur_self_energy(z=2,*,map_kind="conditional"):
    """Exact Sigma(z)=zI-R-G_family(z)^(-1), the complement Schur term."""
    z=s.Rational(z);g=family_resolvent(z,map_kind=map_kind)
    d=_data();r=_clean(d["left"]*_map(map_kind)*d["B"])
    try:
        return _clean(z*s.eye(3)-r-g.inv())
    except s.matrices.exceptions.NonInvertibleMatrixError as exc:
        raise ValueError("the compressed resolvent is singular at this point") from exc


def _components(size,edges):
    neighbours=[set() for _ in range(size)]
    for a,b in edges:
        neighbours[a].add(b);neighbours[b].add(a)
    remaining=set(range(size));groups=[]
    while remaining:
        reached={min(remaining)};pending=list(reached)
        for a in pending:
            new=neighbours[a]-reached;reached|=new;pending.extend(sorted(new))
        groups.append(tuple(sorted(reached)));remaining-=reached
    return tuple(groups)


def _fixed_projector(groups,size):
    entries={}
    for group in groups:
        for i in group:
            for j in group:
                entries[i,j]=s.Rational(1,len(group))
    return s.ImmutableSparseMatrix(size,size,entries)


@lru_cache(maxsize=1)
def conserved_sectors():
    """Exact fixed sectors of symmetric conditional population maps."""
    d=_data();states=native.states(binary=True)
    edges12=[(i,j) for i in range(12) for j in range(i+1,12) if d["C"][i,j]!=0]
    groups12=_components(12,edges12)
    edges66=[(i,j) for (i,j),v in d["C66"].todok().items() if i<j and v!=0]
    groups66=_components(66,edges66)
    if sorted(map(len,groups12))!=[1,1,10] or sorted(map(len,groups66))!=[1,10,10,45]:
        raise ArithmeticError("protected core/outer connected sectors changed")
    fixed12=_fixed_projector(groups12,12);fixed66=_fixed_projector(groups66,66)
    if d["C"].T!=d["C"] or d["C66"].T!=d["C66"]:
        raise ArithmeticError("the full population response lost detailed balance")
    if d["C"]*fixed12!=fixed12 or d["C66"]*fixed66!=fixed66:
        raise ArithmeticError("claimed conserved projectors are not fixed")
    core=local_preparation.CORE
    if any(len({tuple(states[i][p] for p in core) for i in group})!=1 for group in groups66):
        raise ArithmeticError("the66 fixed sectors do not preserve each core occupation")
    residue=_clean(d["left"]*fixed12*d["B"])
    if residue.rank()!=2:
        raise ArithmeticError("the family response's conserved residue must have ranktwo")
    return {"groups12":groups12,"groups66":groups66,"P12":fixed12,"P66":fixed66,
            "family_fixed_residue":residue}


@lru_cache(maxsize=1)
def minimal_realization():
    """Closed exact minimal B-input/B-output realization of C12.

    Column elimination uses the existing source's exact Q(sqrt5) arithmetic.
    Reachability is minimal here because C12 is real symmetric and the
    observation is the transpose of the same input frame, up to its scalar
    Gram. The full12-map and its separately observable charge remain kept.
    """
    d=_data();c=slow_geometry.source();zero=c.q5(0)
    operator=[[F(d["C"][i,j]) for j in range(12)] for i in range(12)]
    initial=[list(column) for column in zip(*slow_geometry.raw_frame(),strict=True)]
    vectors=[];pivots=[]

    def multiply(vector):
        result=[]
        for row in operator:
            value=zero
            for a,x in zip(row,vector,strict=True):
                value=c.q5_add(value,c.q5_mul(c.q5(a),x))
            result.append(value)
        return result

    def reduce(vector):
        remaining=list(vector);coefficients=[]
        for pivot,basis in zip(pivots,vectors,strict=True):
            coefficient=remaining[pivot];coefficients.append(coefficient)
            remaining=[c.q5_sub(x,c.q5_mul(coefficient,y)) for x,y in zip(remaining,basis,strict=True)]
        return remaining,coefficients

    pending=list(initial)
    for vector in pending:
        residual,_=reduce(vector)
        pivot=next((i for i,x in enumerate(residual) if x!=zero),None)
        if pivot is None:
            continue
        inverse=c.q5_inv(residual[pivot])
        added=[c.q5_mul(inverse,x) for x in residual]
        vectors.append(added);pivots.append(pivot);pending.append(multiply(added))
        if len(vectors)>12:
            raise ArithmeticError("source reachability exceeded its full dimension")
    columns=[]
    for vector in vectors:
        residual,coordinate=reduce(multiply(vector))
        if any(x!=zero for x in residual):
            raise ArithmeticError("computed reachable space is not closed")
        columns.append(coordinate)
    inputs=[]
    for vector in initial:
        residual,coordinate=reduce(vector)
        if any(x!=zero for x in residual):
            raise ArithmeticError("source input was omitted from minimal realization")
        inputs.append(coordinate)

    def expression(x):
        return s.Rational(x[0])+s.Rational(x[1])*s.sqrt(5)
    embedding=s.ImmutableMatrix([[expression(vectors[j][i]) for j in range(len(vectors))] for i in range(12)])
    reduced=s.ImmutableMatrix([[expression(columns[j][i]) for j in range(len(columns))] for i in range(len(columns))])
    input_map=s.ImmutableMatrix([[expression(inputs[j][i]) for j in range(3)] for i in range(len(vectors))])
    output_map=_clean(d["left"]*embedding)
    if not _zero(d["C"]*embedding-embedding*reduced) or not _zero(embedding*input_map-d["B"]):
        raise ArithmeticError("closed minimal realization does not intertwine full source")
    if not _zero(s.ones(1,12)*embedding):
        raise ArithmeticError("a total-charge mode became reachable from zero-sum family inputs")
    return {"dimension":len(vectors),"embedding":embedding,"conditional_map":reduced,
            "input_map":input_map,"output_map":output_map,
            "unreachable_total_charge_retained_in_full_source":True}


@lru_cache(maxsize=1)
def incident_clock_maps():
    """Existing CAP0 next-incident stopping law, with terminal failures kept.

    Accepted incident J, accepted nonincident R and timeout tauI give
    S=J(I-R)^(-1), F=tau(I-R)^(-1). Both act on unnormalized live work.
    They retain all possible nonincident waits; a finite wait cut has its
    own retained unfinished block. No physical-clock interpretation follows.
    """
    d=_data();j12=s.zeros(12);r12=s.zeros(12);jentries=defaultdict(lambda:s.Rational(0));rentries=defaultdict(lambda:s.Rational(0))
    clock_cap=set(local_preparation.CAPS[0])
    selected=local_preparation.selectors(32)["completion_accepted_seam_weights"]
    for (a,b),weight in selected.items():
        endpoint=s.eye(12);endpoint[a,a]=endpoint[a,b]=endpoint[b,a]=endpoint[b,b]=s.Rational(1,2)
        if {a,b}&clock_cap:
            j12+=s.Rational(weight)*endpoint
        else:
            r12+=s.Rational(weight)*endpoint
    for before,flag,seam,coin,current,after,probability in source.proposal_branches():
        if flag!="accepted":continue
        target=jentries if set(seam)&clock_cap else rentries
        target[after,before]+=s.Rational(probability)
    j12,r12=s.ImmutableSparseMatrix(j12),s.ImmutableSparseMatrix(r12)
    j66=s.ImmutableSparseMatrix(66,66,dict(jentries));r66=s.ImmutableSparseMatrix(66,66,dict(rentries))
    p=sum(j12[:,0]);r=sum(r12[:,0])
    if p+r+d["tau"]!=1 or j66+r66!=d["A66"] or j12+r12!=d["A"]:
        raise ArithmeticError("incident clock split omitted accepted or terminal outcomes")
    waiting=(s.eye(12)-r12).inv()
    success=s.ImmutableMatrix(j12*waiting)
    failure=s.ImmutableMatrix(d["tau"]*waiting)
    survival=p/(p+d["tau"])
    ones=s.ones(1,12)
    if ones*(success+failure)!=ones or ones*success!=survival*ones:
        raise ArithmeticError("stopped clock success/failure map is not complete")
    loads=s.ImmutableMatrix(12,66,lambda i,k:native.states(binary=True)[k][i])
    if loads*j66!=j12*loads or loads*r66!=r12*loads:
        raise ArithmeticError("clock stopping does not intertwine all66 population states")
    sigma=s.ImmutableMatrix(list(map(s.Rational,local_preparation.reference())))
    if j66*sigma!=p*sigma or r66*sigma!=r*sigma:
        raise ArithmeticError("incident and nonincident branches do not preserve the same reference")
    return {"J12":j12,"R12":r12,"J66":j66,"R66":r66,"success12":success,
            "failure12":failure,"conditional_success12":s.ImmutableMatrix(success/survival),
            "incident_probability":p,"nonincident_probability":r,"tick_survival":survival}


def finite_tick_history(wait_cutoff,initial=None):
    """Full66-state first-incident history with finite waiting and unfinished output."""
    if type(wait_cutoff) is not int or not 0<=wait_cutoff<=32:
        raise ValueError("a finite wait cutoff from0 through32 required")
    data=_data();clock=incident_clock_maps()
    vector=s.ImmutableMatrix(list(map(s.Rational,local_preparation.reference()))) if initial is None else s.ImmutableMatrix(initial)
    if vector.shape!=(66,1) or any(not x.is_Rational or x<0 for x in vector) or sum(vector)!=1:
        raise ValueError("an exact normalized66-state population column required")
    successes=[];failures=[]
    for waits in range(wait_cutoff):
        successes.append({"nonincident_waits":waits,"work":s.ImmutableMatrix(clock["J66"]*vector)})
        failures.append({"nonincident_waits":waits,"work":s.ImmutableMatrix(data["tau"]*vector)})
        vector=clock["R66"]*vector
    if sum(vector)+sum(sum(x["work"]) for x in successes+failures)!=1:
        raise ArithmeticError("clock wait histories lost a terminal or unfinished block")
    return {"incident_successes":tuple(successes),"terminal_timeouts":tuple(failures),
            "unfinished":s.ImmutableMatrix(vector),"total_mass":s.Integer(1)}


def _spectrum(matrix,*,family_residues=False):
    numeric=np.array(matrix,dtype=float)
    values,vectors=np.linalg.eigh(numeric)
    if not np.isfinite(values).all() or np.min(values)<=0:
        raise ArithmeticError("population pole diagnostic lost its positive spectrum")
    d=_data();basis=np.array(d["B"],float);left=np.array(d["left"],float)
    groups=[]
    for i,value in enumerate(values):
        if groups and abs(value-groups[-1]["value"])<1e-11:
            groups[-1]["indices"].append(i)
        else:groups.append({"value":float(value),"indices":[i]})
    rows=[]
    for group in groups:
        ix=group["indices"];row={"pole":group["value"],"multiplicity":len(ix)}
        if family_residues:
            residue=left@vectors[:,ix]@vectors[:,ix].T@basis
            row["family_residue"]=residue.tolist();row["family_residue_trace"]=float(np.trace(residue))
        rows.append(row)
    nonfixed=values[values<1-1e-10]
    gaps=1-nonfixed;rates=-np.log(nonfixed)
    return {"dimension":len(values),"conserved_multiplicity_diagnostic":int(np.count_nonzero(values>1-1e-10)),
            "decaying_multiplicity_diagnostic":len(nonfixed),"poles":rows,
            "positive_gap_spread_diagnostic":float(max(gaps)/min(gaps)),
            "positive_log_decay_rate_spread_diagnostic":float(max(rates)/min(rates)),
            "pole_parameter":"z in (zI-C)^(-1); rate=-log(lambda) per accepted event",
            "floating_diagnostics_are_interval_certificates":False,"rates_are_particle_masses":False}


def _clock_poles(matrix):
    """Record and guard imaginary residuals before displaying real poles."""
    values=np.linalg.eigvals(matrix)
    imaginary=float(np.max(abs(values.imag)))
    if not np.isfinite(values).all() or imaginary>1e-10 or np.min(values.real)<=0:
        raise ArithmeticError("stopped-clock poles violate the declared positive-real diagnostic")
    return np.sort(values.real),imaginary


def report():
    d=_data();fixed=conserved_sectors();minimal=minimal_realization();clock=incident_clock_maps()
    kernels=memory_kernel(3,map_kind="conditional")
    g=family_resolvent(2);sigma=schur_self_energy(2)
    r=kernels["instantaneous"]
    if not _zero(g*(2*s.eye(3)-r-sigma)-s.eye(3)):
        raise ArithmeticError("full-family resolvent and complement Schur term disagree")
    history=proposal_history(3);tick=finite_tick_history(3)
    incident_numeric=np.array(clock["success12"],float)/float(clock["tick_survival"])
    vals,imaginary12=_clock_poles(incident_numeric);decay=vals[vals<1-1e-10]
    j66,r66=np.array(clock["J66"],float),np.array(clock["R66"],float)
    s66=np.linalg.solve(np.eye(66)-r66.T,j66.T).T/float(clock["tick_survival"])
    vals66,imaginary66=_clock_poles(s66);decay66=vals66[vals66<1-1e-10]
    lower12=2*min(d["C"].diagonal())-1
    lower66=2*min(d["C66"].diagonal())-1
    if lower12<=0 or lower66<=0:
        raise ArithmeticError("the exact Gershgorin positivity guard changed")
    bindings=dict(bound_pins(d["packet"]["source_bindings"]))
    for name in SOURCE_FILES:
        bindings[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return {
        "schema":"oph.full_repair_response.memory.v1",
        "source_scope":{"full66_population_map":True,"full4356_quantum_matrix_unit_map":False,
                        "seam_coin_and_current_are_unread_in_this_population_response":True,
                        "continuing_accepted_updates_until_timeout_is_declared_extension":True,
                        "autonomous_physical_matter_law_selected":False,"empirical_mass_or_coupling_inputs":[]},
        "accepted_dynamics":{"A12":_serialize(d["A"]),"conditional_C12":_serialize(d["C"]),
                             "terminal_timeout_probability_exact":str(d["tau"]),"survival_per_proposal_exact":str(d["q"]),
                             "live_after_n":"A66^n initial, mass q^n",
                             "timeout_after_j_accepted":"tau*A66^j initial, stopped forever with event label j",
                             "three_proposal_reference_live_mass":str(history["live_mass"]),
                             "three_proposal_total_retained_mass":str(history["total_mass"]),
                             "conditioning_is_only_a_diagnostic":True,
                             "hazard_pole_is_not_a_mass":"The q pole on conserved work modes is controller termination; changing selector rounds changes it."},
        "conserved":{"port_component_sizes":list(map(len,fixed["groups12"])),
                     "population_component_sizes":list(map(len,fixed["groups66"])),
                     "port_fixed_dimension":3,"population_fixed_dimension":4,
                     "labels":"n0,n1 are separately protected; the66-state map also preserves their joint occupation",
                     "family_conserved_pole_residue":_serialize(fixed["family_fixed_residue"]),
                     "family_conserved_residue_rank":2,
                     "reference_stationary_does_not_select_conserved_sector_weights":True},
        "exact_memory":{"conditional_R3":_serialize(r),
                        "coefficients_M0_M1_M2":list(map(_serialize,kernels["coefficients"])),
                        "complement_map_full12":_serialize(kernels["complement_map"]),
                        "complement_not_discarded":True,
                        "recurrence":"y[n+1]=R*y[n]+sum(j=0..n-1) M[n-1-j]*y[j] plus propagated initial-complement term",
                        "resolvent_formula":"G(z)=[zI-R-left*C*Q*(zI-QCQ)^(-1)*Q*C*B]^(-1)",
                        "exact_resolvent_at_z2":_serialize(g),"exact_Sigma_at_z2":_serialize(sigma),
                        "Schur_identity_at_z2_checked_exact":True},
        "closed_minimal_observable_realization":{"dimension":minimal["dimension"],
                                                "embedding12_by_r":_serialize(minimal["embedding"]),
                                                "conditional_map_r_by_r":_serialize(minimal["conditional_map"]),
                                                "accepted_map_r_by_r":_serialize(d["q"]*minimal["conditional_map"]),
                                                "input_r_by3":_serialize(minimal["input_map"]),
                                                "output3_by_r":_serialize(minimal["output_map"]),
                                                "exact_closure_and_input_intertwiners":True,
                                                "minimality_reason":"For real symmetric C, the reachable span of B equals the observable span of its transpose; independent reachable columns give a positive Gram.",
                                                "full12source_and_unreachable_global_charge_retained":True},
        "pole_diagnostics":{"exact_C12_positive_Gershgorin_lower":str(lower12),
                            "exact_C66_positive_Gershgorin_lower":str(lower66),
                            "accepted_tick_port_response":_spectrum(d["C"],family_residues=True),
                            "accepted_tick_full66_population_response":_spectrum(d["C66"])},
        "existing_incident_clock":{"cap":list(local_preparation.CAPS[0]),
                                   "J12_incident":_serialize(clock["J12"]),"R12_nonincident":_serialize(clock["R12"]),
                                   "next_tick_success12":_serialize(clock["success12"]),
                                   "terminal_failure12":_serialize(clock["failure12"]),
                                   "success_per_tick_exact":str(clock["tick_survival"]),
                                   "formula":"S=J(I-R)^(-1), failure=tau(I-R)^(-1); retained unfinished R^h at finite wait h",
                                   "infinite_wait_resolvent_is_stopping_limit":True,
                                   "original_finite_deadline_protocol_replayed":False,
                                   "same_reference_in_incident_and_nonincident_branches_exact":True,
                                   "first_moment_intertwiner_checked_on_all66_inputs":True,
                                   "finite_wait_cut3_total_retained_mass":str(tick["total_mass"]),
                                   "conditional_port_poles_numeric":vals.tolist(),
                                   "port_pole_maximum_imaginary_residual_diagnostic":imaginary12,
                                   "port_log_decay_rate_spread_diagnostic":float(max(-np.log(decay))/min(-np.log(decay))),
                                   "conditional_full66_poles_numeric":vals66.tolist(),
                                   "population_pole_maximum_imaginary_residual_diagnostic":imaginary66,
                                   "positive_real_pole_checks_are_interval_certificates":False,
                                   "real_similarity_reason":"S=J(I-R)^(-1) is similar to the positive symmetric (I-R)^(-1/2) J (I-R)^(-1/2).",
                                   "population_log_decay_rate_spread_diagnostic":float(max(-np.log(decay66))/min(-np.log(decay66))),
                                   "conditional_tick_count_is_physical_time":False},
        "mass_readout_scope":{"three_positive_mass_candidates_selected":False,
                              "full_response_poles_or_relaxation_rates_identified_as_charged_masses":False,
                              "coherent_propagation_charge_and_kinetic_bridge_supplied":False,
                              "correct_incident_clock_checked_instead_of_assuming_proposal_count":True,
                              "Standard_Model_mass_ratios_derived":False,
                              "failure_of_this_population_rate_route_excludes_all_OPH":False},
        "source_bindings":bound_pins(dict(sorted(bindings.items())))}


if __name__=="__main__":
    print(json.dumps(report(),indent=2,allow_nan=False))
