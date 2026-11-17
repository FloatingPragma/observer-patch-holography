"""Source-selected entropy envelope with two conserved core moments.

The binary two-occupation shell has 66 states. For core I=(0,1), its
charge q has multiplicities g=(45,20,1); c=n0*n1=q(q-1)/2. The actual
mixed complementary-cap C/E channel P in common_caps preserves both.
Its established fixed algebra is A_I, with expectation E_I and gap>=1/600.
One additional unread native repair Q on seam (0,1) satisfies, exactly,

    Q E_I = E_I Q = Pi,   Pi(rho)=sum_q Tr(P_q rho) P_q/g_q.

All complementary seams commute with Q. Consequently Q P^n converges to
Pi with the SAME reference-uniform bound 33(1-1/600)^n, and one extra
native attempt. This removes the retained core matrix constraint: the
remaining input data are the two conserved scalar moments m=<q>, c=<n0 n1>.
It does not select their values. Global and incident-movement preparations
independently select two different, explicit points of this family.

For 0<c, 2c<m, m<1+c, put p=(1-m+c,m-2c,c), a_q=p_q/g_q.
sigma(m,c)=sum_q a_q P_q is faithful and is the unique entropy maximizer
at these two moments, for arbitrary coherent source states:

  S(rho)-S(sigma)=h1 delta m+h2 delta c-D(rho||sigma),
  -log(sigma)=(-log a0)I+h1 q+h2 c,
  h1=log(a0/a1), h2=log(a1^2/(a0 a2)).

Existence of the smooth MaxEnt family is explicit, not assumed. Its
entropy Hessian is -A^T diag(1/p) A, A=((-1,1),(1,-2),(0,1)); its
determinant is 1/(p0 p1 p2)>0. For any full-matrix fixed-moment tangent X,
dS_sigma[X]=0 and d²S_sigma[X,X]=-sum_ij |Xij|² L(ai,aj)<0 when X!=0,
where L(x,y)=(log x-log y)/(x-y), L(x,x)=1/x. Thus the curvature is
between -||X||_2²/min(a) and -||X||_2²/max(a).

Every cap R containing I has sigma in A_R. Its exact direct-sum entropy
is S_bulk(R)+sum_j p_R,j log(binomial(12-|R|,2-j))=S(sigma).
All 1024 cap choices give the SAME envelope; changing the discrete region
cancels its bulk/edge entropy change exactly at sigma. This is an algebraic
region change, not a differentiable metric, physical volume or area law.
For any rho with the same two moments, E_R rho has entropy deficit
D(E_R rho||sigma), and equality iff E_R rho=sigma. E_R is a readout
completion; native convergence is asserted only for the declared connected
mixed complement protocol, not every possible cap complement.

With an untouched reference B, Pi tensor id retains the conditional
reference states Tr_source(P_q rho_SB). It need not equal sigma tensor
rho_B. The entropy-maximization claim above concerns the source marginal.
The coherent unread theorem in this module includes no adaptive feedback or
conditioning on fine repair records. The attached history_entropy theorem
instead supplies a population envelope on complete native records, including
adaptive crossings, bounded stopping costs and an operational core census.
The unread marginal needs no physical tape erasure; the complete-record
accounting in record_thermodynamics applies with the entropy of the declared
independent choice word replacing n log(m) for a uniform choice alphabet.
The explicit first-law/envelope data instantiate the finite population
Lean contracts; the coherent variational proof is analytic, not Lean code.
Neither the conserved moments nor counted edge multiplicities are identified
with stress energy, metric variations, gravitational area or laboratory units.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import lcm
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from tools import strict_json
from . import construction as mixing, cut_completion as cut
from ..native_clock import regional_preparation as incident
from ..native_modular_reference import common_caps as common, local_preparation as local
from ..observer import native, quantum_mixing

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name("variation_results.json")
CORE = common.CORE
DIMENSIONS = (45, 20, 1)
SOURCE_PATHS = tuple(sorted(set((*local.SOURCE_PATHS, *common.SOURCE_PATHS,
    "model/native_equilibration/common_variation.py",
    "model/native_equilibration/record_thermodynamics.py",
    "model/native_equilibration/history_entropy.py",
    "Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean"))))


@lru_cache(maxsize=1)
def charges():
    result = tuple(sum(x[i] for i in CORE) for x in native.states(binary=True))
    if tuple(result.count(q) for q in range(3)) != DIMENSIONS:
        raise ArithmeticError("native core-charge multiplicities changed")
    return result


def _exact(value):
    if type(value) not in (int, F):
        raise ValueError("exact integer or rational moments required")
    return F(value)


def family(mean, pair):
    """An explicit faithful MaxEnt family, without fitting any physical data."""
    m, c = _exact(mean), _exact(pair)
    probabilities = (1-m+c, m-2*c, c)
    if min(probabilities) <= 0:
        raise ValueError("faithful domain requires 0<c, 2c<m, m<1+c")
    return tuple(probabilities[q]/DIMENSIONS[q] for q in charges())


def moments(weights):
    if (type(weights) is not tuple or len(weights) != 66
            or any(type(x) not in (int, F) or x < 0 for x in weights)
            or sum(weights) != 1):
        raise ValueError("66 exact normalized nonnegative populations required")
    return (sum((q*w for q, w in zip(charges(), weights, strict=True)), F()),
            sum((q*(q-1)//2*w for q, w in zip(charges(), weights, strict=True)), F()))


def family_certificate(mean, pair):
    mean, pair = _exact(mean), _exact(pair)
    w = family(mean, pair)
    p = (1-mean+pair, mean-2*pair, pair)
    a = tuple(pq/g for pq, g in zip(map(F, p), DIMENSIONS, strict=True))
    hessian = ((-1/p[0]-1/p[1], 1/p[0]+2/p[1]),
               (1/p[0]+2/p[1], -1/p[0]-4/p[1]-1/p[2]))
    return {"moments": moments(w), "charge_probabilities": p,
            "weights_per_source_state": a,
            "linear_multiplier_log_argument": a[0]/a[1],
            "pair_multiplier_log_argument": a[1]**2/(a[0]*a[2]),
            "entropy_hessian": hessian,
            "entropy_hessian_determinant": 1/(p[0]*p[1]*p[2]),
            "full_matrix_negative_hessian_HS_bounds": (1/max(w), 1/min(w))}


@lru_cache(maxsize=1)
def center_counts():
    """Integer representation of Pi, defined independently of the update."""
    scale = lcm(*DIMENSIONS)
    rows, cols, values = [], [], []
    for i, qi in enumerate(charges()):
        for j, qj in enumerate(charges()):
            if qi == qj:
                rows.append(67*i); cols.append(67*j)
                values.append(scale//DIMENSIONS[qi])
    return coo_matrix((values, (rows, cols)), shape=(4356, 4356), dtype=np.int64).tocsr(), scale


def source_audit():
    e, es = cut.completion_counts(CORE)
    z, zs = center_counts()
    if not cut._zero(z@z-zs*z) or not cut._zero(z-z.T):
        raise ArithmeticError("charge-sector expectation lost idempotence/symmetry")
    qvec = np.zeros(4356, np.int64)
    cvec = np.zeros(4356, np.int64)
    qvec[::67], cvec[::67] = charges(), [q*(q-1)//2 for q in charges()]
    _, weights, schedule = common.cap_schedule()
    result = {}
    for law in ("C", "E"):
        core = incident.seam_counts(CORE, law)
        if (not cut._zero(zs*(core@e)-2*es*z)
                or not cut._zero(zs*(e@core)-2*es*z)):
            raise ArithmeticError("one native core attempt does not realize Pi after E_I")
        for seam in (*weights, CORE):
            channel = incident.seam_counts(seam, law)
            if (not np.array_equal(channel.T@qvec, 2*qvec)
                    or not np.array_equal(channel.T@cvec, 2*cvec)):
                raise ArithmeticError("an actual repair changes a retained moment")
            if (not cut._zero(channel@core-core@channel)
                    or not cut._zero(channel@z-2*z)
                    or not cut._zero(z@channel-2*z)):
                raise ArithmeticError("native core/complement commutation or fixed center failed")
        p, ps = common.mixed_counts(law)
        if not cut._zero(core@p-p@core):
            raise ArithmeticError("source mixed complement channel does not commute with core repair")
        result[law] = {"matrix_units": 4356, "actual_seams": len(weights)+1,
                       "both_moments_conserved": True,
                       "one_core_attempt_completes_charge_expectation": True,
                       "fixed_center_dimension": 3}
    return {"laws": result, "complement_schedule": schedule,
            "fixed_core_algebra_dimension_before_core_attempt": 6}


def preparation_budget(precision_bits=24):
    """Uniform arbitrary-reference convergence to Pi(input), not to a product."""
    _, _, schedule = common.cap_schedule()
    b = mixing.preparation_budget(66, F(schedule["gap_lower_per_attempt"]), precision_bits)
    epsilon = F(b["trace_distance_upper"])
    return {"mixed_complement_attempts": b["attempts"], "final_core_attempts": 1,
            "native_work_attempts": b["attempts"]+1,
            "schedule": "iid fair cap choice, then uniform complementary seam; final fixed core seam",
            "auxiliary_selector_calls_in_work_subtotal": False,
            "reference_uniform_halftrace_upper": epsilon,
            "source_entropy_deficit_upper_nats": epsilon*(1+F(7, 10)*(precision_bits+9)),
            "target": "Pi tensor identity applied to the original input",
            "schedule_randomness_and_fine_records_in_distance": False}


def stationarity_budget(kind="local"):
    """One actual selected reference certifies all admitted cap tangents.

    Both existing preparations sandwich their normalized success marginal
    between (1-v)/(1+v)*sigma and its reciprocal. E_R and the final core
    channel are positive unital and fix sigma, so the same sandwich holds
    after them. Log monotonicity gives |dS[X]|<=2v/(1-v)||X||_1 for every
    normalized fixed-moment tangent, uniformly over all core-containing R.
    This is finite approximate stationarity, with every terminal flag kept
    in the preparer's complete bound; it is not exact finite invariance.
    """
    if kind == "local":
        b = local.budget()
        flagged = b["complete_flagged_reference_halftrace_upper"]
        count = b["total_native_attempts_upper"]+1
        work = b["work_native_attempts_upper"]+1
    elif kind == "global":
        b = common.preparation_budget()
        flagged = b["complete_flagged_reference_half_trace_upper"]
        count = None  # this ideal schedule does not budget auxiliary selectors
        work = b["warmup_native_attempts"]+b["global_movement_herald_attempts"]+b["mixed_cap_native_attempts_after_success"]+1
    else:
        raise ValueError("source reference must be local or global")
    v = b["operator_sandwich_radius"]
    return {"reference": kind, "operator_sandwich_radius": v,
            "entropy_tangent_error_per_trace_norm_upper": 2*v/(1-v),
            "complete_flagged_reference_halftrace_upper": flagged,
            "additional_core_native_attempts_on_success": 1,
            "native_work_attempts_per_trial_upper": work,
            "native_total_attempts_upper": count,
            "simultaneously_admitted_caps": 1024,
            "tangent_constraint": "Tr X=Tr(q X)=Tr(n0*n1 X)=0",
            "bound_parameter": "dimensionless entropy per declared variation parameter"}


def _log_add(form, coefficient, argument):
    """Canonical exact rational linear combination of logarithms of primes."""
    argument = F(argument)
    for integer, sign in ((argument.numerator, 1), (argument.denominator, -1)):
        prime = 2
        while prime*prime <= integer:
            while integer % prime == 0:
                form[prime] += sign*coefficient
                integer //= prime
            prime += 1
        if integer > 1:
            form[integer] += sign*coefficient


def _clean(form):
    return {str(prime): value for prime, value in sorted(form.items()) if value}


def region_entropy(cap, weights):
    """Exact cap-record and counted-complement entropies, in nats.

    Empty/full caps are allowed as readouts; no native-mixing certificate
    for their complements is inferred. Non-core caps are allowed as controls.
    """
    if (type(cap) is not tuple or any(type(i) is not int or not 0 <= i < 12 for i in cap)
            or tuple(sorted(set(cap))) != cap):
        raise ValueError("an increasing tuple of native ports is required")
    if weights != family(*moments(weights)):
        raise ValueError("a faithful core-charge reference from the explicit family is required")
    groups = defaultdict(list)
    for x, w in zip(native.states(binary=True), weights, strict=True):
        groups[tuple(x[i] for i in cap)].append(w)
    bulk, edge, total = defaultdict(F), defaultdict(F), defaultdict(F)
    for row in groups.values():
        probability = sum(row, F())
        _log_add(bulk, -probability, probability)
        _log_add(edge, probability, len(row))
    for w in weights:
        _log_add(total, -w, w)
    return {"bulk": _clean(bulk), "counted_edge": _clean(edge), "source": _clean(total),
            "completion": _clean({p: bulk[p]+edge[p] for p in set(bulk)|set(edge)})}


def all_regions(weights):
    checksum = hashlib.sha256()
    for mask in range(1 << 10):
        cap = (*CORE, *(p for bit, p in enumerate(range(2, 12)) if mask & (1 << bit)))
        row = region_entropy(cap, weights)
        if row["completion"] != row["source"]:
            raise ArithmeticError("bulk/edge cancellation failed on a core-containing region")
        checksum.update(json.dumps([cap, row], default=str, sort_keys=True).encode())
    return {"regions_checked": 1024, "exact_prime_log_ledger_sha256": checksum.hexdigest(),
            "core_region": region_entropy(CORE, weights),
            "four_port_region": region_entropy(common.CAPS[0], weights)}


def run():
    from . import history_entropy

    references = {}
    for kind, weights in (("global", common.reference()), ("local", local.reference())):
        m, c = moments(weights)
        if weights != family(m, c):
            raise ArithmeticError("the native selected reference is outside the explicit family")
        references[kind] = {"family": family_certificate(m, c), "region_envelope": all_regions(weights),
                            "finite_stationarity": stationarity_budget(kind)}
    # Mean charge alone agrees with the maximum-entropy tracial control.
    if moments(common.reference())[0] != moments((F(1, 66),)*66)[0]:
        raise ArithmeticError("mean-charge-only countercontrol changed")
    payload = {"schema": "oph.native_equilibration.common_variation.v1",
        "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        "native_source": source_audit(), "arbitrary_input_preparation": preparation_budget(),
        "selected_references": references,
        "continuing_record_entropy": history_entropy.run(),
        "countercontrols": {"global_and_tracial_mean_charge": F(1, 3),
            "global_pair_moment": F(1, 75), "tracial_pair_moment": F(1, 66),
            "missing_pair_constraint_allows_positive_entropy_derivative": True,
            "unprotected_cap": region_entropy((2, 3, 4, 6), common.reference())},
        "contracts": {"finite_population_first_law": "OPH.Thermodynamics.cap_firstLaw_exact",
            "explicit_unit_speed_moment_family": "OPH.EinsteinBranch.maxEnt_envelope_identity",
            "full_coherent_variational_result_machine_checked_in_Lean": False,
            "single_source_region_decompositions": True,
            "arbitrary_reference_decoupled": False,
            "geometric_area_or_metric_stationarity_identified": False,
            "physical_stress_energy_or_laboratory_units_identified": False}}
    return json.loads(json.dumps(payload, default=str))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    payload = run()
    if args.check:
        if strict_json.load(args.check) != payload:
            raise SystemExit("native common entropy-envelope receipt drift")
    else:
        with OUTPUT.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    print("native common entropy envelope: all checks passed")


if __name__ == "__main__":
    main()
