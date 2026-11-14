"""A declared all-left Weyl/Higgs action using exact public M15 tensors.

The internal representation and Yukawa tensor identities are inherited and
rechecked from retained public artifacts. Native reference log odds and its
covariant family feature supply a chosen trial coupling map. No observed
particle masses or couplings enter; no selected physical vacuum is claimed.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from ..particle_identification import family_attachment
from ...native_modular_reference import local_preparation
from ...observer import native


ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "model/observer_particles/results"
SOURCE_FILES = (
    "model/observer_particles/results/representation_identification_operators.json",
    "model/observer_particles/results/representation_identification.json",
    "model/observer_particles/results/unified_yukawa_operators.json",
    "model/observer_particles/results/unified_yukawa_action.json",
    "model/observer_particles/representation_identification.py",
    "model/observer_particles/unified_yukawa_action.py",
    "Lean/Screen/WeylYukawaConventions.lean",
    "Lean/Screen/FermionSectorAssembly.lean",
)
CHANNELS = ("up", "down", "charged_lepton")
SPIN_EPSILON = s.ImmutableSparseMatrix([[0, 1], [-1, 0]])


def _clean(matrix):
    matrix = s.SparseMatrix(matrix)
    return s.ImmutableSparseMatrix(matrix.rows, matrix.cols,
                                   {key: s.expand(value) for key, value in matrix.todok().items()
                                    if s.expand(value) != 0})


def _zero(matrix):
    return not _clean(matrix).todok()


def _load_sparse(record):
    return _clean(s.SparseMatrix(*record["shape"],
                                 {(i,j):s.sympify(value) for i,j,value in record["entries"]}))


def _serialize(matrix):
    matrix = _clean(matrix)
    return {"shape": list(matrix.shape),
            "entries": [[i,j,str(value)] for (i,j),value in sorted(matrix.todok().items())]}


def _pinned(name, receipt_name):
    path = RESULTS / name
    receipt = json.loads((RESULTS / receipt_name).read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest() != receipt["artifacts"][name]:
        raise ArithmeticError(f"retained public artifact changed: {name}")
    return json.loads(path.read_text())


def _exterior_derivative(generator, subsets):
    """dGamma on the actual retained ordered exterior basis."""
    indices = {tuple(row): j for j,row in enumerate(subsets)}
    entries = {}
    for column, subset in enumerate(subsets):
        for position, old in enumerate(subset):
            rest = [x for x in subset if x != old]
            for new in range(5):
                coefficient = generator[new, old]
                if coefficient == 0 or new in rest:
                    continue
                destination = tuple(sorted([*rest,new]))
                insertion = sum(x < new for x in rest)
                key = indices[destination], column
                entries[key] = entries.get(key,0) + (-1)**(position+insertion)*coefficient
    return _clean(s.SparseMatrix(len(subsets),len(subsets),entries))


def source_coefficients():
    states, probabilities = native.states(binary=True), local_preparation.reference()
    levels = []
    for charge in range(3):
        weights = {p for x,p in zip(states,probabilities,strict=True)
                   if sum(x[i] for i in local_preparation.CORE)==charge}
        if len(weights) != 1:
            raise ArithmeticError("reference no longer has three core-charge levels")
        value = weights.pop()
        levels.append(s.Rational(value.numerator,value.denominator))
    w0,w1,w2 = levels
    return {"a": s.log(w1/w0), "b": s.log(w2/w1),
            "J": s.log(w1*w1/(w0*w2)), "Omega": s.log(w2/w0)}


@lru_cache(maxsize=1)
def _representation():
    packet = _pinned("representation_identification_operators.json", "representation_identification.json")
    subsets = tuple(tuple(packet["basis_subsets"][i]) for i in packet["even_indices"])
    if len(subsets) != 15 or {len(row) for row in subsets} != {2,4}:
        raise ArithmeticError("retained nonvacuum even exterior basis changed")
    g5 = tuple(_load_sparse(row) for row in packet["spin_lifted_one_body_generators"])
    currents = tuple(_exterior_derivative(g,subsets) for g in g5)
    operators = {key:_load_sparse(value) for key,value in packet["even_identification_operators"].items()}
    if len(currents) != 12 or any(not _zero(g.H+g) for g in currents):
        raise ArithmeticError("twelve anti-Hermitian source currents required")
    vertices = _pinned("unified_yukawa_operators.json", "unified_yukawa_action.json")["vertices"]
    channel_for_origin = {(-4,1):"up", (1,2):"down", (-3,6):"charged_lepton"}
    channels = {}
    for row in vertices:
        origin = tuple(map(int,row["constituent_q6Y"]))
        channel = channel_for_origin[origin]
        channels[channel] = {
            "origin": origin,
            "conjugate_Higgs": bool(row["uses_conjugate_of_same_Higgs_doublet"]),
            "components": tuple(_load_sparse(x) for x in row["coefficient_matrices"]),
        }
    if set(channels) != set(CHANNELS):
        raise ArithmeticError("the three exact charged-fermion invariant lines are required")
    return subsets,g5,currents,operators,channels


def generators(*, families=True):
    """Twelve exact anti-Hermitian internal source generators, in source order."""
    currents = _representation()[2]
    return tuple(_clean(s.kronecker_product(g,s.eye(3))) for g in currents) if families else currents


def electric_charge(*, families=True):
    """Q=T3+Y in the all-left charge-conjugated-singlet convention."""
    operators = _representation()[3]
    charge = _clean(operators["T3"] + operators["q"]/6)
    return _clean(s.kronecker_product(charge,s.eye(3))) if families else charge


def identification_operators(*, families=True):
    """Exact charge, Cartan and Casimir operators in the same component order."""
    operators = dict(_representation()[3])
    operators["Y"] = _clean(operators["q"]/6)
    operators["Q"] = electric_charge(families=False)
    return {key:_clean(s.kronecker_product(value,s.eye(3))) if families else value
            for key,value in operators.items()}


def spin_electric_charge():
    """Q90=I_spin tensor Q45, before Nambu doubling."""
    return _clean(s.kronecker_product(s.eye(2),electric_charge()))


def mode_labels():
    """Labels for all90 annihilators, without relabeling family axes as species."""
    components = ordering()["all_left_components"]
    return [{**row,"internal_index":row["index"],"family":family,"spin":spin,
             "index":spin*45+3*row["index"]+family}
            for spin in range(2) for row in components for family in range(3)]


@lru_cache(maxsize=1)
def family_feature():
    """Native F3 in its existing metric-orthogonal source-coordinate frame."""
    packet = family_attachment.report()
    rows = packet["compression"]["matrix_Qsqrt5"]
    return _clean(s.Matrix([[s.Rational(a)+s.Rational(b)*s.sqrt(5) for a,b in row] for row in rows]))


def yukawa_matrices():
    coefficients, feature = source_coefficients(), family_feature()
    return {"up":_clean(coefficients["a"]*feature),
            "down":_clean(coefficients["b"]*feature),
            "charged_lepton":_clean(coefficients["J"]*feature)}


def _vacuum(v):
    value = source_coefficients()["Omega"] if v is None else s.sympify(v)
    if value.free_symbols or value.is_real is not True or value.is_nonnegative is not True:
        raise ValueError("a declared real nonnegative neutral-Higgs component norm is required")
    return value


def mass_kernel(v=None, *, component=1):
    """Symmetric all-left K45 for a real basis component H=v e_component.

    The executed neutral choice is component1, H=(0,v). v defaults to Omega,
    an explicit scalar-vacuum architecture choice. The family trial matrices
    are symmetric, so their tensor product with each symmetric internal
    invariant component preserves the all-left transpose pairing.
    """
    if type(component) is not int or component not in (0,1):
        raise ValueError("Higgs component0 or1 required")
    value, ys = _vacuum(v), yukawa_matrices()
    channels = _representation()[4]
    return _clean(value*sum((s.kronecker_product(channels[name]["components"][component],ys[name])
                             for name in CHANNELS),s.zeros(45)))


def pairing_kernel(v=None):
    """Neutral-vacuum creation Delta=epsilon_spin tensor K45, dimension90.

    H_pair=1/2 sum Delta_ij c_i†c_j†+h.c. The equivalent all-left Lagrangian
    convention is +1/2 psi^T Delta* psi+h.c.; H_mass=-L_mass and CAR reorder
    signs are included. Neutral trial coefficients are real. For a generic
    complex Higgs, creation coefficients require conjugating its action kernel.
    """
    return _clean(s.kronecker_product(SPIN_EPSILON,mass_kernel(v)))


def ordering():
    """Exact one-generation subset/charge order and tensor index dictionary."""
    subsets,_,_,operators,_ = _representation()
    charge = electric_charge(families=False)
    labels = {-4:"u^c",1:"Q",2:"d^c",-3:"L",6:"e^c"}
    return {"all_left_components":[
        {"index":i,"exterior_subset":list(row),"internal_type":labels[int(operators["q"][i,i])],
         "q6Y":str(operators["q"][i,i]),"T3":str(operators["T3"][i,i]),"electric_Q":str(charge[i,i])}
        for i,row in enumerate(subsets)],
        "family_index": "existing native coordinate axes0,1,2; not electron/muon/tau labels",
        "internal_family_index": "k=3*internal_component+family (0<=k<45)",
        "massless_neutrino_internal_indices":[i for i in range(15)
                                             if operators["q"][i,i]==-3 and charge[i,i]==0],
        "massless_neutrino_internal_family_indices":[3*i+j for i in range(15) for j in range(3)
                                                    if operators["q"][i,i]==-3 and charge[i,i]==0],
        "spin_order": "spin outermost: spin*45+k (0<=index<90)",
        "Nambu_order": "(c(p), c†(-p)); this doubling is not an extra particle family"}


def matrices(v=None):
    """Exact immutable matrix API for the propagator/source adapters."""
    return {"generators5":_representation()[1],
            "generators15":generators(families=False),"generators45":generators(),
            "Q15":electric_charge(families=False),"Q45":electric_charge(),
            "Q90":spin_electric_charge(),
            "F3":family_feature(),"Yukawa":yukawa_matrices(),
            "K45":mass_kernel(v),"epsilon2":SPIN_EPSILON,"Delta90":pairing_kernel(v)}


def _identities(v=None):
    _,g5,currents,operators,channels = _representation()
    family_currents = generators()
    ys, count = yukawa_matrices(), 0
    for name in CHANNELS:
        row = channels[name]
        cs = row["components"]
        if any(c.T != c for c in cs):
            raise ArithmeticError("internal all-left coefficients must be symmetric")
        family_cs = tuple(_clean(s.kronecker_product(c,ys[name])) for c in cs)
        for source,g,gf in zip(g5,currents,family_currents,strict=True):
            h = source[3:,3:]
            if row["conjugate_Higgs"]:
                h = -h.T
            for component in (0,1):
                correction = sum((h[j,component]*cs[j] for j in (0,1)),s.zeros(15))
                if not _zero(g.T*cs[component]+cs[component]*g+correction):
                    raise ArithmeticError("retained current/Higgs invariant identity failed")
                family_correction = sum((h[j,component]*family_cs[j] for j in (0,1)),s.zeros(45))
                if not _zero(gf.T*family_cs[component]+family_cs[component]*gf+family_correction):
                    raise ArithmeticError("native family-coupled current/Higgs identity failed")
                count += 1
    y = operators["q"]/6
    h1,h2 = operators["H1"],operators["H2"]
    anomalies = {
        "grav2_Y":s.trace(y),"Y3":s.trace(y**3),
        "SU3_2_Y_Casimir_trace":s.trace(y*operators["color_C2"]),
        "SU2_2_Y_Casimir_trace":s.trace(y*operators["weak_C2"]),
        "SU3_Cartan_111":s.trace(h1**3),"SU3_Cartan_112":s.trace(h1*h1*h2),
        "SU3_Cartan_122":s.trace(h1*h2*h2),"SU3_Cartan_222":s.trace(h2**3)}
    if any(s.expand(value) != 0 for value in anomalies.values()):
        raise ArithmeticError("the all-left internal anomaly traces did not cancel")
    k = mass_kernel(v); delta = pairing_kernel(v); electric = electric_charge()
    if k.T != k or delta.T != -delta or not _zero(electric.T*k+k*electric):
        raise ArithmeticError("neutral-vacuum CAR or electromagnetic pairing failed")
    row_counts = [sum(k[i,j] != 0 for j in range(45)) for i in range(45)]
    if any(count > 1 for count in row_counts):
        raise ArithmeticError("neutral kernel lost its disjoint nonzero pairing blocks")
    rank = sum(row_counts)
    if rank != (0 if _vacuum(v)==0 else 42):
        raise ArithmeticError("neutral kernel rank or massless-neutrino rows changed")
    pauli = (s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    if any(not _zero(j.T*SPIN_EPSILON+SPIN_EPSILON*j)
           for p in pauli for j in (-s.I*p/2,-p/2)):
        raise ArithmeticError("declared left-Weyl epsilon invariant failed")
    return {"one_generation_current_Higgs_identities":count,
            "native_family_current_Higgs_identities":count,
            "symmetric_internal_mass_and_antisymmetric_spin_pairing":True,
            "neutral_vacuum_electric_pairing_exact":True,
            "left_Weyl_rotation_and_boost_epsilon_identities":6,
            "anomaly_traces_one_generation":{key:str(s.expand(value)) for key,value in anomalies.items()},
            "weak_doublets_one_generation":4,"weak_doublets_three_families":12,
            "neutral_mass_kernel_rank":rank,
            "minimal_massless_neutrino_components":3,
            "additional_zero_components_when_v_is_zero":42 if _vacuum(v)==0 else 0}


def report(v=None):
    """JSON-ready declared-action report; tree parameters are not pole masses."""
    value = _vacuum(v)
    coefficients, feature = source_coefficients(), family_feature()
    feature_values = sorted((feature[i,i] for i in range(3)),key=lambda x:float(x))
    if (s.exp(2*coefficients["a"]) >= s.exp(5*coefficients["b"])
            or s.simplify(3*feature_values[0]/2-feature_values[-1]).is_positive is not True):
        raise ArithmeticError("the exact trial hierarchy bound changed")
    rows = []
    for name,key in (("up","a"),("down","b"),("charged_lepton","J")):
        exact = [value*coefficients[key]*x for x in feature_values]
        rows.append({"channel":name,"tree_mass_parameters_exact":list(map(str,exact)),
                     "tree_mass_over_clock_float_diagnostic":[float(x/coefficients["Omega"]) for x in exact],
                     "physical_pole_mass_or_species_assignment":False})
    bindings = dict(family_attachment.report()["source_bindings"])
    for name in (*SOURCE_FILES,Path(__file__).relative_to(ROOT).as_posix()):
        bindings[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return {
        "schema":"oph.chiral_candidate.action.v1",
        "status":"EXACT_DECLARED_NATIVE_FAMILY_ALL_LEFT_HIGGS_ACTION",
        "source_coefficients_exact":{key:str(x) for key,x in coefficients.items()},
        "architecture":{"internal_module":"all-left M15, three declared family-band components",
                        "Yukawa_map":{"up":"a F_sigma","down":"b F_sigma","charged_lepton":"J F_sigma"},
                        "neutral_Higgs_component_norm_exact":str(value),
                        "default_v_equals_Omega_is_declared_choice":v is None,
                        "source_only_trial_defaults_used":v is None,
                        "canonical_kinetic_metric_I45_is_declared":True,
                        "one_same_Higgs_and_its_conjugate":True},
        "ordering":ordering(),"identities":_identities(v),
        "tree_parameters":rows,
        "mass_hierarchy":{"same_family_singular_value_ratios_in_all_three_channels":True,
                          "tree_mass_ratios_require_positive_v":True,
                          "tree_CKM_if_canonical_flavor_interpretation":"identity up to rephasings and ordering permutations",
                          "tree_flavor_CP_phase_generated":False,
                          "physical_CKM_or_PMNS_prediction":False,
                          "family_ratio_exact":str(feature_values[-1]/feature_values[0]),
                          "family_ratio_strictly_below_three_halves":True,
                          "all_charged_tree_parameters_max_over_min":str(coefficients["a"]*feature_values[-1]/(coefficients["b"]*feature_values[0])),
                          "all_charged_tree_parameters_ratio_strictly_below_fifteen_quarters":True,
                          "hierarchy_bound_reason":"a/b<5/2 follows (25/8)^2<(8/5)^5; max(F)/min(F)<3/2 by exact radical comparison",
                          "observed_charged_fermion_hierarchy_reproduced":False},
        "matrix_dimensions":{"M15":15,"internal_times_family":45,"left_spin_times_internal_family":90},
        "interpretation":{"no_empirical_mass_or_coupling_inputs":True,
                          "reference_is_exact_native_preparation_target":True,
                          "finite_native_preparation_error_propagated":False,
                          "physical_chiral_regulator_and_vacuum_selected":False,
                          "Yukawa_architecture_forced_by_repair":False,
                          "physical_Standard_Model_pole_masses_derived":False,
                          "quantum_anomaly_free_regulator_measure_constructed_here":False,
                          "right_handed_mirror_fields_added_to_make_mass":False,
                          "Nambu_doubling_counts_as_extra_matter":False,
                          "time_parameter":"declared coherent parameter; masses reported relative to native-reference clock frequency"},
        "source_bindings":dict(sorted(bindings.items()))}


if __name__ == "__main__":
    print(json.dumps(report(),indent=2,allow_nan=False))
