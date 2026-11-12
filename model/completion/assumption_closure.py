"""Source-bound assumption disclosure; no theorem or physical-gate promotion.

Run with no arguments to validate the retained inventory.  --write rebuilds its
source bindings after a scientific review; it does not alter any claim registry.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import subprocess
from functools import cache
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "model/completion/assumption_closure.json"
BASELINE = "44d42ff4e6f754c67fd6e67d1b3156d76ff3bf34"
SCHEMA = "oph.assumption_closure.v1"
FIELDS = ("state", "algebra", "operation", "reference", "clock", "domain")
CLASSES = {
    "proved_on_this_source", "imported_theorem_with_matched_hypotheses",
    "declared_architectural_primitive", "empirical_calibration",
    "open_implementation", "open_physical_map", "rejected_candidate",
}
PROVED = {"proved_on_this_source", "imported_theorem_with_matched_hypotheses"}
CLAIMS = (
    "OPH-MODEL-NATIVE-OBSERVER-EXECUTION", "OPH-LEPTON-GENERATION-REGISTER",
    "OPH-PARTICLE-DEFECT-INTEGER-LIFT", "OPH-PARTICLE-DEFECT-CENTRE-CLASSES",
    "OPH-GAUGE-PORT-PHASE-GLOBAL-FORM", "OPH-EW-SEAM-CLASS-BALANCE",
    "OPH-EW-LATTICE-CROSSOVER", "OPH-COLOUR-LATTICE-DECONFINEMENT-AND-BINDING",
    "OPH-ALPHA-PIXEL", "OPH-HIER-EW", "OPH-KOIDE-CIRCULANT-IDENTITY",
    "OPH-GR-D5", "OPH-GR-D4C-LOCAL-STRESS", "OPH-GR-D5A-ABSOLUTE-EINSTEIN",
    "OPH-GR-E2E-BRANCH-ENTRY", "OPH-FINITE-LOCALITY-NOSIGNALLING",
)
# Each tuple states a distinct mathematical domain. Equality of identifiers or
# numerical outputs supplies no transfer between these tuples.
DOMAIN_DATA = {
 "lepton_cycle": (
  "Uniform vector on 60 addresses, two fresh independent factors; arbitrary system/reference input.",
  "M_3 system and two 60-dimensional source factors with complete record/environment algebra.",
  "Householder preparation, oriented-current or mass controlled unitary, all four outcomes, SWAP archive.",
  "Arbitrary external quantum reference retained; full failure/environment coherences included.",
  "Declared dimensionless phase t and N fresh rounds; J spectrum {-1,0,1}.",
  "Finite declared instrument and finite histories; physical pole and autonomous service maps absent."),
 "record_observer": (
  "Known binary total-two preparation on twelve sites.",
  "Classical 66-state source, addressed seam/coin/current records and bounded controller.",
  "Uniform seam and fair coin; record-only state reconstruction with halt on full archive.",
  "Classical record distribution; immutable code/geometry and hardware cost excluded.",
  "Native attempted seam updates, including no-operation outcomes.",
  "One icosahedral binary sector and explicit bounded observer extension."),
 "population": (
  "All 1296 states; quoted plateau uses independent uniform loads 0..5 and zero initial flux.",
  "Classical three-port load and opposite Z6 face-register pair.",
  "Uniform seam and fair coin, with all 7776 attempted branches retained.",
  "Exact probability law on finite source states, no quantum thermal reference.",
  "One attempted native update per integer step.",
  "Closed three-port cell; no full-federation or cosmological identification."),
 "conjugation": (
  "Conjugation-even initial probability law.",
  "Finite native state functions with complement loads, negative class/lift, complemented coin.",
  "Conjugation-equivariant transition kernel and iteration.",
  "Same finite probability law; odd observables need a declared physical interpretation.",
  "Integer kernel iterations.",
  "Exact native local and closed-cell kernels, not a baryon asymmetry mechanism."),
 "charged_cube": (
  "Fixed compensated charged frame and retained reference with cost 4.",
  "Orthogonal Wilson cube at g=1/64, photon capacity 3, c=21/2, including native complement.",
  "Fixed finite Hamiltonian spectral enclosure and neutral-reference subtraction.",
  "Literal compensating reference with its energy included.",
  "This cube's Hamiltonian time and model mass m.",
  "Fixed finite source; no transfer to primary Wilson/Gauss family or physical pole."),
 "primary_wilson_gauss": (
  "Specified polynomial-selected finite state for each declared member.",
  "Open Wilson family, sixty CAR modes/cell, product link-word cutoff, Gauss ghost constraints.",
  "L_n=2^(2n+1), a_n=2^(-n), ell_n=n+1; electric105, magnetic1/7.",
  "Retained canonical conjugate reference and original energy zero.",
  "Declared primary Hamiltonian time.",
  "Finite primary members; physical joint refinement, vacuum and field attachment conditional."),
 "conditional_chiral_extension": (
  "Supplied common vacuum/state satisfying the extension's conditional hypotheses.",
  "Declared chiral gauge-Higgs extension with its regulator and matter content.",
  "Specified extension action, detector and refinement; no primary-program equivalence assumed.",
  "Supplied common reference, distinct from primary selected-state proof.",
  "Supplied common Hamiltonian/physical comparison clock.",
  "Conditional interacting chiral field reconstruction; physical realization unproved."),
 "gravity_branch": (
  "Supplied common-domain vacuum and fixed-cap variation family.",
  "Typed modular inclusion/stress/geometry tower with specified continuum domain.",
  "Conditional null-stress tomography and small-ball Einstein implication.",
  "Supplied faithful reference, entropy normalization and vacuum subtraction.",
  "Supplied compatible modular and calibrated physical clock.",
  "Conditional geometric branch; finite cap evidence does not supply branch nonemptiness."),
 "coupling_comparison": (
  "Declared pixel/screen branch and imported comparison inputs.",
  "Finite closure/endpoint maps and compared coupling observables.",
  "Declared transmutation, load and matching maps; empirical transport tracked separately.",
  "Specified closure normalization and empirical comparison scheme.",
  "External Planck/energy calibration where consumed.",
  "Conditional coupling chain, no inferred source-measure or physical endpoint identity."),
 "electroweak": (
  "Stored lattice peak estimates and marginal errors for declared run sets.",
  "Three-dimensional effective lattice theory and scalar comparison estimators.",
  "Finite-size weighted fit and linear covariance propagation under explicit covariance assumptions.",
  "Shared reference peaks; paired replicas/raw series unavailable in inspected exports.",
  "Externally calibrated GeV scale; no native cosmological clock.",
  "Finite numerical comparison; delta-method bounds conditional on marginal errors/linearization."),
 "species_dictionary": (
  "Nonnegative normalized distributions on declared 32 or separately restricted 30 patterns.",
  "Six centre-class records; separate boson-extended 48 or 46 component columns.",
  "Exact class aggregation and row-space identifiability/linear-program witnesses.",
  "Declared component probabilities; no physical species-state preparation.",
  "Static finite census; no thermal or cosmological time.",
  "Class labels and ambiguity; physical decoder, spin and generation attachments unproved."),
 "thermal_benchmark": (
  "External equilibrium ideal gas at 1 TeV, zero chemical potentials and massless species.",
  "90 fermionic and 28 bosonic helicity states; three active neutrino families, no populated RH neutrinos.",
  "Fermi/Bose equilibrium number-density comparison with declared degeneracy convention.",
  "External thermal state, not a native absorbing-class distribution.",
  "Temperature is an external parameter; no cooling history supplied.",
  "Ideal-gas comparison excludes masses and interaction precision corrections."),
 "generic_local_channel": (
  "Arbitrary finite bipartite state subject to the stated split.",
  "Named tensor factors and local matrix/record algebras.",
  "Normalized local Kraus map or row-stochastic map.",
  "Arbitrary retained spectator reference in its theorem domain.",
  "One specified local operation or fixed exogenous word.",
  "Conditional locality theorem; source region/factor and spacetime interpretations require maps."),
 "conditional_circulant": (
  "Positive-chamber C3 circulant with supplied balance and phase.",
  "Three-dimensional square-root carrier and declared tracial GNS balance.",
  "Algebraic eigenvalue identity and chosen squared-eigenvalue interpretation.",
  "Supplied tracial/balance reference, distinct from physical vacuum.",
  "No physical mass clock supplied by the identity.",
  "Exact conditional identity; target-informed phase diagnostic and physical mass map separated."),
 "quenched_lattice": (
  "Declared pure-gauge and quenched Wilson-quark lattice ensembles.",
  "Continuous port phases/links with declared cutoff and quark determinant omitted.",
  "Monte Carlo lattice action and measurement protocol.",
  "Declared lattice measure and comparison scale.",
  "Finite simulation updates and externally interpreted energy units.",
  "Finite ensemble measurements; source link-writing law and controlled physical limit unproved."),
 "native_claim_family": (
  "Separate states specified by each native headline evidence package; no joint preparation asserted.",
  "Disjoint family of finite native/classical/coherent/charged comparison carriers.",
  "Each declared package operation retains its own services and boundaries.",
  "Package-specific references, never pooled by matching labels or hashes.",
  "Package-specific update/phase/Hamiltonian clocks.",
  "Assumption disclosure for a multi-package claim, not one integrated executable source."),
}
DOMAINS = {name: dict(zip(FIELDS, values)) for name, values in DOMAIN_DATA.items()}
CLAIM_DOMAINS = {
 CLAIMS[0]: "native_claim_family", CLAIMS[1]: "lepton_cycle",
 CLAIMS[2]: "population", CLAIMS[3]: "species_dictionary",
 CLAIMS[4]: "species_dictionary", CLAIMS[5]: "coupling_comparison",
 CLAIMS[6]: "electroweak", CLAIMS[7]: "quenched_lattice",
 CLAIMS[8]: "coupling_comparison", CLAIMS[9]: "coupling_comparison",
 CLAIMS[10]: "conditional_circulant", CLAIMS[15]: "generic_local_channel",
 **{key: "gravity_branch" for key in CLAIMS[11:15]},
}
DOMAIN_OWNERS = {
 "lepton_cycle": "CP-02", "record_observer": "CP-01", "population": "CP-07",
 "conjugation": "CP-07", "charged_cube": "CP-04", "primary_wilson_gauss": "CP-06",
 "conditional_chiral_extension": "CP-06", "gravity_branch": "CP-09",
 "coupling_comparison": "CP-07", "electroweak": "CP-11", "species_dictionary": "CP-05",
 "thermal_benchmark": "CP-02", "generic_local_channel": "CP-13",
 "conditional_circulant": "CP-07", "quenched_lattice": "CP-07",
 "native_claim_family": "CP-01",
}
# These are input classifications. Conditional mathematical statements in a
# claim do not discharge the physical supply of their hypotheses.
CALIBRATIONS = {"externally_calibrated_Planck_energy_unit", "historically_target_informed_MCPR_diagnostic_coordinate", "same_scheme_hadronic_transport_boundary"}
OPEN_MAPS = {
 "target_independent_pixel_map_selection", "closure_same_invariant_bridge", "Maxwell_endpoint_lane",
 "COMMON_EW_LOAD_CARRIER", "source_endpoint_separation", "family_load_normalization",
 "direct_public_record_fixed_point", "finite_source_operator_and_region_factor_attachment_deferred",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def read_json_bytes(data: bytes) -> Any:
    return json.loads("\n".join(line for line in data.decode().splitlines() if not line.lstrip().startswith("#")))


@cache
def frozen_bytes(path: str) -> bytes:
    return subprocess.check_output(["git", "show", f"{BASELINE}:{path}"], cwd=ROOT)


@cache
def frozen(path: str) -> Any:
    return read_json_bytes(frozen_bytes(path))


def local_file(path: str) -> Path:
    if not isinstance(path, str) or not path or Path(path).is_absolute():
        raise ValueError("invalid relative source path")
    resolved = (ROOT / path).resolve()
    if not resolved.is_relative_to(ROOT) or not resolved.is_file():
        raise ValueError(f"source file missing or outside repository: {path}")
    return resolved


def obligation(owner: str, witness: str, counterexample: str, kind: str = "open_physical_map") -> dict:
    return {"classification": kind, "owner": owner, "smallest_witness": witness, "counterexample": counterexample}


def input_row(identifier: str, statement: str, classification: str, domain: str, witness: str,
              counterexample: str, *, owner: str | None = None, interpretation: bool = False) -> dict:
    owner = owner or DOMAIN_OWNERS[domain]
    return {
        "id": identifier, "statement": statement, "classification": classification,
        "source_domain": domain, "owner": owner, "source_specific_discharge": False,
        "independent_prediction": False, "physical_interpretation": interpretation,
        "disclosure_closed": True,
        "obligation": obligation(owner, witness, counterexample,
            "open_implementation" if classification in {"declared_architectural_primitive", "open_implementation"} and not interpretation else "open_physical_map"),
    }


# Services refine broad claim names without pretending they are canonical PR IDs.
# Every row supplies a bounded witness and an explicit way the proposed transfer can fail.
SERVICE_DATA = [
 ("fresh_joint_preparation", "Fresh sources are jointly independent of system, reference and earlier archives; correct one-source marginals are insufficient.", "declared_architectural_primitive", "lepton_cycle", "Implement one fresh source on the selected continuing source and certify joint trace distance to the product preparation.", "A correlated two-round preparation with correct marginals changes the complete channel."),
 ("preparation_analyzer", "The uniform-address Householder preparation and chosen analyzer are declared services.", "declared_architectural_primitive", "lepton_cycle", "Compile the specified preparation and analyzer into an allowed finite service with complete outcomes and bounded implementation error.", "Equal address marginals with changed phases produce different analyzer statistics."),
 ("controlled_kick", "The oriented controlled kick and mass-round coherent control are specified finite operations.", "declared_architectural_primitive", "lepton_cycle", "Implement one kick on the chosen source and bound arbitrary-reference channel distance, including every failure outcome.", "A classical current sampler cannot reproduce the off-diagonal controlled phase channel."),
 ("reset_archive_service", "A fresh-factor SWAP exports old correlations to an archive; indefinite archive production is a separate resource service.", "open_implementation", "lepton_cycle", "Supply one complete bounded cycle with work, controller and archive accounting and explicit full-capacity behavior.", "Discarding the correlated source or an uncounted fresh factor invalidates the cycle budget."),
 ("fixed_analyzer_error", "Reference trace-distance bounds require a common fixed analyzer and jointly controlled preparation; changing readout with the input is excluded.", "declared_architectural_primitive", "lepton_cycle", "Certify the joint source/reference error and identical readout channel in one comparison experiment.", "Input-dependent analyzer choice violates the claimed common-channel trace-distance comparison."),
 ("common_phase_clock", "Dimensionless current phase and mass-read phase need one implemented clock before their rates are compared.", "open_physical_map", "lepton_cycle", "Exhibit one intervention/readout experiment identifying both generators against the same retained clock.", "Independent rescaling of the two phase clocks changes their physical rate comparison."),
 ("squared_eigenvalue_mass", "Assigning squared finite eigenvalues to charged-lepton pole masses is a physical interpretation.", "open_physical_map", "lepton_cycle", "Derive one gauge-invariant spectral/detector observable with a controlled pole or threshold map and calibrated clock.", "An equally admissible nonlinear eigenvalue map changes the ratios while preserving the finite instrument."),
 ("lepton_calibration_ancestry", "Electron/Planck energy units and exposure to target mass ratios are recorded comparison inputs.", "empirical_calibration", "lepton_cycle", "Freeze a target-independent architecture and measurement protocol before an independent observable comparison.", "Fitting the phase or balance to a comparison mass invalidates independence."),
 ("coupling_measure", "Declared equal seam weights, sector Gibbs laws and continuous-link measures are distinct probability laws.", "open_physical_map", "coupling_comparison", "Produce the intended invariant measure from its actual proposal multiplicities and connect its readout to the compared coupling.", "One-sided proposal multiplicity or sector conservation produces a different stationary law."),
 ("continuous_link_supply", "Continuous port links and their action need a declared source service or a proved narrower implementation.", "open_implementation", "quenched_lattice", "Specify one link-writing operation with complete records and verify its claimed law.", "A Z6 centre label alone cannot distinguish two continuous links with different action."),
 ("species_decoder", "Six centre classes do not identify physical species, generations or Lorentz spin.", "open_physical_map", "species_dictionary", "Add the minimal separating record and prove one species observable is a function of the augmented record on the chosen source.", "Two nonnegative populations in the boson extension have identical six-class census and different lepton/Higgs/gauge totals."),
 ("component_mask", "The full 32-pattern carrier and restricted 30-pattern grammar are separate declarations.", "declared_architectural_primitive", "species_dictionary", "Declare the mask for one census and prove that its source preparation and operation preserve the selected carrier.", "Applying the restricted dictionary to the full carrier omits two allowed states."),
 ("population_preparation", "The 1007/13824 occupied limiting fraction consumes independent uniform loads 0..5 and zero initial flux.", "declared_architectural_primitive", "population", "Prepare that finite law with its independent factors or recompute the absorption law for a specified alternative.", "An absorbing initial state has a different occupied fraction."),
 ("population_federation", "Closed-cell absorption does not establish the same population law on an interacting federation.", "open_physical_map", "population", "Prove a finite federation projection intertwines its complete kernel with the closed-cell kernel or give a controlled approximation.", "A boundary exchange changes the load/flux transition law of the projected cell."),
 ("thermal_density_statistics", "External ideal-gas number densities use 1 TeV, zero chemical potentials, massless states, 90 fermionic and 28 bosonic helicities, active neutrinos and no populated right-handed neutrinos.", "declared_architectural_primitive", "thermal_benchmark", "Define one native energy, equilibrium measure and species decoder; compare its density observable with the stated ideal-gas convention.", "Using energy-density weighting instead of number-density weighting changes the ratio."),
 ("cosmological_history", "A thermal comparison needs a source clock, energy exchange and cooling/chemical history.", "open_physical_map", "thermal_benchmark", "Construct one finite energy-exchange/clock experiment with a declared expansion or cooling readout and controlled evolution.", "The same class census admits distinct temperatures and chemical potentials."),
 ("conjugation_preparation", "Zero odd mean requires conjugation-even preparation, fair paired coins and an equivariant scheduler.", "declared_architectural_primitive", "conjugation", "Verify these three inputs on the exact proposed source law, or quantify one explicit symmetry-breaking modification.", "An asymmetric initial law or unpaired coin changes the odd mean."),
 ("baryon_interpretation", "A native odd class observable is not supplied with a physical baryon-number assignment.", "open_physical_map", "conjugation", "Define the baryon observable on the physical species record algebra and intertwine it with the native update.", "Two species dictionaries assign different baryon numbers to the same native class."),
 ("record_known_preparation", "Record-only reconstruction assumes a known preparation and addressed seam, coin and current packets.", "declared_architectural_primitive", "record_observer", "Implement one complete packet stream from the declared source and verify the bounded controller without load access.", "Omitting an address or coin yields indistinguishable packets with different next states."),
 ("bounded_resources", "298 retained logical bits and 311 including the incoming packet exclude immutable code, geometry and physical hardware costs.", "open_implementation", "record_observer", "Account for one implemented finite controller/work/archive cycle including its full-capacity transition.", "An uncounted external address or archival channel breaks the stated implemented resource budget."),
 ("observer_agreement", "Compatible records need an operational inter-observer comparison and self-reading map on one construction.", "open_implementation", "native_claim_family", "Exhibit two bounded observers comparing one retained event and recovering the declared operational specification.", "A hidden external service changes a record without a recoverable shared specification."),
 ("primary_vacuum", "The selected polynomial primary state is not identified with the physical interacting vacuum.", "open_physical_map", "primary_wilson_gauss", "Prove one uniform selected-state approximation bound for the same-source vacuum and one observable norm.", "A fixed observable separates the selected state from every proposed vacuum approximation."),
 ("joint_chiral_field", "Joint noncollapse, chiral measure and mirror control must hold for the same declared interacting field family.", "open_physical_map", "conditional_chiral_extension", "Prove one interacting chiral regulator estimate with nonzero observable normalization and matched state/reference/clock.", "A mirror mode or vanishing connected correlation violates the proposed field interpretation."),
 ("gravity_kernel", "Conditional modular/null-stress identities need the same realized geometric tower and physical clock.", "open_physical_map", "gravity_branch", "Construct one tower member with all typed modular/stress inputs and a controlled common-domain limit.", "Separate valid cap and stress receipts use incompatible references or clocks."),
 ("gravity_entropy_scale", "Proper area, entropy stationarity and Newton normalization require matched physical scale and geometry inputs.", "open_physical_map", "gravity_branch", "Verify one small-ball area/entropy/stress comparison with independent scale and uniform remainder bound.", "An allowed entropy or area rescaling changes inferred Newton normalization."),
 ("full_mission_cosmology", "The audited finite particle set does not cover all cosmological, dark-sector or internal-closure observables.", "open_physical_map", "gravity_branch", "Choose one cosmological or dark-sector observable and specify its source record, preparation, scale and falsifier.", "Two compatible finite particle constructions disagree on that observable."),
 ("pure_yang_mills_scope", "A pure-gauge gap requires its own source/sector; a massless-photon theory has no positive global vacuum gap.", "open_physical_map", "primary_wilson_gauss", "Specify the pure-gauge instance and one matched sector/decay estimate with noncollapse.", "A photon excitation below every positive threshold refutes a full-theory positive gap claim."),
 ("independent_run_sets", "Electroweak covariance bounds with block structure assume different run sets are independent; marginal errors do not establish independence.", "declared_architectural_primitive", "electroweak", "Retain paired replicas and sampling provenance sufficient to estimate one joint covariance, or use the unrestricted PSD bound.", "A shared random or reweighting source gives nonzero cross-set covariance."),
 ("linear_covariance_budget", "The sharp PSD range applies to supplied marginal-error estimates and local linear influence; nonlinear remainder and coverage are separate.", "open_physical_map", "electroweak", "Retain one paired block sample and bound the nonlinear remainder of the declared estimator on its uncertainty domain.", "A curved estimator has a bias/remainder exceeding its linearized spread."),
 ("EW_scale_calibration", "Planck units, top/Higgs input maps, reduction allowance and finite fit-design sensitivity retain separate ancestry.", "empirical_calibration", "electroweak", "Specify one independent comparison with units, input covariances and controlled reduction error.", "Combining a fixed 1% allowance with statistical SD as calibrated coverage gives an unsupported confidence claim."),
]

SERVICE_DATA += [
 ("charged_cube_architecture", "The charged min-max theorem fixes g=1/64, photon capacity 3, compensated reference cost 4 and reducing native complement c=21/2.", "declared_architectural_primitive", "charged_cube", "Reconstruct the literal fixed Hamiltonian, frame and reducing complement before applying its enclosure.", "Changing c or the compensating reference changes the admissible spectral bound."),
 ("supplied_marginal_errors", "The covariance theorem consumes supplied marginal standard-error estimates; their sampling uncertainty is outside its algebraic variance range.", "empirical_calibration", "electroweak", "Retain one paired sample/error estimator and quantify uncertainty in its supplied marginals.", "A misestimated marginal variance invalidates interpretation of the plug-in range as calibrated uncertainty."),
 ("linear_influence_coefficients", "The covariance calculation uses fixed rational fit weights and numerical local reduction derivatives.", "declared_architectural_primitive", "electroweak", "Enclose the derivative remainder on one declared comparison interval and propagate its effect separately.", "Step-halving agreement alone does not bound every derivative truncation error."),
]

# Scope-limited mathematical witnesses; none is a canonical premise discharge.
PROOF_DATA = [
 ("complete_lepton_instrument", "proved_on_this_source", "lepton_cycle",
  "All four outcomes form a complete finite instrument; finite-cycle bounds retain arbitrary references and failures.",
  "model/lepton_cycle/PROOF.md", "## Repeated cycles, references and controlled limits",
  "model/lepton_cycle/test_cycle.py", ["fresh_joint_preparation", "preparation_analyzer", "controlled_kick", "reset_archive_service", "common_phase_clock"],
  {"generator": "weak J spectrum {-1,0,1}; full diameter 2", "time": "dimensionless t", "bound": "half diamond <= 14*t^2/(81*N); mass <= 349*t^2/(1296*N)", "references": "arbitrary untouched references and all outcome environments"}),
 ("fixed_analyzer_effect_error", "proved_on_this_source", "lepton_cycle",
  "For independent source states and one fixed analyzer, trace distance eta bounds effect operator norm by eta and each ordered eigenvalue by eta.",
  "model/lepton_cycle/COMPARISON.md", "A separate preparation theorem applies to a fixed analyzer projector",
  "model/lepton_cycle/test_compare.py", ["fixed_analyzer_error", "fresh_joint_preparation", "squared_eigenvalue_mass"],
  {"denominator": "eta < electron effect eigenvalue", "mode_order": "2*eta less than every spectral separation", "output": "limiting effect-rate ratio enclosure only"}),
 ("bounded_record_reconstruction", "proved_on_this_source", "record_observer",
  "Complete addressed packets reconstruct the finite state; exact quotient relaxation reaches worst TV <=1/100 first at 39 attempts.",
  "model/repair_sector_bridge/RECORD_OBSERVER.md", "## Exact relaxation of the recorded observable",
  "model/repair_sector_bridge/test_record_observer.py", ["record_known_preparation", "bounded_resources"],
  {"source": "66 binary total-two states", "records": "3960 addressed seam/coin branches", "memory": "298 retained logical bits, 311 including incoming packet; full archive halts"}),
 ("finite_population_absorption", "proved_on_this_source", "population",
  "The exact load-plus-register cost is nonincreasing; absorption and unresolved-mass bounds hold on the declared closed cell.",
  "model/native_abundance_epochs/PROOF.md", "## Abundance law and its limit",
  "model/native_abundance_epochs/test_epochs.py", ["population_preparation", "population_federation"],
  {"uniform_tail": "P(unresolved after n attempts) <= 25^(-floor(n/32)) for every initial distribution", "occupied_limit": "1007/13824 only for independent uniform loads 0..5 and initial flux zero", "meaning": "class occupancy, not physical species abundance or thermal equilibrium"}),
 ("finite_conjugation_evenness", "imported_theorem_with_matched_hypotheses", "conjugation",
  "Finite kernel equivariance and even preparation propagate evenness; odd bounded observable means vanish.",
  "Lean/Dynamics/ChargeConjugation.lean", "theorem iterate_kernelPush_conjugationEven",
  "model/charge_conjugation_dynamics/test_symmetry.py", ["conjugation_preparation", "baryon_interpretation"],
  {"finite_state": "native local finite states and 1296 closed-cell states", "permutation": "complement loads, negate registers/lifts, complement coin; involution", "kernel_equivariance": "all 2592 local and 7776 closed-cell branches checked in symmetry.py", "even_preparation": "uniform loads and paired fair coins", "odd_observable": "normalized class imbalance f(Cx)=-f(x)", "boundary": "executable premise checks plus generic Lean theorem, no Lean extraction of Python"}),
 ("fixed_charged_minmax", "proved_on_this_source", "charged_cube",
  "Fixed-source min-max bounds include the reducing native complement and neutral/reference subtraction.",
  "model/native_charged_cluster_minmax/PROOF.md", "The native complement is reducing",
  "model/native_charged_cluster_minmax/test_minmax.py", ["primary_vacuum"],
  {"completion": "c=21/2; c=23/2 is a separate helper", "gap": "full rank-two cluster gap >0.17175185 in model units", "transition_over_m": "[1.11134603,1.12278560]", "scope": "fixed compensated orthogonal cube, no primary-source or pole transfer"}),
 ("species_linear_identifiability", "proved_on_this_source", "species_dictionary",
  "A linear aggregate is identified exactly when its coefficients are constant within each observed class; classwise extrema attain sharp nonnegative bounds.",
  "model/native_abundance_epochs/SPECIES_PROOF.md", "## Identifiability and sharp bounds",
  "model/native_abundance_epochs/test_species_identifiability.py", ["species_decoder", "component_mask", "thermal_density_statistics"],
  {"matrices": "rank 6; 32/30/48/46 columns; nullities 26/24/42/40", "minimum_label": "one fermion/boson bit for four extended totals; three label values for all stored aggregates", "domain": "declared code-component census; native face records rejected"}),
 ("linear_covariance_extrema", "proved_on_this_source", "electroweak",
  "For fixed linear coefficients and marginal errors, PSD covariance extrema are attained by aligned or closed-polygon Gram vectors.",
  "model/electroweak_crossover/PROOF.md", "## Covariance of the shared reference peaks",
  "model/electroweak_crossover/test_comparison_budget.py", ["independent_run_sets", "linear_covariance_budget", "EW_scale_calibration"],
  {"formula": "max(0,2*max(l)-sum(l))^2 <= variance <= sum(l)^2, l_i=abs(a_i)*s_i", "conditional_inputs": "supplied estimated marginal errors and local numerical influence coefficients", "scope": "delta-method variance range, not nonlinear coverage or a rigorous derivative remainder"}),
]

# An open physical transfer is not a hypothesis of a finite conditional theorem.
PROOF_HYPOTHESES = {
 "complete_lepton_instrument": ["fresh_joint_preparation", "preparation_analyzer", "controlled_kick"],
 "fixed_analyzer_effect_error": ["fresh_joint_preparation", "fixed_analyzer_error"],
 "bounded_record_reconstruction": ["record_known_preparation"],
 "finite_population_absorption": ["population_preparation"],
 "finite_conjugation_evenness": ["conjugation_preparation"],
 "fixed_charged_minmax": ["charged_cube_architecture"],
 "species_linear_identifiability": ["component_mask"],
 "linear_covariance_extrema": ["independent_run_sets", "supplied_marginal_errors", "linear_influence_coefficients"],
}

JOIN_DATA = [
 ("lepton_to_primary", "lepton_cycle", "primary_wilson_gauss", "CP-06", "Intertwine a complete lepton cycle and its reference/clock with one primary member, with a controlled error.", "Different preparation or clock changes the complete output law."),
 ("cube_to_primary", "charged_cube", "primary_wilson_gauss", "CP-06", "Map the literal charged operator, state, reference and energy subtraction into one primary member.", "A nonreducing complement or changed reference invalidates the transported min-max enclosure."),
 ("population_to_species", "population", "species_dictionary", "CP-05", "Implement a complete count decoder with efficiency and normalization, distinguishing flux faces from code quanta.", "A zero-flux face is compatible with several species or no detected quantum."),
 ("population_to_thermal", "population", "thermal_benchmark", "CP-02", "Identify energy, equilibrium measure and number-count apparatus on the same source.", "A finite absorbing law and a thermal gas have incompatible stationary probabilities."),
 ("primary_to_gravity", "primary_wilson_gauss", "gravity_branch", "CP-10", "Construct the modular/stress/geometry inputs on one primary state/reference/clock with controlled limit.", "Individually valid geometry and field receipts refer to different states."),
 ("extension_to_primary", "conditional_chiral_extension", "primary_wilson_gauss", "CP-06", "Prove operation and complete-record comparison if equivalence to the primary source is claimed.", "A changed chiral regulator or action gives distinct output correlations."),
]

# These scientific artifacts are admitted separately from the baseline claim
# input graph. Their proof boundaries forbid a hidden promotion of that graph.
BATCH = [
 {"id": "lepton_comparison", "path": "model/lepton_cycle/physical_comparison.json", "producer": "model/lepton_cycle/compare.py", "proof": "model/lepton_cycle/COMPARISON.md", "independent_check": "model/lepton_cycle/test_compare.py", "inherited_inputs": ["fresh_joint_preparation", "preparation_analyzer", "fixed_analyzer_error", "squared_eigenvalue_mass", "lepton_calibration_ancestry"], "physical_map_supplied": False, "scope": "Fixed architecture, exposed mass comparison and conditional implementation/effect sensitivity; no physical pole prediction."},
 {"id": "electroweak_covariance", "path": "model/electroweak_crossover/comparison_budget.json", "producer": "model/electroweak_crossover/comparison_budget.py", "proof": "model/electroweak_crossover/PROOF.md", "independent_check": "model/electroweak_crossover/test_comparison_budget.py", "inherited_inputs": ["independent_run_sets", "supplied_marginal_errors", "linear_influence_coefficients", "linear_covariance_budget", "EW_scale_calibration"], "physical_map_supplied": False, "scope": "Sharp conditional linear covariance range; no nonlinear coverage or continuum error theorem."},
 {"id": "species_identifiability", "path": "model/native_abundance_epochs/species_identifiability_results.json", "producer": "model/native_abundance_epochs/species_identifiability.py", "proof": "model/native_abundance_epochs/SPECIES_PROOF.md", "independent_check": "model/native_abundance_epochs/test_species_identifiability.py", "inherited_inputs": ["species_decoder", "component_mask", "thermal_density_statistics", "cosmological_history", "population_federation"], "physical_map_supplied": False, "scope": "Exact inverse-problem obstruction and external ideal-gas control; no native species decoder or cosmological thermal history."},
]
SOURCE_FILES = (
 "model/lepton_cycle/instrument.py", "model/lepton_cycle/composition.py", "model/lepton_cycle/sensitivity.py",
 "model/repair_sector_bridge/record_observer.py", "model/native_abundance_epochs/model.py",
 "model/defect_integer_lift/model.py", "model/charge_conjugation_dynamics/symmetry.py",
 "model/charge_conjugation_dynamics/PROOF.md", "model/native_charged_cluster_minmax/audit.py",
 "model/native_charged_cluster_minmax/input_provenance.json",
 "model/particle_common_contract/PROOF.md", "model/particle_common_contract/contract.py", "model/particle_common_contract/contract.json",
 "model/particle_joint_reconstruction/dependencies.json", "model/particle_joint_reconstruction/CANDIDATE.md",
 "model/defect_particle_classes/PROOF.md", "model/defect_particle_classes/model.py",
 "model/completion/assumption_closure.py", "model/completion/test_assumption_closure.py", "model/completion/ASSUMPTION_CLOSURE.md",
)


def service_rows() -> list[dict]:
    interpretations = {"squared_eigenvalue_mass", "species_decoder", "thermal_density_statistics", "baryon_interpretation", "common_phase_clock"}
    owners = {
        "common_phase_clock": "CP-04", "squared_eigenvalue_mass": "CP-07",
        "lepton_calibration_ancestry": "CP-07", "observer_agreement": "CP-12",
        "population_federation": "CP-02", "population_preparation": "CP-02",
        "full_mission_cosmology": "CP-11", "gravity_kernel": "CP-09",
        "species_decoder": "CP-05", "component_mask": "CP-05",
    }
    return [input_row(*row, owner=owners.get(row[0]), interpretation=row[0] in interpretations) for row in SERVICE_DATA]


def proof_rows() -> list[dict]:
    return [{
        "id": identifier, "classification": category, "source_domain": domain,
        "interface": copy.deepcopy(DOMAINS[domain]), "statement": statement,
        "proof": {"path": proof, "anchor": anchor}, "independent_check": checker,
        "hypothesis_inputs": PROOF_HYPOTHESES[identifier], "downstream_obligations": [name for name in inputs if name not in PROOF_HYPOTHESES[identifier]], "matched_hypotheses": hypotheses,
        "source_specific_discharge": True, "canonical_premises_discharged": [],
        "physical_map_supplied": False,
    } for identifier, category, domain, statement, proof, anchor, checker, inputs, hypotheses in PROOF_DATA]


def classification_for_claim(claim_id: str, name: str) -> str:
    if name in CALIBRATIONS:
        return "empirical_calibration"
    if claim_id in CLAIMS[11:15] or name in OPEN_MAPS:
        return "open_physical_map"
    return "declared_architectural_primitive"


def claim_rows(claims: list[dict]) -> list[dict]:
    rows = []
    for claim in claims:
        domain = CLAIM_DOMAINS[claim["claim_id"]]
        for name in claim["assumptions"]:
            category = classification_for_claim(claim["claim_id"], name)
            witness = (f"Instantiate the exact input '{name}' on one declared {domain} member; check the consuming statement's hypotheses and complete records. "
                       "For a physical interpretation, supply its observable map and controlled error; for a declared service, give an implementation only if claiming a narrower source realizes it.")
            counter = f"One member violating '{name}' while satisfying the other declared inputs invalidates the corresponding conditional application."
            row = input_row(f"{claim['claim_id']}::{name}", name, category, domain, witness, counter)
            row["claim_id"] = claim["claim_id"]
            row["canonical_consumption_inferred"] = False
            row["source_specific_discharge"] = False
            rows.append(row)
    return rows


EMPIRICAL_IMPORTS = {
 "CODATA_2022_electron_and_muon_masses", "PDG_2024_tau_mass",
 "PDG_2024_nucleon_and_rho_masses", "published_SU3_deconfinement_coupling_and_string_tension_ratio",
 "Standard_Model_lattice_crossover_temperature", "endpoint_of_the_first_order_electroweak_transition_line",
}


def imported_rows(claims: list[dict]) -> list[dict]:
    result = []
    for claim in claims:
        domain = CLAIM_DOMAINS[claim["claim_id"]]
        for name in claim["imported_results"]:
            category = "empirical_calibration" if name in EMPIRICAL_IMPORTS else "open_implementation"
            item = input_row(f"{claim['claim_id']}::import::{name}", name, category, domain,
                f"Name the exact statement or data version for '{name}' and check each applicability hypothesis against this claim's actual source or comparison estimator.",
                f"A violated hypothesis, changed source/clock or omitted data covariance prevents this application of '{name}'.")
            item.update({"claim_id": claim["claim_id"], "role": "empirical_comparison_or_scale_input" if category == "empirical_calibration" else "named_import_with_application_unreviewed",
                         "hypotheses_matched_here": False,
                         "boundary": "An import name alone does not discharge its source-specific applicability; reviewed finite applications are separately named in proved_inputs."})
            if category == "empirical_calibration":
                item["empirical_use"] = {
                    "roles": ["comparison_reference", "optional_display_calibration"] if name == "CODATA_2022_electron_and_muon_masses" else ["comparison_reference"],
                    "architecture_fit_inferred": False,
                    "note": "The category tracks empirical ancestry. Frozen comparison targets do not imply numerical fitting of architecture; the CODATA electron energy separately calibrates an optional MeV display. Retrospective target exposure is retained.",
                }
            result.append(item)
    return result


def canonical_rows(premises: dict, ownership: dict, claims: list[dict]) -> list[dict]:
    owners = {row["premise"]: row for row in ownership["canonical_ownership"]}
    result = []
    for premise in premises["rows"]:
        uses = []
        for claim in claims:
            for role in ("consumed", "open", "boundary"):
                if premise["id"] not in claim["premise_dependencies"][role]:
                    continue
                domain = CLAIM_DOMAINS[claim["claim_id"]]
                category = "declared_architectural_primitive" if role == "consumed" and premise["id"] in {"PR-13", "PR-17", "PR-37", "PR-44", "PR-59"} else "open_physical_map"
                uses.append({"claim_id": claim["claim_id"], "role": role, "consumed": role == "consumed",
                             "classification": category if role != "boundary" else None,
                             "source_domain": domain, "owner": DOMAIN_OWNERS[domain],
                             "source_specific_discharge": False,
                             "obligation": obligation(DOMAIN_OWNERS[domain],
                                 f"Supply one {domain} instance of '{premise['name']}' with the exact canonical statement, or prove the consuming theorem without this input.",
                                 f"A same-domain instance failing '{premise['name']}' refutes that source-specific application.")})
        result.append({"id": premise["id"], "canonical_contract": {key: copy.deepcopy(premise[key]) for key in ("name", "type", "disposition", "statement", "discharge")},
                       "canonical_row_sha256_at_revision": sha(canonical(premise)),
                       "owners": owners[premise["id"]]["owners"], "issues": owners[premise["id"]]["issues"],
                       "coverage": "selected_claim_edges" if uses else "not_consumed_in_audited_claims",
                       "uses": uses, "canonical_status_changed": False,
                       "non_consumption_scope": "Absence from this selected claim graph is not global irrelevance or a premise discharge."})
    return result


def joint_rows(graph: dict) -> list[dict]:
    consumers: dict[str, list[str]] = {}
    for node in graph["nodes"]:
        for name in node["conditional_assumptions"]:
            consumers.setdefault(name, []).append(node["id"])
    result = []
    for name, nodes in sorted(consumers.items()):
        domain = "primary_wilson_gauss" if name.startswith("pure_YM") else "conditional_chiral_extension"
        category = "declared_architectural_primitive" if name == "extension_carrier_and_action" else "open_physical_map"
        row = input_row(name, name, category, domain,
                        f"Exhibit one member satisfying '{name}' with common state, algebra, operation, reference, clock and domain, or a counterexample to the proposed transfer.",
                        f"One typed member satisfying the other graph inputs but violating '{name}' rejects this source join.")
        row["consumer_nodes"] = nodes
        result.append(row)
    return result


def build_inventory(*, bindings: bool = True) -> dict:
    registry = frozen("claims/claim_registry.yaml")
    by_id = {claim["claim_id"]: claim for claim in registry["claims"]}
    claims = [by_id[key] for key in CLAIMS]
    premises = frozen("tracking/premise_register.json")
    graph = frozen("model/particle_joint_reconstruction/dependencies.json")
    rows = service_rows()
    proof = proof_rows()
    frozen_paths = ("tracking/premise_register.json", "claims/axiom_registry.yaml", "model/completion/finite_attachment_map.json", "model/completion/finite_record_dynamics.json", "model/particle_joint_reconstruction/dependencies.json", "claims/claim_registry.yaml")
    value = {
        "schema": SCHEMA,
        "audit": {"baseline_revision": BASELINE,
            "scope": f"All {len(premises['rows'])} canonical rows, all conditional assumptions of the named joint graph, the explicit assumptions, imported results and premise edges of {len(claims)} selected baseline claim rows, {len(rows)} refined services/physical interpretations, and {len(proof)} reviewed finite mathematical applications. Wider claim evidence is listed as unreviewed pointers. This is disclosure/coverage on this set, not a consistency or completeness proof for all OPH claims.",
            "claim_inputs_frozen": True,
            "current_batch_policy": "Separately admitted receipts inherit explicit input obligations. Mutable claim prose/evidence and completion-program pointers are excluded from live hash dependencies.",
            "integrated_gate_accepted": False, "canonical_premises_discharged": [],
            "physical_comparison_count_added": 0, "common_source_join_established": False},
        "baseline_manifest": [{"path": path, "sha256_at_revision": sha(frozen_bytes(path))} for path in frozen_paths],
        "normative_premise_reference": {"path": "claims/axiom_registry.yaml", "use": "Consume the registered theorem-premise declarations in their own stated domains; this inventory supplies no derivation or unique-architecture claim."},
        "theorem_premise_references": [{"id": row["id"], "path": "claims/axiom_registry.yaml", "sha256_at_revision": sha(canonical(row)), "classification": "declared_architectural_primitive", "source_specific_discharge": False, "scope": "Registered theorem-premise declaration consumed by reference; no realization or uniqueness inferred."} for row in frozen("claims/axiom_registry.yaml")["axioms"]],
        "domains": copy.deepcopy(DOMAINS),
        "claim_manifest": [{"claim_id": c["claim_id"], "status_at_revision": c["status"], "claim_class": c["claim_class"],
                            "assumptions": c["assumptions"], "imported_results": c["imported_results"], "premise_dependencies": c["premise_dependencies"],
                            "whole_row_sha256_at_revision": sha(canonical(c)), "unreviewed_evidence_pointers": c["evidence"]} for c in claims],
        "canonical_rows": canonical_rows(premises, frozen("model/completion/finite_attachment_map.json"), claims),
        "joint_inputs": joint_rows(graph),
        "joint_node_statuses_at_revision": {n["id"]: n["status"] for n in graph["nodes"]},
        "claim_inputs": claim_rows(claims), "imported_inputs": imported_rows(claims), "service_inputs": rows, "proved_inputs": proof,
        "source_join_obligations": [{"id": identifier, "from": source, "to": target, "join_established": False,
                                     "obligation": obligation(owner, witness, counter)} for identifier, source, target, owner, witness, counter in JOIN_DATA],
        "batch_admissions": copy.deepcopy(BATCH),
        "coverage_counts": {"canonical_rows": len(premises["rows"]), "joint_assumptions": len(joint_rows(graph)), "selected_claims": len(claims), "claim_inputs": len(claim_rows(claims)), "imported_inputs": len(imported_rows(claims)), "service_inputs": len(rows), "proved_applications": len(proof)},
    }
    if bindings:
        paths = set(SOURCE_FILES)
        for row in proof:
            paths.update((row["proof"]["path"], row["independent_check"]))
        for row in BATCH:
            paths.update(row[key] for key in ("path", "producer", "proof", "independent_check"))
        value["source_manifest"] = [{"path": path, "sha256": sha(local_file(path).read_bytes()),
                                     "role": "reviewed_finite_application_or_inventory_contract"} for path in sorted(paths)]
    return value


def validate(value: dict, *, verify_bindings: bool = True) -> dict:
    """Reject missing inputs, altered classifications, false joins and stale witnesses.

    The reviewed classification policy is this module's explicit input data.
    Equality with it is a disclosure invariant, not an automated proof judgment.
    """
    if not isinstance(value, dict) or value.get("schema") != SCHEMA:
        raise ValueError("assumption closure schema")
    expected = build_inventory(bindings=False)
    for key, wanted in expected.items():
        if value.get(key) != wanted:
            raise ValueError(f"assumption coverage/classification mismatch: {key}")
    if set(value) != set(expected) | {"source_manifest"}:
        raise ValueError("unexpected/missing inventory section")
    service_ids = {row["id"] for row in value["service_inputs"]}
    owners = {row["id"] for row in frozen("claims/completion_program.json")["stages"]}
    for row in value["claim_inputs"] + value["imported_inputs"] + value["service_inputs"] + value["joint_inputs"]:
        if row["owner"] not in owners or row["obligation"]["owner"] not in owners:
            raise ValueError("unknown completion owner")
        if row["classification"] not in CLASSES or row["classification"] in PROVED:
            raise ValueError("unproved input promoted or unclassified")
        if row["source_domain"] not in DOMAINS or row["source_specific_discharge"]:
            raise ValueError("input lacks domain or claims false discharge")
        if row["classification"] == "empirical_calibration" and row["independent_prediction"]:
            raise ValueError("calibration presented as independent prediction")
        if not all(row["obligation"].get(key) for key in ("owner", "smallest_witness", "counterexample")):
            raise ValueError("unresolved input lacks owner/witness/counterexample")
    for row in value["proved_inputs"]:
        if row["interface"] != DOMAINS[row["source_domain"]] or set(row["interface"]) != set(FIELDS):
            raise ValueError("source-specific proof lacks exact six-field interface")
        if not set(row["hypothesis_inputs"]) <= service_ids or not row["matched_hypotheses"]:
            raise ValueError("unmatched theorem hypotheses")
        proof = local_file(row["proof"]["path"])
        checker = local_file(row["independent_check"])
        if proof.suffix not in {".md", ".lean", ".tex"} or proof.samefile(checker):
            raise ValueError("receipt/checker substituted for mathematical proof")
        if row["proof"]["anchor"] not in proof.read_text():
            raise ValueError("mathematical witness anchor absent")
    paths = set(SOURCE_FILES)
    for row in value["proved_inputs"]:
        paths.update((row["proof"]["path"], row["independent_check"]))
    for row in BATCH:
        paths.update(row[key] for key in ("path", "producer", "proof", "independent_check"))
    manifest = value.get("source_manifest", [])
    if len(manifest) != len(paths) or {row.get("path") for row in manifest} != paths:
        raise ValueError("reviewed source/proof/checker manifest coverage")
    for row in manifest:
        path = local_file(row["path"])
        if set(row) != {"path", "sha256", "role"} or row["role"] != "reviewed_finite_application_or_inventory_contract":
            raise ValueError("source binding shape")
        if verify_bindings and sha(path.read_bytes()) != row["sha256"]:
            raise ValueError(f"stale reviewed source binding: {row['path']}")
    # Canonical statements/statuses are immutable under this inventory. Claim
    # prose and evidence may evolve without creating a circular hash dependency.
    # The generated model/theorem navigation (model_supply) follows the Lean and
    # claim inventories through tools/build_model_theorem_map.py.
    current = read_json_bytes((ROOT / "tracking/premise_register.json").read_bytes())
    retained = frozen("tracking/premise_register.json")
    if ({key: value for key, value in current.items() if key != "model_supply"}
            != {key: value for key, value in retained.items() if key != "model_supply"}):
        raise ValueError("canonical register changed: inventory requires explicit review")
    return copy.deepcopy(value["coverage_counts"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="write reviewed input inventory and source bindings")
    args = parser.parse_args()
    if args.write:
        value = build_inventory()
        validate(value)
        PATH.write_text(json.dumps(value, indent=2) + "\n")
    counts = validate(json.loads(PATH.read_text()))
    print(json.dumps({"schema": SCHEMA, "valid": True, "coverage": counts, "canonical_premises_discharged": [], "integrated_gate_accepted": False}, sort_keys=True))


if __name__ == "__main__":
    main()
