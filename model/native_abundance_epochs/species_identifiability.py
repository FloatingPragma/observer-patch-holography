"""Exact inverse problem for centre-class code counts and a declared thermal comparison.

Flux-face records have no admitted species decoder. The thermal model is an
external ideal-gas comparison, independent of the native attempted-event law.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from model.defect_particle_classes.model import fields, patterns

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / "species_identifiability_results.json"
GRAMMARS = ("full32", "restricted30")
AGGREGATES = ("quark", "lepton", "up", "down", "charged_lepton", "neutrino",
              "right_neutrino", "fermion", "higgs", "gauge", "boson")


def exact_vector(values, length, *, nonnegative=True):
    if len(values) != length or any(type(x) is not int and not isinstance(x, F) for x in values):
        raise ValueError("an exact integer/rational vector of the declared length is required")
    answer = tuple(F(x) for x in values)
    if nonnegative and any(x < 0 for x in answer):
        raise ValueError("negative populations are inadmissible")
    return answer


def components(grammar="full32", *, bosons=False):
    """Internal component menu; its columns do not supply physical helicities."""
    if grammar not in GRAMMARS or type(bosons) is not bool:
        raise ValueError("declare full32 or restricted30 and a Boolean boson extension")
    rows = []
    for p in patterns():
        if grammar == "restricted30" and p["mask"] in (0, 31):
            continue
        name = p["name"].removeprefix("anti-")
        family = "quark" if name[0] in ("u", "d") else "lepton"
        species = ("up" if name[0] == "u" else "down") if family == "quark" else (
            "neutrino" if name.startswith("nu") else "charged_lepton")
        rows.append(dict(key=f"mask:{p['mask']:02d}", mask=p["mask"],
                         degree=p["n3"] + p["n2"], k=p["k"], name=p["name"],
                         multiplet=p["multiplet"], kind="fermion", family=family,
                         species=species, right_neutrino=name == "nu_R"))
    if bosons:
        dimensions = {"H": 2, "anti-H": 2, "G": 8, "W": 3, "B": 1}
        for field in fields():
            if field["name"] not in dimensions:
                continue
            family = "higgs" if field["lorentz"] == "S" else "gauge"
            for index in range(dimensions[field["name"]]):
                rows.append(dict(key=f"{field['name']}:{index}", mask=None, degree=None,
                                 k=field["k"], name=field["name"], multiplet=field["name"],
                                 kind="boson", family=family, species=family, right_neutrino=False))
    return tuple(rows)


def groups(rows):
    if not rows or len({r["key"] for r in rows}) != len(rows):
        raise ValueError("distinct component identifiers are required")
    if any(type(r["k"]) is not int or not 0 <= r["k"] < 6 for r in rows):
        raise ValueError("component classes must be integers zero through five")
    return tuple(tuple(j for j, r in enumerate(rows) if r["k"] == k) for k in range(6))


def matrix(rows):
    grouped = groups(rows)
    return tuple(tuple(int(j in group) for j in range(len(rows))) for group in grouped)


def pushforward(rows, population):
    population = exact_vector(population, len(rows))
    return tuple(sum((population[j] for j in group), F(0)) for group in groups(rows))


def kernel_basis(rows):
    """e_j-e_p for each nonpivot column in each nonempty class, over Q."""
    answer = []
    for group in groups(rows):
        for j in group[1:]:
            vector = [0] * len(rows)
            vector[j], vector[group[0]] = 1, -1
            answer.append(tuple(vector))
    return tuple(answer)


def coefficients(rows, aggregate):
    if aggregate not in AGGREGATES:
        raise ValueError("unknown aggregate")
    return tuple(int(r["right_neutrino"] if aggregate == "right_neutrino" else
                     aggregate in (r["kind"], r["family"], r["species"])) for r in rows)


def sharp_bounds(rows, census, weights):
    """Attained extrema of w.x subject to x>=0 and A.x=c; no uniform prior."""
    census = exact_vector(census, 6)
    weights = exact_vector(weights, len(rows), nonnegative=False)
    low, high = [F(0)] * len(rows), [F(0)] * len(rows)
    lo_coeff, hi_coeff = [], []
    for k, group in enumerate(groups(rows)):
        if not group:
            if census[k]:
                raise ValueError("positive census in a class absent from the component menu")
            lo_coeff.append(None); hi_coeff.append(None)
            continue
        imin, imax = min(group, key=lambda j: weights[j]), max(group, key=lambda j: weights[j])
        low[imin], high[imax] = census[k], census[k]
        lo_coeff.append(weights[imin]); hi_coeff.append(weights[imax])
    lower = sum((a * b for a, b in zip(low, weights)), F(0))
    upper = sum((a * b for a, b in zip(high, weights)), F(0))
    return dict(lower=lower, upper=upper, lower_witness=low, upper_witness=high,
                lower_class_coefficients=lo_coeff, upper_class_coefficients=hi_coeff,
                identifiable_for_this_census=lower == upper,
                identifiable_for_every_census=lo_coeff == hi_coeff)


def infer_aggregates(rows, census, *, observable_kind):
    if observable_kind != "code_quantum_component_counts":
        raise ValueError("flux-face and thermal records require an independently proved species decoder")
    return {name: sharp_bounds(rows, census, coefficients(rows, name)) for name in AGGREGATES}


def minimal_extra_label(rows, aggregates):
    """Smallest deterministic per-quantum label alphabet for these linear aggregates.

    Labels may depend on the class, which the observer retains. This proves
    an information requirement; it supplies no physical readout instrument.
    """
    if not aggregates or len(set(aggregates)) != len(aggregates):
        raise ValueError("distinct nonempty aggregate names are required")
    weights = tuple(coefficients(rows, name) for name in aggregates)
    labels = [None] * len(rows)
    counts = []
    for group in groups(rows):
        signatures = []
        for j in group:
            signature = tuple(w[j] for w in weights)
            if signature not in signatures:
                signatures.append(signature)
            labels[j] = signatures.index(signature)
        counts.append(len(signatures))
    alphabet = max(counts)
    return dict(aggregates=tuple(aggregates), labels=labels,
                distinct_signatures_per_class=counts, minimum_alphabet_size=alphabet,
                minimum_fixed_binary_bits=(alphabet - 1).bit_length(),
                joint_class_label_rank=sum(counts), physical_instrument_implemented=False)


def thermal_benchmark():
    """One declared massless equilibrium SM comparison, with no thermal nu_R."""
    matter = components("restricted30")
    bosons = tuple(r for r in components("restricted30", bosons=True) if r["kind"] == "boson")
    physical = []
    # Each Weyl particle/antiparticle entry has one helicity; do not double it.
    for generation in range(3):
        for row in matter:
            physical.append(dict(row, key=f"generation:{generation}:{row['key']}", generation=generation,
                                 polarization="one Weyl particle or antiparticle helicity"))
    for row in bosons:
        polarizations = ("scalar",) if row["family"] == "higgs" else ("transverse:+", "transverse:-")
        for polarization in polarizations:
            physical.append(dict(row, key=f"{row['key']}:{polarization}", generation=None,
                                 polarization=polarization))
    quantities = {}
    for name, fermion_factor, unit in (("number", F(3, 4), "zeta(3)*T^3/pi^2"),
                                       ("energy", F(7, 8), "pi^2*T^4/30")):
        population = tuple(fermion_factor if row["kind"] == "fermion" else F(1) for row in physical)
        total = sum(population)
        law = pushforward(physical, population)
        aggregates = {a: sum((p * w for p, w in zip(population, coefficients(physical, a))), F(0))
                      for a in AGGREGATES}
        quantities[name] = dict(fermion_weight=fermion_factor, boson_weight=F(1), density_unit=unit,
                               effective_degrees=total, class_weights=law,
                               class_fractions=tuple(x / total for x in law),
                               aggregate_weights=aggregates,
                               aggregate_fractions={a: x / total for a, x in aggregates.items()},
                               quark_fraction_among_fermions=aggregates["quark"] / aggregates["fermion"])
    return dict(status="comparison_only_no_native_decoder", temperature_GeV=1000,
                phase="electroweak symmetric ideal-gas comparison",
                equilibrium="one common equilibrium temperature for all included species",
                chemical_potentials="all zero; no baryon or lepton asymmetry",
                mass_approximation="all masses set to zero; interaction and thermal-mass corrections omitted",
                neutrinos="three active Weyl neutrinos and their antiparticles; no populated right-handed neutrinos",
                helicity_convention="one helicity per Weyl particle/antiparticle, two transverse gauge polarizations, one per scalar particle/antiparticle component; no extra spin doublet",
                generation_count=3, physical_components=physical,
                degree_counts=dict(Counter(r["kind"] for r in physical)), quantities=quantities,
                precision_claim=False, simulated_temperature=False,
                reference=dict(title="Saikawa and Shirai, Primordial gravitational waves, precisely: The role of thermodynamics in the Standard Model",
                               url="https://arxiv.org/abs/1803.01038", doi="10.1088/1475-7516/2018/05/035",
                               use="Equilibrium Standard Model equation of state and the distinction between the ideal 106.75 energy count and interacting thermodynamics; T=1 TeV is a declared comparison setting."))


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def receipt():
    checks, cases = {}, {}
    for grammar in GRAMMARS:
        for bosons in (False, True):
            label = grammar + ("_bosons" if bosons else "")
            rows = components(grammar, bosons=bosons)
            a, basis = matrix(rows), kernel_basis(rows)
            multiplicities = pushforward(rows, [1] * len(rows))
            census = tuple(x / len(rows) for x in multiplicities)
            bounds = infer_aggregates(rows, census, observable_kind="code_quantum_component_counts")
            rank = sum(bool(g) for g in groups(rows))
            checks[f"{label}.partition_and_rank"] = all(sum(a[k][j] for k in range(6)) == 1 for j in range(len(rows))) and rank == 6
            checks[f"{label}.nullspace"] = len(basis) == len(rows) - rank and all(
                all(sum(row[j] * v[j] for j in range(len(rows))) == 0 for row in a) for v in basis)
            checks[f"{label}.sharp_bound_witnesses"] = all(
                pushforward(rows, b[side + "_witness"]) == census and
                sum(x * y for x, y in zip(b[side + "_witness"], coefficients(rows, name))) == b[side]
                for name, b in bounds.items() for side in ("lower", "upper"))
            uniform = [F(1, len(rows))] * len(rows)
            uniform_aggregates = {name: sum(x * y for x, y in zip(uniform, coefficients(rows, name)))
                                  for name in AGGREGATES}
            labels = {"family_totals": minimal_extra_label(rows, ("quark", "lepton", "fermion", "boson", "higgs", "gauge")),
                      "all_selected_aggregates": minimal_extra_label(rows, AGGREGATES)}
            checks[f"{label}.extra_labels_determine_aggregates"] = all(
                all(tuple(coefficients(rows, name)[i] for name in item["aggregates"]) ==
                    tuple(coefficients(rows, name)[j] for name in item["aggregates"])
                    for i in range(len(rows)) for j in range(len(rows))
                    if rows[i]["k"] == rows[j]["k"] and item["labels"][i] == item["labels"][j])
                for item in labels.values())
            cases[label] = dict(components=rows, pushforward_matrix=a, rank=rank, nullity=len(basis),
                               nullspace_basis=basis, class_multiplicities=multiplicities,
                               uniform_component_control=dict(census=census, aggregate_fractions=uniform_aggregates,
                                                              interpretation="declared equal component weights; no dynamical population law"),
                               aggregate_bounds_at_control_census=bounds, minimal_extra_labels=labels)

    examples = []
    for grammar, bosons, first, second, aggregate in (
            ("full32", False, "u_R", "d_R", "up"),
            ("full32", False, "e_R", "nu_R", "neutrino"),
            ("restricted30", False, "e_L", "nu_L", "neutrino"),
            ("full32", True, "e_R", "B", "fermion"),
            ("restricted30", True, "e_L", "H", "lepton")):
        rows = components(grammar, bosons=bosons)
        populations = []
        for name in (first, second):
            chosen = next(j for j, row in enumerate(rows) if row["name"] == name)
            populations.append(tuple(F(j == chosen) for j in range(len(rows))))
        laws = [pushforward(rows, p) for p in populations]
        values = [sum(x * y for x, y in zip(p, coefficients(rows, aggregate))) for p in populations]
        examples.append(dict(grammar=grammar, bosons=bosons, names=[first, second], populations=populations,
                             common_census=laws[0], aggregate=aggregate, aggregate_values=values))
        checks[f"indistinguishable.{first}.{second}"] = laws[0] == laws[1] and values[0] != values[1]

    thermal = thermal_benchmark()
    checks["restricted_grammar_excludes_exactly_neutral_degree_zero_and_five"] = {
        r["name"] for r in components() if r["mask"] in (0, 31)} == {"nu_R", "anti-nu_R"}
    checks["uniform_full32_count_is_24_to_8"] = cases["full32"]["uniform_component_control"]["aggregate_fractions"]["quark"] == F(3, 4)
    checks["uniform_restricted30_count_is_24_to_6"] = cases["restricted30"]["uniform_component_control"]["aggregate_fractions"]["quark"] == F(4, 5)
    checks["thermal_helicities_are_90_fermion_28_boson"] = thermal["degree_counts"] == {"fermion": 90, "boson": 28}
    checks["number_and_energy_weights_are_distinct"] = thermal["quantities"]["number"]["effective_degrees"] == F(191, 2) and thermal["quantities"]["energy"]["effective_degrees"] == F(427, 4)
    checks["thermal_quark_fraction_within_fermions_is_four_fifths"] = all(q["quark_fraction_among_fermions"] == F(4, 5) for q in thermal["quantities"].values())
    try:
        infer_aggregates(components(), (1, 0, 0, 0, 0, 0), observable_kind="native_flux_faces")
    except ValueError:
        checks["native_flux_faces_rejected_as_species_counts"] = True
    else:
        checks["native_flux_faces_rejected_as_species_counts"] = False

    inputs = ("model/native_abundance_epochs/species_identifiability.py",
              "model/native_abundance_epochs/test_species_identifiability.py",
              "model/native_abundance_epochs/SPECIES_PROOF.md",
              "model/native_abundance_epochs/SPECIES_README.md",
              "model/defect_particle_classes/model.py", "model/defect_particle_classes/exact.py",
              "model/defect_particle_classes/PROOF.md", "model/native_abundance_epochs/model.py",
              "model/native_abundance_epochs/PROOF.md", "Lean/Screen/ExteriorSelection.lean",
              "model/defect_integer_lift/model.py")
    return encode(dict(schema="oph.species_identifiability.v1", cases=cases,
        scope="Exact inverse problem for the declared code component menus. No species decoder for native flux-face records is supplied.",
        observable_domain=dict(native="oriented face-flux class, including zero flux; not a code-quantum count",
            conditional="one counted code component per record, correct class, common efficiency/normalization, retained zero-class events and all outcomes",
            absent_transfer="a source-bound preparation/readout instrument must implement those component labels and probabilities on actual run events; flux class zero cannot be equated with a gauge particle"),
        identified_linear_combinations="precisely coefficients constant within each populated class; over all censuses, the row space of the six-class matrix",
        indistinguishable_populations=examples, thermal_benchmark=thermal,
        checks=checks, all_checks_pass=all(checks.values()), physical_species_claim=False,
        input_sha256={name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in inputs}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = receipt()
    if not result["all_checks_pass"]:
        raise SystemExit("species-identifiability control failed")
    if args.write:
        RESULT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif args.check:
        if json.loads(RESULT.read_text()) != result:
            raise SystemExit("species-identifiability receipt differs from exact replay")
    else:
        print(json.dumps(result, indent=2, sort_keys=True))
    if args.write or args.check:
        print("Species-identifiability receipt verified; no physical species decoder admitted")


if __name__ == "__main__":
    main()
