"""Additive topological-charge constraints from declared SM vertices.

Massless fields have zero defect charge in this hypothesis; antiparticles
carry the inverse. Each nonzero vertex conserves the charge. These premises
constrain a dictionary without constructing a defect. Neutrinos are massless.
"""
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

SPECIES = ("up", "charm", "top", "down", "strange", "bottom",
           "electron", "muon", "tau", "W-", "Z", "Higgs")


def relation_matrix(*, electroweak_instanton=False):
    """Present the universal abelian charge group by integer relations.

    All nine quark charged-current entries are nonzero. The optional row is
    the three-generation (udd nu)^3 electroweak 't Hooft vertex. Realizing
    this process in OPH is an additional premise.
    """
    if type(electroweak_instanton) is not bool:
        raise ValueError("explicit Boolean instanton hypothesis required")
    rows, labels = [], []

    def add(label, coefficients):
        rows.append([coefficients.get(name, 0) for name in SPECIES])
        labels.append(label)

    for up in SPECIES[:3]:
        for down in SPECIES[3:6]:
            add(f"charged-current vertex: {up}, {down}, W+", {up: 1, down: -1, "W-": 1})
    for lepton in SPECIES[6:9]:
        add(f"charged-current vertex: W-, {lepton}, anti-neutrino", {"W-": 1, lepton: -1})
    add("Z coupling to a fermion-antifermion pair", {"Z": 1})
    add("Higgs coupling to a charged conjugate pair", {"Higgs": 1})
    if electroweak_instanton:
        add("three-generation electroweak (udd nu)^3 vertex",
            {"up": 1, "charm": 1, "top": 1,
             "down": 2, "strange": 2, "bottom": 2})
    return s.ImmutableMatrix(rows), labels


def allowed_assignment(values, *, modulus=None, electroweak_instanton=False):
    """Check integer or Z/n charges against every declared vertex."""
    if type(values) is not dict or set(values) != set(SPECIES):
        raise ValueError("a complete species-to-integer charge dictionary required")
    if any(type(value) is not int for value in values.values()):
        raise ValueError("exact integer charge representatives required")
    if modulus is not None and (type(modulus) is not int or modulus < 2):
        raise ValueError("modulus must be an integer at least two")
    matrix, _ = relation_matrix(electroweak_instanton=electroweak_instanton)
    residual = matrix*s.Matrix([values[name] for name in SPECIES])
    return all(value == 0 if modulus is None else value % modulus == 0 for value in residual)


def _group_receipt(include_instanton):
    matrix, labels = relation_matrix(electroweak_instanton=include_instanton)
    smith = smith_normal_form(matrix, domain=s.ZZ)
    nonzero = [abs(int(smith[i, i])) for i in range(min(smith.shape)) if smith[i, i] != 0]
    return {
        "relations": labels,
        "integer_relation_matrix": [list(map(int, matrix.row(i))) for i in range(matrix.rows)],
        "rank": int(matrix.rank()),
        "universal_additive_charge_group_free_rank": len(SPECIES)-len(nonzero),
        "universal_additive_charge_group_torsion_orders": [d for d in nonzero if d > 1],
        "Smith_nonzero_diagonal_absolute_values": nonzero,
    }


def constraint_certificate():
    d, e = s.symbols("D E", integer=True)
    general = [d-e]*3+[d]*3+[e]*4+[s.S.Zero]*2
    matrix, _ = relation_matrix()
    if matrix*s.Matrix(general) != s.zeros(matrix.rows, 1) or matrix.rank() != 10:
        raise ArithmeticError("perturbative charge-kernel derivation failed")
    instanton, _ = relation_matrix(electroweak_instanton=True)
    extra = s.expand((instanton*s.Matrix(general))[-1])
    if extra != 9*d-3*e:
        raise ArithmeticError("three-generation instanton relation changed")
    integer_vector = [-2]*3+[1]*3+[3]*4+[0]*2
    if instanton*s.Matrix(integer_vector) != s.zeros(instanton.rows, 1):
        raise ArithmeticError("integer topological-charge normalization failed")
    return {
        "species_order": list(SPECIES),
        "hypotheses": {
            "charge_group_is_abelian_and_additive": True,
            "antiparticles_carry_inverse_topological_charge": True,
            "massless_photon_gluons_and_minimal_neutrinos_carry_zero": True,
            "all_nine_quark_charged_current_entries_nonzero": True,
            "every_declared_interaction_conserves_the_charge": True,
            "these_hypotheses_are_derived_for_OPH_defects": False,
        },
        "perturbative_vertices": _group_receipt(False),
        "general_integer_solution": dict(zip(SPECIES, map(str, general), strict=True)),
        "optional_electroweak_instanton": {
            **_group_receipt(True),
            "extra_relation_on_D_E": str(extra)+" = 0",
            "integer_target_group_solution_basis": integer_vector,
            "integer_basis_is_minus_three_times_electric_charge": True,
            "torsion_is_discarded": False,
            "native_instanton_constructed": False,
        },
        "conclusions": {
            "charged_lepton_generations_share_charge": True,
            "up_generations_share_charge": True,
            "down_generations_share_charge": True,
            "Z_and_Higgs_additive_topological_charge": 0,
            "distinct_generation_windings_compatible_with_hypotheses": False,
            "neutral_composite_defect_states_excluded": False,
            "nonadditive_or_approximate_topology_excluded": False,
            "neutrino_mass_extension_included": False,
            "physical_particle_defect_dictionary_constructed": False,
        },
    }
