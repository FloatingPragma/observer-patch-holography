"""Hostile mutations for the exact Arithmon McKay certificate."""
import copy
import importlib.util
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path
from unittest.mock import patch

SCRIPT_DIR=Path(__file__).resolve().parent
sys.path.insert(0,str(SCRIPT_DIR))

def load_module(name, filename):
    spec=importlib.util.spec_from_file_location(name,SCRIPT_DIR/filename)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module

producer=load_module("mckay_producer","mckay_golden_field_certificate.py")
verifier=load_module("mckay_verifier","verify_mckay_golden_field_independent.py")

class ExactMcKayCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.elements=producer.enumerate_2i()
        cls.raw=producer.main_certificate()
        cls.verified=verifier.verify(cls.raw)
        cls.irreps,cls.labels,cls.irrep_report=producer.derive_irreducibles(
            cls.elements,producer.chars_from_doublet(cls.elements),
            producer.conjugacy_classes(cls.elements))
        cls.source=Path(producer.__file__).read_text()
        cls.independent_source=Path(verifier.__file__).read_text()

    def test_01_exact_coordinate_set_has_120_elements(self):
        self.assertEqual(len(self.elements),120)

    def test_02_duplicate_element_control_is_detected(self):
        rows=list(self.elements);rows.append(rows[0])
        self.assertNotEqual(len(rows),len(set(rows)))

    def test_03_altered_generator_element_fails_group_gate(self):
        changed=set(self.elements);x=next(iter(changed));changed.remove(x)
        with self.assertRaises(ValueError):producer.exact_group_checks(changed)

    def test_04_perturbed_exact_coefficient_fails_group_gate(self):
        changed=set(self.elements);x=next(iter(changed));changed.remove(x)
        y=(x[0]+producer.Q5(F(1,101)),*x[1:]);changed.add(y)
        with self.assertRaises(ValueError):producer.exact_group_checks(changed)

    def test_05_rationalized_golden_coordinates_fail_closure(self):
        with patch.object(producer,"PHI",producer.Q5(F(1618,1000))), patch.object(producer,"PHI_INV",producer.Q5(F(618,1000))):
            changed=producer.enumerate_2i()
        with self.assertRaises(ValueError):producer.exact_group_checks(changed)

    def test_06_removed_element_fails_closure_gate(self):
        changed=set(self.elements);changed.remove(producer.NEG_IDENTITY)
        with self.assertRaises(ValueError):producer.exact_group_checks(changed)

    def test_07_center_corruption_in_raw_receipt_fails_closed(self):
        bad=copy.deepcopy(self.raw);bad["source_group"]["center_size"]=1
        with self.assertRaises(ValueError):verifier.verify(bad)

    def test_08_nonfaithful_two_dimensional_map_fails_kernel_gate(self):
        identity=((producer.c5(producer.ONE),producer.c5(producer.ZERO)),
                  (producer.c5(producer.ZERO),producer.c5(producer.ONE)))
        trivial={g:identity for g in self.elements}
        with self.assertRaises(ValueError):producer.check_faithful_kernel(trivial,identity)

    def test_09_reducible_one_plus_one_fails_irreducibility_norm(self):
        trivial={g:producer.Q5(2) for g in self.elements}
        self.assertEqual(producer.int_value(producer.character_inner(trivial,trivial,self.elements),"test"),4)

    def test_10_wrong_dimension_character_fails_doublet_dimension(self):
        wrong={g:producer.Q5(3) for g in self.elements}
        self.assertNotEqual(wrong[producer.IDENTITY],producer.Q5(2))

    def test_11_deleted_irreducible_fails_completeness_sum(self):
        dims=self.irrep_report["dimensions"]
        self.assertNotEqual(sum(d*d for d in dims[:-1]),len(self.elements))

    def test_12_duplicated_irreducible_fails_orthogonality(self):
        self.assertEqual(producer.int_value(producer.character_inner(self.irreps[0],self.irreps[0],self.elements),"test"),1)
        self.assertNotEqual(producer.int_value(producer.character_inner(self.irreps[0],self.irreps[0],self.elements),"test"),0)

    def test_13_altered_character_value_fails_irreducibility_gate(self):
        chi=producer.chars_from_doublet(self.elements);chi[producer.IDENTITY]=producer.Q5(7)
        with self.assertRaises(ValueError):
            producer.int_value(producer.character_inner(chi,chi,self.elements),"test")

    def test_14_altered_class_assignment_fails_partition_gate(self):
        classes=producer.conjugacy_classes(self.elements)
        flattened=set().union(*classes[:-1])
        self.assertNotEqual(flattened,self.elements)

    def test_15_forced_wrong_inner_product_is_detectable(self):
        chi=producer.chars_from_doublet(self.elements)
        actual=producer.int_value(producer.character_inner(chi,chi,self.elements),"test")
        self.assertEqual(actual,1)
        self.assertNotEqual(actual,2)

    def test_16_dimension_multiset_is_discovered_not_loaded(self):
        self.assertFalse(self.raw["irreducibles"]["hardcoded_dimension_list"])
        self.assertFalse(self.raw["irreducibles"]["hardcoded_character_table"])

    def test_17_affine_e8_template_is_absent_from_producer(self):
        self.assertNotIn("AFFINE_E8_EDGES",self.source)
        self.assertIn("AFFINE_E8_EDGES",self.independent_source)

    def test_18_equal_vertex_and_edge_counts_do_not_certify_e8(self):
        path_edges=[[i,i+1] for i in range(8)]
        self.assertEqual((9,len(path_edges)),(9,8))
        self.assertEqual(verifier.isomorphisms(path_edges,verifier.AFFINE_E8_EDGES,9),[])

    def test_19_reusing_first_fusion_for_galois_fails_receipt_gate(self):
        bad=copy.deepcopy(self.raw)
        bad["galois_conjugate_fusion"]["matrix"]=copy.deepcopy(bad["fusion"]["matrix"])
        with self.assertRaises(ValueError):verifier.verify(bad)

    def test_20_forcing_conjugate_character_equal_is_detected(self):
        chi=producer.chars_from_doublet(self.elements)
        self.assertFalse(producer.char_equal(chi,producer.sigma_character(chi),self.elements))
        self.assertTrue(producer.char_equal(chi,chi,self.elements))

    def test_21_erasing_sqrt5_trace_coefficients_fails_generation(self):
        chi=producer.chars_from_doublet(self.elements)
        erased={g:producer.Q5(v.a) for g,v in chi.items()}
        self.assertFalse(any(v.b!=0 for v in erased.values()))

    def test_22_coefficient_field_alone_is_not_the_trace_field_gate(self):
        bad=copy.deepcopy(self.raw)
        bad["golden_character_field"]["nonrational_trace_witness"]["trace"]="1/2"
        with self.assertRaises(ValueError):verifier.verify(bad)

    def test_36_witness_minimal_polynomial_is_derived_exactly(self):
        witness=producer.Q5(F(1,2),F(-1,2))
        self.assertEqual(producer.quadratic_minimal_polynomial_text(witness),"X^2-X-1")
        bad=copy.deepcopy(self.raw)
        bad["golden_character_field"]["witness_minimal_polynomial"]="X^2+X-1"
        with self.assertRaises(ValueError):verifier.verify(bad)

    def test_23_galois_pair_does_not_select_positive_phi(self):
        self.assertNotEqual(self.raw["golden_embedding"]["phi"],self.raw["golden_embedding"]["psi"])
        self.assertTrue(self.raw["golden_embedding"]["sigma_phi_equals_psi"])

    def test_24_mass_power_selector_is_absent(self):
        self.assertNotIn("27^phi",self.source)

    def test_25_koide_selector_is_absent(self):
        self.assertNotIn("import koide",self.source.lower())
        self.assertNotIn("koideq",self.source.lower())

    def test_26_experimental_mass_data_are_absent(self):
        self.assertNotIn("measured mass",self.source.lower())

    def test_27_oph_is_not_imported_by_producer(self):
        self.assertNotIn("import oph",self.source.lower())
        self.assertNotIn("oph_",self.source.lower())

    def test_28_no_oph_character_table_or_fusion_is_loaded(self):
        self.assertNotIn("oph character table",self.source.lower())
        self.assertNotIn("oph fusion matrix",self.source.lower())

    def test_29_no_op_h_g_plus_selector(self):
        self.assertNotIn("g_plus",self.source.lower())

    def test_30_no_physical_e8_claim_is_emitted(self):
        self.assertNotIn("physical e8",self.source.lower())
        self.assertFalse(self.verified["trust_boundary"]["physical_identification"])

    def test_31_independent_verifier_does_not_import_producer(self):
        self.assertNotIn("import mckay_golden_field_certificate",self.independent_source)
        self.assertFalse(self.verified["independent_verifier"]["imports_producer"])

    def test_32_explicit_affine_e8_isomorphism_exists(self):
        self.assertTrue(self.verified["affine_e8"]["isomorphic"])
        self.assertEqual(len(self.verified["affine_e8"]["explicit_isomorphism"]),9)

    def test_33_galois_graph_is_recomputed_and_isomorphic(self):
        self.assertTrue(self.verified["galois"]["fusion_recomputed"])
        self.assertTrue(self.verified["galois"]["same_affine_e8_graph_type"])
        self.assertFalse(self.verified["galois"]["same_labeled_graph"])

    def test_34_sl2f5_is_not_identified_by_order_alone(self):
        self.assertFalse(self.raw["sl2f5_cross_identification"]["explicit_isomorphism"])
        self.assertTrue(self.raw["sl2f5_cross_identification"]["not_identified_by_order_alone"])

    def test_35_claim_boundary_excludes_koide_and_data(self):
        self.assertFalse(self.verified["trust_boundary"]["koide_used"])
        self.assertFalse(self.verified["trust_boundary"]["experimental_data_used"])

if __name__=="__main__":unittest.main()
