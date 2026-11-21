"""Exact attainment, premise controls and finite incidence mutation tests."""
import unittest

import sympy as s

from . import certificate as c


class CarrierSelectionTest(unittest.TestCase):
    def test_weighted_inverse_identity_for_an_arbitrary_complete_frame(self):
        # This is a nonicosahedral basis with unequal weights. It checks the
        # square-frame identity used in the universal weighted argument.
        vectors = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
                   (1, 1, 0), (1, 0, 1), (0, 1, 1)]
        projectors = [c.projector(v) for v in vectors]
        weights = [s.Rational(k, 21) for k in range(1, 7)]
        m = c.information(projectors, weights)
        q = s.Matrix.hstack(*(c.coordinates(p) for p in projectors))
        self.assertEqual(c.clean(q.T*m.inv()*q), s.diag(*(1/w for w in weights)))
        e0 = s.Matrix([1, 1, 1, 0, 0, 0])/s.sqrt(3)
        self.assertEqual(c.clean((e0.T*m*e0)[0]), s.Rational(1, 3))
        self.assertEqual(c.clean(s.trace(m)-(e0.T*m*e0)[0]), s.Rational(2, 3))

    def test_optimum_and_the_derived_incidence(self):
        result = c.witness_certificate()
        self.assertEqual(result['channels'], 6)
        self.assertEqual(result['weights'], ['1/6']*6)
        incidence = result['paired_positive_overlap_incidence']
        self.assertEqual((incidence['ports'], incidence['edges'], incidence['triangle_faces']), (12, 30, 20))
        self.assertEqual(incidence['degrees'], [5]*12)
        self.assertEqual(incidence['faces_per_edge'], [2])
        self.assertTrue(incidence['all_vertex_links_are_five_cycles'])

    def test_sign_classification_and_one_edge_mutation(self):
        result = c.sign_classification()
        self.assertEqual(result['normalized_sign_matrices_checked'], 2**10)
        self.assertEqual(result['conference_solutions'], 12)
        matrix = c.gram_signs()
        matrix[0][1] *= -1
        matrix[1][0] *= -1
        self.assertFalse(c.is_conference(matrix))
        with self.assertRaises(ValueError):
            c.paired_incidence(matrix)

    def test_switching_and_permutation_preserve_incidence(self):
        matrix = c.gram_signs()
        signs = [-1, 1, 1, -1, 1, -1]
        permutation = [4, 2, 0, 5, 1, 3]
        transformed = [[signs[i]*signs[j]*matrix[permutation[i]][permutation[j]]
                        for j in range(6)] for i in range(6)]
        self.assertTrue(c.is_conference(transformed))
        self.assertEqual(c.paired_incidence(transformed), c.paired_incidence(matrix))

    def test_essential_premise_controls(self):
        controls = c.control_certificates()
        self.assertEqual(controls['tetrahedral_four_axis_rank'], 4)
        self.assertEqual(controls['orthogonal_three_axis_rank'], 3)
        self.assertEqual(controls['nonoptimal_minimal_six_axis_rank'], 6)
        self.assertLess(s.Rational(controls['nonoptimal_minimal_unit_response_information']), s.Rational(2, 15))
        self.assertLess(s.Rational(controls['unequal_icosahedral_minimum_information']), s.Rational(2, 15))
        self.assertEqual(controls['nonminimal_optimum_paired_ports'], 24)
        self.assertEqual(controls['two_dimensional_optimum_paired_ports'], 6)
        self.assertFalse(controls['perturbed_optimum']['metric_is_equiangular'])
        self.assertTrue(controls['perturbed_optimum']['exact_icosahedral_sign_incidence'])

    def test_exact_stability_constants_have_a_strict_integer_gap(self):
        bounds = {k: s.Rational(v) for k, v in c.stability_arithmetic().items()}
        self.assertLess(bounds['inverse_norm_difference_upper'], s.Rational(1, 50))
        self.assertLess(bounds['weight_lower'], s.Rational(1, 6))
        self.assertGreater(bounds['weight_upper'], s.Rational(1, 6))
        self.assertEqual(bounds['conference_integer_gap_bound'], s.Rational(57729, 50000))
        self.assertLess(bounds['conference_integer_gap_bound'], bounds['nonzero_even_integer_gap'])

    def test_invalid_readout_inputs_are_rejected(self):
        axes = c.icosahedral_axes()
        bad_weights = [[s.Rational(1, 5)]*6,
                       [0]+[s.Rational(1, 5)]*5,
                       [-1]+[s.Rational(2, 5)]*5]
        for weights in bad_weights:
            with self.assertRaises(ValueError):
                c.information(axes, weights)
        with self.assertRaises(ValueError):
            c.information([s.eye(3)], [1])
        with self.assertRaises(ValueError):
            c.information([s.eye(3)/3], [1])
        with self.assertRaises(ValueError):
            c.projector([0, 0, 0])


if __name__ == '__main__':
    unittest.main()
