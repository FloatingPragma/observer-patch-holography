"""Exact counterexamples to capacity and response-domain identifications."""
import unittest

import sympy as sp

from . import gauge_capacity as g


class GaugeCapacityTests(unittest.TestCase):
    def test_exact_export_capacity_and_full_even_algebra_kernel(self):
        capacity = g.compiler_capacity()
        self.assertEqual((capacity["total_modes"], capacity["edge_qubits"],
                          capacity["local_penalty_cycles"], capacity["independent_cycles"]),
                         (174, 432, 260, 259))
        self.assertEqual(capacity["code_dimension"], 2 ** 173)
        fock = g.fock_capacity(174)
        self.assertEqual(fock["full_even_algebra_dimension"], 2 * fock["one_parity_image_dimension"])
        witness = g.parity_kernel_witness()
        odd, embedding = witness["odd_projection"], witness["even_embedding"]
        self.assertEqual(odd ** 2, odd)
        self.assertEqual(sp.trace(odd), 4)
        self.assertEqual(embedding.H * embedding, sp.eye(4))
        self.assertEqual(witness["restricted_odd_projection"], sp.zeros(4))
        operators = g.fock_operators(3)
        for a in operators:
            for b in operators:
                bilinear = a.H * b
                self.assertEqual(odd * bilinear, bilinear * odd)

    def test_seven_species_require_full_128_state_capacity(self):
        capacity = g.fock_capacity(7)
        self.assertEqual(capacity["full_dimension"], 128)
        self.assertEqual(capacity["parity_block_dimension"], 64)
        self.assertGreater(64, 7)
        self.assertGreater(64, 2 * 7)
        self.assertEqual(capacity["full_fock_qubits"], 7)
        for owners, singlets in ((10, 16796), (15, 9694845)):
            result = g.singlet_capacity(owners)
            self.assertEqual(result["collective_singlet_dimension"], singlets)
            self.assertTrue(result["binary_code_fits"])
        self.assertFalse(g.singlet_capacity(4)["binary_code_fits"])

    def test_inherited_current_algebra_is_faithful_on_raw_owners(self):
        self.assertEqual(g.current_algebra_dimensions(),
                         {"current_rank": 12, "bracket_rank": 11,
                          "center_rank": 1, "closed_span_rank": 12})
        for current in g.owner_currents():
            self.assertEqual(current.H, current)
            self.assertEqual(current[:, 0], sp.zeros(7, 1))

    def test_collective_singlet_code_annihilates_all_twelve_currents(self):
        witness = g.two_owner_witness()
        code, number, charged = (witness[key] for key in ("code", "scalar_number", "charged"))
        self.assertEqual(code.H * code, sp.eye(2))
        for current in witness["currents"]:
            self.assertEqual(current * code, sp.zeros(49, 2))
            self.assertEqual(current * number, number * current)
        self.assertEqual(code.H * number * code, sp.diag(0, 1))
        self.assertEqual(witness["currents"][-1] * charged, charged)
        self.assertEqual(number * charged, sp.zeros(49, 1))
        self.assertEqual(witness["independent_owner_charge_squared"], sp.Rational(2, 3))
        raw_x = witness["raw_x"]
        self.assertEqual(raw_x.H * raw_x, sp.eye(49))
        self.assertEqual(raw_x * charged, charged)
        self.assertEqual(raw_x * code, code * sp.Matrix([[0, 1], [1, 0]]))
        # A nonconstant scalar phase is independent of every old current:
        # all old currents kill code; number has a nonzero code compression.
        self.assertNotEqual(number * code, sp.zeros(49, 2))

    def test_relative_phase_is_observable_and_global_number_is_not(self):
        witness = g.phase_response_witness()
        splitter, total = witness["splitter"], witness["total_number"]
        self.assertEqual(sp.simplify(splitter.H * splitter), sp.eye(8))
        self.assertEqual(splitter * total, total * splitter)
        self.assertEqual(total * witness["initial"], 2 * witness["initial"])
        self.assertNotEqual(witness["hop"] * witness["numbers"][0],
                            witness["numbers"][0] * witness["hop"])
        records = witness["records"]
        self.assertEqual(records["none"]["fuel_one"], 1)
        self.assertEqual(records["source_pi"]["fuel_one"], 0)
        self.assertEqual(records["none"], records["global_pi_over_two"])
        for result in records.values():
            self.assertEqual(sum(result["complete_basis_law"]), 1)
            self.assertEqual(result["reference_one"], 1)
            for word, probability in enumerate(result["complete_basis_law"]):
                if word.bit_count() != 2:
                    self.assertEqual(probability, 0)


if __name__ == "__main__":
    unittest.main()
