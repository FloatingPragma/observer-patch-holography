"""Admission and custody checks for downstream use of the common clock."""
from pathlib import Path
import hashlib
import json
import unittest

import numpy as np

from common_clock import CommonClock, HERE, ROOT


class ClockContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.clock = CommonClock()
        cls.rho, _ = cls.clock.preparation()

    def test_number_block_compression_rejects_cross_number_coherence(self):
        # N=0/N=2 superposition respects fermion parity. Silently deleting its
        # off-diagonal number block would change the admitted preparation.
        indices0 = self.clock.model.spectral[0][0]
        indices2 = self.clock.model.spectral[2][0]
        psi = np.zeros(384, complex)
        psi[indices0[0]] = 1/np.sqrt(2)
        psi[indices2[0]] = 1j/np.sqrt(2)
        state = np.outer(psi, psi.conj())
        with self.assertRaisesRegex(ValueError, "number-block-diagonal"):
            self.clock.semigroup(state, .1, .8)
        full = self.clock.evolve(state, .1)
        self.assertAlmostEqual(float(np.trace(full).real), 1)
        self.assertGreater(np.linalg.norm(full[np.ix_(indices0,indices2)]), .49)

    def test_no_event_path_retains_coherent_duration_and_full_state(self):
        state, record = self.clock.trajectory(self.rho, .37, .8, 1.6, [])
        self.assertEqual(record["native_attempts"], 0)
        self.assertEqual(record["native_commit_pulses"], 0)
        self.assertEqual(record["coherent_service_calls"], 1)
        self.assertEqual(record["supplied_coherent_duration"], .37)
        self.assertGreater(np.linalg.norm(state-self.rho), .1)

    def test_invalid_uniformization_cannot_masquerade_as_complete_instrument(self):
        with self.assertRaises(ValueError):
            self.clock.trajectory(self.rho, .37, .8, .4, [])
        with self.assertRaises(ValueError):
            self.clock.semigroup(self.rho, .37, .8, attempt_mark=-1)

    def test_every_reservoir_and_count_label_reconstructs_same_output(self):
        receipt = json.loads((HERE/"COMMON_CLOCK.json").read_text())
        archive = np.load(HERE/receipt["state_archive"])
        branches = [archive[key] for row in receipt["preparation_branches"]
                    for key in row["complete_final_count_state_archive_keys"]]
        self.assertEqual(len(branches), 12)
        self.assertLess(np.max(abs(sum(branches)-archive["final_state"])), 1e-12)
        self.assertTrue(all(np.linalg.eigvalsh(state).min()>-1e-12 for state in branches))

    def test_frozen_receipts_and_typed_clock_contract(self):
        contract = json.loads((HERE/"CLOCK_CONTRACT.json").read_text())
        receipt = json.loads((HERE/"COMMON_CLOCK.json").read_text())
        review = json.loads((HERE/"REFERENCE_REVIEW.json").read_text())
        self.assertTrue(all(receipt["checks"].values()))
        self.assertTrue(all(review["checks"].values()))
        for document in (contract, receipt, review):
            for path, expected in document["source_sha256"].items():
                self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(), expected, path)
        self.assertEqual(hashlib.sha256((HERE/contract["receipt"]).read_bytes()).hexdigest(), contract["receipt_sha256"])
        self.assertEqual(hashlib.sha256((HERE/contract["state_archive"]).read_bytes()).hexdigest(), contract["state_archive_sha256"])
        self.assertEqual(hashlib.sha256(np.asarray(self.clock.h, dtype="<f8").tobytes()).hexdigest(), contract["source_H_array_sha256_le_float64"])
        self.assertEqual(hashlib.sha256(np.asarray(self.rho, dtype="<c16").tobytes()).hexdigest(), contract["preparation_array_sha256_le_complex128"])
        self.assertFalse(any(contract["transfer_boundaries"].values()))
        self.assertEqual(contract["official_rows_changed"], [])


if __name__ == "__main__":
    unittest.main()
