"""Source custody for scoped observable claims, independent of physical awards."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from tools import build_openwave_joint_claims as claims


class ObservableCustodyTests(unittest.TestCase):
    def fixture(self, root):
        source=root/'producer.py';source.write_text('declared source')
        receipt=root/'result.json'
        value={'official_rows_promoted':0,'source_sha256':{'producer.py':hashlib.sha256(source.read_bytes()).hexdigest()}}
        receipt.write_text(json.dumps(value))
        return source,receipt,value

    def test_source_mutation_cannot_reuse_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source,_,_=self.fixture(root)
            with patch.object(claims,'BUNDLE',root):
                self.assertEqual(claims.observable_receipt('result.json')['official_rows_promoted'],0)
                source.write_text('different construction')
                with self.assertRaisesRegex(ValueError,'stale observable source'):
                    claims.observable_receipt('result.json')

    def test_experiment_cannot_award_a_row(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);_,receipt,value=self.fixture(root)
            value['official_rows_promoted']=1;receipt.write_text(json.dumps(value))
            with patch.object(claims,'BUNDLE',root):
                with self.assertRaisesRegex(ValueError,'cannot award'):
                    claims.observable_receipt('result.json')

    def test_external_source_or_missing_custody_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);_,receipt,value=self.fixture(root)
            with patch.object(claims,'BUNDLE',root):
                for pins in ({'../external.py':'hash'},{}):
                    value['source_sha256']=pins;receipt.write_text(json.dumps(value))
                    with self.assertRaises(ValueError):claims.observable_receipt('result.json')

    def test_component_authority_and_custody(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);_,receipt,value=self.fixture(root)
            del value['official_rows_promoted'];value.update(official_rows_changed=[],full_criterion_met=False)
            receipt.write_text(json.dumps(value))
            with patch.object(claims,'BUNDLE',root):
                self.assertFalse(claims.observable_receipt('result.json')['full_criterion_met'])
                value['official_rows_changed']=[31];receipt.write_text(json.dumps(value))
                with self.assertRaisesRegex(ValueError,'cannot award'):claims.observable_receipt('result.json')

    def test_boolean_award_is_not_integer_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);_,receipt,value=self.fixture(root)
            value['official_rows_promoted']=False;receipt.write_text(json.dumps(value))
            with patch.object(claims,'BUNDLE',root):
                with self.assertRaisesRegex(ValueError,'cannot award'):claims.observable_receipt('result.json')

    def test_C8_preserves_physical_and_preparation_domains(self):
        values={r['id']:r for r in claims.c8_finite_claims()}
        self.assertEqual(len(values),6)
        self.assertIn('fails under',values['native_pinned_charge_spectrum']['domain'])
        self.assertIn('No signed sea',values['shared_many_load_native']['domain'])
        self.assertIn('not native-only',values['declared_uniform_fiber_preparation']['domain'])
        self.assertIn('rather than an extra dephasing',values['continuing_record_apparatus']['domain'])
