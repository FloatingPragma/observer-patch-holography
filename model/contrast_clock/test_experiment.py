"""Adversarial receipt and command-line tests, including optimized Python."""
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from . import experiment as e


class ReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt=e.read_json(e.HERE/'results.json')

    def test_receipt_passes_independent_scientific_relations(self):
        e.validate(self.receipt)

    def test_every_mutation_is_rejected_by_science_or_fresh_replay(self):
        def mutate(path,value):
            data=deepcopy(self.receipt)
            obj=data
            for key in path[:-1]: obj=obj[key]
            obj[path[-1]]=value
            return data
        attacks=[
            (['cases',0,'complete_instrument_error'],1e-8),
            (['cases',0,'complete_instrument_error'],True),
            (['cases',0,'final_no_pair'],0.),
            (['cases',0,'branches',0,'public'],[[0,0],[0,0]]),
            (['cases',0,'branches',0,'all_input_weights',0],-1e-14),
            (['cases',0,'prefix_records',1,'archive',0,'weight'],-1e-14),
            (['cases',0,'prefix_records',1,'counts',0],False),
            (['cases',0,'prefix_records',1,'counts',0],0),
            (['cases',0,'prefix_records',0,'step'],True),
            (['clock_records',0,'position_weights',0],float('nan')),
            (['clock_records',2,'unfinished_weight'],1e-12),
            (['clock_records',0,'public_records',0,'weight'],-1e-14),
            (['clock_records',6,'unfinished_weight'],0.),
            (['resources','clock_operator_norm'],0.),
            (['resources','fresh_reader_bits'],0),
            (['resources','window_unfinished_bound'],0.),
            (['resources','window_unfinished_bound'],-1e-12),
            (['clock_records',2,'unfinished_weight'],0.),
            (['competing_controls',0,'public_total_variation'],0.),
            (['literal_clock_check','operator_errors',0],3e-10),
            (['claims','physical_clock_or_particle_identification'],True),
        ]
        for path,value in attacks:
            with self.subTest(path=path):
                bad=mutate(path,value)
                with self.assertRaises(ValueError):
                    e.validate(bad)
                    e.compare(self.receipt,bad)

    def test_complete_prefix_and_branch_sets_cannot_be_pruned(self):
        for parent,key in [('prefix_records',50),('branches',0)]:
            bad=deepcopy(self.receipt)
            bad['cases'][0][parent].pop(key)
            with self.assertRaises(ValueError):
                e.validate(bad)
                e.compare(self.receipt,bad)
        bad=deepcopy(self.receipt)
        bad['clock_records'].pop()
        with self.assertRaises(ValueError): e.validate(bad)

    def test_wrong_agreement_is_detectably_wrong_and_reverse_sign_is_retained(self):
        rows=self.receipt['cases']
        self.assertLess(rows[2]['final_no_pair'],rows[0]['final_no_pair'])
        self.assertGreater(rows[1]['final_no_pair'],.4)
        self.assertGreater(self.receipt['competing_controls'][-1]['public_total_variation'],.8)
        self.assertEqual(len(rows[0]['branches']),16)
        # The complete doubled complement contributes deterministic empty-pair
        # outputs; it cannot be discarded from the input instrument.
        self.assertEqual(len(rows[0]['branches'][0]['all_input_weights']),78)

    def test_strict_json_rejects_duplicate_nonfinite_and_overflow(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.json'
            for text in ('{"x":1,"x":2}','[NaN]','[Infinity]','[1e999]'):
                path.write_text(text)
                with self.assertRaises(ValueError): e.read_json(path)

    def test_source_pin_coverage_cannot_be_reduced(self):
        with patch.object(e,'read_json',return_value={}):
            with self.assertRaisesRegex(ValueError,'pin coverage'): e.source_pins()

    def test_refrozen_protocol_cannot_change_control_or_precision_gate(self):
        protocol=e.read_json(e.HERE/'protocol.json')
        for key,value in [('angle',.25),('laws',[8.,0.]),('empirical_selection',True)]:
            bad=deepcopy(protocol); bad[key]=value
            with patch.object(e,'read_json',return_value=bad):
                with self.assertRaises(ValueError): e.protocol()

    def test_failed_record_is_retained_and_existing_file_is_preserved(self):
        bad=deepcopy(self.receipt)
        bad['literal_clock_check']['operator_errors'][0]=1.
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'failure.json'
            with patch.object(e,'run',return_value=bad),patch.object(sys,'argv',['clock','--record',str(path)]):
                with self.assertRaisesRegex(ValueError,'literal clock'): e.main()
                self.assertEqual(e.read_json(path),bad)
                before=path.read_bytes()
                with self.assertRaises(FileExistsError): e.main()
                self.assertEqual(path.read_bytes(),before)

    def test_optimized_verifier_rejects_changed_clock_and_false_completion(self):
        script='''
from model.contrast_clock import experiment as e
r=e.read_json(e.HERE/'results.json')
r['clock_records'][2]['unfinished_weight']=1e-12
e.validate(r)
'''
        proc=subprocess.run([sys.executable,'-B','-O','-c',script],cwd=e.ROOT,capture_output=True,text=True)
        self.assertNotEqual(proc.returncode,0)
        self.assertIn('unfinished recomputation',proc.stderr)


if __name__=='__main__': unittest.main()
