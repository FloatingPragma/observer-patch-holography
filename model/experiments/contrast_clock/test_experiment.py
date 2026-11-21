"""Adversarial receipt and command-line tests, including optimized Python."""
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from . import experiment as e


class ReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt=e.read_json(e.HERE/'results.json')

    def test_receipt_passes_independent_scientific_relations(self):
        e.validate(self.receipt)

    def test_coherent_nonunitarity_cannot_hide_below_entrywise_tolerance(self):
        gate=np.eye(32)+5e-12*np.ones((32,32))
        self.assertGreater(np.linalg.norm(gate.T@gate-np.eye(32),2),3e-10)
        with self.assertRaisesRegex(ValueError,'nonunitary'):
            e.c.controller_matrix([gate])

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
            (['cases',0,'prefix_records',0,'no_pair'],-1e-12),
            (['cases',0,'prefix_records',0,'no_pair'],1e-12),
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

    def test_tiny_false_unfinished_or_early_ready_record_is_not_tolerated(self):
        bad=deepcopy(self.receipt)
        rows=bad['clock_records'][3]['public_records']
        blank=next(row for row in rows if row['readings']==[[0,0],[0,0]])
        largest=max(rows,key=lambda row:row['weight'])
        blank['weight']+=5e-12
        largest['weight']-=5e-12
        with self.assertRaisesRegex(ValueError,'impossible clock public support'): e.validate(bad)
        bad=deepcopy(self.receipt)
        bad['cases'][0]['prefix_records'][0]['archive'][0]['words'][0]=64
        with self.assertRaisesRegex(ValueError,'readiness'): e.validate(bad)

    def test_coordinated_erasure_of_resolved_unfinished_tails_is_rejected(self):
        bad=deepcopy(self.receipt)
        for index in (2,4):
            row=bad['clock_records'][index]
            self.assertGreater(row['unfinished_weight'],1e-171)
            self.assertLess(row['unfinished_weight'],1e-169)
            row['position_weights'][:1336]=[0.]*1336
            row['unfinished_weight']=0.
            recomputed={}
            for k,weight in enumerate(row['position_weights']):
                for branch in bad['cases'][0]['prefix_records'][min(k,1336)]['archive']:
                    key=tuple(tuple(p) for p in e.public(branch['words']))
                    recomputed[key]=recomputed.get(key,0.)+weight*branch['weight']
            row['public_records']=[{'readings':[list(p) for p in key],'weight':weight}
                                   for key,weight in sorted(recomputed.items())]
        # A coordinated mutation satisfies normalization and every saved
        # marginal identity. Comparison must preserve the analytic clock tail.
        e.validate(bad)
        with self.assertRaisesRegex(ValueError,'relative probability replay'):
            e.compare(self.receipt,bad)

    def test_impossible_nominal_mismatch_and_doubled_load_readings_are_rejected(self):
        bad=deepcopy(self.receipt)
        for i in range(78):
            bad['cases'][0]['branches'][1]['all_input_weights'][i]=1e-11
            bad['cases'][0]['branches'][0]['all_input_weights'][i]-=1e-11
        e.validate(bad)
        with self.assertRaisesRegex(ValueError,'impossible probability support'):
            e.compare(self.receipt,bad)
        bad=deepcopy(self.receipt)
        doubled=next(i for i,row in enumerate(e.c.source_sector()) if 2 in row)
        bad['cases'][0]['branches'][3]['all_input_weights'][doubled]=1e-11
        bad['cases'][0]['branches'][0]['all_input_weights'][doubled]-=1e-11
        e.validate(bad)
        with self.assertRaisesRegex(ValueError,'impossible probability support'):
            e.compare(self.receipt,bad)

    def test_zero_prepared_complex_amplitudes_and_prefix_support_are_exact(self):
        for part in ('prepared_output_real','prepared_output_imag'):
            bad=deepcopy(self.receipt)
            bad['cases'][0]['branches'][1][part][0]=1e-12
            e.validate(bad)
            with self.assertRaisesRegex(ValueError,'impossible prepared-output support'):
                e.compare(self.receipt,bad)
        bad=deepcopy(self.receipt)
        row=next(row for row in bad['cases'][0]['prefix_records'] if row['step']>=668 and
                 any(branch['weight']==0. for branch in row['archive']))
        zero=next(branch for branch in row['archive'] if branch['weight']==0.)
        largest=max(row['archive'],key=lambda branch:branch['weight'])
        zero['weight']=1e-12
        largest['weight']-=1e-12
        with self.assertRaisesRegex(ValueError,'impossible probability support'):
            e.compare(self.receipt,bad)

    def test_engineered_probability_tolerance_resolves_tiny_mass_and_subnormal_rounding(self):
        expected={'kind':'engineered','position_weights':[1e-170,0.]}
        e.compare(expected,{'kind':'engineered','position_weights':[1e-170*(1+1e-12),0.]})
        for weights in ([0.,0.],[1e-170,1e-300]):
            with self.assertRaises(ValueError):
                e.compare(expected,{'kind':'engineered','position_weights':weights})
        e.compare({'kind':'engineered','position_weights':[1e-323]},
                  {'kind':'engineered','position_weights':[5e-324]})
        for key in ('unfinished_weight','public_records'):
            expected={'kind':'engineered','position_weights':[1.],key:
                      1e-170 if key=='unfinished_weight' else [{'weight':1e-170}]}
            bad=deepcopy(expected)
            if key=='unfinished_weight': bad[key]=0.
            else: bad[key][0]['weight']=0.
            with self.assertRaisesRegex(ValueError,'relative probability replay'):
                e.compare(expected,bad)

    def test_nonzero_amplitudes_and_numerical_clock_controls_keep_backend_tolerance(self):
        bad=deepcopy(self.receipt)
        vector=bad['cases'][0]['branches'][0]['prepared_output_real']
        index=next(i for i,value in enumerate(vector) if abs(value)>1e-4)
        vector[index]+=1e-12
        # A numerical tridiagonal control can have cancellation-scale weights
        # with large relative variation. It has no analytic tiny-tail claim.
        bad['clock_records'][7]['position_weights'][1075]=2.4711587660189603e-31
        e.validate(bad)
        e.compare(self.receipt,bad)

    def test_optimized_cli_rejects_structural_weight_and_amplitude_mutations(self):
        for mutation in ('weight','amplitude'):
            bad=deepcopy(self.receipt)
            if mutation=='weight':
                bad['cases'][0]['branches'][1]['all_input_weights'][0]=1e-11
            else:
                bad['cases'][0]['branches'][1]['prepared_output_imag'][0]=1e-12
            with tempfile.TemporaryDirectory() as directory:
                path=Path(directory)/'mutant.json'
                path.write_text(json.dumps(bad),encoding='utf8')
                # Isolate receipt acceptance from source freeze changes while
                # exercising the production CLI under optimized Python.
                script=("import sys\nfrom model.experiments.contrast_clock import experiment as e\n"
                        "baseline=e.read_json(e.HERE/'results.json')\n"
                        "e.run=lambda:baseline\n"
                        f"sys.argv=['clock','--verify',{str(path)!r}]\ne.main()\n")
                proc=subprocess.run([sys.executable,'-B','-O','-c',script],cwd=e.ROOT,
                                    capture_output=True,text=True)
                self.assertNotEqual(proc.returncode,0,proc.stdout+proc.stderr)
                self.assertIn('impossible',proc.stderr)

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
from model.experiments.contrast_clock import experiment as e
r=e.read_json(e.HERE/'results.json')
r['clock_records'][2]['unfinished_weight']=1e-12
e.validate(r)
'''
        proc=subprocess.run([sys.executable,'-B','-O','-c',script],cwd=e.ROOT,capture_output=True,text=True)
        self.assertNotEqual(proc.returncode,0)
        self.assertIn('unfinished recomputation',proc.stderr)


if __name__=='__main__': unittest.main()
