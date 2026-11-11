import unittest
from fractions import Fraction as F
from .calibration import Anchor,common_scale,convert_energy,ratio
from .catalogue import build

class CalibrationTests(unittest.TestCase):
    def test_exact_intersection_and_exclusion(self):
        a=Anchor('a',(1,2),(4,6),'control')
        b=Anchor('a',(2,3),(9,12),'control')
        self.assertEqual(common_scale([a,b]),(F(3),F(6)))
        self.assertIsNone(common_scale([a,Anchor('a',(1,1),(7,8),'control')]))
        with self.assertRaises(ValueError):common_scale([a,Anchor('b',(1,1),(4,4),'control')])
        with self.assertRaises(ValueError):common_scale([])
        with self.assertRaises(ValueError):common_scale([a,Anchor('a',(1,1),(4,4),'control','different scheme')])

    def test_complete_signed_energy_conversion(self):
        self.assertEqual(convert_energy('a','a',(-3,2),(2,4)),(-12,8))
        self.assertEqual(convert_energy('a','a',(0,0),(2,4)),(0,0))
        with self.assertRaises(ValueError):convert_energy('a','b',(1,2),(3,4))
        with self.assertRaises(ValueError):ratio((1,2),(0,1))

    def test_native_dictionary_and_all_species(self):
        r=build();rows=r['target_catalogue'];q=r['single_clock_quark_test']
        self.assertEqual(len(rows),17)
        self.assertTrue(all(x['physical_source_rest_mass_GeV'] is None for x in rows))
        self.assertIsNone(q['common_scale_interval'])
        lo,hi=map(F,q['source_ratio_interval']);target=F(q['target_fixed_point_ratio'])
        self.assertLess(target,lo);self.assertLess(lo,1);self.assertGreater(hi,1)
        self.assertFalse(q['probabilities_changed_by_unit_conversion'])
        from .catalogue import read
        original=read('model/native_quark_response/results.json')['full_interacting_prepared_records']
        for f in ['uc','dc','Q']:
            self.assertEqual(q['original_probability_accounting'][f]['prepared_pump_failure_interval'],original[f]['prepared_pump_failure_interval'])
            self.assertEqual(q['u_anchor_calibrated_all_bins'][f][0]['source_energy_interval'],original[f]['records'][0]['source_excitation_relative_actual_ground'])
        self.assertIsNone(q['failure_total_source_energy_GeV'])
        self.assertEqual(len(q['u_anchor_calibrated_all_bins']['uc']),len(q['u_anchor_calibrated_all_bins']['dc']))
        self.assertEqual(r['conditional_electron_anchored_leptons'][1]['value_GeV'][:13],'0.10565941449')

if __name__=='__main__':unittest.main()
