import unittest
from fractions import Fraction as F
from .catalogue import UnitMap,anchored_interval,transport_complete_measure,exact_unit_execution,base_catalogue,strong_and_weak_catalogue,read_parent,PARENTS
class UnitTransportTests(unittest.TestCase):
    def test_anchor_uncertainty_contains_all_corners(self):
        got=list(map(F,anchored_interval(['2','3'],['4','5'],['6','7'])))
        for e in (2,3):
            for a in (4,5):
                for g in (6,7):self.assertTrue(got[0]<=F(e*g,a)<=got[1])
        self.assertEqual(got,[F(12,5),F(21,4)])
    def test_failure_cannot_be_conditioned_away(self):
        u=UnitMap('a','r','t',F(2,3))
        with self.assertRaises(ValueError):transport_complete_measure(u,['a'],['2'],['1/2'],'0','3')
        v=transport_complete_measure(u,['a'],['2'],['1/2'],'1/2','3')
        self.assertEqual(v['energy_GeV_exact'],['4/3']);self.assertEqual(v['time_GeV_inverse_exact'],'9/2')
        self.assertEqual(v['failure_probability_exact'],'1/2')
    def test_source_change_not_hidden_by_equal_mass(self):
        a=UnitMap('a','r','t',F(1));b=UnitMap('b','r','t',F(1))
        with self.assertRaises(ValueError):a.assert_same_source(b)
        with self.assertRaises(ValueError):UnitMap('a','r','t',F(0))
        with self.assertRaises(ValueError):a.interval(['2','1'])
    def test_noncommuting_complete_instrument_and_reference(self):
        r=exact_unit_execution();self.assertLess(r['four_level_probability_error'],2e-14)
        self.assertAlmostEqual(sum(r['complete_probabilities']),1)
    def test_complete_catalogue_preserves_conventions(self):
        p={name:read_parent(name)[0]for name in PARENTS}
        rows,ref,extra,baryon=base_catalogue(p)
        more=strong_and_weak_catalogue(p)
        self.assertEqual(len(more['weak_partial_widths']),3)
        self.assertFalse(more['tau_leptonic_width_ratio_comparison']['used_in_parameter_selection'])
        self.assertFalse(more['nuclear_mass_calibration']['independent_prediction'])
        self.assertEqual(len(rows),6);self.assertEqual(len(ref),17);self.assertEqual(len(extra),6);self.assertEqual(len(baryon),2)
        self.assertTrue(all(not x['native_prediction']for x in ref))
        self.assertTrue(all(x['sterile_right_chirality_is_same_mass_species']for x in rows if x['particle'].startswith('nu_')))
        self.assertTrue(all(x['original_window_probability_lower']<x['conditional_window_probability_lower']for x in baryon))
        self.assertTrue(all('MSbar' in x['observable']for x in ref if x['particle']in ('up','down','strange','charm','bottom','top')))
if __name__=='__main__':unittest.main()
