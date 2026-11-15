"""Physical-domain and official-status publication controls."""
from copy import deepcopy
import json
from pathlib import Path
import unittest
from tools.native_bridge_results import NAMES,project,render

ROOT=Path(__file__).resolve().parents[1]
class BridgeResultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data={n:json.loads((ROOT/('model/openwave_'+n)/'results.json').read_text())for n in NAMES}
    def test_complete_source_inventory_has_no_physical_promotion(self):
        view=project(self.data);self.assertEqual(len(view['source_packages']),10)
        self.assertEqual(view['physical_comparison_count_added'],0)
        self.assertEqual(view['official_openwave_untested_rows'],31)
    def test_missing_source_and_promoted_quark_pole_fail(self):
        d=deepcopy(self.data);d.pop('local_gauge_memory')
        with self.assertRaises(KeyError):project(d)
        d=deepcopy(self.data);d['quark_string_transport']['physical_quark_pole_or_GeV_rest_mass_proved']=True
        with self.assertRaises(ValueError):project(d)
    def test_late_preparation_cannot_be_renamed_primordial(self):
        d=deepcopy(self.data);d['capacity_finite_wavelength']['scope']['late_preparation_is_early_primordial_solution']=True
        with self.assertRaises(ValueError):project(d)
    def test_arithmetic_certification_cannot_award_a_cell(self):
        d=deepcopy(self.data);d['m10_review_gate']['official_status']['upstream_counts']['M10']['validated']=1
        with self.assertRaises(ValueError):project(d)

    def test_optional_counterterm_cannot_transfer_decay_or_claim_physical_masses(self):
        for key in ('physical_mass_prediction_or_native_process_closure','unchanged_current_Ward_work_or_many_body_decay_transferred'):
            d=deepcopy(self.data);d['matter_flux_source']['optional_charged_counterterm'][key]=True
            with self.assertRaises(ValueError):project(d)
    def test_ideal_clock_theorem_cannot_hide_controller_bias(self):
        d=deepcopy(self.data);d['repeated_native_clock']['scope']['static_controller_exact_strength_independent_frequency']=True
        with self.assertRaises(ValueError):project(d)

if __name__=='__main__':unittest.main()
