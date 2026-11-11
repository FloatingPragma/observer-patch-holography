import copy
import json
import unittest

from .bridge import HERE
from .verify_live import certificate, verify_snapshot


class LiveVerificationTests(unittest.TestCase):
    def test_adopted_statuses_do_not_claim_complete_wording(self):
        result = certificate()
        self.assertEqual((result['live_completed'], result['live_total']), (58, 235))
        self.assertTrue(result['all_eight_status_updates_applied'])
        self.assertEqual(result['status_changes_needed'], 0)
        self.assertFalse(result['all_required_wording_corrections_applied'])
        self.assertEqual({row['id'] for row in result['remaining_text_corrections']},
                         {'exec-coherent-screen', 'exec-integer-lift'})
        self.assertFalse(result['website_written_by_this_package'])
        self.assertEqual(result['score_increase_from_remaining_text_corrections'], 0)

    def test_changed_membership_status_or_wording_requires_review(self):
        live = json.loads((HERE / 'live_status.json').read_text())
        proposed = json.loads((HERE / 'results.json').read_text())
        duplicate = copy.deepcopy(live)
        duplicate['observables'][0] = duplicate['observables'][1]
        with self.assertRaises(ValueError):
            verify_snapshot(duplicate, proposed)
        status = copy.deepcopy(live)
        next(row for row in status['observables'] if row['id'] == 'measurement')['status'] = 'theorem'
        with self.assertRaises(ValueError):
            verify_snapshot(status, proposed)
        wording = copy.deepcopy(live)
        next(row for row in wording['selected_source_rows']
             if row['id'] == 'exec-integer-lift')['why'] = 'Corrected physical scope.'
        with self.assertRaises(ValueError):
            verify_snapshot(wording, proposed)


if __name__ == '__main__':
    unittest.main()
