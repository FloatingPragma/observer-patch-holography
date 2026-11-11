import copy
import json
import unittest
from .bridge import HERE,build,validate_updates


class BridgeTests(unittest.TestCase):
    def test_exact_public_score_rule(self):
        r=build()
        self.assertEqual((r['baseline']['completed'],r['baseline']['total']),(50,235))
        self.assertEqual(r['proposed']['completed'],58)
        self.assertEqual(r['proposed']['status_counts']['conditional'],75)
        self.assertFalse(r['website_modified'])
        self.assertEqual(r['new_empirical_observations'],0)

    def test_no_duplicate_or_physical_promotion(self):
        r=build();snapshot=json.loads((HERE/'dashboard_inputs.json').read_text())
        updates=copy.deepcopy(r['updates']);updates[0]['physical_premise_promotions']=['physical clock']
        with self.assertRaises(ValueError):validate_updates(snapshot,updates)
        updates=copy.deepcopy(r['updates']);updates[0]=updates[1]
        with self.assertRaises(ValueError):validate_updates(snapshot,updates)
        updates=copy.deepcopy(r['updates']);updates[0]['status']='theorem'
        with self.assertRaises(ValueError):validate_updates(snapshot,updates)


if __name__=='__main__':unittest.main()
