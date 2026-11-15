"""Complete matrix and tally corruption controls."""
import unittest
from .source import MODELS,parse

def fixture():
    lines=['### Summary Count','| ✅ validated in-platform | 9 | 0 | 0 | 3 | 0 | 0 | 0 |',
        '| ⚠️ partial | 11 | 10 | 10 | 2 | 1 | 1 | 0 |','| ❌ negative | 2 | 0 | 3 | 3 | 0 | 0 | 0 |',
        '| 🚧 untested | 9 | 21 | 18 | 23 | 30 | 30 | 31 |','### Summary Status']
    columns=[['✅']*9+['⚠️']*11+['❌']*2+['🚧']*9,['⚠️']*10+['🚧']*21,
       ['⚠️']*10+['❌']*3+['🚧']*18,['✅']*3+['⚠️']*2+['❌']*3+['🚧']*23,
       ['⚠️']+['🚧']*30,['⚠️']+['🚧']*30,['🚧']*31]
    for n,row in enumerate(zip(*columns)):lines.append('| criterion'+str(n)+' | '+' | '.join(row)+' | static |')
    return '\n'.join(lines+['### Simplest Test per Criterion'])

class StatusTests(unittest.TestCase):
    def test_complete_official_column(self):
        x=parse(fixture());self.assertEqual(x['oph_official_rank'],7)
        self.assertEqual(x['upstream_counts']['M10']['untested'],31)
    def test_local_promotion_cannot_leave_old_tally(self):
        s=fixture().replace('| 🚧 | static |','| ✅ | static |',1)
        with self.assertRaises(ValueError):parse(s)
    def test_truncation_unknown_icon_and_false_total(self):
        for s in [fixture().replace('criterion0 |','criterion1 |',1),fixture().replace('| ✅ |','| ✔ |',1),fixture().replace('| 9 | 0 |','| 10 | 0 |',1)]:
            with self.assertRaises(ValueError):parse(s)
if __name__=='__main__':unittest.main()
