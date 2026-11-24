"""Test exact terms and first-parent rules on an unrelated synthetic history."""
from pathlib import Path
import subprocess
import tempfile
from unittest import TestCase
from unittest.mock import patch

from model.candidate_laws.interval_echo.openwave.row_30_collective.registration import verify


class Registration(TestCase):
    def test_first_landing_not_restoration_and_wrong_bytes(self):
        with tempfile.TemporaryDirectory() as path:
            repo=Path(path)
            def command(*args):
                return subprocess.check_output(['git','-C',path,*args]).decode().strip()
            command('init','-q')
            command('config','user.email','fixture@example.invalid')
            command('config','user.name','Fixture')
            def commit(data):
                (repo/'unrelated.txt').write_bytes(data)
                command('add','unrelated.txt')
                command('commit','-qm','unrelated synthetic fixture')
                return command('rev-parse','HEAD')
            base=commit(b'initial\n')
            terms=b'fictional task for a pytest assertion\n'
            first=commit(terms)
            commit(b'other\n')
            head=commit(terms)
            self.assertEqual(verify.first_exact(repo,base,head,'unrelated.txt',terms),first)
            self.assertIsNone(verify.first_exact(repo,base,head,'unrelated.txt',terms.rstrip()))
            with self.assertRaises(ValueError):
                verify.first_exact(repo,first,head,'unrelated.txt',terms)

    def test_absence_stops_before_selection_and_metrology(self):
        with patch.object(verify,'official_anchor',return_value=None),patch.object(verify,'settings') as selection:
            argv=['verify','acquire','--openwave','.','--release','unused']
            with patch('sys.argv',argv),self.assertRaisesRegex(SystemExit,'No exact public registration'):
                verify.main()
            selection.assert_not_called()
