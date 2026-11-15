"""Freeze the opposing photon audit and strict immutable-data adapter."""
from pathlib import Path
import hashlib
import json
import sys

sys.dont_write_bytecode = True
import frozen

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    result = json.loads((HERE/'AUDIT.json').read_text())
    controls = json.loads((HERE/'FACADE_CONTROL.json').read_text())
    spatial = json.loads((HERE/'SPATIAL_DOMAIN.json').read_text())
    assert all(result['checks'].values()) and all(controls['checks'].values())
    assert all(spatial['checks'].values())
    assert set(result['checks']) == {'continuous', 'spectral_frame', 'photon_symbol', 'native', 'cap_boundary', 'spatial_route'}
    assert result['source_sha256'] == frozen.verify_inputs()
    names = ('audit.py', 'frozen.py', 'facade_control.py', 'spatial_domain.py', 'freeze.py', 'PROOF.md',
             'AUDIT.json', 'FACADE_CONTROL.json', 'SPATIAL_DOMAIN.json', 'IMPORTS.json')
    manifest = {'schema':'oph.quadratic_photon.opposing_freeze.v1', 'worker':'DYN20',
        'reviewed_manifest':frozen.MANIFEST_SHA256,
        'files':{str((HERE/name).relative_to(ROOT)):sha(HERE/name) for name in names},
        'source_pins':result['source_sha256'], 'checks':result['checks'],
        'facade_checks':controls['checks'],
        'spatial_domain_checks':spatial['checks'],
        'interfaces':{'certificate':'frozen.spectral_certificate',
                      'analytic_phase_targets':'frozen.phase_targets',
                      'frozen_diagonal':'frozen.photon_phase', 'frozen_experiment':'frozen.word',
                      'custody_check':'frozen.verify_inputs'},
        'native_default_domain':'Exact Fraction time1, X/Y clock, emitted eight-bit spectral code and ninth clock bit; all padding, resets and complete outcomes retained.',
        'scope':'Independent free-field spectral/cap/symbol/native-record audit and immutable-data replay facade. No local-edge or harmonic native encoder, charged interaction or physical species is certified.',
        'producer_bytes_unchanged':True, 'reserved_acceptance_runs':[], 'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(sha(HERE/'MANIFEST.json'))


if __name__ == '__main__': run()
