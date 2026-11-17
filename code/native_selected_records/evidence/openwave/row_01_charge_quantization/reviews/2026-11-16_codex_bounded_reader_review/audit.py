"""Independent spectral interval audit; no producer or proposed reader import."""
import gzip
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROW = HERE.parents[1]
SOURCE = ROW/'work/2026-11-15_codex_common_selected_source'
READER = ROW/'work/2026-11-16_codex_bounded_cosine_reader'

def unpack(value):
    if isinstance(value, dict):
        return np.array(value['real']) + 1j*np.array(value['imag'])
    return np.array(value)

for folder, expected in [(SOURCE, '89bd1d876a6c428e20f1f0477f3b539c84b1e27292bb32f05e209f97f85887fe'),
                         (READER, '9315e2f41cb3029b26e8e6ee12482f22c6353d05ab95d5f83d1a2ed4c2815259')]:
    manifest = folder/'MANIFEST.sha256'
    assert hashlib.sha256(manifest.read_bytes()).hexdigest() == expected
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ', 1)
        assert hashlib.sha256((folder/name).read_bytes()).hexdigest() == digest

base = json.loads(gzip.decompress((SOURCE/'v2_run/EVIDENCE.json.gz').read_bytes()))['base_execution']
result = json.loads((READER/'OUTPUT.json').read_text())
word = np.linalg.matrix_power(unpack(base['clock']['W']), 4)
hermitian = (word+word.conj().T)/2
values, vectors = np.linalg.eigh(hermitian)
interior = sorted({round(float(x), 10) for x in values if abs(x) < 1-1e-8})
assert len(interior) == 2
cosine = interior[-1]  # lowest nonconstant unsigned frequency, no source f supplied
mask = abs(values-cosine) < 1e-8
assert sum(mask) == 2
projector = vectors[:, mask] @ vectors[:, mask].conj().T
frequency = math.acos(float(np.mean(values[mask])))
checks = []
for source, row in zip(base['cases'], result['rows'], strict=True):
    density = unpack(source['clock_density'])
    assert np.max(abs(density-density.conj().T)) < 1e-12
    assert abs(np.trace(density)-1) < 1e-11
    assert min(np.linalg.eigvalsh(density)) > -1e-12
    eta = float(1-np.trace(density@projector).real)
    assert -1e-11 <= eta <= .01
    real = row['record']['x_plus']-row['record']['x_minus']
    assert abs(real-np.trace(density@hermitian).real) < 1e-11
    low, high = row['reader']['frequency_interval']
    assert low <= frequency <= high
    for value in [low, high]:
        # Endpoint admissibility in the exact mixture interval, checked by inversion.
        other = (real-(1-.01)*math.cos(value))/.01
        assert abs(other) <= 1+1e-10
    assert abs(eta-row['separate_numerical_source_audit']['actual_off_band_population']) < 1e-10
    checks.append({'source_case': row['case'], 'independent_off_band': eta,
                   'frequency_in_interval': True})

# A falsely exact reader would reject a representable target/zero-band mixture.
planted_eta = .01
planted_real = (1-planted_eta)*cosine+planted_eta
wrong_exact_frequency = math.acos(planted_real)
assert abs(wrong_exact_frequency-frequency) > .004
print(json.dumps({'passed': True, 'producer_imports': False,
                  'independent_spectral_frequency': frequency,
                  'spectral_target_rank': int(sum(mask)),
                  'checked_densities': checks,
                  'false_zero_impurity_reader_error': abs(wrong_exact_frequency-frequency),
                  'scope': 'Independent numerical source-band certificate and interval verification, with exact mixture proof reviewed; no finite-shot or external G6 claim.'}, indent=2))
