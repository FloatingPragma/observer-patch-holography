"""Source, record, consumer and receiver-custody controls."""
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from model.native_binary_source_refinement import certificate as cert
from model.native_binary_source_refinement import native_validation as kernel
from model.native_binary_source_refinement import source_profile as profile
from tools import strict_json

HAS_LAKE = shutil.which('lake') is not None
if os.environ.get('OPH_REQUIRE_KERNEL_TESTS') == '1' and not HAS_LAKE:
    raise RuntimeError('The declared kernel test job must install Lean')
needs_kernel = pytest.mark.skipif(not HAS_LAKE, reason='executed by native-binary-source-refinement.yml with Lean')


@pytest.fixture(scope='module')
def python_replay():
    data = cert.carrier()
    function = profile.compile_native(profile.native_text(cert.ROOT))
    index = {state: i for i, state in enumerate(data['states'])}
    rows, counts = [], [[0]*66 for _ in range(66)]
    for i, state in enumerate(data['states']):
        for e, edge in enumerate(data['seams']):
            for coin in (0, 1):
                final, record = function(state, edge, coin)
                j = index[final]
                rows.append([i, e, coin, j, *record])
                counts[i][j] += 1
    return {'native_text': profile.native_text(cert.ROOT), 'branches': rows, 'counts': counts}


@pytest.fixture(scope='module')
def semantic_data(python_replay):
    return cert._receipt_math(python_replay, cert.carrier())


def test_actual_source_count_moments_and_regional_input(semantic_data):
    assert len(semantic_data['branches']) == 3960
    assert len(semantic_data['moments']) == 4356
    assert semantic_data['scene']['moment'] == '13/21465'
    assert semantic_data['scene']['coefficient_omission'] == '13/1287900'
    assert semantic_data['scene']['indices'] == [0, 6]
    assert semantic_data['scene']['tail_count'] == 104
    assert sum(row[-1] < 0 for row in semantic_data['branches']) == 300
    assert sum(row[-1] == 0 for row in semantic_data['branches']) == 3360


@pytest.mark.parametrize('name', tuple(profile.grammar_cases(profile.native_text(cert.ROOT))))
def test_unsupported_grammar_is_rejected(name):
    text = profile.grammar_cases(profile.native_text(cert.ROOT))[name]
    with pytest.raises(ValueError, match='grammar:'):
        profile.parse_source(text)


def test_arbitrary_tail_and_signed_record():
    run = profile.compile_native(profile.native_text(cert.ROOT))
    tail = (7, -11, 1000000000000000000000000)
    for x in (0, 1):
        for y in (0, 1):
            for coin in (0, 1):
                final, record = run((x, y, *tail), (0, 1), coin)
                assert final[2:] == tail and record == (coin, x-final[0])
    assert run((0, 1, *tail), (0, 1), 1) == ((1, 0, *tail), (1, -1))


def test_source_slice_preserves_raw_line_endings(tmp_path):
    path = tmp_path/'model/observer/native.py'
    path.parent.mkdir(parents=True)
    raw = profile.native_text(cert.ROOT).replace('\n', '\r\n')
    path.write_bytes(raw.encode('utf-8'))
    extracted = profile.native_text(tmp_path)
    assert '\r\n' in extracted
    with pytest.raises(ValueError, match='grammar:source-character-profile'):
        profile.parse_source(extracted)


def test_valid_diagonal_control_rejected_by_actual_consumer():
    source = cert.reconstruction
    original = source.reconstructed_density
    actual = cert.gravity.source_data()[2]
    diagonal = np.diag(np.diag(actual))
    assert cert.gravity.trace(diagonal) == 1
    for port in range(12):
        assert cert.blocks(cert.gravity.restrict(actual, (port,))) == cert.blocks(cert.gravity.restrict(diagonal, (port,)))
    try:
        source.reconstructed_density = lambda steps: tuple(map(tuple, diagonal))
        cert.gravity.source_data.cache_clear()
        with pytest.raises(ValueError, match='regional state differs'):
            cert.gravity.source_data()
    finally:
        source.reconstructed_density = original
        cert.gravity.source_data.cache_clear()


def test_coefficient_omission_rejected_by_actual_reconstruction():
    source = cert.reconstruction
    original = source.matrix_unit_moment
    try:
        source.matrix_unit_moment = lambda steps, a, b: source.closed_path_probability(
            steps, source.matrix_unit_paths()[a, b][0])/source.return_probability(steps)
        with pytest.raises(ValueError, match='closed-history moment'):
            source.reconstructed_density(2)
    finally:
        source.matrix_unit_moment = original


def test_actual_consumer_floor_tamper_rejected(python_replay, monkeypatch):
    counts, total, state, floor = cert.gravity.source_data()
    monkeypatch.setattr(cert.gravity, 'source_data', lambda: (counts, total, state, floor+1))
    with pytest.raises(ValueError, match='consumer:faithful-floor-value'):
        cert._receipt_math(python_replay, cert.carrier())


def test_equal_but_wrong_one_port_controls_rejected(python_replay, monkeypatch):
    original = cert.gravity.restrict
    def altered(matrix, region):
        if len(region) == 1:
            return {0: np.array([[Q(1, 6)]], dtype=object),
                    1: np.array([[Q(5, 6)]], dtype=object)}
        return original(matrix, region)
    monkeypatch.setattr(cert.gravity, 'restrict', altered)
    with pytest.raises(ValueError, match='control:one-port-values'):
        cert._receipt_math(python_replay, cert.carrier())


def test_lean_absence_fails_instead_of_skipping(monkeypatch, tmp_path):
    monkeypatch.setattr(kernel.shutil, 'which', lambda name: None)
    with pytest.raises(ValueError, match='kernel:lean-required'):
        kernel.NativeValidationSession(cert.ROOT, tmp_path)


def test_ambiguous_json_is_rejected():
    with pytest.raises(ValueError, match='duplicate JSON key'):
        strict_json.loads('{"schema":"a","schema":"b"}')


def test_transitive_source_graph_is_live(tmp_path):
    bindings = profile.source_bindings(cert.ROOT)
    assert 'model/observer/native.py' in bindings
    assert 'model/native_gravity_response/certificate.py' in bindings
    assert 'Lean/EventAlgebra/NativeBinarySourceSemantics.lean' in bindings
    assert 'Lean/EventAlgebra/NativeL2CountMoment.lean' in bindings
    for relative in bindings:
        target = tmp_path/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((cert.ROOT/relative).read_bytes())
    assert profile.source_bindings(tmp_path) == bindings
    consumer = tmp_path/'model/native_gravity_response/certificate.py'
    consumer.write_bytes(consumer.read_bytes()+b'\n# transitive source control\n')
    assert profile.source_bindings(tmp_path) != bindings


@pytest.fixture(scope='module')
def live_session():
    session = kernel.current_session(cert.ROOT)
    report = cert.verify(strict_json.load(cert.HERE/'results.json'), session=session)
    return session, report


@needs_kernel
def test_current_kernel_and_complete_receipt(live_session):
    session, report = live_session
    cert.compare_receipt(strict_json.load(cert.HERE/'results.json'), report)
    assert any(row['event'] == 'native-check' for row in session.events)
    assert any(row['event'] == 'validated-reuse' for row in session.events)


@needs_kernel
@pytest.mark.parametrize('field', ('counts', 'branches', 'moments', 'normalizer', 'denominator',
    'density', 'faithful_floor', 'region', 'regional_input', 'entropy', 'one_port_states', 'scene', 'grammar_controls'))
@pytest.mark.parametrize('mode', ('empty', 'changed'))
def test_all_load_bearing_receipt_fields_rejected(live_session, field, mode):
    session, report = live_session
    altered = deepcopy(report)
    altered['data'][field] = [] if mode == 'empty' else {'altered': True}
    with pytest.raises(ValueError, match='receipt:complete-semantic-equality'):
        cert.compare_receipt(altered, report)


@needs_kernel
@pytest.mark.parametrize('mutation', ('record', 'divisor', 'update'))
def test_source_mutations_reach_semantic_rejection(live_session, mutation):
    session, report = live_session
    text = report['source_text']
    changed = {'record': text.replace('(coin, current)', '(coin, 0)'),
               'divisor': text.replace(')) // 2', ')) // 3'),
               'update': text.replace('out[a] -= current', 'out[a] += current')}[mutation]
    profile.parse_source(changed)
    with pytest.raises(ValueError, match='source:signed-record-law|source:binary-endpoint-law|kernel:execution-rejected'):
        cert.check_source_case(changed, cert.carrier(), session)


@needs_kernel
def test_grammar_case_input_output_binding(live_session):
    session, report = live_session
    cases = profile.grammar_cases(report['source_text'])
    path = session.output/'grammar.json'
    path.write_text(json.dumps([[name, text] for name, text in cases.items()]), encoding='utf-8')
    result = strict_json.loads(session.execute(kernel.PREFIX+'GrammarCheck', [path]).decode())
    assert result == [{'name': name, 'source_text': text, 'grammar_admitted': False}
                      for name, text in cases.items()]


@needs_kernel
@pytest.mark.parametrize('member', ('object', 'stdout', 'companion'))
def test_private_artifact_tamper_is_rejected(live_session, member):
    session, _ = live_session
    record = session.unit(kernel.PREFIX+'Endpoints')
    path = (record['object'] if member == 'object' else record['object'].parents[1]/'stdout.log'
            if member == 'stdout' else Path(str(record['object'])+'.private'))
    original = path.read_bytes() if path.is_file() else None
    try:
        path.write_bytes((original or b'')+b'altered')
        with pytest.raises(ValueError, match='kernel:private-artifact-changed'):
            session.unit(kernel.PREFIX+'Endpoints')
    finally:
        if original is None:
            path.unlink()
        else:
            path.write_bytes(original)


@needs_kernel
@pytest.mark.parametrize('cache_layer', ('snapshot', 'module'))
@pytest.mark.parametrize('mutation', ('changed', 'missing', 'metadata', 'metadata_missing',
    'transitive', 'transitive_missing', 'transitive_metadata_missing', 'shadow', 'companion'))
def test_vendor_cache_layers_revalidate_every_member(tmp_path, live_session, cache_layer, mutation):
    # Real SDK copies enter the receiver's actual import-custody method.
    # This tests freshness, not a theorem assertion about altered object bytes.
    session, _ = live_session
    original_roots = session.libraries
    lo, hi = tmp_path/'lo', tmp_path/'hi'
    for base in (lo, hi):
        base.mkdir()
    copied = {}
    for module in ('Init.Notation', 'Init.Coe'):
        canonical = kernel._resolve_import(module, original_roots)
        relative = Path(*module.split('.')).with_suffix('.olean')
        target = lo/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        for source in (canonical, canonical.with_suffix('.ilean'), Path(str(canonical)+'.private'), Path(str(canonical)+'.server')):
            if source.is_file():
                destination = target.with_suffix('.ilean') if source.suffix == '.ilean' else Path(str(target)+str(source)[len(str(canonical)):])
                destination.write_bytes(source.read_bytes())
        copied[module] = target
    session.libraries = [hi.resolve(), lo.resolve(), *original_roots]
    runtime = session.settings()
    try:
        identity = session.library_snapshot(['Init.Notation'], runtime)
        assert session.library_snapshot(['Init.Notation'], runtime) == identity
        if cache_layer == 'module':
            for key in list(kernel._LIBRARY_SNAPSHOTS):
                if str(lo.resolve()) in key[2]:
                    del kernel._LIBRARY_SNAPSHOTS[key]
        child = copied['Init.Coe' if mutation.startswith('transitive') else 'Init.Notation']
        if mutation in ('missing', 'transitive_missing'):
            child.unlink()
        elif mutation == 'metadata':
            path = copied['Init.Coe'].with_suffix('.ilean')
            value = json.loads(path.read_bytes())
            value['directImports'] = []
            path.write_text(json.dumps(value), encoding='utf-8')
        elif mutation == 'shadow':
            path = hi/'Init/Coe.olean'
            path.parent.mkdir(parents=True)
            path.write_bytes(copied['Init.Coe'].read_bytes())
        elif mutation in ('metadata_missing', 'transitive_metadata_missing'):
            child.with_suffix('.ilean').unlink()
        elif mutation == 'companion':
            Path(str(child)+'.private').write_bytes(b'changed-companion')
        else:
            child.write_bytes(child.read_bytes()+b'altered')
        with pytest.raises(ValueError, match='kernel:library-'):
            session.library_snapshot(['Init.Notation'], runtime)
    finally:
        session.libraries = original_roots


@needs_kernel
@pytest.mark.parametrize('field', ('native_text', 'states', 'seams', 'branches', 'counts'))
def test_case_output_substitution_is_rejected(live_session, monkeypatch, field):
    session, report = live_session
    original = session.execute
    def substituted(module, args):
        raw = original(module, args)
        result = strict_json.loads(raw.decode('utf-8'))
        if field == 'native_text':
            result[field] += ' '
        elif field == 'branches':
            result[field][0], result[field][1] = result[field][1], result[field][0]
        else:
            result[field][0][0] += 1
        return json.dumps(result).encode('utf-8')
    monkeypatch.setattr(session, 'execute', substituted)
    with pytest.raises(ValueError, match='source:replay-input-binding|source:full-branch-record-equality|source:count-table-equality'):
        cert.check_source_case(report['source_text'], cert.carrier(), session)


@needs_kernel
def test_rehashed_false_source_proof_rejected_by_actual_entry(live_session):
    _, report = live_session
    # Keep Windows compiler output paths below the native file-name ceiling.
    # Preserve the isolated case and its raw kernel logs for inspection.
    root = Path(tempfile.mkdtemp(prefix='nf-', dir=cert.ROOT/'temp'))
    for relative in profile.source_bindings(cert.ROOT):
        target = root/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((cert.ROOT/relative).read_bytes())
    source = root/'model/observer/native.py'
    endpoints = root/'Lean/EventAlgebra/NativeBinarySourceEndpoints.lean'
    for path in (source, endpoints):
        text = path.read_text(encoding='utf-8')
        path.write_text(text.replace('(coin, current)', '(coin, current + 1)'), encoding='utf-8', newline='\n')
    altered = deepcopy(report)
    altered['source_text'] = profile.native_text(root)
    altered['source_bindings'] = profile.source_bindings(root)
    receipt = root/'false-receipt.json'
    receipt.write_text(json.dumps(altered), encoding='utf-8')
    packages = root/'Lean/.lake/packages'
    packages.parent.mkdir(parents=True)
    provided = (cert.ROOT/'Lean/.lake/packages').resolve()
    if os.name == 'nt':
        subprocess.run(['cmd', '/c', 'mklink', '/J', str(packages), str(provided)],
                       check=True, capture_output=True)
    else:
        packages.symlink_to(provided, target_is_directory=True)
    (root/'old-proof.log').write_text('all source theorems accepted\n', encoding='utf-8')
    result = subprocess.run([sys.executable, '-B', '-m', 'model.native_binary_source_refinement',
                             '--check', str(receipt)], cwd=root, capture_output=True, text=True,
                            encoding='utf-8', env=dict(os.environ, PYTHONUTF8='1'))
    assert result.returncode != 0
    assert 'kernel:current-source-rejected:EventAlgebra.NativeBinarySourceEndpoints' in result.stderr
    events = list((root/'temp/native-binary-source').glob('session-*/events.json'))
    assert len(events) == 1
    rows = json.loads(events[0].read_text(encoding='utf-8'))
    assert any(row['module'] == kernel.PREFIX+'Endpoints' and row['exit_code'] != 0 for row in rows)


@needs_kernel
def test_changed_runtime_input_is_rejected(live_session, monkeypatch):
    session, report = live_session
    inputs = session.output/'changing-inputs'
    inputs.mkdir(exist_ok=True)
    source, data = inputs/'source.py', inputs/'carrier.json'
    source.write_text(report['source_text']+'\n', encoding='utf-8', newline='\n')
    data.write_text(cert.canonical(cert.carrier()), encoding='utf-8')
    original = kernel.subprocess.run
    def changing(*args, **kwargs):
        result = original(*args, **kwargs)
        data.write_bytes(data.read_bytes()+b' ')
        return result
    monkeypatch.setattr(kernel.subprocess, 'run', changing)
    with pytest.raises(ValueError, match='kernel:runtime-input-changed'):
        session.execute(kernel.PREFIX+'Replay', [source, data])
