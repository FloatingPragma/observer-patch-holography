"""Actual entry controls for receiving Python source epochs."""
from pathlib import Path
import os
import shutil
import subprocess
import sys
import tempfile

import pytest

from model.native_binary_source_refinement import certificate as cert
from model.native_binary_source_refinement import source_profile as profile

HAS_LAKE = shutil.which('lake') is not None
if os.environ.get('OPH_REQUIRE_KERNEL_TESTS') == '1' and not HAS_LAKE:
    raise RuntimeError('The declared kernel test job must install Lean')
needs_kernel = pytest.mark.skipif(not HAS_LAKE,
    reason='executed by native-binary-source-refinement.yml with Lean')


DRIVER = r'''
import copy
from functools import lru_cache, partial
from types import ModuleType
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys

root = Path.cwd()
mode = sys.argv[1]
donor = Path(sys.argv[2]) if len(sys.argv) > 2 else None


def purge_modules(*roots):
    for name, module in tuple(sys.modules.items()):
        locations = [getattr(module, '__file__', None),
                     *getattr(module, '__path__', ())]
        if any(location and any(Path(location).resolve().is_relative_to(base) for base in roots)
               for location in locations):
            sys.modules.pop(name, None)


def old_seam_module():
    purge_modules(root, donor)
    sys.path.insert(0, str(donor))
    try:
        module = importlib.import_module('model.observer.quantum_mixing')
        module.seams()
    finally:
        purge_modules(root, donor)
        sys.path.remove(str(donor))
    module.__name__ = 'external_display_only_seam_cache'
    module.__file__ = str(root/'retained.json')
    sys.modules['model.observer.quantum_mixing'] = module
    return module


if mode == 'transitive_import':
    held_seams = old_seam_module()
if mode == 'shared_preloaded':
    from model.observer import quantum_mixing as preloaded_mixing
if mode in ('preloaded', 'preloaded_missing'):
    # Ordinary earlier model imports must remain collectable, while this
    # already-loaded producer has no enrollment in the receiving epoch.
    from model.native_gravity_response import certificate as preloaded
    if mode == 'preloaded_missing':
        # Package attributes still hold the actual modules after the canonical
        # registry names disappear, before the receiving adapter even exists.
        for name, module in tuple(sys.modules.items()):
            location = getattr(module, '__file__', None)
            if (location and Path(location).resolve().is_relative_to(root)
                    and Path(location).name != '__init__.py'):
                sys.modules.pop(name, None)
from model.native_binary_source_refinement import certificate as cert
from model.native_binary_source_refinement.native_validation import NativeValidationSession

receipt = json.loads((root/'retained.json').read_bytes())
receipt['source_bindings'] = cert.source_bindings(root)
consumer = root/'model/native_gravity_response/certificate.py'
original = consumer.read_bytes()
needle = b'    square = C@C\n'
assert original.count(needle) == 1


def mutate():
    consumer.write_bytes(original.replace(needle,
        b'    square = C@C + np.diag([1] + [0]*(DIM-1))\n'))
    value = copy.deepcopy(receipt)
    value['source_bindings'] = cert.source_bindings(root)
    assert value['source_bindings'] != receipt['source_bindings']
    return value


def fresh_adapter():
    purge_modules(root)
    return importlib.import_module('model.native_binary_source_refinement.certificate')


def rejected(value, layer, session=None):
    try:
        cert.verify(value, session=session)
    except ValueError as error:
        assert str(error) == layer, str(error)
    else:
        raise AssertionError('the actual public verifier accepted source drift')


if mode == 'lean_required':
    os.environ['PATH'] = ''
    rejected(receipt, 'kernel:lean-required')
elif mode == 'first':
    rejected(mutate(), 'source:python-source-epoch')
    assert not list((root/'temp/native-binary-source').glob('session-*/events.json'))
elif mode == 'data':
    spec = root/'evidence/native_geometric_source_20260925/spec.json'
    spec.write_bytes(spec.read_bytes()+b'\n')
    value = copy.deepcopy(receipt)
    value['source_bindings'] = cert.source_bindings(root)
    rejected(value, 'source:python-source-epoch')
elif mode in ('preloaded', 'preloaded_missing'):
    # Pure controls are allowed even when the strict accepting context is
    # unavailable. The accepting entry requires a fresh receiving process.
    assert len(cert.carrier()['states']) == 66
    rejected(receipt, 'source:unbound-python-context')
elif mode == 'shared_preloaded':
    from model.native_born_bridge import certificate as born
    from model.native_thermodynamic_bridge import certificate as thermo
    assert born.carrier_seams() == preloaded_mixing.seams()
    assert thermo.heat_bath_certificate()['source_coin_branches'] == 3960
    rejected(receipt, 'source:unbound-python-context')
elif mode == 'shared_cached_acquisition':
    from model.native_born_bridge import certificate as born
    cert._PYTHON_EPOCH.check(root)
    name = 'model.observer.quantum_mixing'
    original_module = sys.modules[name]
    alternate = ModuleType(name)
    alternate.seams = lambda: ((0, 1),)
    born.carrier_seams.cache_clear()
    try:
        sys.modules[name] = alternate
        assert born.carrier_seams() == ((0, 1),)
    finally:
        sys.modules[name] = original_module
    assert born.carrier_seams() == ((0, 1),)
    rejected(receipt, 'source:unbound-python-context')
elif mode == 'enrolled_missing':
    sys.modules.pop('model.native_gravity_response.certificate')
    rejected(receipt, 'source:python-module-replaced')
elif mode == 'module_display':
    cert.gravity.__file__ = str(root/'display-only.py')
    rejected(receipt, 'source:python-module-replaced')
elif mode == 'package_context':
    cert.reconstruction.quantum_mixing.__package__ = 'donor_namespace.observer'
    rejected(receipt, 'source:python-import-context')
elif mode == 'function_reload':
    old_function = cert.gravity.source_data
    module = cert.gravity
    importlib.reload(module)
    assert cert.gravity is module and module.source_data is not old_function
    rejected(receipt, 'source:python-reference-replaced')
elif mode in ('cache_reference', 'partial_reference'):
    original_function = cert.gravity.source_data
    if mode == 'cache_reference':
        cert.gravity.source_data = lru_cache()(original_function.__wrapped__)
    else:
        cert.gravity.source_data = partial(original_function)
    rejected(receipt, 'source:python-reference-replaced')
elif mode == 'transitive_import':
    assert cert.reconstruction.quantum_mixing is held_seams
    rejected(receipt, 'source:unbound-python-context')
elif mode == 'warm':
    session = NativeValidationSession(root, root/'e')
    fresh = cert.verify(receipt, session=session)
    first_hits = cert.gravity.source_data.cache_info().hits
    before_native_checks = sum(row['event'] == 'native-check' for row in session.events)
    repeated = cert.verify(fresh, session=session)
    assert repeated == fresh
    assert cert.gravity.source_data.cache_info().hits > first_hits
    assert sum(row['event'] == 'native-check' for row in session.events) == before_native_checks
    assert any(row['event'] == 'validated-reuse' for row in session.events)
    events = len(session.events)
    value = mutate()
    rejected(value, 'source:python-source-epoch', session=session)
    assert len(session.events) == events
    # Clearing an old function's cache cannot authorize a different source.
    cert.gravity.source_data.cache_clear()
    rejected(value, 'source:python-source-epoch', session=session)
    cold = subprocess.run([sys.executable, '-B', '-c',
        'from model.native_gravity_response.certificate import source_data; source_data()'],
        cwd=root, capture_output=True, text=True, encoding='utf-8',
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))
    assert cold.returncode != 0
    assert 'regional state differs from the native count reconstruction' in cold.stderr
    assert 'python-source-epoch' not in cold.stderr
    consumer.write_bytes(original)
    assert cert.verify(fresh, session=session) == fresh
    # Keep the actual producer and receiver-owned kernel session. Reload only
    # the unchanged adapter after changing the consumer's current source.
    old_consumer = cert.gravity
    value = mutate()
    importlib.reload(cert)
    assert cert.gravity is old_consumer
    keep = {'model', 'model.native_binary_source_refinement',
            'model.native_binary_source_refinement.certificate',
            'model.native_binary_source_refinement.python_epoch'}
    removed = {}
    for name, module in tuple(sys.modules.items()):
        location = getattr(module, '__file__', None)
        if location and Path(location).resolve().is_relative_to(root) and name not in keep:
            removed[name] = sys.modules.pop(name)
    assert 'model.native_gravity_response.certificate' not in sys.modules
    events = len(session.events)
    rejected(value, 'source:unbound-python-context', session=session)
    assert len(session.events) == events
    # Re-registering the held objects and clearing their caches cannot turn
    # that invalid incarnation into a source witness either.
    sys.modules.update(removed)
    for function in (old_consumer.source_data, cert.reconstruction.source_counts,
                     cert.reconstruction.path_counts, cert.reconstruction.matrix_unit_paths):
        function.cache_clear()
    rejected(value, 'source:unbound-python-context', session=session)
    assert len(session.events) == events
    consumer.write_bytes(original)
    # A fresh receiving context must fix actual callable bindings as well as
    # module names. Preserve this process's genuine baseline cache and native
    # session, then load the current source through a fully fresh namespace.
    old_wrapper = old_consumer.source_data
    original_wrapped = old_wrapper.__wrapped__
    old_wrapper()
    assert old_wrapper.cache_info().currsize == 1
    consumer.write_bytes(original.replace(needle,
        b'    square = C@C + np.diag([1] + [0]*(DIM-1))\n'))
    cert = fresh_adapter()
    cert._PYTHON_EPOCH.check(root)
    bound_wrapper = cert.gravity.source_data
    value = copy.deepcopy(fresh)
    value['source_bindings'] = cert.source_bindings(root)
    assert bound_wrapper is not old_wrapper
    old_wrapper.__wrapped__ = bound_wrapper.__wrapped__
    events = len(session.events)
    try:
        cert.gravity.source_data = old_wrapper
        rejected(value, 'source:python-reference-replaced', session=session)
        old_wrapper.cache_clear()
        rejected(value, 'source:python-reference-replaced', session=session)
        cert.gravity.source_data = partial(old_wrapper)
        rejected(value, 'source:python-reference-replaced', session=session)
        assert len(session.events) == events
        cert.gravity.source_data = bound_wrapper
        cert._PYTHON_EPOCH.check(root)
        rejected(value, 'regional state differs from the native count reconstruction',
                 session=session)
    finally:
        old_wrapper.__wrapped__ = original_wrapped
        cert.gravity.source_data = bound_wrapper
        consumer.write_bytes(original)
    # Display paths cannot turn a held local producer into an external runtime
    # dependency. Keep the actual code/cache and genuine receiving session.
    cert = fresh_adapter()
    assert cert.verify(fresh, session=session) == fresh
    held_module, held_wrapper = cert.gravity, cert.gravity.source_data
    held_code = held_wrapper.__wrapped__.__code__
    consumer.write_bytes(original.replace(needle,
        b'    square = C@C + np.diag([0]*(DIM-1) + [1])\n'))
    keep = {'model', 'model.native_binary_source_refinement',
            'model.native_binary_source_refinement.certificate',
            'model.native_binary_source_refinement.python_epoch'}
    labels = []
    try:
        for name, module in tuple(sys.modules.items()):
            location = getattr(module, '__file__', None)
            if name not in keep and location and Path(location).resolve().is_relative_to(root):
                labels.append((module, location))
                module.__file__ = str(root/'display-only'/(name.replace('.', '_')+'.py'))
        importlib.reload(cert)
        assert cert.gravity is held_module and cert.gravity.source_data is held_wrapper
        assert held_wrapper.__wrapped__.__code__ is held_code
        for name, module in tuple(sys.modules.items()):
            locations = [getattr(module, '__file__', None),
                         *getattr(module, '__path__', ())]
            if name not in keep and any(location and Path(location).resolve().is_relative_to(root)
                                        for location in locations):
                sys.modules.pop(name, None)
        value = copy.deepcopy(fresh)
        value['source_bindings'] = cert.source_bindings(root)
        assert value['source_bindings'] != fresh['source_bindings']
        events = len(session.events)
        rejected(value, 'source:unbound-python-context', session=session)
        held_wrapper.cache_clear()
        rejected(value, 'source:unbound-python-context', session=session)
        assert len(session.events) == events
    finally:
        for module, location in labels:
            module.__file__ = location
        consumer.write_bytes(original)
    # Every local import edge must bind its actual object, including a
    # transitive producer acquired from a different source copy.
    held_qm = old_seam_module()
    old_seams = held_qm.seams()
    old_code = held_qm.seams.__wrapped__.__code__
    qm_path = root/'model/observer/quantum_mixing.py'
    qm_original = qm_path.read_bytes()
    qm_needle = b'    result = tuple(tuple(edge) for edge in source().seams())\n'
    assert qm_original.count(qm_needle) == 1
    try:
        qm_path.write_bytes(qm_original.replace(qm_needle,
            qm_needle+b'    result = result[1:] + result[:1]\n'))
        cert = importlib.import_module('model.native_binary_source_refinement.certificate')
        assert cert.reconstruction.quantum_mixing is held_qm
        assert held_qm.seams.__wrapped__.__code__ is old_code
        assert not Path(old_code.co_filename).resolve().is_relative_to(root)
        sys.modules.pop('model.observer.quantum_mixing', None)
        value = copy.deepcopy(fresh)
        value['source_bindings'] = cert.source_bindings(root)
        events = len(session.events)
        rejected(value, 'source:unbound-python-context', session=session)
        held_qm.seams.cache_clear()
        rejected(value, 'source:unbound-python-context', session=session)
        assert len(session.events) == events
        purge_modules(root, donor)
        cert = importlib.import_module('model.native_binary_source_refinement.certificate')
        cert._PYTHON_EPOCH.check(root)
        assert cert.carrier()['seams'] == old_seams[1:]+old_seams[:1]
        rejected(value, 'receipt:complete-semantic-equality', session=session)
    finally:
        qm_path.write_bytes(qm_original)
    # Relative requests retain their loader-recorded caller context even when
    # mutable package metadata names another valid namespace.
    carrier_path = donor/'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py'
    carrier_raw = carrier_path.read_bytes()
    carrier_needle = b'    return exact_reference_edges()\n'
    assert carrier_raw.count(carrier_needle) == 1
    carrier_path.write_bytes(carrier_raw.replace(carrier_needle,
        b'    edges = exact_reference_edges()\n    return edges[1:] + edges[:1]\n'))
    import hashlib
    spec_path = donor/'evidence/native_geometric_source_20260925/spec.json'
    spec_data = json.loads(spec_path.read_bytes())
    spec_data['source_sha256']['oph_exact/carrier.py'] = hashlib.sha256(carrier_path.read_bytes()).hexdigest()
    spec_path.write_text(json.dumps(spec_data)+'\n', encoding='utf-8')
    purge_modules(root, donor)
    sys.path.insert(0, str(donor))
    try:
        donor_geometry = importlib.import_module('model.common_source.slow_geometry')
        donor_edges = donor_geometry.source().seams()
    finally:
        purge_modules(root, donor)
        sys.path.remove(str(donor))
    cert = importlib.import_module('model.native_binary_source_refinement.certificate')
    assert cert.verify(fresh, session=session) == fresh
    qm = cert.reconstruction.quantum_mixing
    package = qm.__package__
    alias = 'donor_namespace'
    top = ModuleType(alias)
    top.__path__ = [str(donor/'model')]
    common = ModuleType(alias+'.common_source')
    common.__path__ = [str(donor/'model/common_source')]
    sys.modules[alias] = top
    sys.modules[alias+'.common_source'] = common
    sys.modules[alias+'.common_source.slow_geometry'] = donor_geometry
    value = copy.deepcopy(fresh)
    value['carrier']['seams'] = list(map(list, donor_edges))
    assert value['carrier']['seams'] == fresh['carrier']['seams'][1:]+fresh['carrier']['seams'][:1]
    value['data']['branches'] = []
    for state in range(66):
        for edge in range(30):
            for coin in (0, 1):
                row = fresh['data']['branches'][state*60+((edge+1)%30)*2+coin].copy()
                row[1] = edge
                value['data']['branches'].append(row)
    try:
        qm.__package__ = alias+'.observer'
        qm.seams.cache_clear()
        events = len(session.events)
        rejected(value, 'source:python-import-context', session=session)
        assert len(session.events) == events
        qm.__package__ = package
        qm.seams.cache_clear()
        cert._PYTHON_EPOCH.check(root)
        rejected(value, 'receipt:complete-semantic-equality', session=session)
    finally:
        qm.__package__ = package
        qm.seams.cache_clear()
        for name in (alias+'.common_source.slow_geometry', alias+'.common_source', alias):
            sys.modules.pop(name, None)
else:
    raise AssertionError('unknown test mode')
print(json.dumps({'mode': mode, 'accepted': False if mode != 'warm' else 'fresh-and-warm-only'}))
'''


def _isolated_source(with_libraries=False):
    # Keep native output paths short; never alter the canonical source tree.
    (cert.ROOT/'temp').mkdir(exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix='pe-', dir=cert.ROOT/'temp'))
    for relative in profile.source_bindings(cert.ROOT):
        source, target = (cert.ROOT/relative).resolve(), root/relative
        assert source.is_relative_to(cert.ROOT.resolve())
        assert target.resolve().is_relative_to(root.resolve())
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as stream:
            stream.write(source.read_bytes())
    (root/'retained.json').write_bytes((cert.HERE/'results.json').read_bytes())
    if with_libraries:
        packages = root/'Lean/.lake/packages'
        packages.parent.mkdir(parents=True)
        provided = (cert.ROOT/'Lean/.lake/packages').resolve()
        if os.name == 'nt':
            subprocess.run(['cmd', '/c', 'mklink', '/J', str(packages), str(provided)],
                           check=True, capture_output=True)
        else:
            packages.symlink_to(provided, target_is_directory=True)
    return root


def _run(mode, with_libraries=False):
    root = _isolated_source(with_libraries)
    donor = _isolated_source() if mode in ('warm', 'transitive_import') else None
    command = [sys.executable, '-B', '-c', DRIVER, mode]
    if donor is not None:
        command.append(str(donor))
    result = subprocess.run(command, cwd=root,
        capture_output=True, text=True, encoding='utf-8',
        env=dict(os.environ, PYTHONUTF8='1', PYTHONDONTWRITEBYTECODE='1'))
    assert result.returncode == 0, result.stdout+'\n'+result.stderr


@pytest.mark.parametrize('mode', ('first', 'data', 'preloaded', 'preloaded_missing',
                                  'enrolled_missing', 'module_display', 'package_context', 'function_reload', 'cache_reference',
                                  'partial_reference', 'transitive_import', 'lean_required'))
def test_actual_verify_rejects_unbound_python_epoch_without_kernel(mode):
    _run(mode)


@pytest.mark.parametrize('mode', ('shared_preloaded', 'shared_cached_acquisition'))
def test_shared_producers_keep_ordinary_calls_and_receiving_admission(mode):
    _run(mode)


@needs_kernel
def test_actual_verify_reuses_bound_epoch_and_rejects_rehashed_consumer_drift():
    # One isolated receiving process owns baseline and warm reuse. No fixture
    # substitutes a retained output for the native evaluator or its custody.
    _run('warm', with_libraries=True)
