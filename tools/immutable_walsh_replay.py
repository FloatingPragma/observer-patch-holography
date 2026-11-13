"""Run the complete original12 tests with an explicit exact-arithmetic adapter.

This is NOT a metadata-only or unchanged-algorithm replay.  Three native
control compiler functions use a proved exact radix-two Walsh algorithm.
All emitted native pulses, including zero/no-op requests, inverse checks,
resource counts and strict receipt comparisons are unchanged.  Only host
compilation arithmetic is accelerated; no native observer clock is sped up.

Run against a dedicated immutable checkout with its genuinely computed
current receipt.  Actual original12-test success and an original-receipt
comparison are required before this route is admissible in CI.  The source
implementation, adapter, exact-arithmetic helper and metadata guards are
all authenticated.  No scientific result or expected receipt is cached.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import tempfile
import time

SIBLING_SHA256 = {
    "immutable_custody_replay.py": "19f030d5c246efeb8a4556154a3ac0789ea2974d1483813f52f8bfb6d9c2ad61",
    "immutable_custody_pytest.py": "6867da943ec188b82eca69c60aba018f42ac2ed480068395ae7c8376a7f85ac5",
    "exact_walsh.py": "913ef7ab6d7bd9e567cfab5cc78dffbff36a5ea98f9e49f25e152252982b41f1",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_siblings():
    modules, paths = {}, {}
    for filename, expected in SIBLING_SHA256.items():
        path = Path(__file__).with_name(filename).resolve()
        if sha(path) != expected:
            raise ValueError("arithmetic replay dependency requires renewed review: " + filename)
        name = "_oph_arithmetic_replay_" + path.stem
        spec = importlib.util.spec_from_file_location(name,path)
        if spec is None or spec.loader is None:
            raise ValueError("cannot load reviewed arithmetic dependency")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        modules[filename], paths[filename] = module, path
    return modules, paths


def verify_helpers(paths):
    for name, path in paths.items():
        if sha(path) != SIBLING_SHA256[name]:
            raise ValueError("arithmetic replay dependency changed during execution: " + name)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--audit",type=Path,required=True)
    args=parser.parse_args(argv)
    root,audit,adapter=args.root.resolve(),args.audit.resolve(),Path(__file__).resolve()
    modules,paths=load_siblings()
    guard=modules["immutable_custody_replay.py"]
    harness=modules["immutable_custody_pytest.py"]
    arithmetic=modules["exact_walsh.py"]
    if (audit.is_relative_to(root) or guard.same_file(audit,adapter)
            or any(guard.same_file(audit,p) for p in paths.values()) or args.audit.is_symlink()
            or (audit.exists() and (not audit.is_file() or audit.stat().st_nlink != 1))):
        parser.error("audit must be an external unaliased report path")
    if any(os.environ.get(name) for name in ("PYTEST_ADDOPTS","PYTEST_PLUGINS")):
        parser.error("PYTEST_ADDOPTS and PYTEST_PLUGINS must be empty")
    harness.authenticate_test(root);harness.authenticate_harness(root)
    snapshot=guard.FrozenSnapshot(root)
    expected=root/guard.OUTPUT
    expected_sha=sha(expected)
    old_cwd,old_path,old_argv=Path.cwd(),list(sys.path),list(sys.argv)
    old_autoload=os.environ.get("PYTEST_DISABLE_PLUGIN_AUTOLOAD")
    started=time.monotonic();slots=[];math_slot=None;plugin=harness.ExactTests(root)
    report={"schema":"oph.adapted_exact_walsh_pytest_replay.v1",
            "mode":"audited_exact_walsh_and_metadata_path_cache_pytest",
            "root":str(root),"adapter_sha256":sha(adapter),
            "helper_sha256":dict(SIBLING_SHA256),"input_sha256":snapshot.bindings,
            "expected_receipt_sha256_before":expected_sha,
            "original_test":harness.TEST,"original_test_sha256":harness.TEST_SHA256,
            "python":sys.version,"platform":platform.platform(),
            "arithmetic_algorithm_changed":True,"metadata_only_replay":False,
            "unmodified_cli_run":False,"original_tests_modified":False,
            "scientific_result_cache":False,"comparison_function_changed":False,
            "native_observer_clock_or_resource_counts_changed":False,
            "optimization_scope":"Host compiler arithmetic only; all native pulses and no-op requests retained.",
            "reviewed_repository_pytest_ini_sha256":harness.PYTEST_INI_SHA256,
            "reviewed_absent_ancestor_conftests":list(harness.ABSENT_CONFTESTS)}
    try:
        os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"]="1"
        os.chdir(root);sys.path.insert(0,str(root))
        if any(guard.same_file(audit,root/name) for name in snapshot.names):
            raise ValueError("audit aliases an immutable source")
        guard.authenticate_owner(root);guard.authenticate_loaded_modules(snapshot)
        guard.install(snapshot,frozenset({guard.OUTPUT}),slots)
        source=importlib.import_module("model.native_effect_synthesis.synthesis")
        math_slot=arithmetic.SynthesisAdapter(source,root);math_slot.install()
        report["arithmetic_adapter"]=math_slot.evidence()
        import pytest
        report["pytest_version"]=pytest.__version__
        with tempfile.TemporaryDirectory(prefix="oph-walsh-pytest-") as temporary:
            temp=Path(temporary);config=temp/"pytest.ini";config.write_text("[pytest]\n")
            pytest_args=[harness.TEST,"-q","-W","error","-p","no:cacheprovider","--noconftest",
                         "--import-mode=importlib","-c",str(config),"--rootdir",str(root),
                         "--basetemp",str(temp/"cases")]
            report["pytest_arguments"]=pytest_args
            code=pytest.main(pytest_args,plugins=[plugin])
            report["pytest_exit_code"]=int(code);plugin.validate(code)
        guard.check_expected(expected,expected_sha)
        math_slot.verify()
        for slot in slots:slot.verify()
        harness.authenticate_test(root);harness.authenticate_harness(root)
        guard.authenticate_loaded_modules(snapshot);snapshot.check();verify_helpers(paths)
        if sha(adapter)!=report["adapter_sha256"]:
            raise ValueError("arithmetic replay adapter changed during execution")
        report.update(verdict="PASS_EXACT_ARITHMETIC_ADAPTED_ORIGINAL_12_TEST_MODULE",
                      expected_receipt_sha256_after=sha(expected),tests_passed=12,
                      original_report_called_with_exact_arithmetic_substitution=True,
                      original_strict_comparison_executed=True,
                      immutable_source_bytes_and_filenames_preserved=True)
    except BaseException as error:
        report.update(verdict="FAIL",error=repr(error));raise
    finally:
        if math_slot is not None:math_slot.restore()
        for slot in slots:slot.module.source_bindings=slot.original
        report["test_execution"]=plugin.evidence()
        report["metadata_caches"]=[slot.evidence() for slot in slots]
        report["elapsed_seconds"]=time.monotonic()-started
        if old_autoload is None:os.environ.pop("PYTEST_DISABLE_PLUGIN_AUTOLOAD",None)
        else:os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"]=old_autoload
        os.chdir(old_cwd);sys.path[:],sys.argv=old_path,old_argv
        if not any(guard.same_file(audit,root/name) for name in snapshot.names):
            audit.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
        print(json.dumps({key:report.get(key) for key in ("verdict","mode","tests_passed","elapsed_seconds")}),flush=True)


if __name__=="__main__":main()
