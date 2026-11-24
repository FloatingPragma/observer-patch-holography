"""Run the OW-30 confirmation from its frozen RER source or verified OpenWave port."""

import argparse
import hashlib
import importlib
from importlib.machinery import PathFinder
from contextlib import contextmanager
import os
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile


PACKAGE = "model.candidate_laws.interval_echo.openwave.row_30_collective.confirmation"
LOCAL = Path(*PACKAGE.split("."))
PORT = Path("openwave/xperiments/m10_oph/research/scripts/m10_5")
HERE = Path(__file__).resolve().parent
PILOT = (927016, 927016, 927016)
REPEAT = (374621, 614357, 801173)
VALIDATION = (785029, 514823, 690271)
PUBLIC = "https://github.com/openwave-labs/openwave.git"
BASE = "d4fc3e4046c472f4ad38468df8dddfa1f97bac97"
TASK = "openwave/xperiments/m10_oph/research/tasks/m10_5_task_details.md"


class NoPublicTask(ValueError):
    """Official first-parent history contains no task at this path."""


class ScientificUnavailable(Exception):
    """A required numerical library cannot be loaded."""


def git(root, *args):
    return subprocess.run(
        ["git", "-C", str(root), *args],
        check=True,
        env=dict(os.environ, GIT_NO_REPLACE_OBJECTS="1", GIT_TERMINAL_PROMPT="0"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def official_head(platform):
    advertised = git(
        platform, "ls-remote", "--exit-code", "--refs", PUBLIC, "refs/heads/main"
    )
    match = re.fullmatch(rb"([0-9a-f]{40})\trefs/heads/main\n", advertised)
    if not match:
        raise ValueError("Cannot authenticate official OpenWave main")
    return match[1].decode()


def authenticate_port(platform):
    """Authenticate the launcher's and vendor's official bytes without importing either package."""
    head = official_head(platform)
    if git(platform, "rev-parse", "--is-shallow-repository").strip() != b"false":
        raise ValueError("Complete OpenWave first-parent history required")
    history = git(platform, "rev-list", "--first-parent", head).decode().splitlines()
    if BASE not in history:
        raise ValueError("Official baseline missing from first-parent history")
    if not any(
        git(platform, "ls-tree", commit, "--", TASK)
        for commit in history[: history.index(BASE) + 1]
    ):
        raise NoPublicTask(
            "Official OpenWave main has no first-parent landing of the task"
        )
    paths = [PORT / "run.py", PORT / "vendor"]
    tree = git(platform, "ls-tree", "-rz", "--full-tree", head, "--", *map(str, paths))
    published = set()
    for entry in tree.split(b"\0"):
        if not entry:
            continue
        metadata, name = entry.split(b"\t", 1)
        mode, kind, oid = metadata.split()
        rel = Path(name.decode())
        path = platform / rel
        if (
            mode not in (b"100644", b"100755")
            or kind != b"blob"
            or path.resolve() != path
        ):
            raise ValueError(f"Indirect or nonregular OpenWave runtime file: {rel}")
        if not path.is_file() or path.read_bytes() != git(
            platform, "cat-file", "blob", oid.decode()
        ):
            raise ValueError(
                f"OpenWave working port differs from the official blob: {rel}"
            )
        published.add(rel)
    vendor = platform / PORT / "vendor"
    present = {
        p.relative_to(platform)
        for p in vendor.rglob("*")
        if (p.is_file() or p.is_symlink()) and "__pycache__" not in p.parts
    }
    if not present or published != present | {PORT / "run.py"}:
        raise ValueError(
            "Official OpenWave launcher or exact vendor inventory is missing"
        )
    return head


def public_settings(anchor):
    if not re.fullmatch(r"[0-9a-f]{40}", anchor):
        raise ValueError("Invalid public registration SHA")
    seeds = []
    for level in range(3):
        key = f"OPH-OW30-collective-confirmation-v3:{anchor}:{level}".encode("ascii")
        seeds.append(
            (1 << 63) | int.from_bytes(hashlib.sha256(key).digest()[:8], "big")
        )
    if len(set(seeds)) != len(seeds):
        raise ValueError("Public selector collision. This draw is invalid.")
    return seeds


def verified_root(source, port_root, revision):
    """Mount only release bytes before importing model code, on either host."""
    terms = (LOCAL / "TASK.md").as_posix()
    manifest = (LOCAL / "files.json").as_posix()
    manifest_blob = git(source, "show", f"{revision}:{manifest}")
    terms_blob = git(source, "show", f"{revision}:{terms}")
    expected = json.loads(manifest_blob)
    if (
        not isinstance(expected, dict)
        or not isinstance(expected.get("working_inputs"), dict)
        or (not isinstance(expected.get("adopted_native_inputs", {}), dict))
    ):
        raise ValueError("Invalid frozen file manifest")
    paths = dict(expected["working_inputs"])
    paths.update(expected.get("adopted_native_inputs", {}))
    paths[manifest] = hashlib.sha256(manifest_blob).hexdigest()
    paths[terms] = hashlib.sha256(terms_blob).hexdigest()
    root = port_root / PORT / "vendor" if port_root else source
    for rel, required in paths.items():
        path = root / rel
        if (
            not path.resolve().is_relative_to(root)
            or path.resolve() != path
            or digest(path) != required
        ):
            raise ValueError(f"Execution tree differs from the source binding: {rel}")
    if digest(Path(__file__)) != paths[(LOCAL / "run.py").as_posix()]:
        raise ValueError("Launcher differs from the source binding")
    if port_root:
        actual = {
            p.relative_to(root).as_posix()
            for p in root.rglob("*")
            if p.is_file() and "__pycache__" not in p.parts
        }
        if actual != set(paths):
            raise ValueError(
                "Incomplete or unexpected contents in the OpenWave runtime"
            )
    return root, set(paths)


def model_module(name):
    return name == "model" or name.startswith("model.")


class FrozenModelImports:
    """Resolve each model import against the manifest before executing its module body."""

    def __init__(self, runtime, inputs):
        self.runtime = runtime
        self.inputs = inputs

    def check(self, name, origin, locations):
        if origin:
            path = Path(origin).resolve()
            if (
                not path.is_relative_to(self.runtime)
                or path.relative_to(self.runtime).as_posix() not in self.inputs
            ):
                raise ValueError(f"Unbound model import: {name}")
        if locations is not None and any(
            not Path(p).resolve().is_relative_to(self.runtime) for p in locations
        ):
            raise ValueError(f"Unbound model namespace: {name}")

    def find_spec(self, fullname, path, target=None):
        if not model_module(fullname):
            return None
        spec = PathFinder.find_spec(fullname, path, target)
        if spec is None:
            raise ValueError(f"Missing frozen model module: {fullname}")
        self.check(fullname, spec.origin, spec.submodule_search_locations)
        return spec


@contextmanager
def model_runtime(runtime, inputs):
    guard = FrozenModelImports(runtime, inputs)
    for name, module in tuple(sys.modules.items()):
        if model_module(name) and module is not None:
            guard.check(
                name,
                getattr(module, "__file__", None),
                getattr(module, "__path__", None),
            )
    sys.path.insert(0, str(runtime))
    sys.meta_path.insert(0, guard)
    try:
        yield
    finally:
        sys.meta_path.remove(guard)
        sys.path.remove(str(runtime))


def scientific_engine():
    try:
        for dep in ("numpy", "scipy.linalg", "scipy.optimize", "scipy.sparse.csgraph"):
            importlib.import_module(dep)
    except (ImportError, OSError) as exc:
        raise ScientificUnavailable(f"Scientific runtime unavailable: {exc}") from exc
    acq = importlib.import_module(PACKAGE + ".acquire")
    importlib.import_module(PACKAGE + ".spatial")
    importlib.import_module(PACKAGE.rsplit(".", 1)[0] + ".coupling_meter")
    return acq


def run(args):
    local_tree = HERE.as_posix().endswith("/" + LOCAL.as_posix())
    base = Path(git(HERE, "rev-parse", "--show-toplevel").decode().strip()).resolve()
    if local_tree:
        source = base
        platform = args.openwave.resolve() if args.openwave else None
        if args.source and args.source.resolve() != source:
            raise ValueError(
                "Run the launcher at the root of the requested frozen source"
            )
    else:
        if HERE != base / PORT or not args.source:
            raise ValueError(
                "The OpenWave port requires its frozen RER checkout via --source"
            )
        platform, source = base, args.source.resolve()
        if args.openwave and args.openwave.resolve() != base:
            raise ValueError("The platform checkout and execution path must agree")
    public_port_head = None
    formal = args.command == "status" or (
        args.command != "check-sources" and args.mode == "public"
    )
    if not local_tree and formal:
        public_port_head = authenticate_port(platform)
    revision = git(source, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if not (source / LOCAL / "files.json").is_file() or not git(
        source, "ls-tree", revision, "--", (LOCAL / "files.json").as_posix()
    ):
        raise ValueError(
            "Use the frozen confirmation checkout named in the package manifest"
        )
    runtime, inputs = verified_root(source, None if local_tree else platform, revision)
    with model_runtime(runtime, inputs):
        return execute(args, source, platform, revision, local_tree, public_port_head)


def execute(args, source, platform, revision, local_tree, public_port_head):
    binding = importlib.import_module(PACKAGE + ".binding")
    candidate = binding.check_local(source)
    task_hash = digest(source / LOCAL / "TASK.md")
    if args.command == "check-sources":
        print(
            json.dumps(
                dict(
                    source_commit=revision,
                    candidate_sha256=candidate,
                    protocol_sha256=task_hash,
                )
            )
        )
        return 0
    mode = args.mode
    if args.command == "status" or mode == "public":
        if not platform:
            raise ValueError(
                "Public verification requires an official OpenWave checkout"
            )
        try:
            decision = binding.check_public(source, platform)
        except binding.Unaccepted as err:
            print(json.dumps({"status": "unaccepted", "reason": str(err)}))
            return 5
        if public_port_head is not None and decision.platform_head != public_port_head:
            raise ValueError(
                "Official OpenWave head changed during runtime authentication"
            )
        if decision.source_commit != revision or decision.candidate_sha256 != candidate:
            raise ValueError("Accepted candidate differs from this source checkout")
        if args.command == "status":
            print(json.dumps({"status": "accepted", **vars(decision)}))
            return 0
        if local_tree:
            raise ValueError(
                "Run the public measurement from the verified OpenWave port"
            )
        git(platform, "merge-base", "--is-ancestor", decision.anchor, "HEAD")
        for path in (binding.TASK, binding.POINTER, binding.DECISION):
            if (platform / path).read_bytes() != git(
                platform, "show", decision.anchor + ":" + path
            ):
                raise ValueError(
                    "Platform execution tree differs from the accepted registration"
                )
        if binding._canonical_pin((platform / binding.CANON).read_bytes()) != revision:
            raise ValueError(
                "The platform working specification names a different candidate"
            )
        seeds = public_settings(decision.anchor)
        anchor = decision.anchor
    else:
        seeds = list({"pilot": PILOT, "repeat": REPEAT, "validation": VALIDATION}[mode])
        anchor = None
    if not args.release:
        raise ValueError("The measurement requires a release directory")
    record = dict(
        mode=mode,
        anchor=anchor,
        source_commit=revision,
        candidate_pin=revision,
        candidate_sha256=candidate,
        protocol_sha256=task_hash,
    )
    acq = scientific_engine()
    action = {
        "acquire": acq.acquire,
        "inspect": acq.inspect_release,
        "replay": acq.replay,
    }[args.command]
    outcome = action(record, seeds, args.release)
    print(
        json.dumps(
            dict(
                status=outcome,
                mode=mode,
                source_commit=revision,
                release=str(args.release),
            )
        )
    )
    return acq.EXIT_CODES[outcome]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=("check-sources", "status", "acquire", "inspect", "replay")
    )
    parser.add_argument(
        "--mode", choices=("public", "pilot", "repeat", "validation"), default="public"
    )
    parser.add_argument("--openwave", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--release", type=Path)
    args = parser.parse_args()
    if not sys.flags.isolated or not sys.dont_write_bytecode or not __debug__:
        parser.error("use python -I -B with assertions enabled")
    with tempfile.TemporaryDirectory() as empty_cache:
        sys.pycache_prefix = empty_cache
        try:
            return run(args)
        except NoPublicTask as err:
            print(json.dumps({"status": "unaccepted", "reason": str(err)}))
            return 5
        except ScientificUnavailable as err:
            print(
                json.dumps({"status": "incomplete", "reason": str(err)}),
                file=sys.stderr,
            )
            return 3
        except (
            OSError,
            ImportError,
            KeyError,
            ValueError,
            TypeError,
            subprocess.CalledProcessError,
        ) as err:
            print(
                json.dumps({"status": "invalid", "reason": str(err)}), file=sys.stderr
            )
            return 4


if __name__ == "__main__":
    raise SystemExit(main())
