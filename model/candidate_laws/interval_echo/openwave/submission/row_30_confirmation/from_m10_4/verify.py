"""Check both review stages against the recorded official base without a public draw."""

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = "c531265eb6085725c6268e156804499c810ba5e4"
PREVIOUS = "09f062f77f7f540f2783f128ae6314c5d8b818ea"
FROZEN = "92370140edfbdb0ac2be647b2b8166d821e5f1a2"
AREA = "openwave/xperiments/m10_oph/research"
READER = "model/candidate_laws/interval_echo/openwave/row_30_collective/confirmation/binding.py"


def git(repo, *args, env=None, data=None):
    return subprocess.run(
        ["git", "-C", str(repo), *args],
        input=data,
        stdout=subprocess.PIPE,
        env=env,
        check=True,
    ).stdout


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def verify(platform, source):
    """Use separate scratch indexes and unreachable commits, leaving all refs and files alone."""
    source = source.resolve()
    legacy = json.loads((HERE.parent / "PACKAGE.json").read_text())
    require(
        git(source, "rev-parse", "HEAD").decode().strip() == FROZEN,
        "wrong source checkout",
    )
    code = git(source, "show", f"{FROZEN}:{READER}")
    require(
        (source / READER).read_bytes() == code, "binding differs from frozen source"
    )
    spec = importlib.util.spec_from_file_location(
        "frozen_ow30_binding", source / READER
    )
    binding = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = binding
    spec.loader.exec_module(binding)
    declaration = binding.check_local(source)
    require(declaration == legacy["candidate_sha256"], "candidate declaration changed")
    terms = (source / binding.TERMS).read_bytes()
    require(digest(terms) == legacy["task_sha256"], "task terms changed")
    require(
        digest((source / binding.MANIFEST).read_bytes()) == legacy["manifest_sha256"],
        "frozen source manifest changed",
    )
    require(
        binding._canonical_pin(git(platform, "show", f"{BASE}:{binding.CANON}"))
        == PREVIOUS,
        "baseline does not select the M10.4 pair",
    )
    require(
        digest((HERE.parent / "OPENWAVE.patch").read_bytes())
        == legacy["artifacts_sha256"]["OPENWAVE.patch"],
        "legacy packet changed",
    )
    with tempfile.TemporaryDirectory() as tmp:
        env = dict(
            os.environ,
            GIT_INDEX_FILE=f"{tmp}/index",
            GIT_NO_REPLACE_OBJECTS="1",
            GIT_AUTHOR_NAME="Offline reviewer",
            GIT_COMMITTER_NAME="Offline reviewer",
            GIT_AUTHOR_EMAIL="offline@invalid.example",
            GIT_COMMITTER_EMAIL="offline@invalid.example",
            GIT_AUTHOR_DATE="2026-11-26T00:00:00Z",
            GIT_COMMITTER_DATE="2026-11-26T00:00:00Z",
        )

        def run(*args, data=None):
            return git(platform, *args, env=env, data=data)

        def stage(patch, parent):
            run("apply", "--cached", "--check", str(patch))
            run("apply", "--cached", str(patch))
            tree = run("write-tree").decode().strip()
            commit = (
                run(
                    "commit-tree",
                    tree,
                    "-p",
                    parent,
                    data=b"Offline packet verification\n",
                )
                .decode()
                .strip()
            )
            return tree, commit

        def absent(commit):
            try:
                binding._landing(platform, commit, terms, declaration)
            except binding.Unaccepted:
                return
            raise ValueError("source checker accepted an unregistered history")

        def objects(tree, *paths):
            return run("ls-tree", "-rz", tree, "--", *paths)

        run("read-tree", legacy["platform_base"])
        original_tree, _ = stage(
            HERE.parent / "OPENWAVE.patch", legacy["platform_base"]
        )
        run("read-tree", BASE)
        request_tree, request = stage(HERE / "01-request.patch", BASE)
        require(
            not objects(request, binding.TASK, binding.POINTER, binding.DECISION),
            "request publishes an operative term or decision",
        )
        require(
            objects(request, binding.CANON) == objects(BASE, binding.CANON),
            "request changes pin",
        )
        absent(BASE)
        absent(request)
        decision_tree, decision = stage(HERE / "02-maintainer-decision.patch", request)
        require(
            binding._landing(platform, decision, terms, declaration)
            == (decision, FROZEN),
            "frozen checker does not recognize the private acceptance fixture",
        )
        # Every supplied scientific, execution and disclosure blob matches the original proposal.
        exact = [
            f"{AREA}/scripts/m10_5",
            f"{AREA}/findings/m10_5_author_validation.md",
            binding.TASK,
            binding.POINTER,
            binding.DECISION,
        ]
        frozen_objects = objects(original_tree, *exact)
        require(
            frozen_objects == objects(decision_tree, *exact),
            "candidate payload diverged",
        )
        for tree in (request_tree, decision_tree):
            require(
                objects(tree, "MODELS.md") == objects(BASE, "MODELS.md"),
                "official grades changed",
            )
            for name in ("m10_4_task_details.md", "m10_4_maintainer_review.md"):
                path = f"{AREA}/{'tasks' if 'task' in name else 'findings'}/{name}"
                require(
                    objects(BASE, path) and objects(BASE, path) == objects(tree, path),
                    "M10.4 evidence changed",
                )
        # Replace only the sole active pin by M10.4 on the positive history: a real negative control.
        run("read-tree", decision_tree)
        run("reset", "-q", BASE, "--", binding.CANON)
        wrong_tree = run("write-tree").decode().strip()
        wrong = (
            run("commit-tree", wrong_tree, "-p", request, data=b"Wrong canonical pin\n")
            .decode()
            .strip()
        )
        try:
            binding._landing(platform, wrong, terms, declaration)
        except binding.InvalidBinding as err:
            require(
                "Canonical pin differs" in str(err),
                "wrong-pin control failed for another reason",
            )
        else:
            raise ValueError("wrong-pin control accepted")
        # Terms posted before the decision permanently invalidate a later
        # "repair" within the authentic first-parent window.
        run("read-tree", request_tree)
        run("reset", "-q", decision_tree, "--", binding.TASK)
        early_tree = run("write-tree").decode().strip()
        early = (
            run(
                "commit-tree",
                early_tree,
                "-p",
                request,
                data=b"Premature operative terms\n",
            )
            .decode()
            .strip()
        )
        belated = (
            run("commit-tree", decision_tree, "-p", early, data=b"Belated decision\n")
            .decode()
            .strip()
        )
        try:
            binding._landing(platform, belated, terms, declaration)
        except binding.InvalidBinding as err:
            require(
                "Missing acceptance evidence" in str(err),
                "early-terms control failed for another reason",
            )
        else:
            raise ValueError(
                "late acceptance incorrectly repaired a prior task landing"
            )
        return dict(
            official_platform_base=BASE,
            source_pin_before=PREVIOUS,
            proposed_source_pin=FROZEN,
            candidate_sha256=declaration,
            task_sha256=digest(terms),
            source_manifest_sha256=legacy["manifest_sha256"],
            request_tree=request_tree,
            decision_tree=decision_tree,
            identical_payload_files=len(frozen_objects.rstrip(b"\0").split(b"\0")),
            request_sha256=digest((HERE / "01-request.patch").read_bytes()),
            decision_sha256=digest(
                (HERE / "02-maintainer-decision.patch").read_bytes()
            ),
            original_guard_before="unaccepted",
            original_guard_request="unaccepted",
            original_guard_decision_fixture="accepted",
            original_guard_wrong_pin="invalid",
            original_guard_terms_before_decision="invalid",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.platform.resolve(), args.source.resolve())
    expected = json.loads((HERE / "receipt.json").read_text())
    require(result == expected["packet"], "reproduced receipt differs")
    print(json.dumps(result, indent=2))
