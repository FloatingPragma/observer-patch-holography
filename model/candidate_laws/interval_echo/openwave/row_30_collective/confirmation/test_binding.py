"""Public-history binding tests on disposable repositories and explicit local heads."""

import hashlib
import json
import subprocess
from types import SimpleNamespace

import pytest

from model.candidate_laws.interval_echo.openwave.row_30_collective.confirmation import (
    binding as b,
)


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args]).decode().strip()


def start(p):
    p.mkdir()
    git(p, "init", "-q", "-b", "main")
    git(p, "config", "user.email", "fixture@example.invalid")
    git(p, "config", "user.name", "Fixture")
    git(p, "config", "commit.gpgsign", "false")
    return p


def put(root, name, data):
    p = root / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(data)


def commit(root):
    git(root, "add", "-A")
    git(root, "commit", "-qm", "fixture")
    return git(root, "rev-parse", "HEAD")


def canonical(revision):
    return (
        "| Pinned candidate (M10.5) | [Candidate at `"
        + revision
        + "`](https://github.com/FloatingPragma/observer-patch-holography/tree/"
        + revision
        + "/model) |\n"
    )


def digest(root, path):
    return hashlib.sha256((root / path).read_bytes()).hexdigest()


@pytest.fixture
def repos(tmp_path, monkeypatch):
    source, platform = [start(tmp_path / s) for s in ("source", "platform")]
    put(source, "carrier.py", "old carrier\n")
    ancestor = commit(source)
    verifier = f"{b.DIRECTORY}/binding.py"
    put(source, verifier, "frozen verifier\n")
    put(source, "runner.py", "work\n")
    impl = {path: digest(source, path) for path in ("runner.py", "carrier.py")}
    put(
        source,
        b.DECLARATION,
        json.dumps(dict(name="coherent pilot", implementation=impl)) + "\n",
    )
    working = {
        path: digest(source, path) for path in (verifier, b.DECLARATION, "runner.py")
    }
    native = {"carrier.py": digest(source, "carrier.py")}
    manifest = dict(
        format=2,
        adopted_source_pin=ancestor,
        adopted_native_inputs=native,
        working_inputs=working,
    )
    put(source, b.MANIFEST, json.dumps(manifest) + "\n")
    terms = (
        f"Protocol v3\nFrozen source manifest SHA-256 `{digest(source, b.MANIFEST)}`\n"
    )
    put(source, b.TERMS, terms)
    revision = commit(source)
    put(platform, b.CANON, canonical(ancestor))
    baseline = commit(platform)
    monkeypatch.setattr(b, "BASE", baseline)

    # Only the remote advertisement is simulated. All ancestry, history, and
    # regular-file content comparisons use the two actual Git repositories.
    def advertised(repo, url):
        assert (repo, url) in {(source, b.RER), (platform, b.OPENWAVE)}
        return git(repo, "rev-parse", "refs/heads/main")

    monkeypatch.setattr(b, "_head", advertised)
    return SimpleNamespace(
        source=source,
        platform=platform,
        revision=revision,
        ancestor=ancestor,
        terms=terms,
        canonical=canonical(revision),
        digest=digest(source, b.DECLARATION),
    )


def acceptance(r, **patch):
    record = dict(
        protocol=b.PROTOCOL,
        row="OW-30",
        decision="accepted",
        source_commit=r.revision,
        candidate_sha256=r.digest,
    )
    return json.dumps(record | patch) + "\n"


def stage(r):
    put(r.platform, b.TASK, r.terms)
    put(r.platform, b.POINTER, r.revision + "\n")
    put(r.platform, b.CANON, r.canonical)
    put(r.platform, b.DECISION, acceptance(r))


def test_acceptance_and_benign_canonical_edit(repos):
    r = repos
    assert b.check_local(r.source) == r.digest
    stage(r)
    anchor = commit(r.platform)
    put(r.platform, b.CANON, "Unrelated clarification\n" + r.canonical)
    current = commit(r.platform)
    bound = b.check_public(r.source, r.platform)
    assert bound == b.Accepted(anchor, r.revision, r.digest, current, r.revision)


def test_public_terms_absent_or_only_on_side_branch(repos):
    r = repos
    with pytest.raises(b.Unaccepted):
        b.check_public(r.source, r.platform)
    git(r.platform, "checkout", "-qb", "unadopted")
    stage(r)
    commit(r.platform)
    git(r.platform, "checkout", "-q", "main")
    git(
        r.platform,
        "merge",
        "--no-ff",
        "-s",
        "ours",
        "-qm",
        "side work excluded",
        "unadopted",
    )
    with pytest.raises(b.Unaccepted):
        b.check_public(r.source, r.platform)


@pytest.mark.parametrize("path", [b.DECISION, b.POINTER, b.CANON])
def test_incomplete_first_landing_cannot_be_repaired_later(repos, path):
    r = repos
    stage(r)
    (r.platform / path).unlink()
    commit(r.platform)
    stage(r)
    commit(r.platform)
    with pytest.raises(b.InvalidBinding, match="Missing acceptance evidence"):
        b.check_public(r.source, r.platform)


@pytest.mark.parametrize("path", [b.DECISION, b.POINTER, b.CANON, b.TASK])
def test_change_and_restoration_mark_the_first_draw_invalid(repos, path):
    r = repos
    stage(r)
    commit(r.platform)
    changed = {
        b.DECISION: acceptance(r, decision="withdrawn"),
        b.POINTER: r.ancestor + "\n",
        b.CANON: canonical(r.ancestor),
        b.TASK: "other terms\n",
    }
    put(r.platform, path, changed[path])
    commit(r.platform)
    stage(r)
    commit(r.platform)
    with pytest.raises(b.InvalidBinding):
        b.check_public(r.source, r.platform)


@pytest.mark.parametrize(
    "patch",
    [
        {"protocol": "oph-ow30-confirmation-v1"},
        {"row": "OW-04"},
        {"decision": "proposed"},
        {"candidate_sha256": "0" * 64},
        {"source_commit": "a" * 40},
    ],
)
def test_acceptance_is_explicit_and_specific(repos, patch):
    r = repos
    stage(r)
    put(r.platform, b.DECISION, acceptance(r, **patch))
    commit(r.platform)
    with pytest.raises(b.InvalidBinding, match="exact affirmative"):
        b.check_public(r.source, r.platform)


def test_rejects_unaccepted_or_ambiguous_canonical_pin(repos):
    r = repos
    stage(r)
    put(r.platform, b.CANON, canonical(r.ancestor) + "\nProposed source " + r.revision)
    commit(r.platform)
    with pytest.raises(b.InvalidBinding, match="Canonical pin"):
        b.check_public(r.source, r.platform)
    with pytest.raises(b.InvalidBinding):
        b._canonical_pin((r.canonical + r.canonical).encode())
    with pytest.raises(b.InvalidBinding):
        b._canonical_pin(r.canonical.replace(r.revision, r.ancestor, 1).encode())


def test_public_source_hashes_are_verified(repos):
    r = repos
    # An actual public-main child carrying a changed runner. Leave the local
    # tree on the coherent inputs so only the committed byte check can catch it.
    correct = r.revision
    put(r.source, "runner.py", "changed public runner\n")
    r.revision = commit(r.source)
    git(r.source, "checkout", "-q", "--detach", correct)
    r.canonical = canonical(r.revision)
    stage(r)
    commit(r.platform)
    with pytest.raises(b.InvalidBinding, match="Committed input differs"):
        b.check_public(r.source, r.platform)


def test_public_source_must_be_reachable_from_official_main(repos):
    r = repos
    git(r.source, "checkout", "-qb", "private")
    put(r.source, "only-here", "private commit\n")
    r.revision = commit(r.source)
    r.canonical = canonical(r.revision)
    stage(r)
    commit(r.platform)
    with pytest.raises(b.InvalidBinding, match="merge-base"):
        b.check_public(r.source, r.platform)


def test_every_worktree_input_is_verified(repos):
    r = repos
    put(r.source, "runner.py", "changed\n")
    with pytest.raises(b.InvalidBinding, match="Local file differs"):
        b.check_local(r.source)
    git(r.source, "restore", "runner.py")
    put(r.source, b.TERMS, r.terms.replace(digest(r.source, b.MANIFEST), "1" * 64))
    with pytest.raises(b.InvalidBinding, match="matching source manifest"):
        b.check_local(r.source)


@pytest.mark.parametrize("impl", [{}, {"other.py": "1" * 64}, {"runner.py": "2" * 64}])
def test_candidate_implementation_must_match_frozen_files(repos, impl):
    r = repos
    put(r.source, b.DECLARATION, json.dumps(dict(implementation=impl)) + "\n")
    manifest = json.loads((r.source / b.MANIFEST).read_text())
    manifest["working_inputs"][b.DECLARATION] = digest(r.source, b.DECLARATION)
    put(r.source, b.MANIFEST, json.dumps(manifest) + "\n")
    put(
        r.source,
        b.TERMS,
        f"Frozen source manifest SHA-256 `{digest(r.source, b.MANIFEST)}`\n",
    )
    with pytest.raises(b.InvalidBinding, match="implementation differs"):
        b.check_local(r.source)
