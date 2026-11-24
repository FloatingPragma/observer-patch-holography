"""Verify the official v3 OW-30 acceptance and its executable inputs before drawing."""

import hashlib
import json
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

BASE = "d4fc3e4046c472f4ad38468df8dddfa1f97bac97"
OPENWAVE = "https://github.com/openwave-labs/openwave.git"
RER = "https://github.com/FloatingPragma/observer-patch-holography.git"
DIRECTORY = "model/candidate_laws/interval_echo/openwave/row_30_collective/confirmation"
MANIFEST = f"{DIRECTORY}/files.json"
DECLARATION = f"{DIRECTORY}/candidate.json"
TERMS = f"{DIRECTORY}/TASK.md"
TASK = "openwave/xperiments/m10_oph/research/tasks/m10_5_task_details.md"
CANON = "openwave/xperiments/m10_oph/research/m10_theory_canonical.md"
POINTER = TASK.removesuffix(".md") + ".source"
DECISION = TASK.removesuffix(".md") + ".acceptance.json"
PROTOCOL = "oph-ow30-confirmation-v3"


class Unaccepted(ValueError):
    """No first-parent landing of these exact terms on official OpenWave main."""


class InvalidBinding(ValueError):
    """A required source, acceptance or continuity check failed."""


@dataclass(frozen=True)
class Accepted:
    anchor: str
    source_commit: str
    candidate_sha256: str
    platform_head: str
    source_head: str


def _git(repo, *args):
    env = dict(os.environ, GIT_NO_REPLACE_OBJECTS="1", GIT_TERMINAL_PROMPT="0")
    result = subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, env=env
    )
    if result.returncode:
        raise InvalidBinding(f"git {args[0]} failed: {result.stderr.decode().strip()}")
    return result.stdout


def _head(repo, url):
    advertised = _git(
        repo, "ls-remote", "--exit-code", "--refs", url, "refs/heads/main"
    )
    match = re.fullmatch(rb"([0-9a-f]{40})\trefs/heads/main\n", advertised)
    if not match:
        raise InvalidBinding(f"Cannot establish official main at {url}")
    return match[1].decode()


def _digest(value):
    return hashlib.sha256(value).hexdigest()


def _json(raw):
    def unique(entries):
        obj = dict(entries)
        if len(obj) != len(entries):
            raise InvalidBinding("Duplicate JSON key")
        return obj

    try:
        obj = json.loads(raw, object_pairs_hook=unique)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise InvalidBinding("Invalid JSON document") from exc
    if not isinstance(obj, dict):
        raise InvalidBinding("Expected a JSON object")
    return obj


def _local(root, name):
    path = root / name
    if path.resolve() != root / name or not path.is_file():
        raise InvalidBinding(f"Missing or indirect local file: {name}")
    return path.read_bytes()


def _tree(repo, revision, paths):
    # Inspect the Git mode too. A symlink is a blob of its link text.
    entries = _git(repo, "ls-tree", "-z", "--full-tree", revision, "--", *paths)
    found = {}
    for entry in entries.rstrip(b"\0").split(b"\0") if entries else ():
        meta, name = entry.split(b"\t", 1)
        mode, kind, oid = meta.split()
        if mode not in (b"100644", b"100755") or kind != b"blob":
            raise InvalidBinding(f"Expected a regular Git file at {revision}:{name!r}")
        found[name.decode()] = oid.decode()
    return found


def _blob(repo, oid):
    return _git(repo, "cat-file", "blob", oid)


def _files(manifest, key):
    entries = manifest.get(key, {})
    if not isinstance(entries, dict):
        raise InvalidBinding(f"Expected a file map in {key}")
    for name, digest in entries.items():
        path = PurePosixPath(name)
        if not re.fullmatch(r"[a-zA-Z0-9_./-]+", name) or (
            path.is_absolute() or ".." in path.parts or path.as_posix() != name
        ):
            raise InvalidBinding(f"Invalid manifest path: {name}")
        if not isinstance(digest, str) or not re.fullmatch("[0-9a-f]{64}", digest):
            raise InvalidBinding(f"Invalid SHA-256 for {name}")
    return entries


def _verify_tree(repo, revision, expected):
    entries = _tree(repo, revision, list(expected))
    if entries.keys() != expected.keys():
        raise InvalidBinding(
            f"Missing frozen input at {revision}: {expected.keys() - entries.keys()}"
        )
    for name, oid in entries.items():
        if _digest(_blob(repo, oid)) != expected[name]:
            raise InvalidBinding(f"Committed input differs at {revision}:{name}")


def _inventory(root):
    terms = _local(root, TERMS)
    raw = _local(root, MANIFEST)
    bound = re.findall(rb"Frozen source manifest SHA-256 `([0-9a-f]{64})`", terms)
    if bound != [_digest(raw).encode()]:
        raise InvalidBinding("Task must bind exactly one matching source manifest")
    manifest = _json(raw)
    if manifest.get("format") != 2:
        raise InvalidBinding("Expected file manifest format 2")
    working = _files(manifest, "working_inputs")
    native = _files(manifest, "adopted_native_inputs")
    if (
        working.keys() & native.keys()
        or not {DECLARATION, f"{DIRECTORY}/binding.py"} <= working.keys()
    ):
        raise InvalidBinding(
            "Overlapping inputs or missing candidate declaration and verifier"
        )
    if (working.keys() | native.keys()) & {TERMS, MANIFEST}:
        raise InvalidBinding("Task and manifest are bound externally")
    if native:
        pin = manifest.get("adopted_source_pin")
        if not isinstance(pin, str) or not re.fullmatch("[0-9a-f]{40}", pin):
            raise InvalidBinding("Missing the exact ancestor for preserved inputs")
        _verify_tree(root, pin, native)
    expected = {**working, **native}
    for name, digest in expected.items():
        if _digest(_local(root, name)) != digest:
            raise InvalidBinding(f"Local file differs: {name}")
    implementation = _files(_json(_local(root, DECLARATION)), "implementation")
    if not implementation or any(
        expected.get(name) != digest for name, digest in implementation.items()
    ):
        raise InvalidBinding(
            "Candidate implementation differs from the frozen manifest"
        )
    return {**expected, MANIFEST: _digest(raw), TERMS: _digest(terms)}


def check_local(root):
    """Return the candidate declaration SHA-256 after checking every local input."""
    return _inventory(Path(root).resolve())[DECLARATION]


def _canonical_pin(raw):
    rows = re.findall(rb"(?m)^\|\s*Pinned candidate(?=[\s|])[^\r\n]*", raw)
    if len(rows) != 1:
        raise InvalidBinding(
            "Canonical specification must have one Pinned candidate row"
        )
    row = rows[0]
    links = re.findall(
        rb"\]\(https://github\.com/FloatingPragma/observer-patch-holography/"
        rb"(?:tree|commit)/([a-f0-9]{40})(?=[/#)])",
        row,
    )
    hashes = set(re.findall(rb"(?<![0-9A-Za-z])[0-9a-f]{40}(?![0-9A-Za-z])", row))
    if len(set(links)) != 1 or hashes != set(links):
        raise InvalidBinding(
            "Ambiguous or missing official RER candidate link in canonical pin"
        )
    return links[0].decode()


def _first_parent_window(repo, head):
    if _git(repo, "rev-parse", "--is-shallow-repository").strip() != b"false":
        raise InvalidBinding(
            "A complete public OpenWave first-parent history is required"
        )
    lineage = _git(repo, "rev-list", "--first-parent", head).decode().splitlines()
    if BASE not in lineage:
        raise InvalidBinding("Official baseline absent from first-parent history")
    return reversed(lineage[: lineage.index(BASE) + 1])


def _landing(platform, head, terms, candidate_digest):
    # Walking each snapshot also detects a change followed by a restoration.
    # Cache committed blobs so unrelated commits do not reread canonical prose.
    anchor = None
    frozen = {}
    revision = None
    cache = {}

    def content(oid):
        if oid is None:
            return None
        if oid not in cache:
            cache[oid] = _blob(platform, oid)
        return cache[oid]

    for commit in _first_parent_window(platform, head):
        if anchor is None:
            at = _tree(platform, commit, [TASK])
            if content(at.get(TASK)) != terms:
                continue
            if commit == BASE:
                raise InvalidBinding("Exact terms already present at the baseline")
            anchor = commit

        tree = _tree(platform, commit, [TASK, POINTER, DECISION, CANON])
        if set(tree) != {TASK, POINTER, DECISION, CANON}:
            raise InvalidBinding(f"Missing acceptance evidence at {commit}")
        if commit == anchor:
            frozen = {path: tree[path] for path in (TASK, POINTER, DECISION)}
            pointer = content(tree[POINTER])
            if not re.fullmatch(rb"[0-9a-f]{40}\n", pointer):
                raise InvalidBinding("Invalid .source pointer at first landing")
            revision = pointer[:-1].decode()
            expected = dict(
                protocol=PROTOCOL,
                row="OW-30",
                decision="accepted",
                source_commit=revision,
                candidate_sha256=candidate_digest,
            )
            if _json(content(tree[DECISION])) != expected:
                raise InvalidBinding(
                    "First landing lacks exact affirmative row candidate acceptance"
                )
        elif any(tree[path] != oid for path, oid in frozen.items()):
            raise InvalidBinding(f"Terms, source or acceptance amended at {commit}")
        if _canonical_pin(content(tree[CANON])) != revision:
            raise InvalidBinding(
                f"Canonical pin differs from the accepted candidate at {commit}"
            )
    if not anchor:
        raise Unaccepted("Exact v3 terms have no landing on official OpenWave main")
    return anchor, revision


def check_public(source, platform):
    """Verify both official main histories and all local/released bytes, without selecting seeds."""
    source, platform = Path(source).resolve(), Path(platform).resolve()
    inventory = _inventory(source)
    platform_head = _head(platform, OPENWAVE)
    anchor, revision = _landing(
        platform, platform_head, _local(source, TERMS), inventory[DECLARATION]
    )
    source_head = _head(source, RER)
    if _git(source, "rev-parse", "--is-shallow-repository").strip() != b"false":
        raise InvalidBinding("Complete public RER ancestry is required")
    _git(source, "merge-base", "--is-ancestor", revision, source_head)
    _verify_tree(source, revision, inventory)
    if _head(platform, OPENWAVE) != platform_head or _head(source, RER) != source_head:
        raise InvalidBinding(
            "Official heads advanced during verification; recheck the same histories"
        )
    return Accepted(
        anchor, revision, inventory[DECLARATION], platform_head, source_head
    )
