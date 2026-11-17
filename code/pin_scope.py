"""Which repository paths a custody pin may bind.

Receipts and checks bind the code that produces or verifies a result and the
data it reads. Prose and repository configuration (papers, documentation,
READMEs, workflow files, Git attribute files) change for editorial reasons,
and test modules check a receipt without producing it, so no receipt pins
them and no verifier compares them. This rule file itself (pin_scope.py)
governs custody, not results, and is exempt as well. A verifier filters a recorded pin map
with ``bound_pins`` before comparing it, so receipts written before this
rule still verify.

The implementation lives here so that code/ lanes import it as ``pin_scope``;
tools/ and model/ load the same file through ``tools/pin_scope.py``.
"""

from __future__ import annotations

from pathlib import PurePosixPath

EXEMPT_PREFIXES = ("docs/", "paper/", "book/", "flagship/", "extra/", ".github/")
EXEMPT_SUFFIXES = (".md", ".tex", ".pdf", ".bib", ".rst")
EXEMPT_NAMES = (".gitattributes", ".gitignore", "conftest.py", "pin_scope.py")


def exempt(path) -> bool:
    """True for prose, configuration and test modules, which pins never bind."""
    text = str(path).replace("\\", "/")
    while text.startswith("./"):
        text = text[2:]
    name = PurePosixPath(text).name
    return (
        text.startswith(EXEMPT_PREFIXES)
        or text.lower().endswith(EXEMPT_SUFFIXES)
        or name in EXEMPT_NAMES
        or name.upper().startswith("README")
        or (name.startswith("test_") and name.endswith(".py"))
    )


def bound(paths):
    """The bindable subset of ``paths``, in the given order."""
    return [path for path in paths if not exempt(path)]


def bound_pins(pins):
    """The bindable entries of a recorded ``{path: digest}`` map.

    A value that is not a map is returned unchanged, so the caller's own
    comparison rejects it.
    """
    if not isinstance(pins, dict):
        return pins
    return {path: digest for path, digest in pins.items() if not exempt(path)}


def repo_sha256(path) -> str:
    """SHA-256 of a repository file, by its root-relative path.

    Code that records the digest of a parent receipt reads it here instead of
    keeping the value as a constant, so a regenerated parent never needs a
    hand-edited hash.
    """
    import hashlib
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    return hashlib.sha256((root / str(path)).read_bytes()).hexdigest()


def scoped(value):
    """``value`` without exempt paths, keeping its type.

    A comparison with a filtered pin map applies this to the other side: a
    map loses its exempt keys, a set, list or tuple of paths its exempt
    members. Any other value is returned unchanged.
    """
    if isinstance(value, dict):
        return bound_pins(value)
    if isinstance(value, (set, frozenset, list, tuple)) and all(isinstance(v, str) for v in value):
        return type(value)(bound(value))
    return value
