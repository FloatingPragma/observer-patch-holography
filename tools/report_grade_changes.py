"""Show register and candidate grade changes from a Git revision to the current checkout."""

import argparse
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
LEDGER = Path("docs/registers/OBSERVABLES.json")
CANDIDATES = Path("model/candidate_laws")


def git(root, *args):
    return subprocess.check_output(["git", *args], cwd=root)


def is_grade_file(path):
    return path == LEDGER or (path.name == "EVALUATION.json" and path.parent.parent == CANDIDATES)


def grades(path, raw):
    data = json.loads(raw)
    result = {}

    def add(scope, key, status):
        if (scope, key) in result or not isinstance(status, str) or not status:
            raise ValueError(f"{path}: invalid or duplicate status for {scope} {key}")
        result[scope, key] = status

    if path == LEDGER:
        for row in data["rows"]:
            add("register", row["id"], row["status"])
    else:
        law = data["candidate"]
        if law != path.parent.name:
            raise ValueError(f"{path}: candidate identifier {law!r} differs from directory")
        for section in data.get("requirements", ()):
            for row in section["criteria"]:
                add(f"{law}/criterion", row["id"], row["status"])
        for row in data.get("premises", ()):
            add(f"{law}/premise", row["id"], row["status"])
        for key, row in data.get("joint", {}).items():
            add(f"{law}/joint", key, row["status"])
    return result


def grade_changes(root, base):
    revision = git(root, "rev-parse", "--verify", "--end-of-options", f"{base}^{{commit}}").decode().strip()
    # List old deposits even for laws deleted from the working tree. Only read
    # their grade files, never the simulations or large physical acquisitions.
    tree = git(root, "ls-tree", "-r", "-z", "--name-only", revision, "--", str(LEDGER), str(CANDIDATES))
    old_files = {p for raw in tree.split(b"\0") if raw and is_grade_file(p := Path(raw.decode()))}
    current_files = {p.relative_to(root) for p in (root / CANDIDATES).glob("*/EVALUATION.json")}
    if (root / LEDGER).is_file():
        current_files.add(LEDGER)
    changes = []
    for path in sorted(old_files | current_files):
        old = grades(path, git(root, "show", f"{revision}:{path.as_posix()}")) if path in old_files else {}
        new = grades(path, (root / path).read_bytes()) if path in current_files else {}
        for scope, key in sorted(old.keys() | new.keys()):
            before, after = old.get((scope, key)), new.get((scope, key))
            if before != after:
                changes.append((scope, key, before, after))
    return changes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, help="Git commit, tag or branch to compare against")
    args = parser.parse_args()
    changes = grade_changes(ROOT, args.base)
    print(f"{len(changes)} grade changes from {args.base} to the current checkout.")
    for scope, key, before, after in changes:
        print(f"{scope} {key}: {before or '[absent]'} → {after or '[absent]'}")


if __name__ == "__main__":
    main()
