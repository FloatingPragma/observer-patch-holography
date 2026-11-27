"""Type-check the experiment-local Lean evidence against the pinned project toolchain."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEAN = ROOT / "Lean"


def sources() -> list[Path]:
    # The GIFT snapshot is external comparison material with its own package.
    external = ROOT / "code/angular_sprint/external"
    files = {p for area in ("model", "code", "openwave", "evidence")
             for p in (ROOT / area).rglob("*.lean")
             if ".lake" not in p.parts and not p.is_relative_to(external)}
    for row in json.loads((ROOT / "docs/registers/OBSERVABLES.json").read_text())["rows"]:
        files.update(ROOT / p for p in row["lean"] if not p.startswith(("Lean/", "http")))
    return sorted(files)


def check(path: Path) -> bool:
    result = subprocess.run(["lean", str(path)], cwd=LEAN,
                            capture_output=True, text=True)
    log = result.stdout + result.stderr
    if result.returncode or "declaration uses `sorry`" in log:
        print(f"{path.relative_to(ROOT)}\n{log}", flush=True)
        return False
    return True


if __name__ == "__main__":
    if sys.argv[1:] != ["--in-lake-env"]:
        os.chdir(LEAN)
        os.execvp("lake", ["lake", "env", sys.executable, str(Path(__file__).resolve()), "--in-lake-env"])
    index = int(os.environ.get("EXPERIMENT_LEAN_SHARD", "0"))
    count = int(os.environ.get("EXPERIMENT_LEAN_SHARDS", "1"))
    if not 0 <= index < count:
        sys.exit("Invalid experiment proof shard")
    files = sources()[index::count]
    with ThreadPoolExecutor(max_workers=2) as pool:
        passed = sum(pool.map(check, files))
    print(f"Experiment Lean files checked: {passed}/{len(files)}")
    sys.exit(int(passed != len(files)))
