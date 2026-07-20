#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


def run(command: list[str], cwd: Path) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=cwd, check=True)


def main() -> int:
    root = Path(__file__).resolve().parent
    fixture = root / "data" / "conditional_smdr_order1_fixture.json"
    receipt = root / "outputs" / "conditional_strict_1l_pole_map_receipt.json"
    run([sys.executable, str(root / "src" / "wz_pole_map.py"), str(fixture), "--output", str(receipt)], root)
    run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], root)
    run([sys.executable, str(root / "checker" / "check_receipt.py")], root)
    manifest_tool = root / "src" / "make_integration_manifest.py"
    run([sys.executable, str(manifest_tool)], root)
    run([sys.executable, str(manifest_tool), "--check"], root)
    result = json.loads(receipt.read_text(encoding="utf-8"))
    print(json.dumps({
        "status": "PASS",
        "classification": result["status"],
        "map_evaluated": result["map_evaluated"],
        "physical_promotion_allowed": result["physical_promotion_allowed"],
        "W": result["complex_poles"]["W"],
        "Z": result["complex_poles"]["Z"]
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
