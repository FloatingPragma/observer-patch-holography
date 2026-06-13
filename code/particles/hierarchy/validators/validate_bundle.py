#!/usr/bin/env python3
import subprocess, sys, pathlib, json
root = pathlib.Path(__file__).resolve().parents[1]
checks = [
    [sys.executable, str(root/"validators/validate_dag.py"), str(root/"certificates/DAG_U.json")],
    [sys.executable, str(root/"validators/validate_ru_interval_certificate.py"), str(root/"certificates/R_U_krawczyk_certificate.json")],
    [sys.executable, str(root/"validators/validate_global_repair_tick_certificate.py"), str(root/"certificates/R_N_global_repair_tick_certificate.json")],
    [sys.executable, str(root/"validators/validate_issue_332_rg_naturality.py"), str(root/"issue_332_rg_naturality_certificate.json")],
    [sys.executable, str(root/"validators/validate_joint_fixed_point_certificate.py"), str(root/"certificates/R_PN_joint_fixed_point_certificate_report.json")],
    [sys.executable, str(root/"validators/validate_manifest.py"), str(root/"manifest.json")]
]
results=[]
for cmd in checks:
    p=subprocess.run(cmd, capture_output=True, text=True)
    results.append({"cmd": cmd, "returncode": p.returncode, "stdout": p.stdout.strip(), "stderr": p.stderr.strip()})
print(json.dumps(results, indent=2))
sys.exit(0 if all(r["returncode"]==0 for r in results) else 1)
