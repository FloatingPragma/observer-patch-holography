#!/usr/bin/env python3
"""Replay the actual native atomic-controller production source."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PACKAGE=ROOT/'evidence/openwave/row_31_orbital_levels/work/2026-11-14_codex_production_field'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--custody-only',action='store_true');args=ap.parse_args()
    pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())
    for name,value in pins.items():assert sha(ROOT/name)==value,name
    if args.custody_only:print(json.dumps({'custody':'pass','files':len(pins)}));return
    t=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='oph_orbital_replay_') as td:
        work=Path(td)/'source';shutil.copytree(PACKAGE,work)
        subprocess.run([sys.executable,str(work/'probe.py')],cwd=work,check=True,stdout=subprocess.PIPE,text=True)
        hashes={n:sha(work/n)==sha(PACKAGE/n) for n in ['RECEIPT.json','RECORDS.json.gz']}
        assert all(hashes.values())
        r=json.loads((work/'RECEIPT.json').read_text());assert r['all_registered_checks']=='pass'
        assert r['native_calls']==347533
        print(json.dumps(dict(candidate=r['candidate'],actual_native_calls=r['native_calls'],
            actual_inner_calls=r['actual_inner_calls'],auxiliary_holds=r['auxiliary_holds'],
            byte_replay=hashes,source_standing_ground_energies=[b['ground']['slice_energy'] for b in r['population'][0]['bands']],
            disabled_standing_ground_energies=[b['ground']['slice_energy'] for b in r['population'][1]['bands']],
            max_covariance_relative_error=r['max_covariance_relative_error'],max_clock_error=r['max_clock_error'],
            unchanged_massless_sector_error=r['massless_sector_change'],actual_clock_ticks=r['actual_clock_ticks'],
            final_open_intervals=r['actual_final_open_intervals'],actual_censors=r['actual_censors'],
            population_is_actual_long_tape=False,hydrogenic_shells='bound written controlled-limit theorem; not measured by this small run',
            elapsed_seconds=time.monotonic()-t,official='untested'),indent=2,sort_keys=True))

if __name__=='__main__':main()
