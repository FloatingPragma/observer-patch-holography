#!/usr/bin/env python3
"""Standalone replay of the frozen common member, with row-specific summaries."""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PACKAGE=ROOT/'evidence/openwave/row_30_klein_gordon/work/2026-11-14_codex_common_cubic'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--row',choices=['29','30','both'],default='both')
    ap.add_argument('--custody-only',action='store_true')
    args=ap.parse_args()
    pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())
    for relative,expected in pins.items():
        assert digest(ROOT/relative)==expected,relative
    if args.custody_only:
        print(json.dumps(dict(custody='pass',files=len(pins))))
        return
    started=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='oph_common_replay_') as td:
        work=Path(td)/'source'
        shutil.copytree(PACKAGE,work)
        subprocess.run([sys.executable,str(work/'probe.py')],check=True,
                       stdout=subprocess.PIPE,text=True,cwd=work)
        r=json.loads((work/'RECEIPT.json').read_text())
        byte_replay={name:digest(work/name)==digest(PACKAGE/name)
                     for name in ['RECEIPT.json','RECORDS.json.gz']}
        assert all(byte_replay.values()),byte_replay
        assert r['all_registered_checks']=='pass'
        assert r['actual_native_calls']==169364
        summary=dict(candidate=r['candidate'],row=args.row,
            actual_native_calls=r['actual_native_calls'],
            population_block_length=r['population_block_length'],
            population_is_actual_tape=False,byte_replay=byte_replay,
            unfinished_clock_intervals=r['final_open_clock_intervals_retained'],
            clock_ticks=r['same_tape_completed_clock_ticks'],
            first_tick_censor_tail=r['first_clock_tail_at_tape_cap'],
            interpretation='finite source replay plus separately bound full-limit theorem',
            official='untested',elapsed_seconds=time.monotonic()-started)
        if args.row in ('29','both'):
            rows=[x for x in r['population'] if x['internal_lambda']<1e-9]
            summary['row29']=dict(
                electric_magnetic_relative_error=max(x['electric_magnetic_relative_error_including_contact'] for x in rows),
                negative_finite_contact_eigenvalues_retained=[x['reflected_electric_eigenvalues_with_contact'][0] for x in rows],
                source_changed_iid_separated_covariance=[x['changed_iid_separated_covariance'] for x in rows],
                full_two_polarization_causal_limit='bound written proof; not inferred from curl null eigenvalue')
        if args.row in ('30','both'):
            summary['row30']=dict(native_internal_bands=r['native_internal_bands'],
                max_source_eigenvector_residual=r['max_source_eigenvector_residual'],
                max_covariance_relative_error=r['max_covariance_relative_error'],
                max_first_hit_population_error=r['max_first_hit_population_error'],
                full_common_speed_KG_limit='bound written proof; finite N=3 is not the continuum')
        print(json.dumps(summary,indent=2,sort_keys=True))

if __name__=='__main__':main()
