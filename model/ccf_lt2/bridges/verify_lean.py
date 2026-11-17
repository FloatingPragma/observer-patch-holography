"""Check the two abstract finite-source Lean interfaces on the pinned toolchain."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
MODULES=('kinematics/FiniteSource.lean','matter/SectorBridge.lean')
OUTPUT=HERE/'LEAN_RESULTS.json'
ALLOWED_AXIOMS={'propext','Classical.choice','Quot.sound'}
EXPECTED={
 'kinematics/FiniteSource.lean': {'CCF2.Kinematics.'+name for name in (
    'finite_source_gibbs','stationary_density_fixed','no_half_electron_in_native_charge_lattice','finite_source_unitary')},
 'matter/SectorBridge.lean': {'CCF2MatterBridge.'+name for name in (
    'support_commutator_entry','support_preserves_sector','support_lapseJoin','support_sum','gibbs_contrast_offset_cancels')} |
    {'OPH.EinsteinBranch.'+name for name in ('geometry_responds_to_all_charge','exists_dark_source_of_geometric_excess',
     'anomalousStress_tendsto_zero','deepProfile_characterization','perCut_a0_dictionary')},
}


def inputs():
    paths={HERE/name for name in MODULES}
    todo=list(paths)
    while todo:
        path=todo.pop()
        for module in re.findall(r'^import\s+(\S+)',path.read_text(),re.M):
            dependency=ROOT/'Lean'/Path(*module.split('.')).with_suffix('.lean')
            if dependency.is_file() and dependency not in paths:
                paths.add(dependency);todo.append(dependency)
    paths.update(ROOT/'Lean'/name for name in ('lean-toolchain','lakefile.lean','lake-manifest.json'))
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def axiom_reports(output):
    reports={}
    for theorem,axioms in re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]",output):
        names=[name.strip() for name in axioms.split(',') if name.strip()]
        if not set(names)<=ALLOWED_AXIOMS:raise ValueError('Unaccepted proof axiom: '+theorem)
        reports[theorem]=names
    if not reports:raise ValueError('Missing Lean axiom reports')
    return reports


def validate_saved():
    data=json.loads(OUTPUT.read_text())
    if data.get('all_passed') is not True or data.get('source_sha256')!=inputs():
        raise ValueError('Lean source proof receipt is stale')
    if set(data.get('modules',{}))!=set(MODULES):raise ValueError('Lean module omitted')
    for module,reports in data['modules'].items():
        if set(reports)!=EXPECTED[module]:raise ValueError('Lean theorem omitted')
        if not reports or any(not set(names)<=ALLOWED_AXIOMS for names in reports.values()):
            raise ValueError('Unaccepted proof axiom')
    return data


def run():
    before=inputs();modules={}
    for name in MODULES:
        command=['lake','env','lean','../model/ccf_lt2/bridges/'+name]
        result=subprocess.run(command,cwd=ROOT/'Lean',capture_output=True,text=True)
        if result.returncode:raise RuntimeError(result.stdout+result.stderr)
        modules[name]=axiom_reports(result.stdout+result.stderr)
        if set(modules[name])!=EXPECTED[name]:raise ValueError('Lean theorem omitted')
    if inputs()!=before:raise ValueError('Lean sources changed during checking')
    return {'all_passed':True,'scope':'Abstract finite matrices; executable extraction and physical interpretation are explicit hypotheses.',
            'modules':modules,'source_sha256':before}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true');args=parser.parse_args()
    actual=run()
    if args.write:OUTPUT.write_text(json.dumps(actual,indent=2,sort_keys=True)+'\n')
    elif validate_saved()!=actual:raise ValueError('Lean proof receipt differs')
    print('PASS_CCF_LT2_FINITE_LEAN_BRIDGES')


if __name__=='__main__':main()
