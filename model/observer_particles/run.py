#!/usr/bin/env python3
"""Replay or verify source-linked matter and conditional field constructions."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
sys.dont_write_bytecode=True

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PRODUCERS=(
    'gauge_particle_identification','representation_identification',
    'fermion_species_identification','spin_matter_join','scalar_bilinear_route',
    'family_naturality','unified_yukawa_action','neutrino_action_extension',
    'entangled_alpha','coefficient_selection','electroweak_mass_form',
    'spacetime_correspondence','shared_spatial_matter','geometric_matter_join',
    'cochain_completion','geometric_stationarity',
)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def premise_witnesses():
    """Restricted finite supply; no global physical premise is retired."""
    entries=(
        ('matter_decoder','representation_identification',['PR-59'],
         'Exact conditional fifteen-state charge/Casimir decoder and all source-current brackets.'),
        ('one_higgs_vertices','unified_yukawa_action',['PR-50','PR-67'],
         'One-doublet invariant lines and all 648 family/current/component identities; coefficients are unselected.'),
        ('spin_matter_controls','spin_matter_join',['PR-47'],
         'Independent finite Spin and matter controls with complete readout; physical chirality is unselected.'),
        ('shared_spatial_operator','shared_spatial_matter',['PR-53','PR-71'],
         'Complete conditional CAR field shares the ideal scalar spatial generator and incident coordinate; physical field and refinement are separate.'),
        ('full_cochain_control','cochain_completion',['PR-53','PR-71'],
         'Complete cone cochains preserve the scalar block and remove truncated edge cycles; this cone is distinct from the radial cylinder.'),
        ('metric_and_null_dictionary','geometric_matter_join',['PR-25','PR-53'],
         'Common pointwise tetrad, field symbols and rank-nine null dictionary; independent Einstein balances are not supplied.'),
        ('native_entropy_and_area','geometric_stationarity',['PR-26','PR-27'],
         'Finite constrained entropy and geometric area/volume identities; physical vacuum fixed-volume stationarity is not supplied.'),
        ('entangled_pair_scalar','entangled_alpha',['PR-17','PR-54'],
         "Wissner-Gross's finite spectral scalar, independently reconstructed with complete outcome controls; physical Thomson and same-cell identifications are separate."),
    )
    return {'schema':'oph.observer_particles.finite_premise_witnesses.v1',
        'classification':'constructed_finite_instances_and_restricted_partial_supply',
        'physical_premise_removal_claimed':False,
        'rows':[{'id':name,'premise_ids':ids,'producer':producer,
                 'receipt':f'model/observer_particles/results/{producer}.json',
                 'finite_supply':scope,'physical_premise_removed':False}
                for name,producer,ids,scope in entries]}


def verify():
    statuses={};source_count=artifact_count=0
    for name in PRODUCERS:
        record=json.loads((HERE/'results'/f'{name}.json').read_text())
        status=record.get('status','')
        if 'PASS' not in status.upper() and name!='gauge_particle_identification':
            raise ValueError(f'Unexpected result status: {name}: {status}')
        if name=='gauge_particle_identification' and status!='EXACT_SOURCE_RESPONSE_INTERNAL_GAUGE_MULTIPLETS_IDENTIFIED':
            raise ValueError('Unexpected gauge response result')
        assert record['producer_sha256']==sha(HERE/f'{name}.py'),name
        for relative,expected in bound_pins(record['source_sha256']).items():
            path=(ROOT/relative).resolve()
            assert not Path(relative).is_absolute() and path.is_relative_to(ROOT)
            assert sha(path)==expected,(name,relative)
            source_count+=1
        for relative,expected in record.get('artifacts',{}).items():
            path=(HERE/'results'/relative).resolve()
            assert path.is_relative_to(HERE/'results') and sha(path)==expected,(name,relative)
            artifact_count+=1
        assert record['evidence_contract']['physical_particle_masses_source_selected'] is False
        statuses[name]=status
    return {'schema':'oph.observer_particles.verification.v1',
        'status':'PASS_REPOSITORY_LOCAL_MATTER_AND_FIELD_CONTRACTS',
        'producers':statuses,'computational_source_bindings':source_count,
        'artifact_bindings':artifact_count,
        'physical_mass_predictions_promoted':False,
        'physical_electromagnetic_identity_established':False,
        'geometric_generalized_entropy_stationarity_established':False,
        'scope':'The results instantiate finite hypotheses on declared branches. Physical selection and continuum premises retain their independent requirements.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-only',action='store_true')
    parser.add_argument('--write-views',action='store_true',
                        help='Verify retained results and render summary/premise views without replaying producers.')
    parser.add_argument('--only',nargs='+',choices=PRODUCERS)
    args=parser.parse_args()
    if not args.verify_only and not args.write_views:
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        for name in PRODUCERS:
            if args.only and name not in args.only:continue
            run=subprocess.run([sys.executable,'-B',str(HERE/f'{name}.py')],cwd=ROOT,
                               env=env,capture_output=True,text=True)
            if run.returncode:
                print(run.stdout+run.stderr,file=sys.stderr)
                raise RuntimeError(f'Producer failed: {name}')
            print(name+': PASS',flush=True)
    result=verify()
    if not args.verify_only:
        (HERE/'results/verification.json').write_text(json.dumps(result,indent=2)+'\n', encoding='utf-8', newline='\n')
        (HERE/'results/premise_witnesses.json').write_text(json.dumps(premise_witnesses(),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
