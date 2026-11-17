"""Replay both finite source calculations without comparison data."""
from pathlib import Path
import argparse
import hashlib
import json
from . import complete_c3, schur, relational

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCES=(
    'model/common_source/hierarchy_diagnosis/native_family_dynamics.py',
    'model/common_source/particle_identification/family_attachment.py',
    'model/common_source/native_family_record_memory.py',
    'model/common_source/native_repair_mass.py',
    'model/common_source/slow_geometry.py',
    'model/standard_model_source_map/koide_source.py',
    'model/charged_lepton_response/phase_selection.py',
    'model/native_modular_reference/local_preparation.py',
    'model/yukawa_source_dynamics/schur.py',
    'model/yukawa_source_dynamics/complete_c3.py',
    'model/yukawa_source_dynamics/certificate.py',
    'model/yukawa_source_dynamics/relational.py',
)


def report():
    return dict(schema='oph.yukawa_source_dynamics.v1',
                empirical_or_calibrated_inputs=[],
                native_completion=schur.report(),complete_cyclic_instrument=complete_c3.report(),
                retained_relational_records=relational.report(),
                result='Complete source response and canonical effective matrices computed; no three charged-fermion Yukawa spectra follow from these two operators.',
                domain='Specified finite completion response and native-record C3 routing; no universal exclusion of source dynamics.',
                source_bindings={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES})


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true'); args=parser.parse_args()
    data=report(); text=json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write:
        OUTPUT.write_text(text)
        print(str(OUTPUT))
    else: print(text,end='')


if __name__=='__main__': main()
