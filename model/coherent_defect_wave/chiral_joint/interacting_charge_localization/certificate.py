"""Custody of full-source charge localization and screening controls."""
import hashlib
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from .geometry import certificate as geometry_certificate
from .screening import certificate as screening_certificate
from .localization import certificate as localization_certificate

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).parent
DEPENDENCIES=tuple('model/coherent_defect_wave/chiral_joint/'+name+'/'+file
    for name in ('interacting_charged_transport','interacting_charged_response',
        'interacting_ground_preparation','atomic_reference_control','local_reference','matter_induced_plaquette')
    for file in ('PROOF.md','results.json'))+(
        'code/electromagnetism/discrete_coulomb_green.py',
        'code/electromagnetism/whitney_interacting_quantum.py',
        'paper/tex_fragments/WHITNEY_MAXWELL_DYNAMICS.tex')


def input_hashes():
    seeds={str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(DEPENDENCIES);seeds.add('model/requirements.txt')
    for name in seeds:
        if not (ROOT/name).is_file():raise FileNotFoundError('missing declared charge localization input: '+name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(seeds))}


def certificate():
    return {'schema':'oph.same_H_interacting_charge_localization.v1',
        'actual_graph_and_full_link_charge_control':geometry_certificate(),
        'screened_remote_pair_counterexample':screening_certificate(),
        'full_source_charged_bins':{family:localization_certificate(family) for family in ('ec','L')},
        'input_sha256':input_hashes()}
