"""Retain the same-action pair findings, complete records and their limits."""
from __future__ import annotations
import hashlib
from pathlib import Path
from . import mass, protection, identity, spectrum

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


def source_bindings():
    paths=list(HERE.glob('*.py'))+[HERE/name for name in (
        'README.md','PROTECTION.md','IDENTITY.md','SPECTRUM.md')]
    paths.extend(ROOT/name for name in (
        'model/interacting_field/lattice.py','model/field_stress/energy.py',
        'model/particle_clock/clock.py','model/completion/particle_premise_campaign.json'))
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def run():
    return {'schema':'oph.bound_pair.v1',
            'campaign_id':'OPH-SAME-ACTION-PARTICLE-PREMISES-2026-10-29-v1',
            'scope':'SUPPLIED_WEYL_CONTACT_ACTION_COMPLETE_TWO_FERMION_SECTOR',
            'source_bindings':source_bindings(),
            'mass':mass.certificate(),
            'protection_and_transport':protection.report(),
            'identity_and_currents':identity.report(),
            'spectrum_and_records':spectrum.report(),
            'integrated_gate_accepted':False,'physical_comparison_count_added':0,
            'physical_species_identified':False,'nontrivial_topological_defect_assigned':False,
            'full_native_field_encoding_executed':False}
