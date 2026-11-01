"""Replay the complete-algebra certificate and prepared-reference energy read."""
from __future__ import annotations

import json

from ..field_thermal.thermal import ThermalField, failure_probability
from ..full_fock_field.field import FieldGeometry
from .energy import certify, charged_probe_certificate, grand_canonical
from .instrument import execute


def build():
    certificates = [certify(FieldGeometry(sides))[0] for sides in
                    ((2, 1, 1), (1, 2, 1), (1, 1, 2), (3, 1, 1), (2, 2, 1))]
    case = ThermalField()
    prepared = case.prepare()
    columns = prepared['state'].reshape(case.dimension, case.dimension).T
    accounting, _, _ = execute(case.geometry, columns)
    preparation_abort = failure_probability(prepared['first_export_failures'])
    accounting['source_preparation_executed'] = True
    accounting['preparation_resources'] = prepared['resources']
    accounting['preparation_export_abort_probability'] = preparation_abort
    accounting['preparation_complete_purification_error'] = prepared['purification_error']
    accounting['complete_preparation_evolution_read_probability'] = (
        accounting['complete_probability']+preparation_abort)
    accounting['reference'] = 'complete sixteen-state Gibbs reference, beta=1, mu=0'
    accounting['time_parameter'] = 'Hamiltonian model coordinate, hbar_model=1'
    return {'schema': 'oph.full_fock_energy.v1', 'certificates': certificates,
            'chemical_potential_accounting': charged_probe_certificate(grand_canonical()),
            'prepared_reference_energy_read': accounting,
            'physical_mass_or_gravity_premises_discharged': False}


if __name__ == '__main__':
    print(json.dumps(build(), indent=2))
