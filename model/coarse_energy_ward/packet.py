"""Replay the exact same-action cell Ward/coarse-graining attachment."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..full_fock_field.field import FieldGeometry
from . import ward as w

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_PATHS=(
    'model/coarse_energy_ward/ward.py','model/coarse_energy_ward/packet.py',
    'model/coarse_energy_ward/test_ward.py','model/coarse_energy_ward/README.md',
    'model/full_fock_energy/energy.py','model/field_stress/energy.py',
    'model/full_fock_field/field.py','model/interacting_field/lattice.py',
    'model/full_fock_field/local_simulation.py','model/full_fock_field/bounds.py',
    'Lean/QFT/CoarseEnergyWard.lean',
)


def encoded(polynomial):
    return [[str(mask),str(real),str(imag)] for mask,(real,imag) in sorted(polynomial.terms.items())]


def time_cell_diagnostic():
    """Independent complete 64-state spectral integral, with no quadrature."""
    geometry=FieldGeometry((3,1,1))
    data=w.certify(geometry)
    h=data['H'].matrix(geometry.modes)
    weights,basis=np.linalg.eigh(h)
    gap=weights[:,None]-weights[None,:]
    def evolve(operator,time):
        return basis@(np.exp(1j*gap*time)*(basis.conj().T@operator@basis))@basis.conj().T
    def integral(operator,start,stop):
        width=stop-start
        factor=width*np.exp(1j*gap*(start+stop)/2)*np.sinc(gap*width/(2*np.pi))
        return basis@(factor*(basis.conj().T@operator@basis))@basis.conj().T
    duration=.75
    currents={edge:value.matrix(geometry.modes) for edge,value in data['routed'].items()}
    maximum,wrong,composition=0.,0.,0.
    for x,energy in enumerate(data['local']):
        hx=energy.matrix(geometry.modes)
        flux=np.zeros_like(hx)
        sample=np.zeros_like(hx)
        for (a,b),current in currents.items():
            sign=int(x==a)-int(x==b)
            flux += sign*integral(current,0,duration)
            sample += sign*duration*current
            composition=max(composition,float(np.linalg.norm(integral(current,0,duration)-
                                       integral(current,0,.25)-integral(current,.25,duration),2)))
        delta=evolve(hx,duration)-hx
        maximum=max(maximum,float(np.linalg.norm(delta+flux,2)))
        wrong=max(wrong,float(np.linalg.norm(delta+sample,2)))
    if maximum > 1e-12 or composition > 1e-12 or wrong < 1e-4:
        raise ArithmeticError('integrated Ward diagnostic or point-sample negative control failed')
    return {'full_Fock_dimension':64,'duration':duration,
            'cell_identity_operator_residual':maximum,
            'time_composition_operator_residual':composition,
            'endpoint_sample_counterexample_norm':wrong,
            'floating_admission':1e-12,'numerics_are_not_the_exact_certificate':True}


def run():
    geometry=FieldGeometry((3,3,3))
    data=w.certify(geometry)
    encoded_action=encoded(data['H'])
    current_degrees=sorted(set(mask.bit_count() for value in data['routed'].values() for mask in value.terms))
    budget=w.readout_error(18,F(1,8),10000000,samples=10**12)
    return {
        'schema':'oph.coarse_energy_ward.v1','issues':[1145],'gates':['CP-09'],'all_pass':True,
        'domain':{'action':'unchanged complete-Fock Weyl/contact H; c=1/8, U=1',
                  'geometry':[3,3,3],'spatial_sites':27,'CAR_modes':54,
                  'full_Fock_dimension':2**54,'basis':'exact rational faithful Clifford monomials',
                  'number_projection':False,'state_or_Gibbs_assumption':False,
                  'Hamiltonian_model_time_only':True},
        'exact':{'unordered_pairs_checked':351,'nonzero_energy_current_pairs':len(data['currents']),
                 'nearest_neighbor_faces':len(data['routed']),
                 'local_Ward_equations':27,'rectangular_block_partitions':data['rectangular_partitions'],
                 'routed_current_Majorana_degrees':current_degrees,
                 'action_terms':len(encoded_action),
                 'action_sha256':hashlib.sha256(canonical(encoded_action).encode()).hexdigest(),
                 'all_shortest_paths_equally_weighted':True,
                 'finite_cell_flux_is_time_integrated':True,
                 'block_divergence_commutes_exactly':True},
        'source_readout_bound':{'radius':18,'analytic_patch_sides':[37,37,37],
                               'bound_is_separate_from_the_27_site_symbolic_instance':True,'duration':'1/8','split_steps':10000000,'midpoint_samples':10**12,
                               'errors_exact':{k:str(v) for k,v in budget.items()},
                               'ambient_volume_independent':True,
                               'patch_reference_restriction_supplied':True,
                               'numerical_hardware_state_preparation_and_readout_errors_not_included':True},
        'independent_numerical_replay':time_cell_diagnostic(),
        'scope':{'closed':'nearest-face energy Ward column and its exact finite spatial/time block composition on the supplied action',
                 'remaining':['The other three momentum columns and a symmetric Lorentz stress tensor.',
                              'Geometric modular cap weights, proper clock/length and same-action continuum scaling.',
                              'Area normalization, generalized-entropy stationarity, positive Newton scale and all Einstein G1-G6 premises.',
                              'Actual source preparation and readout of the stated arbitrary local input; controller/failure-space energy.'],
                 'no_full_Einstein_Ward_or_Einstein_gate_accepted':True},
        'source_sha256':{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--write',action='store_true')
    modes.add_argument('--check',nargs='?',const=OUTPUT,type=Path)
    args=parser.parse_args()
    result=run()
    if args.write:
        OUTPUT.write_text(canonical(result))
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text()),result)
    print(json.dumps({key:result[key] for key in ('all_pass','exact','source_readout_bound','scope')},indent=2))


if __name__=='__main__':
    main()
