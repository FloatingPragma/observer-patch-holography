"""Replay same-action kinetic momentum, stress and force attachments."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer.scalar_response.experiment import canonical,compare_receipt
from ..full_fock_field.field import FieldGeometry
from ..coarse_energy_ward import ward as w
from . import momentum as m

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_PATHS=(
    'model/lattice_momentum_balance/momentum.py','model/lattice_momentum_balance/packet.py',
    'model/lattice_momentum_balance/test_momentum.py','model/lattice_momentum_balance/README.md',
    'model/coarse_energy_ward/ward.py','model/full_fock_field/field.py',
    'model/full_fock_field/local_simulation.py','model/full_fock_field/bounds.py',
    'model/interacting_field/lattice.py','Lean/QFT/CoarseEnergyWard.lean',
    'Lean/QFT/LatticeMomentumBalance.lean',
)


def encoded(operator):
    return [[str(mask),str(a),str(b)] for mask,(a,b) in sorted(operator.terms.items())]


def time_cell_diagnostic():
    """Independent complete 64-state spectral time integrals; no quadrature."""
    geometry=FieldGeometry((3,1,1))
    data=m.construct(geometry)
    h=data['H'].matrix(6)
    energies,vectors=np.linalg.eigh(h)
    differences=energies[:,None]-energies[None,:]
    duration=.75
    def evolve(operator):
        return vectors@(np.exp(1j*differences*duration)*(vectors.conj().T@operator@vectors))@vectors.conj().T
    def integral(operator):
        factor=duration*np.exp(1j*differences*duration/2)*np.sinc(differences*duration/(2*np.pi))
        return vectors@(factor*(vectors.conj().T@operator@vectors))@vectors.conj().T
    column=data['columns'][0]
    divergence=w.divergence(column['stress'],3)
    residual=omitted=sampled=0.
    for x,p in enumerate(column['local']):
        pm=p.matrix(6)
        fm=column['force'][x].matrix(6)
        balance=evolve(pm)-pm+integral(divergence[x].matrix(6))
        residual=max(residual,float(np.linalg.norm(balance-integral(fm),2)))
        omitted=max(omitted,float(np.linalg.norm(balance,2)))
        sampled=max(sampled,float(np.linalg.norm(balance-duration*fm,2)))
    if residual>1e-12 or omitted<=.01 or sampled<=1e-4:
        raise ArithmeticError('finite-cell force identity or negative control failed')
    return {'full_Fock_dimension':64,'duration':duration,'forced_balance_residual':residual,
            'omitted_force_residual':omitted,'point_force_sample_residual':sampled,
            'floating_admission':1e-12,'numerics_are_not_the_operator_certificate':True}


def run():
    geometry=FieldGeometry((3,3,3))
    data=m.construct(geometry)
    gram=m.force_gram(data['columns'])
    expected=tuple(tuple(F(9,4) if i==j else F() for j in range(3)) for i in range(3))
    if gram!=expected:
        raise ArithmeticError('all-three-direction source force separation failed')
    witnesses=m.no_go_matrix_elements(27)
    basis=m.conserved_quadratics(geometry)
    payload={'action':encoded(data['H']),
             'columns':[{'momentum':[encoded(z) for z in col['local']],
                         'stress':[[list(edge),encoded(z)] for edge,z in sorted(col['stress'].items())],
                         'force':[encoded(z) for z in col['force']]} for col in data['columns']]}
    error=m.cell_error(18,F(1,8),10000000,samples=10**12)
    return {'schema':'oph.lattice_momentum_balance.v1','issues':[1145],'gates':['CP-09'],'all_pass':True,
            'domain':{'action':'unchanged open complete-Fock Weyl/contact H',
                      'geometry':[3,3,3],'CAR_modes':54,'full_Fock_dimension':2**54,
                      'speed':'1/8','interaction':'1','momentum_convention':'dGamma(D_j), central difference at model spacing one',
                      'arbitrary_states_and_untouched_references':True,'number_projection':False,
                      'physical_momentum_or_time_identified':False},
            'exact':{'nonzero_momentum_columns':3,'local_forced_Ward_equations':81,
                     'kinetic_stress_faces_per_column':[len(c['stress']) for c in data['columns']],
                     'rectangular_partitions_per_column':data['partitions'],
                     'force_Gram_normalized_HS':[[str(x) for x in row] for row in gram],
                     'total_force_operator_norm_lower_each_axis':'3/2',
                     'general_force_Gram':'delta_ij U^2 E_i/8, E_i the number of axis-i spatial edges',
                     'force_is_quartic_Majorana':True,'kinetic_stress_is_quadratic':True,
                     'source_operator_payload_sha256':hashlib.sha256(canonical(payload).encode()).hexdigest()},
            'quadratic_classification':{'fixed_positive_U_no_go':True,
                     'proof':'N=1 restriction forces [h0,q]=0; N=2 doublon transitions force all intersite q coefficients zero; connected Pauli bonds transport one root 2x2 matrix',
                     'offsite_coefficients_isolated':len(witnesses),
                     'exact_doublon_witness_sha256':hashlib.sha256(canonical(witnesses).encode()).hexdigest(),
                     'conserved_homogeneous_number_preserving_quadratic_real_dimension':len(basis),
                     'basis':'total number and the three existing onsite-frame hidden-spin charges',
                     'additive_identity_allowed_separately':True,
                     'higher_body_or_anomalous_dressings_excluded':False},
            'source_bound':{'radius':18,'analytic_patch_sides':[37,37,37],'duration':'1/8',
                     'split_steps':10000000,'midpoint_samples':10**12,
                     'errors_exact':{k:str(v) for k,v in error.items()},
                     'resource_bounds_are_not_source_execution':True,
                     'input_readout_hardware_and_failure_energy_remain_separate':True},
            'independent_numerical_replay':time_cell_diagnostic(),
            'scope':{'closed':'three exact kinetic momentum/face-stress/contact-force columns, their finite cell and block composition, and the complete quadratic conservation class',
                     'remaining':['An interacting conserved momentum with any needed higher-body dressing and a symmetric Lorentz stress tensor.',
                                  'Physical chart/time and a controlled continuum mapping of source energy and all stress columns.',
                                  'Geometric modular kernels, entropy/area stationarity, Newton normalization and an Einstein-register inhabitant.'],
                     'no_conserved_interacting_Ward_or_Einstein_gate_accepted':True},
            'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--write',action='store_true')
    mode.add_argument('--check',nargs='?',const=OUTPUT,type=Path)
    args=parser.parse_args()
    result=run()
    if args.write:OUTPUT.write_text(canonical(result))
    if args.check:compare_receipt(strict_json.loads(args.check.read_text()),result)
    print(json.dumps({k:result[k] for k in ('all_pass','exact','quadratic_classification','source_bound','scope')},indent=2))


if __name__=='__main__':main()
