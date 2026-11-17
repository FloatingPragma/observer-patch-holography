#!/usr/bin/env python3
"""Consume the conditional gauge solution in the existing EW mass-form theorem.

This computes a mass-form shape per unspecified vacuum norm, not physical
pole masses. It neither selects a Higgs vacuum nor fixes the kinetic/pole
matching assumptions retained by the owning theorem.
"""
from evidence_contract import public_receipt
import hashlib
import json
import math
from pathlib import Path
import numpy as np
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    source = HERE/'results/coefficient_selection.json'
    receipt = json.loads(source.read_text())
    branch = receipt['conditional_matter_trace_Chain_P_entropy_gauge_join']
    g = math.sqrt(4*math.pi*branch['alpha2'])
    gp = math.sqrt(4*math.pi*branch['alphaY'])
    # The exact theorem's quadratic form, divided by the unspecified vev^2.
    # Parameter ordering is (a1,a2,a3,b); gauge-field kinetic normalization
    # remains the owning theorem's declared Euclidean convention.
    matrix = np.array([[g*g,0,0,0],[0,g*g,0,0],
                       [0,0,g*g,-g*gp],[0,0,-g*gp,gp*gp]])/4
    photon = np.array([0,0,gp,g])
    broken = np.array([0,0,g,-gp])
    w = g*g/4
    z = (g*g+gp*gp)/4
    eig = np.linalg.eigvalsh(matrix)
    assert np.linalg.norm(matrix@photon) < 1e-14
    assert np.max(abs(eig-np.array([0,w,w,z]))) < 1e-14
    assert abs(broken@matrix@broken/(broken@broken)-z) < 1e-14
    ratio = w/z
    assert abs(ratio-(1-branch['sin2_theta_same_scale_tree'])) < 1e-14
    theorem = ROOT/'Lean/Screen/ElectroweakBreakingComposition.lean'
    result = {
        'status':'PASS_CONDITIONAL_GAUGE_SOLUTION_TO_ELECTROWEAK_MASS_FORM',
        'producer_sha256':sha(__file__),
        'source_sha256':bound_pins({str(p.relative_to(ROOT)):sha(p) for p in (source,theorem)}),
        'empirical_mass_or_coupling_inputs':[],
        'mass_form_divided_by_unspecified_vev_squared':matrix.tolist(),
        'eigenvalues_per_vev_squared':eig.tolist(),
        'charged_to_neutral_mass_form_ratio':ratio,
        'square_root_ratio_if_canonical_mass_interpretation':math.sqrt(ratio),
        'exact_photon_null_direction':'(0,0,gprime,g)',
        'remaining_continuous_inputs':['Higgs vacuum norm','scalar potential parameters',
            'Yukawa coefficients/matrices','physical kinetic and pole matching'],
        'scope':{'same_conditional_gauge_branch':True,'physical_Higgs_vacuum_selected':False,
            'absolute_mass_scale_selected':False,'quantum_pole_masses_computed':False,
            'charged_neutral_mass_ratio_is_observed_W_Z_prediction':False},
    }
    (HERE/'results/electroweak_mass_form.json').write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'mass_form_ratio':ratio,
        'conditional_square_root_ratio':math.sqrt(ratio)},indent=2))


if __name__ == '__main__':
    main()
