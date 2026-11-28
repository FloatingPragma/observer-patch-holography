"""Compute standard adiabatic radiation transfer for a finite stock of cold matter."""

import numpy as np

RATIOS = {
    "zero": 0,
    "sharp_regular": 17 / 260,
    "federation_pooled_classical": 17 / 82,
    "federation_pooled_coherent": 17 / 80,
    "ordered_sharp": 1 / 4,
    "federation_matching": 3 / 7,
    "sharp_classical": 0.5,
    "sharp_coherent": 17 / 32,
    "max_service": 1.0,
    "double_service_control": 2.0,
    "standard_control": 0.12 / 0.02237,
}
BARYONS = (0.02237, 0.0215, 0.0235)
LMAX = 2800


def compute(baryon, ratio):
    import camb

    p = camb.CAMBparams()
    helium = float(np.asarray(camb.bbn.get_predictor().Y_He(baryon, 0.0)).item())
    p.set_cosmology(
        H0=70,
        ombh2=baryon,
        omch2=baryon * ratio,
        mnu=0.06,
        omk=0,
        tau=0.055,
        YHe=helium,
    )
    p.InitPower.set_params(As=2.1e-9, ns=1.0)
    p.set_for_lmax(LMAX, lens_potential_accuracy=1)
    p.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(p)
    angular = results.get_total_cls(lmax=LMAX, CMB_unit="muK")
    return angular[:, [0, 3, 1]], float(results.z_eq)
