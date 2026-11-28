"""Compare a single written layer mapped to cold matter with measured Planck TT bandpowers."""

import argparse
import hashlib
import json
from pathlib import Path

import camb
import numpy as np
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
PILOT = HERE.parents[1] / "tilt_posterior" / "planck_peaks.json"
DATA = HERE / "planck_tt.dat"
NATIVE = Path(
    "model/experiments/interval_echo_shell_flux/product_address/green_response/material_recoil/overlap_source/native_boundary/budget.json"
)


def spectrum(background, ns, amplitude):
    p = camb.CAMBparams()
    helium = float(
        np.asarray(camb.bbn.get_predictor().Y_He(background["ombh2"], 0.0)).item()
    )
    p.set_cosmology(**background, YHe=helium)
    p.InitPower.set_params(As=amplitude, ns=ns)
    p.set_for_lmax(2600, lens_potential_accuracy=1)
    # Use linear lensing potentials for every cosmology. HMCode's calibrated
    # nonlinear correction fails to integrate in the very low cold-density run.
    p.NonLinear = camb.model.NonLinear_none
    result = camb.get_results(p)
    return result.get_total_cls(lmax=2600, CMB_unit="muK")[:, 0], float(result.z_eq)


def run():
    from fractions import Fraction

    pilot = json.loads(PILOT.read_text())
    assert (
        hashlib.sha256(DATA.read_bytes()).hexdigest() == pilot["planck"]["TT"]["sha256"]
    )
    fraction = float(
        Fraction(
            json.loads((HERE.parents[3] / NATIVE).read_text())[
                "perfect_delivery_max_sharp_mass_fraction"
            ]
        )
    )
    ell, y, err = np.loadtxt(DATA, usecols=(0, 1, 2), unpack=True)
    ix = ell <= 1800
    ell, y, err = ell[ix], y[ix], err[ix]
    points = [
        int(np.abs(ell - v["measured"][0]).argmin())
        for v in pilot["features"]["TT_max"][:3]
    ]
    observed = y[points]
    data = dict(
        ell=ell.tolist(),
        measured=y.tolist(),
        marginal_error=err.tolist(),
        measured_peaks=observed.tolist(),
        measured_peak_ell=ell[points].tolist(),
        third_over_first=float(observed[2] / observed[0]),
        third_over_second=float(observed[2] / observed[1]),
    )
    variants = {}
    for label, cold in (
        ("ordinary_cold_control", pilot["background"]["omch2"]),
        ("maximal_one_layer", fraction * pilot["background"]["ombh2"]),
        ("zero", 0.0),
    ):
        background = pilot["background"] | {"omch2": cold}
        tt, zeq = spectrum(background, pilot["control_ns"], pilot["As"])
        peaks = find_peaks(tt[:1601], prominence=40)[0]
        peaks = peaks[peaks > 100][:3]
        assert len(peaks) == 3
        x = np.arange(len(tt))
        dilation = ell[points[0]] / peaks[0]
        gain = observed[0] / tt[peaks[0]]
        raw = np.interp(ell, x, tt)
        aligned = gain * np.interp(ell / dilation, x, tt)
        keep = (ell >= 350) & (ell <= 1500)
        variants[label] = dict(
            omega_c=cold,
            equality_redshift=zeq,
            peak_ell=peaks.tolist(),
            peaks=tt[peaks].tolist(),
            third_over_first=tt[peaks[2]] / tt[peaks[0]],
            third_over_second=tt[peaks[2]] / tt[peaks[1]],
            first_peak_angular_dilation=dilation,
            first_peak_vertical_gain=gain,
            at_observed_crests_aligned=aligned[points].tolist(),
            raw_at_effective_ell=raw.tolist(),
            aligned_at_effective_ell=aligned.tolist(),
            raw_fractional_rms=float(np.sqrt(np.mean((raw[keep] / y[keep] - 1) ** 2))),
            aligned_fractional_rms=float(
                np.sqrt(np.mean((aligned[keep] / y[keep] - 1) ** 2))
            ),
        )
    return dict(
        camb=camb.__version__,
        nonlinear_lensing=False,
        native_fraction=fraction,
        measured=data,
        variants=variants,
        planck_tt=pilot["planck"]["TT"],
        control_background=pilot["background"],
        spectral_index=pilot["control_ns"],
        initial_amplitude=pilot["As"],
    )


def verify(old, new):
    if isinstance(old, dict):
        assert old.keys() == new.keys()
        for key in old:
            verify(old[key], new[key])
    elif isinstance(old, (float, list)):
        np.testing.assert_allclose(old, new, rtol=1e-6, atol=1e-8)
    else:
        assert old == new


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    result = run()
    dest = HERE / "results.json"
    if p.parse_args().write:
        dest.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        verify(json.loads(dest.read_text()), result)
