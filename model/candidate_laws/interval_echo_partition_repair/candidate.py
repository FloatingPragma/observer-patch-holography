"""Native contrast propagation and a local repair law for partition records."""
import argparse
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate as echo
from model.candidate_laws.contrast_memory_spherical_v1 import law as material
from .regional import PartitionRepair

NAME = "interval_echo_partition_repair"
TITLE = "Interval Echo with partition repair"
ROLE = "experimental"
SUMMARY = "Native contrast writes switch a shared aperture. Reversible repair selects the joint record distribution and fixes its entropy response on overlapping regions."
DOCUMENTS = {"readme": "README.md", "specification": "SPEC.md"}
LABELLING = echo.LABELLING
FEDERATION_LEVELS = echo.FEDERATION_LEVELS
# The free one-load interrogation can be addressed separately from the repair bath.
basis = echo.basis
generator = echo.generator
response = echo.response
represent = echo.represent
PARAMETERS = {
    "free_contrast": echo.PARAMETERS["kappa"],
    "interface_alphabet": {"value": 5, "provenance": "designed",
                           "quantity": "one closed pair and four open endpoint words"},
    "repair_commands": {"value": 4, "provenance": "designed",
                        "quantity": "two unbiased command bits for each addressed interface or volume update"},
    "visibility": {"value": "inside endpoints", "provenance": "designed",
                   "quantity": "local observers retain their endpoint letters and their own material registers"},
    "selection": {"value": "local mixing", "provenance": "designed",
                  "quantity": "positive visitation of each repair, native reservoir and contrast gate"},
}
PREMISES = {"PR-44": {"status": "fulfilled",
    "note": "Each addressed record repair and reservoir acts on a declared local factor by a complete Kraus family. Their complete channels and the native contrast unitary preserve an arbitrary untouched reference marginal.",
    "evidence": ["candidate.py", "SPEC.md", "test_candidate.py"]},
    "PR-26": {"status": "partial",
    "note": "Exact first-order stationarity for all overlapping endpoint regions at fixed recorded volume. A separate cubical experiment has independently surveyed boundaries. Area-dependent partner counts also enforce the fourth-order accuracy requirement on a supplied moving-ball survey. Native proper geometry and independent regional stress remain open.",
    "evidence": ["SPEC.md", "REGIONAL.md", "surveyed_cells/README.md", "proper_balls/README.md"]}}
JOINT = {}

HERE = Path(__file__).resolve().parent
WRITE = material.echo()
WRITTEN_CONTRAST = np.repeat([False, *(d != 0 for _, d in material.LABELS)], material.N_SRC)


def contrast_aperture(amplitudes):
    """Write the native contrast, toggle closed/first-open, and unwrite.

    Each column beyond the first two axes can carry an arbitrary reference.
    Linearity thus specifies the operation on entangled states as well.
    """
    amplitudes = np.asarray(amplitudes)
    if amplitudes.shape[:2] != (material.DIM, len(PartitionRepair.bond_words)):
        raise ValueError("expected native seam and legal aperture axes of size 99 and 5")
    written = (WRITE @ amplitudes.reshape(material.DIM, -1)).reshape(amplitudes.shape)
    toggled = written.copy()
    toggled[WRITTEN_CONTRAST, 0] = written[WRITTEN_CONTRAST, 1]
    toggled[WRITTEN_CONTRAST, 1] = written[WRITTEN_CONTRAST, 0]
    return (WRITE.T @ toggled.reshape(material.DIM, -1)).reshape(amplitudes.shape)


def prepare_reference_payload(local_density):
    """One local reservoir service: replace the native word by a tracial word.

    Input axes are local input ket, untouched ket, local input bra, untouched bra.
    The reference/discarded register is traced correctly even when entangled.
    """
    local_density = np.asarray(local_density)
    if (local_density.ndim != 4 or local_density.shape[:2] != local_density.shape[2:]
            or not local_density.shape[0]):
        raise ValueError("expected paired local and reference density axes")
    n = local_density.shape[0]
    rest = np.trace(local_density, axis1=0, axis2=2)
    return np.einsum("ij,ab->iajb", np.eye(n)/n, rest)


def repair_record(local_density, *, volume=False):
    """Measure the local legal word, apply its two-bit transition, retain the new word.

    Each diagonal block is the unnormalized state of an arbitrary untouched
    reference. Keeping the blocks preserves classical correlations with it.
    """
    density = np.asarray(local_density)
    n = 2 if volume else len(PartitionRepair.bond_words)
    if (density.ndim != 4 or density.shape[:2] != density.shape[2:]
            or density.shape[0] != n):
        raise ValueError("expected paired record and reference density axes")
    out = np.zeros(density.shape, dtype=np.result_type(density.dtype, float))
    advance = PartitionRepair.advance_volume if volume else PartitionRepair.advance_bond
    for before in range(n):
        for command in range(4):
            after = advance(before, command)
            out[after, :, after, :] += density[before, :, before, :]/4
    return out


def native_readout():
    rows = []
    for a, b in material.SRC:
        inputs = np.zeros((material.DIM, 5))
        inputs[material.SRC.index((a, b)), 0] = 1
        final = contrast_aperture(inputs)
        recovered = contrast_aperture(final)
        np.testing.assert_allclose(recovered, inputs, rtol=0, atol=2e-15)
        probability = float(np.sum(abs(final[:, 1:])**2))
        if abs(probability - int(a != b)) > 2e-15:
            raise ArithmeticError("contrast read does not agree with the prepared ordinary loads")
        rows.append({"loads": [a, b], "open": int(round(probability))})
    return {"native_seam_dimension": material.DIM, "legal_words": 5,
            "contrast_aperture": rows, "area_gain_nats_per_open_face": float(np.log(2))}


def run(argv=None):
    from .regional_check import check_saved, region_receipt

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify the native seams and regional receipt (default)")
    parser.add_argument("--full", action="store_true", help="also verify the cubical and smooth-ball comparisons")
    args = parser.parse_args(argv)
    report = native_readout()
    regional = region_receipt()
    check_saved(regional, json.loads((HERE / "regional_results.json").read_text()))
    report["verified_restrictions"] = regional["strict_restrictions_verified"]
    if args.full:
        from .proper_balls.measure import measure, verify_saved
        from .surveyed_cells.measure import replay as replay_survey

        cubical_receipt = json.loads((HERE / "surveyed_cells" / "results.json").read_text())
        cube = replay_survey(cubical_receipt)
        check_saved(cube, cubical_receipt)
        report["cubical_held_residual_99"] = cube["held_residual_99"]
        balls = measure()
        verify_saved(balls)
        report["proper_ball_surveys"] = len(balls["frames"]) * len(balls["inverse_scales"])
    print(json.dumps(report, indent=2))
    return 0
