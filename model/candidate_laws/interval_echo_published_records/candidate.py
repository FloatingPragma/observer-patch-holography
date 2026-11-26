"""A one-way published-record interface to a parallel-sheet comparison sector."""
import argparse

from model.candidate_laws.interval_echo import candidate as echo

NAME = "interval_echo_published_records"
TITLE = "Interval Echo with published records"
ROLE = "experimental"
SUMMARY = ("Published traveller and home records prepare a downstream parallel-sheet experiment. "
           "Its rods, material stress and endpoint entropy respond to the same archived coincidence.")
APPARATUS = "../interval_echo_partition_repair/published_strings/"
DOCUMENTS = {"readme": "README.md", "specification": APPARATUS+"README.md"}
LABELLING = echo.LABELLING
FEDERATION_LEVELS = echo.FEDERATION_LEVELS
basis = echo.basis
generator = echo.generator
response = echo.response
represent = echo.represent

PARAMETERS = {
    **echo.PARAMETERS,
    "services": {**echo.PARAMETERS["services"], "value": "../interval_echo/SPEC.md"},
    "archive_interface": {"value": "after publication", "provenance": "designed",
                          "quantity": "one-way access to already committed same-service bit pairs"},
    "sheet_transport": {"value": "2π", "provenance": "designed",
                       "quantity": "transverse holonomy in radians per ordinary or coincidence sheet"},
    "sheet_tension": {"value": 1, "provenance": "designed",
                      "quantity": "energy per longitudinal rod length of each homogeneous sheet family"},
    "comparison_geometry": {"value": "flat time and longitudinal axis with a homogeneous transverse surface",
                            "provenance": "designed",
                            "quantity": "product structure, rods and common observer frames of the added sector"},
    "reference_pair": {"value": "(5,1,1,1)/8", "provenance": "designed",
                       "quantity": "two-endpoint vacuum spectrum with a π/6 rotation of the even subspace"},
    "pulse_reading": {"value": "weighted sheet work", "provenance": "designed",
                      "quantity": "diamond-weighted work drives the endpoint rotation at its vacuum gap and susceptibility"},
    "material_equilibration": {"value": "three independent fair-bit exchanges per pair",
                               "provenance": "designed",
                               "quantity": "refresh each private ticket bit once before the material work reading"},
}

PREMISES = {
    "PR-23": {"status": "partial", "evidence": [APPARATUS+"README.md", APPARATUS+"memories.py"],
              "note": "The finite endpoint density gives separate entropy, vacuum modular energy and relative entropy on the sheet menu. A central overlap decomposition for freely evolving native fields is open."},
    "PR-24": {"status": "partial", "evidence": [APPARATUS+"GEOMETRY.md"],
              "note": "The two parallel sheet families have independently conserved label measures and Nambu–Goto action stress on the declared product geometry. Identifying the stress of native propagating loads and fields needs an additional measurement."},
    "PR-26": {"status": "partial", "evidence": [APPARATUS+"README.md", APPARATUS+"results.json",
                                              APPARATUS+"nested_stress/selection/README.md"],
              "note": "Three local ticket exchanges per memory select the matched source on every recorded history and from arbitrary incoming memories. The rod transport and sheet work law supply its geometry and stress at completed appointments."},
    "PR-27": {"status": "partial", "evidence": [APPARATUS+"GEOMETRY.md", APPARATUS+"results.json"],
              "note": "Direct fixed-volume areas and the ten independent timelike contractions give the small-ball curvature expansion on the product geometry. The active carrier does not select that smooth chart."},
    "PR-28": {"status": "partial", "evidence": [APPARATUS+"README.md", APPARATUS+"ARCHIVE.md"],
              "note": "A positive ordinary calibration transfers to held sheet coincidences on one finite menu and its specified refinement. A source-derived common geometry and coupling for native matter remain open."},
}
JOINT = {
    "einstein-together": {"status": "partial", "evidence": [APPARATUS+"README.md", APPARATUS+"GEOMETRY.md"],
                          "note": "The same archived density supplies two independently conserving sheet families, a rod survey and boundary memories. Their agreement relies on the added sheet, transport and material pulse laws. The active traveller or home excitation itself has no derived common proper-ball stress in this comparison."}
}


def run(argv=None):
    from model.candidate_laws.interval_echo_partition_repair.published_strings import experiment

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="replay the finite observer survey")
    parser.add_argument("--full", action="store_true", help="also replay the original native law")
    parser.add_argument("--selection", action="store_true", help="check all material bath outcomes and the missed-address control")
    args = parser.parse_args(argv)
    if args.full:
        echo.run(["--full"])
    if args.selection or args.full:
        from model.candidate_laws.interval_echo_partition_repair.published_strings.nested_stress.selection import experiment as bath
        bath.verify()
        print("Verified all 64 two-memory displacement records and the skipped-address control.")
    return experiment.main(["--check"])
