"""Material surveys fed after publication by the original clock instruments."""
import argparse
import json

from model.candidate_laws.interval_echo import candidate as echo

NAME = "interval_echo_published_records"
TITLE = "Interval Echo with published records"
ROLE = "experimental"
SUMMARY = ("Published traveller and home records prepare parallel strings or a finite spherical dust body. "
           "On each material the same stock supplies proper-ball geometry, regional entropy and action stress.")
APPARATUS = "../interval_echo_partition_repair/published_strings/"
DUST = "spherical_dust/"
DOCUMENTS = {"readme": "README.md", "specification": "README.md"}
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
    "pulse_reading": {"value": "weighted material work", "provenance": "designed",
                      "quantity": "diamond-weighted sheet or dust work drives the endpoint rotation at its vacuum gap and susceptibility"},
    "material_equilibration": {"value": "three independent fair-bit exchanges per pair",
                               "provenance": "designed",
                               "quantity": "refresh each private ticket bit once before the material work reading"},
    "dust_flux": {"value": "G=1/4, enclosed mass divided by area radius squared",
                  "provenance": "designed", "quantity": "spherical pull per ordinary or coincidence mass unit"},
    "dust_clocks": {"value": "local Lorentz clocks and signals relative to the incoming zero-energy dust flow",
                    "provenance": "designed", "quantity": "clock and null propagation rule inside and outside the unit-radius stock"},
    "dust_pair_clocks": {"value": "permanent separation bands and two local lag tables with response 1/(1-z²)",
                         "provenance": "designed", "quantity": "stock-dependent change in the joint regional sign distribution at two locally recorded clock times"},
}

PREMISES = {
    "PR-23": {"status": "partial", "evidence": [APPARATUS+"README.md", DUST+"PROOF.md", DUST+"clock_pairs/REFINEMENT.md", DUST+"endpoint_counts.json"],
              "note": "The sheet and spherical dust memories give the stress coefficient from regional entropy on their enrolled states. The pulse couples to independently weighted action work."},
    "PR-24": {"status": "partial", "evidence": [APPARATUS+"GEOMETRY.md", DUST+"PROOF.md"],
              "note": "Worldsheet tension and pressureless worldline mass give separate conserved action stresses in their respective material modes. The native propagating-field identification is open."},
    "PR-26": {"status": "partial", "evidence": [APPARATUS+"README.md", APPARATUS+"results.json",
                                              APPARATUS+"nested_stress/selection/README.md", DUST+"PROOF.md", DUST+"clock_pairs/regions.json"],
              "note": "Two-time sign pairs on conserved spherical dust give fifth-order entropy–area balance throughout the regular core on the specified regional instrument. The single line of memories gives sixth-order balance only on its enrolled radii."},
    "PR-27": {"status": "partial", "evidence": [APPARATUS+"GEOMETRY.md", DUST+"survey.py", DUST+"results.json"],
              "note": "Six-axis finite shell surveys and five-stamp dust ranging independently recover fixed-volume areas through moving frames of the two declared geometries."},
    "PR-28": {"status": "partial", "evidence": [APPARATUS+"counted_tower/README.md", DUST+"PROOF.md", DUST+"EXTERIOR_READER.md"],
              "note": "An ordinary dust release sets the coupling across the proper-ball refinements and predicts the held exterior clock and release. Extending the source to freely formed native mass requires an independent identification."},
}
JOINT = {
    "einstein-together": {"status": "fulfilled", "evidence": [DUST+"PROOF.md", DUST+"clock_pairs/REFINEMENT.md", DUST+"clock_pairs/regions.json", DUST+"null_beam.json"],
                          "note": "On the homogeneous pressureless core, recorded clock pairs, proper balls and conserved material mass supply all three regional estimates at one ordinary coupling. The radial force, clock rules and material lag response are designed."}
}


def run(argv=None):
    from model.candidate_laws.interval_echo_partition_repair.published_strings import experiment

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="replay the finite observer survey")
    parser.add_argument("--full", action="store_true", help="also replay the original native law")
    parser.add_argument("--selection", action="store_true", help="check all material bath outcomes and the missed-address control")
    parser.add_argument("--dust", action="store_true", help="replay the spherical body, moving proper balls and exterior clocks")
    parser.add_argument("--counted", action="store_true", help="replay the parallel-string endpoint banks and tensor bounds")
    args = parser.parse_args(argv)
    if args.full:
        echo.run(["--full"])
    if args.selection or args.full:
        from model.candidate_laws.interval_echo_partition_repair.published_strings.nested_stress.selection import experiment as bath
        bath.verify()
        print("Verified all 64 two-memory displacement records and the skipped-address control.")
    if args.dust or args.check or args.full:
        from .spherical_dust import endpoint_counts, exterior_reader, null_beam
        from .spherical_dust import experiment as dust
        from .spherical_dust.clock_pairs import census, experiment as pair_clocks, modular_reader, regions, registered_reference

        data = dust.main([])
        endpoint_counts.replay(data)
        exterior_reader.check_archive(exterior_reader.run(), json.loads((dust.HERE/"exterior_reader.json").read_text()))
        null_beam.verify()
        pair_clocks.verify()
        registered_reference.verify()
        modular_reader.verify()
        census.verify()
        regions.verify()
        print("Verified the spherical dust surveys, paired-clock apertures, light focusing, and exterior clocks and release.")
    if args.counted or args.check or args.full:
        from model.candidate_laws.interval_echo_partition_repair.published_strings.counted_tower import experiment as counted
        counted.verify()
        print("Verified 105 parallel-string endpoint banks and the signed tensor bounds.")
    if args.dust and not (args.check or args.full):
        return 0
    return experiment.main(["--check"])
