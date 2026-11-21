"""Interval Echo uses the native contrast interaction and an addressed patch-memory rule."""
from model.candidate_laws.contrast_memory_spherical_v1 import candidate as contrast

NAME = "interval_echo"
TITLE = "Interval Echo"
ROLE = "main"
SUMMARY = ("A finite patch carries the signed contrast repair, its hand and its interval memory. "
           "A completed handoff filters the fresh closing mark through a bounded first-return walk.")
DOCUMENTS = {"readme": "README.md", "specification": "SPEC.md"}
PARTICLE_IDENTIFICATIONS = {
    "e": {
        "status": "open",
        "note": "The freely released bare pair responds to the centre of the tested local SU(3) current. "
                "Its terminal records also distinguish all five alternative port-label permutations. "
                "A physical colour-singlet assignment, signed electromagnetic response, spatial spin half "
                "and an electron–positron identification require further evidence on the same excitation.",
        "evidence": [
            "model/candidate_laws/interval_echo/openwave/row_04_identity/README.md",
            "model/candidate_laws/interval_echo/openwave/row_04_conjugates/README.md",
        ],
    },
}
LABELLING = contrast.LABELLING
FEDERATION_LEVELS = contrast.FEDERATION_LEVELS
KAPPA = contrast.KAPPA
native_law = contrast.law

# Both candidates call the same implementation in the one-load sector.
basis = contrast.basis
represent = contrast.represent


def generator(car):
    return contrast.generator(car, KAPPA)


def response(car):
    return generator(car).toarray()

PARAMETERS = {
    "kappa": contrast.PARAMETERS["kappa"],
    "age_bits": {"value": 3, "provenance": "designed", "quantity": "moving interval register width",
                 "note": "the reference replay uses three bits, the baseline comparison four, and the clock-trip precision profile five"},
    "precision_age_bits": {"value": 5, "provenance": "designed", "quantity": "interval width on clock-trip entry",
                           "note": "30 work bits across the moving patch, including five-bit age, length and address"},
    "hand": {"value": 5, "provenance": "designed", "quantity": "load positions on the prepared clock hand"},
    "elapsed_bits": {"value": 9, "provenance": "designed", "quantity": "elapsed count capacity in bits"},
    "walk_budget": {"value": [12, 12], "provenance": "designed",
                    "quantity": "maximum signed steps and address trials per step"},
    "closing_filter": {"value": "first return after a compulsory step", "provenance": "calibrated",
                       "observation": "the duration-weighted reunion ratio on the supplied reciprocal relay",
                       "quantity": "closing-mark suppression rule"},
    "services": {"value": "SPEC.md", "provenance": "designed",
                 "quantity": "local gate durations, preparations, route, relative opportunities and read windows"},
    "strict_repair_window": {"value": {"absolute_written_contrast": 2, "requests": 2},
                             "provenance": "designed", "quantity": "optional two-cell priority service",
                             "note": "the stored priority is served first on blank collars as specified in CONSENSUS.md"},
    "port_current_service": {"value": {"gain": 1, "pulse": 0.41}, "provenance": "designed",
                             "quantity": "gain and pulse setting of the optional rotation-selected port instrument",
                             "note": "source and protected response space are specified in CURRENTS.md"},
}

PREMISES = {
    "PR-01": {"status": "partial", "evidence": ["CONSENSUS.md", "consensus.py", "consensus.json",
                                                "../../experiments/interval_echo_consensus/phase.json"],
              "note": "Two freely ordered strict writes can end in distinguishable terminal seam records. "
                      "The strict contrast guard and stored priority give identical material endpoint channels on the "
                      "six-dimensional intercell two-unit input. The priority service preserves relative phase within "
                      "one accepted address. General schedules and coupled clock transactions remain open."},
    "PR-43": {"status": "fulfilled", "evidence": ["SPEC.md", "../contrast_memory_spherical_v1/law.py"],
              "note": "The supplied free Hamiltonian generates a continuous group of conjugations on "
                      "the finite matrix algebra. Clock measurement and conditional services are separate "
                      "instruments. The flow parameter has no independent physical-time calibration."},
    "PR-44": {"status": "partial", "evidence": ["SPEC.md", "clock.py"],
              "note": "Port and seam registers and the addressed patch give finite subsystems. "
                      "The routing and cadence are supplied."},
    "PR-02": {"status": "partial", "evidence": ["SPEC.md", "clock.py"],
              "note": "Local pointer branches have normalized weights on the declared preparations. "
                      "The quantum state and measurement representation is supplied."},
    "PR-54": {"status": "partial", "evidence": ["CURRENTS.md", "currents.json", "prior_calibration.json"],
              "note": "A declared family of twelve closed current operators acts measurably on the native bound "
                      "pair. Its completeness and a source-selected physical connection, field strength and "
                      "local gauge action in a common spacetime remain open."},
}
JOINT = {
    **contrast.JOINT,
    "matter-together": {
        "status": "partial",
        "evidence": ["currents.json", "CURRENTS.md", "prior_calibration.json"],
        "note": "A compiled signed current pulse acts on the bound pair used for free flight and seam detection. "
                "The source and detector are calibrated on supplied windows. Timing that flight by a hand "
                "on the same object and establishing complete physical gauge currents remain open.",
    },
    "patch-measurement": {
        "status": "partial",
        "evidence": ["currents.json", "prior_calibration.json"],
        "note": "The earlier four-cell finite meter reads the thermal-reference modular generator of the "
                "current-pulsed pair using 4,423 settings with exact probabilities. Its 31-setting sampled "
                "meter did not resolve the during-flight contrast. The exposure is externally timed and "
                "a directional spacetime stress observable remains open.",
    },
    "einstein-together": {
        "status": "partial",
        "evidence": ["prior_calibration.json"],
        "note": "The earlier matched pair experiment diagnoses how regional modular energy changes during "
                "propagation under opposite current pulses. The reference is a supplied thermal state. "
                "Geometric cap normalization, directional matter stress and independent conservation remain open.",
    },
}


def run(argv=None):
    from .replay import main
    return main(argv)
