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
        "comparison": {
            "law": "Written-field Laplacian extension",
            "label": "Strong candidate for the electron’s charge and field response",
            "finding": "An independent boundary-flux read tracks the mobile signed load, and a probe detects "
                       "the confined source 22 cells away. This makes the load a strong candidate for the "
                       "electron’s charge and field response. Its three-input counts "
                       "exclude spin-½ fermions with matched spin-independent transfers. A four-face orbital "
                       "doublet can transfer coherently with the charge. During an addressed material exchange "
                       "the two different orbital states stay in their original cells. In the complete finite "
                       "one-cell sectors tested, all 5,790 vacancy and 1,954 excess energies are simple. Two "
                       "vacancies with resolved field endpoints give either exchange sign under different "
                       "intermediate reads of the same source. A freely transported spin doublet and an "
                       "independently calibrated rest energy in keV are unresolved. A separate first-order mass "
                       "construction realizes the geometric family ratios on a guided seven-cell charge, "
                       "with strong dependence on the guide and a supplied register mass operator. Each rest level "
                       "has one state per charge sign. A spinless orbital realization gives identical charge and flux histories.",
            "source": "model/candidate_laws/interval_echo/openwave/electron_field/README.md",
        },
    },
    "mu": {
        "status": "open",
        "note": "A distinct muon rest branch and a weak decay channel require measurements on Interval Echo.",
        "evidence": [],
        "comparison": {
            "law": "Written-field extension and first-order mass comparison",
            "label": "Charged-mode conversion and conserved guided characters",
            "finding": "The 48-cell route resolves three propagating bands with minimum-energy ratios "
                       "1 : 4.56155 : 10.40388. On a fully written pair of cells, three energy-heralded "
                       "inputs cross the same bridge. Opening the bridge does work and produces lower charged "
                       "modes plus neutral excitation. All six separated downward channels miss the isolated "
                       "parent energy. A separate first-order register comparison yields the geometric ratios on a "
                       "seven-cell charge path. Its free law preserves every family character exactly. Five sequential distal "
                       "reads distinguish all three resting branches with 34.46% optimal accuracy at equal priors.",
            "source": "model/candidate_laws/interval_echo/openwave/electron_field/lepton_family/identity_audit/README.md",
        },
    },
    "tau": {
        "status": "open",
        "note": "A third charged-lepton mass and weak daughter tags require measurements on Interval Echo.",
        "evidence": [],
        "comparison": {
            "law": "Written-field extension and first-order mass comparison",
            "label": "Guided rest ratios and absolute energy control",
            "finding": "Exact spectral bounds exclude a Koide-balanced triplet among the 5,790 vacancy "
                       "and 1,954 excess levels on the tested cell backgrounds. All connected pinned "
                       "single-configuration-per-port guides of the three stored federations have energy "
                       "span below 1,273. A first-order generator with the cyclic read installed as a mass operator reaches a ratio "
                       "of 3477.47284 on an unpinned seven-cell path, falling to 1.32433 upon adding one pinned port. A common "
                       "charged-sector energy offset preserves all charge and flux histories while changing the formation ratios.",
            "source": "model/candidate_laws/interval_echo/openwave/electron_field/lepton_family/identity_audit/README.md",
        },
    },
    "photon": {
        "status": "open",
        "note": "A scalar material load carries a message to an independent receiver. Its orbital angular "
                "patterns and single-load clicks give concrete alternatives to intrinsic vector polarization. A charged "
                "dipole writes a neutral mark. Uniform material makes all free memory terms independent, so "
                "a calibrated distant receiver cannot see that mark. Preparing a neutral excess–vacancy pair on "
                "blank seams costs at least 80 native units above the unit-filled reference.",
        "evidence": ["model/experiments/interval_echo_message/README.md",
                     "model/candidate_laws/interval_echo/openwave/mediator_controls/README.md",
                     "model/candidate_laws/interval_echo/openwave/photon_identity/native_emission/README.md",
                     "model/candidate_laws/interval_echo/openwave/photon_identity/native_emission/NEUTRAL_PAIR.md"],
        "comparison": {
            "law": "Written-field Laplacian extension",
            "label": "Signed probes and remote field-energy detections",
            "finding": "Four source seams and six independent probe seams detect field transport on 80- and "
                       "320-cell spheres. Every tested nonzero frequency has one resolved response pattern. "
                       "At the lowest 80-cell frequency, counts bound a second measured component below 0.97% of the first "
                       "at 99% confidence on these probes. The smaller spectrum is complete. The guided packet has "
                       "quadratic long-wave dispersion. Two independent ground-silent energy probes detect the parked source "
                       "on the 320-cell sphere.",
            "source": "model/candidate_laws/interval_echo/openwave/photon_identity/README.md",
        },
        "additional_comparisons": [{
            "law": "Finite quantized Wilson–Maxwell field with controlled charge circuits",
            "label": "Identified model photon in local emission and reception records",
            "finding": "A local charged-circuit herald prepares over 99.99999% one-quantum occupation. "
                       "A disjoint receiver distinguishes opposite quarter-turn windings by interfering "
                       "two interaction times, with the global field count unread. The transverse field, "
                       "quantization, clock and charged path operations belong to this comparison.",
            "source": "model/candidate_laws/interval_echo/openwave/photon_identity/quantum_transfer/README.md",
        }],
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
    "kappa": {**contrast.PARAMETERS["kappa"],
              "source": "../../experiments/interval_echo_consensus/modular_standard.py",
              "note": "The historical specification supplies 8 in exchange units. Its design search "
                      "compared pair persistence across strengths. A later experimental law, developed "
                      "with those results known, uses the same signed record to compare contrast and "
                      "exchange parity. Its eight incident witnesses on a shielded reference cell give "
                      "magnitude-one energy 8 exactly. Calibrating the source on one-port stationarity "
                      "gives magnitude-two energy 32 and recovers the whole seam operator. A fair parity comparison approaches 4. The "
                      "comparator and the choice to deploy one standard across the screen are design inputs."},
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
                                                "../../experiments/interval_echo_consensus/phase.json",
                                                "../../experiments/interval_echo_consensus/ASYNCHRONOUS.md",
                                                "../../experiments/interval_echo_consensus/coherent_sweep.json"],
              "note": "Two freely ordered strict writes can end in distinguishable terminal seam records. "
                      "The strict contrast guard and stored priority give identical material endpoint channels on the "
                      "six-dimensional intercell two-unit input. The priority service preserves relative phase within "
                      "one accepted address. A protected sweep with stored local precedences preserves the "
                      "full addressed coherent instrument under every legal schedule, including untouched "
                      "entangled references. Interleaved free and clock transactions remain open."},
    "PR-43": {"status": "fulfilled", "evidence": ["SPEC.md", "../contrast_memory_spherical_v1/law.py"],
              "note": "The supplied free Hamiltonian generates a continuous group of conjugations on "
                      "the finite matrix algebra. Clock measurement and conditional services are separate "
                      "instruments. The flow parameter has no independent physical-time calibration."},
    "PR-44": {"status": "partial", "evidence": ["SPEC.md", "clock.py"],
              "note": "Port and seam registers and the addressed patch give finite subsystems. "
                      "The routing and cadence are supplied."},
    "PR-02": {"status": "partial", "evidence": ["SPEC.md", "clock.py",
                                                "../../experiments/interval_echo_consensus/coherent_sweep.json",
                                                "../../experiments/interval_echo_consensus/collar_reference.json"],
              "note": "Local pointer branches give normalized public outcome partitions. Addressed joint "
                      "collar parity separates coherent preparations with identical single-collar margins. "
                      "The collar comparison retains its bypass and rejected records. The quantum state "
                      "and measurement representation are supplied."},
    "PR-23": {"status": "partial", "evidence": ["../../experiments/interval_echo_einstein_station/PAIR.md",
                                                "../../experiments/interval_echo_einstein_station/cap.json"],
              "note": "One parent Gibbs operator predicts the modular response of compatible nested and overlapping regions. "
                      "A prepared finite cap supplies a faithful state and a readable modular charge. "
                      "Six weighted pair-energy variations agree with an independent meter. "
                      "Geometrically normalized modular flow on compatible nested caps remains open."},
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
        "evidence": ["prior_calibration.json", "../../experiments/interval_echo_einstein_station/PAIR.md",
                     "../../experiments/interval_echo_einstein_station/pair.json",
                     "../../experiments/interval_echo_null_balance/interval.json",
                     "../../experiments/interval_echo_null_balance/regional.json",
                     "../../experiments/interval_echo_einstein_station/cap.json",
                     "../../experiments/interval_echo_null_balance/directional.json",
                     "../../experiments/interval_echo_einstein_station/backreaction.json",
                     "../../experiments/interval_echo_einstein_station/volume.json"],
        "note": "The native free generator obeys additive patch energy continuity. On the earlier pair, "
                "a prepared cap reads weighted seam energy at a held port where the flat reference failed. "
                "An interior exchange changes this cap's entropy with its complete gateway state fixed. "
                "Time-reversal-odd one-load seam currents exhibit torque. Load probes are exactly blind to "
                "an antisymmetric collar's stored energy. A parent Hamiltonian predicts held modular changes on "
                "nested patches. Two interior modes with identical exterior single-load histories yield different "
                "entropy slopes. Physical area, a conserved directional stress "
                "and the common curvature response to matter remain open.",
    },
}


def run(argv=None):
    from .replay import main
    return main(argv)
