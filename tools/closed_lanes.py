"""Historical lane issues, their superseding lanes, and the open issues that own open work.

Every register gate, consuming lane, owning issue and blocking issue names a
GitHub issue. The lane issues of the retired V3 and V4 programmes are all
closed; a register keeps such a number as labelled history and names the open
issue that owns the remaining work through ``live_owner``. Generated pages
render a historical number as "historical lane #n" beside its live owner.

``SUPERSEDED_LANE_SUCCESSORS`` records lanes that a later lane replaced while
the programme still ran; a claim gate names the successor lane, never the
superseded one. ``LIVE_OWNERS`` maps every historical lane that a register
still cites to the open issue that owns its remaining work. Add a lane here when
it closes, and update ``OPEN_ISSUES`` when the open physics issues change.
"""

from __future__ import annotations

# The open physics issues, by number, with their titles.
OPEN_ISSUES: dict[int, str] = {1033: 'Ongoing audit and repair of existing evidence',
 1078: 'Search for programme-wide OPH obstructions and stopping criteria',
 1134: 'One native process for observer clocks, particles and gravity',
 1135: 'Freeze the executable program class and decision contract',
 1136: 'Execute one continuing native observer with accounted resources',
 1137: 'Integrate finite theorems and curate the active proof chain',
 1138: 'Construct the continuing regime, faithful reference and local modular clocks',
 1139: 'Exhibit a localized native defect and its protecting invariant',
 1140: 'Derive common defect maintenance energy, kinetics and conserved stress',
 1141: 'Transport the same defects and derive particle identity and currents',
 1142: 'Obtain one interacting chiral 3+1 field description from native records',
 1143: 'Freeze the physical dictionary and derive the matched particle spectrum',
 1144: 'Construct responsive geometry and causal clocks on the same source',
 1145: 'Realize the joint modular, stress and entropy premises G1--G6',
 1146: 'Compose the realized Einstein response with matter and clock readouts',
 1147: 'Run joint independent physical tests and decide the completion milestone',
 1148: 'Audit obstructions and apply candidate, class and program stopping rules'}


# Lanes replaced by later lanes while the programme ran (historical record).
SUPERSEDED_LANE_SUCCESSORS: dict[int, tuple[int, ...]] = {
    # Spacetime adequacy, superseded in the 2026-09-10 V4 review.
    728: (740, 777),
    # Mechanics and thermodynamics, closed as bounded conditional milestones.
    731: (739,),
    732: (739,),
    # Electromagnetism, continued by the same-history Maxwell episode.
    733: (754,),
    # Standard Model structure, continued by common-world integration.
    734: (740,),
    # Simulation instruments, continued by the support-wiring readout.
    737: (776,),
    # Conserved current on a declared flat family.
    780: (735, 740, 745),
}

# Backward-compatible name used by the validators.
CLOSED_LANE_SUCCESSORS = SUPERSEDED_LANE_SUCCESSORS

# Every closed lane or issue that a register cites, mapped to the open issue
# that owns its remaining work. Physical attachments, sector completions and
# the constants belong to the completion program (#1134); discriminating tests and
# freezes to #1147; repair-law, scheduler and provenance derivations to #1135;
# source-coded histories and their readouts to #1136.
LIVE_OWNERS: dict[int, int] = {
    1025: 1147,  # superseded physical-test coordination
    1026: 1134,  # superseded common-model coordination
    1076: 1135,  # superseded repair-law coordination
    1077: 1136,  # original scalar-history benchmark retained, not marked executed
    314: 1134,  # physical matter and Spin action proof packet
    425: 1134,  # hadron backend
    506: 1147,  # independent same-scheme alpha/HVP falsification
    545: 1134,  # fine-structure endpoint reduction
    546: 1134,  # charged-lepton Yukawa landing
    554: 1134,  # discrete structural selectors
    565: 1134,  # twelve-port unit splitting and pairing
    566: 1134,  # twelve-port current from source tomography
    607: 1147,  # frozen-prediction register establishment
    618: 1134,  # gravity premise-elimination ladder
    625: 1134,  # integer load fibre and quadratic readback
    628: 1134,  # operational integer-counting readback
    643: 1147,  # A5 angular test
    649: 1134,  # N reserve attachment
    655: 1147,  # FZ-11 physical bridge
    664: 1147,  # carrier scale and positive lower bound
    666: 1147,  # seam-current discriminator
    667: 1147,  # FZ-12 comparison
    688: 1134,  # thermodynamic physical receipts
    692: 1134,  # finite causal observer net
    708: 1147,  # H0 preflight
    712: 1134,  # source-linked observer quantum-event tower
    726: 1134,  # observation ledger establishment
    727: 1134,  # premise register establishment
    728: 1134,  # spacetime adequacy
    729: 1134,  # gravitation
    730: 1134,  # quantum law and physical readout
    731: 1134,  # mechanics
    732: 1134,  # thermodynamics
    733: 1134,  # electromagnetism
    734: 1134,  # Standard Model structure
    735: 1134,  # matter-action coefficients
    736: 1134,  # masses and constants: P, N and the particle spectrum
    737: 1136,  # simulation instruments
    739: 1134,  # premise discharge and axiomatization queue
    740: 1134,  # one fixed model and cross-observable comparisons
    742: 1134,  # cosmology and astrophysics adequacy
    743: 1134,  # interactions, running and scattering
    744: 1134,  # QCD, hadrons, nuclear and atomic physics
    745: 1134,  # EWSB, weak interactions, flavour and neutrinos
    750: 1135,  # scheduler and source hypotheses from the repair process
    751: 1134,  # dark sector
    754: 1136,  # finite observer Maxwell episode on one history
    763: 1135,  # event precedence from source provenance
    776: 1136,  # support-wiring provenance readout
    777: 1136,  # the read law as a production
    778: 1136,  # population selection of committed record states
    779: 1135,  # one source law and its geometry response
    780: 1134,  # current readout on the sourced causal poset
    782: 1134,  # count-volume and ordering-fraction limits
    1087: 1078,  # audit sub-issues of the obstruction search
    1088: 1078,
    1089: 1078,
    1090: 1078,
}


def successors(lane: int) -> tuple[int, ...]:
    """Return the successor lanes of a superseded lane, or the lane itself."""
    return SUPERSEDED_LANE_SUCCESSORS.get(lane, (lane,))


def is_open(issue: int) -> bool:
    return issue in OPEN_ISSUES


def live_owner(issue: int) -> int:
    """Return the open issue that owns the work a cited issue scoped."""
    if issue in OPEN_ISSUES:
        return issue
    if issue in LIVE_OWNERS:
        return LIVE_OWNERS[issue]
    raise KeyError(f"issue #{issue} is neither open nor a registered historical lane")


def live_owners(issues) -> tuple[int, ...]:
    """Return the distinct live owners of several cited issues, in ascending order."""
    return tuple(sorted({live_owner(issue) for issue in issues}))


def is_known(issue: int) -> bool:
    return issue in OPEN_ISSUES or issue in LIVE_OWNERS


ISSUE_URL = "https://github.com/FloatingPragma/observer-patch-holography/issues"


def issue_link(issue: int) -> str:
    return f"[#{issue}]({ISSUE_URL}/{issue})"


def owner_text(issues) -> str:
    """Render live owners as links, with cited historical lanes as labelled history.

    A historical number is printed as plain text, never as a link, so a reader
    never mistakes a closed lane for a live owner.
    """
    issues = list(issues)
    if not issues:
        return ""
    rendered = ", ".join(issue_link(owner) for owner in live_owners(issues))
    historical = [issue for issue in issues if not is_open(issue)]
    if historical:
        noun = "lanes" if len(historical) > 1 else "lane"
        listed = ", ".join(f"#{issue}" for issue in historical)
        rendered += f" (historical {noun} {listed}, closed)"
    return rendered
