"""Check direct M1 consumer coverage against the current canonical registry."""

import hashlib
from pathlib import Path

from .check import canonical, require, strict_load

ROOT = Path(__file__).resolve().parents[2]
KEYS = {"declared_golden_source_population_read_law_clock_and_count_measure",
        "M1_full_family_read_feedback_and_routing",
        "declared_metric_coarsening_voronoi_kernel_scalar_action",
        "declared_lossless_independent_interventions_and_dense_observation"}
EXPECTED = {
    "OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM": "named_limits_and_sparse_scalar_action",
    "OPH-COSMO-FLRW-RECORD-DENSITY": "same_weight_identities_and_profile_boundary",
    "OPH-BH-CROSSING-READ-AREA-LAW": "cut_separation_and_entropy_preserving_replacement",
    "OPH-GOLDEN-SOURCE-QUADRATURE": "analogous_limits_not_golden_arithmetic",
    "OPH-RG-INFORMATION-RESOURCE-BOUND": "general_information_bound_retained_dense_counts_not_transferred"}
OWN = {"OPH-M1-SPARSE-CAUSAL-COUNT-REPLACEMENT", "OPH-M1-CUT-OBSERVABLE-NONUNIVERSALITY",
       "OPH-M1-SPARSE-ACTION-AND-TICK-OBSTRUCTION", "OPH-M1-ENTROPY-PRESERVING-SPARSIFICATION"}


def consumer_audit(registry=None, inventory=None):
    if registry is None:
        registry = strict_load(ROOT/"claims/claim_registry.yaml")
    if inventory is None:
        inventory = strict_load(Path(__file__).with_name("consumers.json"))
    require(canonical(inventory) == canonical(EXPECTED), "consumer disposition mismatch")
    selected = [row for row in registry["claims"] if KEYS & set(row["assumptions"])]
    ids = [row["claim_id"] for row in selected]
    require(len(ids) == len(set(ids)) and set(ids) == set(inventory), "M1 consumer coverage changed")
    # These commitments identify the complete original claims; the disposition
    # is deliberately not a statement that an entire mixed-scope claim transfers.
    return {row["claim_id"]: {"disposition":inventory[row["claim_id"]],
                              "original_claim_sha256":hashlib.sha256(canonical(row)).hexdigest()}
            for row in selected}


def downstream_audit(registry=None, graph=None, expected=None):
    if registry is None: registry=strict_load(ROOT/'claims/claim_registry.yaml')
    if graph is None: graph=strict_load(ROOT/'claims/dependency_graph.json')
    if expected is None: expected=strict_load(Path(__file__).with_name('downstream.json'))
    reached=set(EXPECTED)
    while True:
        extra={row['to'] for row in graph['edges'] if row['from'] in reached}
        if extra<=reached: break
        reached.update(extra)
    reached-=OWN
    require(type(expected) is list and all(type(x) is str for x in expected) and
            len(expected)==len(set(expected)) and set(expected)==reached,
            'registered downstream coverage changed')
    claims={row['claim_id']:row for row in registry['claims']}
    require(reached<=set(claims),'unknown downstream claim')
    # Reachability can include contextual/boundary edges: retain every edge's
    # role and every claim's own assumptions, without claiming premise transfer.
    links=sorted([row for row in graph['edges'] if row['from'] in reached and row['to'] in reached],
                 key=lambda row:(row['from'],row['to'],row['role']))
    return {'claims':{name:{'assumptions':claims[name]['assumptions'],
                            'original_claim_sha256':hashlib.sha256(canonical(claims[name])).hexdigest(),
                            'disposition':EXPECTED.get(name,'existing_contract_retained_without_automatic_transfer')}
                      for name in sorted(reached)},'links':links}
