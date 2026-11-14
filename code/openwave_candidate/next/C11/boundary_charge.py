"""Complete all-number boundary/cut instrument on the selected C10 source.

Known total occupation N and Pauli occupancy bound the regional charge by
max(0,N-(d-r)) <= Q_R <= min(N,r). The interval width is at most floor(d/2).
A modulus exceeding that width suffices. Reading boundary B=3-N is part of
the instrument; it must not be forgotten when comparing continued states.
"""
from pathlib import Path
import argparse
import hashlib
import json
import operator
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
sys.path.insert(0, str(BASE / 'next/C10'))
from charged_reservoir import BoundarySource


def integer(value):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError('Boolean is not an integer apparatus label')
    try:
        return operator.index(value)
    except TypeError as error:
        raise ValueError('Exact integer apparatus label required') from error


def regional_interval(total, region_size, positions=6):
    total, region_size, positions = map(integer, (total, region_size, positions))
    if not 0 <= total <= positions or not 0 <= region_size <= positions:
        raise ValueError('Occupation or region outside the declared carrier')
    return max(0, total - positions + region_size), min(total, region_size)


def decode(residue, total, region_size, modulus=6, positions=6):
    residue, modulus = map(integer, (residue, modulus))
    low, high = regional_interval(total, region_size, positions)
    if modulus < 2 or modulus <= high - low:
        raise ValueError('Modulus must exceed the complete admitted charge interval width')
    charge = low + (residue - low) % modulus
    if charge > high:
        raise ValueError('Residue has no admitted regional occupation')
    return charge


def cut_record(source, state, region):
    value = 0
    for edge, (a, b) in enumerate(source.construction.seams):
        transit = source.construction.ports + edge
        value += (int(a in region) - int(transit in region)) * state[1][2 * edge]
        value += (int(b in region) - int(transit in region)) * state[1][2 * edge + 1]
    return value % source.modulus


def execute():
    model = BoundarySource(6)
    numbers = np.array([mask.bit_count() for mask, _ in model.states])
    comparisons = 0
    all_records = []
    max_width = 0
    for bits in range(64):
        region = {z for z in range(6) if bits >> z & 1}
        gauss = sum(z < 3 for z in region)
        charge = np.array([sum((mask >> z) & 1 for z in region) for mask, _ in model.states])
        cuts = np.array([cut_record(model, state, region) for state in model.states])
        records = []
        for total in range(7):
            low, high = regional_interval(total, len(region))
            max_width = max(max_width, high - low)
            boundary = 3 - total
            for cut in range(6):
                indices = np.flatnonzero((numbers == total) & (cuts == cut))
                residue = (gauss - cut - boundary * int(0 in region)) % 6
                try:
                    value = decode(residue, total, len(region))
                except ValueError:
                    value = None
                expected = np.array([], dtype=int) if value is None else np.flatnonzero((numbers == total) & (charge == value))
                if not np.array_equal(indices, expected):
                    raise ValueError('Complete boundary/cut and number/charge projectors differ')
                records.append({'boundary': boundary, 'total_number': total, 'cut_residue': cut,
                                'regional_charge': value, 'basis_indices': indices.tolist(),
                                'zero_outcome_retained': len(indices) == 0})
                comparisons += len(indices)
        all_records.append({'region': sorted(region), 'outcomes': records})
    # The control uses N=0/N=2, the same even fermion parity.
    # A Q_R-only apparatus preserves cross-number coherence in the same local
    # charge branch. Reading B destroys it. Use an actual pair of source states.
    region = {0}
    i = next(j for j, (mask, _) in enumerate(model.states) if mask == 0)
    j = next(j for j, (mask, _) in enumerate(model.states) if mask == 6)
    local = np.array([[.5, .5j], [-.5j, .5]], dtype=complex)
    with_boundary = np.diag(np.diag(local))
    difference = float(np.linalg.svd(local - with_boundary, compute_uv=False).sum() / 2)
    # Actual full-region cut is zero for both vacuum and full occupation. A
    # cut-only decoder cannot reconstruct all N=0..6 without B or a bound.
    alias = [next(k for k, (mask, _) in enumerate(model.states) if mask == m) for m in (0, 63)]
    full_region_alias = [cut_record(model, model.states[k], set(range(6))) for k in alias]
    invalid = [(0, 3, 3, 3), (0, True, 1, 6), (0, 7, 1, 6), (0, 1, 7, 6), (1, 0, 0, 6)]
    rejections = []
    for args in invalid:
        try:
            decode(*args)
        except ValueError:
            rejections.append(True)
        else:
            rejections.append(False)
    sources = [Path(__file__), BASE / 'next/C10/charged_reservoir.py',
               BASE / 'next/C9/fermionic_shared_records.py',
               BASE / 'candidate/oph_openwave/construction.py',
               BASE / 'candidate/oph_openwave/fiber_preparation.py']
    return {
        'schema': 'oph.c11.boundary_charge.v1',
        'source': 'C10 BoundarySource(6), same complete384-state Gauss carrier; no Hamiltonian or reference change',
        'checks': {'every_regional_complete_projector_equal': comparisons == 384 * 64,
                   'maximum_interval_width_three': max_width == 3,
                   'all_2688_boundary_cut_outcomes_retained': sum(len(x['outcomes']) for x in all_records) == 2688,
                   'Q_only_continuation_transfer_rejected': abs(difference - .5) < 1e-14,
                   'cut_only_all_number_alias_retained': full_region_alias == [0, 0],
                   'unsafe_or_inconsistent_intervals_rejected': all(rejections)},
        'proof': 'For d Pauli orbitals and r regional orbitals, Q lies in [max(0,N-d+r),min(N,r)]. Its width is at most min(r,d-r,N,d-N) <= floor(d/2). Equal residues in a narrower-than-modulus interval imply equal integers. Summing dressed Gauss gives Q_R=G_R-cut_F-B*1(0 in R) mod6. Equality of each diagonal projector gives equality of the complete instrument branch, including any external reference.',
        'instrument': 'First read integer boundary B, then the commuting modular cut flux; retain both labels. Each branch is P_(B,cut) rho P_(B,cut). No normalization or removal of zero outcomes.',
        'complete_regions': all_records,
        'basis_comparisons': comparisons,
        'cross_number_control': {'region': [0], 'source_basis_indices': [i, j], 'occupation_masks': [0, 6],
                                 'both_regional_charges': [0, 0], 'Q_only_vs_boundary_then_Q_trace_distance': difference,
                                 'external_reference': 'Projector equality tensors with the identity on every external reference; the cross-number control explicitly rejects forgetting the added boundary measurement.'},
        'full_region_cut_alias': {'source_basis_indices': alias, 'numbers': [0, 6], 'cut_residues': full_region_alias},
        'scope': 'Exact finite Gauss and Pauli occupancy instrument. Integer B is a declared additional apparatus record. It is redundant only when total N is already sharp and known. No physical charge species, electron winding, fractional charge or all31 criterion is inferred.',
        'full_criterion_met': False, 'official_rows_changed': [],
        'source_sha256': {str(p.resolve().relative_to(BASE.resolve())): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true'); mode.add_argument('--check', action='store_true')
    args = parser.parse_args(); result = execute()
    if not all(value is True for value in result['checks'].values()):
        raise ValueError(result['checks'])
    path = HERE / 'BOUNDARY_CHARGE.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2) + '\n')
    elif json.loads(path.read_text()) != result:
        raise ValueError('Complete boundary charge receipt differs')
    print(json.dumps({'checks': result['checks'], 'basis_comparisons': result['basis_comparisons']}))
