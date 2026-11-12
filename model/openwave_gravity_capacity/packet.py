"""Source-bound replay of the gravity and capacity consequences."""
import argparse
import hashlib
import json
from pathlib import Path

from model.finite_area_response import packet as source_packet
from model.matched_einstein_limit import packet as tensor_packet
from . import capacity, response

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
LEAN_SOURCES = (
    'Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean',
    'Lean/ObserverPatchHolography/DeSitterCapacityShock.lean',
    'Lean/ObserverPatchHolography/CapacityNonidentifiability.lean',
    'Lean/ObserverPatchHolography/GlobalCapacityAttachment.lean',
)


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def source_bindings():
    paths = set(source_packet.source_bindings()) | set(tensor_packet.source_bindings())
    paths.update(LEAN_SOURCES)
    paths.update(str(p.relative_to(ROOT)) for p in Path(__file__).parent.iterdir()
                 if p.suffix in ('.py', '.md', '.lean', '.tex'))
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def report():
    pins = source_bindings()
    result = {'schema': 'oph.finite_gravity_capacity_consequences.v1', 'all_pass': True,
              'source_bindings': pins, 'inherited_source_custody': source_packet.inherited_custody(),
              'metric_response': response.source_response(),
              'weak_field_curvature': response.linearized_curvature(),
              'regional_force_transfer': response.newton_error_bound(2, '1/1000', '1/2000',
                                                                    '1/4000', 1, 2),
              'capacity': capacity.certificate(),
              'scope': {
                  'new_source_local_metric_readout_response': True,
                  'finite_theorem_reused_for_Newton_error_budget': True,
                  'capacity_amplitude_cancels_from_normalized_history': True,
                  'native_source_supplies_Newton_limit': False,
                  'native_source_supplies_spatial_metric_phenomena': False,
                  'numeric_cosmological_Lambda_selected': False,
                  'observed_cosmology_compared': False,
                  'source_premises_or_Lambda_selection_proved_by_hashes': False,
                  'Lean_kernel_reexecuted_by_this_receipt': False}}
    if source_bindings() != pins:
        raise ValueError('source changed during replay')
    return result


def validate(retained, fresh):
    if canonical(retained) != canonical(fresh):
        raise ValueError('gravity/capacity receipt differs from complete replay')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = ap.parse_args()
    value = report()
    if args.write:
        args.write.write_text(canonical(value))
    else:
        validate(json.loads(args.check.read_text()), value)
    print('PASS: native finite metric response; conditional Newton and capacity consequences')


if __name__ == '__main__':
    main()
