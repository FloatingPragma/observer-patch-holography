"""Verify retained scientific inputs without executing their producer text."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[2]
OWNERS = (
    'colour_deconfinement/inputs', 'defect_integer_lift/inputs',
    'defect_particle_classes/inputs', 'electroweak_crossover/inputs',
    'finite_field_comparison_inputs', 'finite_repair_comparison_inputs',
    'gauge_higgs_comparison_inputs', 'lepton_register_round/inputs',
    'particle_construction_scope/comparison_inputs', 'quenched_hadron_binding/inputs',
)


def inspect_array(stream):
    import numpy.lib.format as fmt
    version = fmt.read_magic(stream)
    if version == (1, 0):
        shape, _, dtype = fmt.read_array_header_1_0(stream)
    elif version == (2, 0):
        shape, _, dtype = fmt.read_array_header_2_0(stream)
    else:
        raise ValueError('Unsupported array header version: ' + repr(version))
    return {'shape': list(shape), 'dtype': str(dtype),
            'contains_python_objects': bool(dtype.hasobject)}


def verify_manifest(path):
    path = Path(path)
    value = json.loads(path.read_text())
    if value.get('schema') != 'oph.scientific_comparison_inputs.v1':
        raise ValueError('Unexpected comparison-input schema')
    seen, size, containers, object_arrays, nonfinite_json = set(), 0, 0, [], []
    for row in value['files']:
        name = row['path']
        rel = Path(name)
        if rel.is_absolute() or '..' in rel.parts or name in seen:
            raise ValueError('Duplicate or escaping input path: ' + name)
        seen.add(name)
        target = path.parent / rel
        if target.is_symlink() or not target.resolve().is_relative_to(path.parent.resolve()):
            raise ValueError('Input path is a symlink or escapes its owner: ' + name)
        raw = target.read_bytes()
        if len(raw) != row['bytes'] or hashlib.sha256(raw).hexdigest() != row['sha256']:
            raise ValueError('Scientific input bytes differ: ' + name)
        size += len(raw)
        if target.suffix in ('.json', '.jsonl'):
            constants = []
            def retain_constant(text):
                constants.append(text)
                return text
            lines = raw.decode().splitlines() if target.suffix == '.jsonl' else [raw.decode()]
            for line in lines:
                if line.strip():
                    json.loads(line, parse_constant=retain_constant)
            if constants:
                if (row.get('nonfinite_constant_count') != len(constants)
                        or row.get('numerical_status') != 'raw_nonfinite_values'
                        or row.get('eligible_for_automatic_prediction_or_passing_numeric_receipt') is not False):
                    raise ValueError('Nonfinite input lacks its ineligible numerical classification: ' + name)
                nonfinite_json.append({'path': name, 'nonstandard_constants': len(constants),
                                       'eligible_for_prediction_or_passing_numeric_receipt': False})
        elif target.suffix in ('.npz', '.npy'):
            containers += 1
            if target.suffix == '.npz':
                with zipfile.ZipFile(target) as archive:
                    if archive.testzip() is not None:
                        raise ValueError('Corrupt array container: ' + name)
                    for member in archive.namelist():
                        with archive.open(member) as stream:
                            info = inspect_array(stream)
                        if info['contains_python_objects']:
                            object_arrays.append({'path': name, 'member': member})
            else:
                with target.open('rb') as stream:
                    info = inspect_array(stream)
                if info['contains_python_objects']:
                    object_arrays.append({'path': name})
    return {'files': len(seen), 'bytes': size, 'numpy_containers': containers,
            'object_arrays_not_deserialized': object_arrays,
            'nonstandard_json_constants_retained': nonfinite_json}


def report():
    values = {owner: verify_manifest(ROOT / 'model' / owner / 'input_manifest.json')
              for owner in OWNERS}
    return {'schema': 'oph.scientific_comparison_input_integrity.v1',
            'scope': 'Byte integrity and container syntax only; no producer execution or physical claim admission',
            'numerical_claims_admitted': False,
            'owners': values, 'files': sum(v['files'] for v in values.values()),
            'bytes': sum(v['bytes'] for v in values.values())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    value = report()
    encoded = json.dumps(value, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(f"Verified {value['files']} scientific inputs, {value['bytes']} bytes")


if __name__ == '__main__':
    main()
