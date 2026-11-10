"""Portable dyadic proposals for independently checked rational spectra.

An archive contains integer arrays and a JSON manifest, without pickle.
Its hashes bind proposals to exact generalized-energy matrices and labels;
they do not replace the consuming Gram, residual or projector certificates.
"""
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from zipfile import BadZipFile

import numpy as np


BITS = 44
SCALE = 1 << BITS
SCHEMA = 'oph.rational-eigenbasis-witness.v1'


class WitnessError(ValueError):
    """A stored proposal does not match its declared exact source block."""


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def block_key(boundary, L, twoSz, sign=0, component=None):
    key = {'boundary': list(boundary), 'L': L, 'twoSz': twoSz,
           'translation_sign': sign}
    if component is not None:
        key['component_coordinates'] = list(component)
    return _json(key)


def matrix_sha256(k, g):
    """Hash exact sparse entries, shapes and their K/G order."""
    matrices = [{'shape': list(m.shape),
                 'entries': [[i, j, str(value)]
                             for (i, j), value in sorted(m.todok().items())]}
                for m in (k, g)]
    return hashlib.sha256(_json(matrices).encode()).hexdigest()


def _witness_sha256(d, z):
    digest = hashlib.sha256()
    for array in (d, z):
        digest.update(np.asarray(array, dtype='<i8').tobytes(order='C'))
    return digest.hexdigest()


def _validate_arrays(d, z, dimension, key):
    if d.shape != (dimension,) or z.shape != (dimension, dimension):
        raise WitnessError(f'eigenbasis witness dimension mismatch: {key}')
    if d.dtype.kind != 'i' or z.dtype.kind != 'i' or d.dtype.itemsize != 8 or z.dtype.itemsize != 8:
        raise WitnessError(f'eigenbasis witness requires int64 arrays: {key}')
    # The public numerical API reconstructs these dyadics exactly in float64.
    if any(np.any((array <= -(1 << 53)) | (array >= (1 << 53))) for array in (d, z)):
        raise WitnessError(f'eigenbasis witness exceeds exact float64 conversion: {key}')
    if np.any(d[1:] < d[:-1]):
        raise WitnessError(f'eigenbasis witness energies are not ordered: {key}')


@lru_cache(maxsize=None)
def _archive(path, stamp, size):
    try:
        with np.load(path, allow_pickle=False) as archive:
            raw = archive['manifest']
            if raw.dtype != np.dtype('uint8') or raw.ndim != 1:
                raise WitnessError('eigenbasis witness manifest must be UTF-8 bytes')
            manifest = json.loads(raw.tobytes().decode('utf-8'))
            arrays = {name: archive[name] for name in archive.files if name != 'manifest'}
    except (OSError, ValueError, KeyError, UnicodeError, BadZipFile) as exc:
        raise WitnessError(f'unreadable eigenbasis witness archive {path}: {exc}') from exc
    if not isinstance(manifest, dict) or manifest.get('schema') != SCHEMA or manifest.get('bits') != BITS:
        raise WitnessError(f'eigenbasis witness schema or dyadic precision mismatch: {path}')
    if not isinstance(manifest.get('blocks'), dict):
        raise WitnessError(f'eigenbasis witness block manifest missing: {path}')
    for array in arrays.values():
        array.flags.writeable = False
    return manifest, arrays


def load_witness(path, domain, key, k, g):
    path = Path(path)
    try:
        status = path.stat()
    except FileNotFoundError as exc:
        raise WitnessError(f'missing eigenbasis witnesses {path}; run --generate-witnesses explicitly') from exc
    manifest, arrays = _archive(str(path), status.st_mtime_ns, status.st_size)
    if manifest.get('domain') != domain:
        raise WitnessError(f'eigenbasis witness source domain mismatch: {domain}')
    entry = manifest['blocks'].get(key)
    if entry is None:
        raise WitnessError(f'missing eigenbasis witness block: {key}; run --generate-witnesses explicitly')
    if (not isinstance(entry,dict) or k.shape != g.shape or k.rows != k.cols or entry.get('dimension') != k.rows
            or entry.get('matrix_sha256') != matrix_sha256(k, g)):
        raise WitnessError(f'eigenbasis witness exact K/G mismatch: {key}')
    prefix = hashlib.sha256(key.encode()).hexdigest()
    try:
        d, z = arrays[prefix + '_d'], arrays[prefix + '_z']
    except KeyError as exc:
        raise WitnessError(f'missing eigenbasis witness integer array: {key}') from exc
    _validate_arrays(d, z, k.rows, key)
    if entry.get('witness_sha256') != _witness_sha256(d, z):
        raise WitnessError(f'eigenbasis witness integer digest mismatch: {key}')
    d,z=d.astype(object),z.astype(object)
    d.flags.writeable=False;z.flags.writeable=False
    return d,z


def numerical_proposal(k, g):
    """Run the floating solver only during explicit witness generation."""
    import scipy.linalg as la
    values, vectors = la.eigh(np.asarray(k.tolist(), dtype=float),
                             np.asarray(g.tolist(), dtype=float), driver='gvd')
    return (np.rint(values * SCALE).astype(np.int64),
            np.rint(vectors * SCALE).astype(np.int64))


def write_witnesses(path, domain, blocks):
    """Persist explicitly supplied (key,K,G,d,Z) proposals atomically."""
    manifest = {'schema': SCHEMA, 'bits': BITS, 'domain': domain, 'blocks': {}}
    arrays = {}
    for key, k, g, d, z in blocks:
        if key in manifest['blocks']:
            raise WitnessError(f'duplicate eigenbasis witness block: {key}')
        _validate_arrays(d, z, k.rows, key)
        prefix = hashlib.sha256(key.encode()).hexdigest()
        arrays[prefix + '_d'] = np.asarray(d, dtype='<i8')
        arrays[prefix + '_z'] = np.asarray(z, dtype='<i8')
        manifest['blocks'][key] = {'dimension': k.rows,
                                 'matrix_sha256': matrix_sha256(k, g),
                                 'witness_sha256': _witness_sha256(d, z)}
    arrays['manifest'] = np.frombuffer(_json(manifest).encode(), dtype=np.uint8)
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    try:
        with temporary.open('wb') as output:
            np.savez_compressed(output, **arrays)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    _archive.cache_clear()
    return len(manifest['blocks'])
