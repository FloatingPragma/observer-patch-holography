"""Repository-local provenance for finite observer-field calculations."""
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[2]


def public_receipt(result):
    """Bind computational inputs without depending on prose renderings."""
    kept = {}
    for relative, expected in result.get('source_sha256', {}).items():
        relative = relative.replace('\\', '/')
        path = (ROOT / relative).resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError('Evidence input escapes the repository')
        relative = path.relative_to(ROOT).as_posix()
        if path.suffix.lower() in {'.md', '.tex', '.pdf'}:
            continue
        if relative.startswith(('model/', 'code/', 'Lean/')):
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            if actual != expected:
                raise ValueError(f'Changed computational input: {relative}')
            kept[relative] = actual
    helper = Path(__file__).resolve()
    kept[helper.relative_to(ROOT).as_posix()] = hashlib.sha256(helper.read_bytes()).hexdigest()
    result['source_sha256'] = kept
    result['evidence_contract'] = {
        'domain': 'Explicit finite state, representation, action and observation contracts in the named producer.',
        'numerical_residuals_are_outward_certificates': False,
        'physical_particle_masses_source_selected': False,
        'physical_electromagnetic_normalization_established': False,
    }
    return result
