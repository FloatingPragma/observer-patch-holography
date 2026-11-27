"""Read the original descriptive sources bound by the message and local-interface deposits."""
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile

BUNDLE = Path(__file__).with_name('archived_record_sources.zip')
DESCRIPTIONS = frozenset({
    'model/candidate_laws/interval_echo/SPEC.md',
    'model/candidate_laws/interval_echo/EVALUATION.json',
    'model/candidate_laws/contrast_memory_spherical_v1/SPEC.md',
})


def verify_deposit(bindings, *, bundle=BUNDLE):
    """Require the original bytes for each named path and digest, including on shallow checkouts."""
    with ZipFile(bundle) as sources:
        for name, digest in bindings.items():
            key = f'{digest}/{name}'
            if key not in sources.namelist() or sha256(sources.read(key)).hexdigest() != digest:
                raise ValueError(f'archived source differs from the deposit: {name}')
