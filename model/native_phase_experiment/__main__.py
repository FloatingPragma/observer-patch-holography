"""Execute, retain or replay the source one-port phase experiment."""
import argparse
from contextlib import ExitStack
import json
from pathlib import Path
import time

import numpy as np

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from .experiment import ROOT, report
from .independent import compute
from .verify import verify


def _persist(result, arrays, receipt_path, archive_path):
    created = []
    try:
        with ExitStack() as files:
            receipt = archive = None
            if receipt_path:
                receipt = files.enter_context(receipt_path.open('x', encoding='utf-8', newline='\n'))
                created.append(receipt_path)
            if archive_path:
                archive = files.enter_context(archive_path.open('xb'))
                created.append(archive_path)
            if archive:
                np.savez_compressed(archive, **arrays)
            if receipt:
                receipt.write(canonical(result))
    except BaseException:
        for path in created:
            path.unlink(missing_ok=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', type=Path)
    actions.add_argument('--check', type=Path)
    parser.add_argument('--archive', type=Path)
    args = parser.parse_args()
    if args.archive and not str(args.archive).endswith('.npz'):
        args.archive = Path(str(args.archive) + '.npz')
    if args.write and args.write.exists():
        parser.error('output exists; choose a distinct receipt path')
    if args.archive and args.archive.exists():
        parser.error('archive exists; choose a distinct archive path')
    start = time.perf_counter()
    actual, arrays = report(keep_arrays=True)
    independent = compute(ROOT, compact=False)
    validation = verify(actual, independent)
    if args.check:
        retained = strict_json.load(args.check)
        verify(retained, independent)
        compare_receipt(retained, actual)
    _persist(actual, arrays, args.write, args.archive)
    print(json.dumps({**validation, 'seconds': time.perf_counter() - start,
                      'trials': [{key: trial[key] for key in ('setting', 'law', 'score')}
                                 for trial in actual['trials']]}, sort_keys=True))


if __name__ == '__main__':
    main()
