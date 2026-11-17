"""Run the official CCF-LT-2 model and verify its complete result set."""
from __future__ import annotations

import argparse


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--headline', action='store_true',
                        help='Run only the finite component checks; omit full native source checks')
    output = parser.add_mutually_exclusive_group()
    output.add_argument('--output', help='Regenerate a current result file at this path')
    output.add_argument('--verify', help='Strictly compare an explicit result file')
    args = parser.parse_args(argv)
    from .ccf_lt2.verify import main as verify
    command = [] if args.headline else ['--full']
    for name in ('output', 'verify'):
        value = getattr(args, name)
        if value is not None:
            command += ['--'+name, value]
    return verify(command)


if __name__ == '__main__':
    main()
