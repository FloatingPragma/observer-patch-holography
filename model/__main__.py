"""Run the executable model with a candidate law, list the candidates, evaluate them, or check the carrier.

    python -m model                         # the main candidate against its stored results
    python -m model --full                  # the main candidate, recomputing every result
    python -m model run --law NAME [ARGS]   # any candidate law; ARGS go to its replay
    python -m model list                    # the candidate laws and their roles
    python -m model evaluate [--write | --check | --law NAME]
    python -m model carrier                 # the checks of the fixed carrier
"""
from __future__ import annotations

import json
import sys

COMMANDS = ('run', 'list', 'evaluate', 'carrier')


def main(argv=None):
    from .executable_model import laws
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] not in COMMANDS + ('-h', '--help'):
        return laws.load(laws.main_candidate()).run(argv)
    command, rest = (argv[0], argv[1:]) if argv[0] in COMMANDS else ('help', [])
    if command == 'help':
        print(__doc__)
        return 0
    if command == 'run':
        name = laws.main_candidate()
        if rest[:1] == ['--law']:
            if len(rest) < 2:
                raise SystemExit('run --law needs a candidate name')
            name, rest = rest[1], rest[2:]
        return laws.load(name).run(rest)
    if command == 'list':
        for name in laws.names():
            cand = laws.load(name)
            print(f'{name:32s} {cand.ROLE:8s} {cand.TITLE}')
        return 0
    if command == 'evaluate':
        from .evaluation.evaluate import main as evaluate
        return evaluate(rest)
    from .executable_model import carrier
    result = carrier.run(levels=(0, 1, 2))
    problems = carrier.check(result)
    summary = {key: {k: v for k, v in entry.items() if not isinstance(v, (list, dict)) or k == 'degrees'}
               for key, entry in result.items()}
    print(json.dumps({'checks': summary, 'problems': problems}, indent=1))
    return 1 if problems else 0


if __name__ == '__main__':
    code = main()
    sys.exit(code if isinstance(code, int) else 0)
