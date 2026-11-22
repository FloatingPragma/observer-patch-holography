"""Reconstruct the incident priority of a sharp repair service from addressed hits."""

import argparse
from pathlib import Path

from model.executable_model import federation

from .cells import cell_inputs, read_cell
from .equipment import install
from .observed import Screen, need, separator
from .priority import authenticate, compare, observe, recover_document
from .replay import replay
from .source import PIN, bind_source, board_hash, canonical, read


def compare_candidates(args):
    screen, left = install(args.level, args.left)
    _, right = install(args.level, args.right)
    need(0 <= args.cell < screen.n // 12, "invalid cell")
    relations = []
    for source, path in ((left, args.left_hits), (right, args.right_hits)):
        document = read(path) if path else observe(args.level, screen, source)
        relations.append(authenticate(args.level, screen, source, document))
    answer = compare(screen, (left, right), relations)
    _, profile = bind_source()
    answer.update(source_commit=PIN, matched_source_profile=profile, level=args.level, cell=args.cell,
                  authenticated_wedges=len(screen.wedges))
    if not answer["equivalent"]:
        return answer

    # The cell theorem uses the observed class of the installed native priority.
    if "native" not in (args.left, args.right):
        _, native = install(args.level, "native")
        reference = authenticate(args.level, screen, native, observe(args.level, screen, native))
        if separator(screen, reference, relations[0]) is not None:
            return dict(answer, native_cell_class=False)
    local, boundary = cell_inputs(screen, relations[1].order)[args.cell]
    result = read_cell(screen, right, relations[1], board_hash(args.level), args.cell, local, boundary)
    return dict(answer, native_cell_class=True, **result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    get = sub.add_parser("observe", help="measure every incident pair")
    get.add_argument("--candidate", default="native", help="native, reverse, relabeled, or a JSON rank array")
    get.add_argument("--output", type=Path, help="write the JSON transcript to this file")
    recover = sub.add_parser("recover", help="reconstruct a relation from a recorded transcript")
    recover.add_argument("--hits", required=True, type=Path)
    pair = sub.add_parser("compare", help="compare two installed candidates and read the selected cell")
    pair.add_argument("--cell", type=int, default=0)
    for side in ("left", "right"):
        pair.add_argument("--" + side, required=True, help="native, reverse, relabeled, or a JSON rank array")
        pair.add_argument("--" + side + "-hits", type=Path, help="authenticate a recorded transcript")
    for command in (get, recover, pair):
        command.add_argument("--level", required=True, type=int, choices=range(3))
    sub.add_parser("replay", help="check all 420 cells and 3,360 preparations against the expected results")
    args = parser.parse_args()
    try:
        bind_source()
        if args.action == "replay":
            result = replay()
        elif args.action == "compare":
            result = compare_candidates(args)
        elif args.action == "recover":
            carrier = federation(args.level)
            relation = recover_document(args.level, Screen(carrier.n, carrier.seams), read(args.hits))
            result = dict(evidence="structural transcript only", predecessors=relation.higher, order=relation.order)
        else:
            screen, source = install(args.level, args.candidate)
            result = observe(args.level, screen, source)
            if args.output:
                need(args.output.suffix != ".gz", "use a plain JSON output file")
                with args.output.open("x", encoding="ascii", newline="\n") as stream:
                    stream.write(canonical(result) + "\n")
                result = dict(recorded_wedges=len(result["rows"]))
        print(canonical(result))
    except (OSError, ValueError, EOFError) as exc:
        parser.exit(2, f"{exc}\n")


if __name__ == "__main__":
    main()
