"""Exploratory geometry test of the equal-tension three-arm opening angle."""
import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_30_local_routes.experiment import (
    HERE, close_to_trimer, grow, native_routing, portable, test_route,
)
from model.executable_model import federation


def run():
    port_cosines = {1: 1/np.sqrt(5), 2: -1/np.sqrt(5), 3: -1.}
    out = {}
    for level in (0,1,2):
        car = federation(level)
        options, across = native_routing(car)
        geometric = {}
        for entry, paths in options.items():
            score = [abs(port_cosines[len(p)-1] + .5) for p in paths]
            geometric[entry] = [p for p,s in zip(paths,score) if s == min(score)]
        rows = {}
        for high in (False,True):
            lengths, trimer, longest, longest_pure = Counter(),Counter(),[],[]
            for start in sorted(geometric):
                # The angular comparison has already resolved paths of different
                # lengths. Only equal lengths can tie, so the previous walker applies.
                route,_ = grow(geometric,across,start,'farthest',high)
                lengths[len(route)] += 1
                if close_to_trimer(car,route):
                    trimer[len(route)] += 1
                    if len(route)>len(longest_pure):
                        longest_pure=route
                if len(route)>len(longest):
                    longest=route
            rows['high' if high else 'low'] = dict(length_histogram=lengths,
                all_trimer_histogram=trimer, longest=test_route(car,longest,4200+level),
                longest_whole_trimer=test_route(car,longest_pure,7161+level))
        out[car.cells]=rows
    return portable(out)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args=parser.parse_args()
    result=run()
    if args.output:
        with args.output.open('x') as handle:
            json.dump(result,handle,separators=(',',':'),allow_nan=False)
            handle.write('\n')
    elif result != json.loads((HERE/'junction.json').read_text()):
        raise SystemExit('Three-arm receipt differs from archive')
    for size,rows in result.items():
        for name,r in rows.items():
            print(size,name,max(map(int,r['length_histogram'])),
                  sum(r['all_trimer_histogram'].values()),
                  max(map(int,r['all_trimer_histogram']),default=0),
                  [i['inside'] for i in r['longest']['native']['interior_gaps']])
