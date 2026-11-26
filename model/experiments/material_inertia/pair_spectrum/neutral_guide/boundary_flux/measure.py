"""Frozen nested-work and independent mobile-probe comparison of written and ordinary material."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model.carrier import federation
from . import field, scatter
from .readback import intervals


def passive(local):
    echo = law.echo()
    entries = []
    for a, b in [(0, 2), (2, 0), (1, 1)]:
        written = echo[:, 3*a + b]
        ports = np.sum(abs(written.reshape(-1, 9))**2, axis=0)
        sharp = np.nonzero(ports > 1e-12)[0].tolist()
        entries.append(dict(before=[a, b], after_port_words=sharp,
                            written_probability=float(np.sum(abs(written[9:])**2)),
                            full_seam_generator_residual=float(np.linalg.norm(local @ written))))
    return dict(blank_to_written_operator_norm=float(np.linalg.norm(local[9:, :9], 2)),
                sharp_source_columns=entries)


def summarize(work, mobile):
    summary = []
    for bank in work['summaries']:
        own = next(r for r in work['summaries'] if r['seam'] == bank['seam']
                   and r['radius'] == 1 and r['dwell'] == bank['dwell'])
        first = next(r for r in work['summaries'] if r['radius'] == 1 and r['dwell'] == bank['dwell'])
        def gain(r):
            return (r['flux'][1] - r['flux'][2]) / 2
        summary.append(dict(seam=bank['seam'], radius=bank['radius'], dwell=bank['dwell'],
                            own_inner_gain=gain(own), first_inner_gain=gain(first),
                            own_normalized=[v / gain(own) for v in bank['flux']],
                            first_normalized=[v / gain(first) for v in bank['flux']]))
    tv, endpoint = [], []
    for bank in mobile['banks']:
        a, b = [bank['receiver_records'][s] for s in ('write02', 'write11')]
        value = sum(abs(np.asarray(a['owner_counts']) - b['owner_counts'])) / (2 * mobile['shots'])
        key = {k: bank[k] for k in ('seam', 'sign', 'time', 'launch')}
        tv.append((float(value), key))
        endpoint.append((abs(a['endpoint_hits'] - b['endpoint_hits']) / mobile['shots'], key))
    return dict(flux=summary, max_written_vs_ordinary_owner_tv=max(tv, key=lambda row: row[0]),
                max_written_vs_ordinary_endpoint_difference=max(endpoint, key=lambda row: row[0]))


def report():
    board = federation(0)
    seams = [board.external[0], board.external[-1]]
    h0, v = law.local_term_parts()
    h = h0 + 8 * v
    work = field.run(board, seams)
    mobile = scatter.measure(board, seams, h)
    return dict(carrier=board.name, mounts=[list(e) for e in seams], passive_memory=passive(h),
                boundary_work=work, mobile_probe=mobile, summary=summarize(work, mobile),
                count_only_readback=intervals(work, mobile))


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--output', type=Path, required=True)
    args = cli.parse_args()
    args.output.write_text(json.dumps(report(), indent=1, allow_nan=False) + '\n')
