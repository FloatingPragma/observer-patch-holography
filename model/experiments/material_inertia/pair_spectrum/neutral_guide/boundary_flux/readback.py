"""Read flux and discrimination intervals from terminal counts alone."""

import json
from pathlib import Path

import numpy as np
from scipy.special import betaincinv


def intervals(work, mobile):
    bank_count = 4 * len(work['source_order']) * sum(len(row['edges']) for row in work['banks'])
    shots = work['shots']
    delta = np.sqrt(2 * np.log(2 * bank_count / .01) / shots)
    records = {}
    for row in work['banks']:
        tau = row['dwell']
        total = np.zeros(len(work['source_order']))
        for edge in row['edges']:
            c = 2 * np.array(edge['plus_counts']) / shots - 1
            z = c[:, :, 0] + 1j * c[:, :, 1]
            shifts = -np.angle(z / z[0]) / tau
            total += (shifts[:, 0] - shifts[:, 1]) / 2
        err = len(row['edges']) * 2 * np.arcsin(np.sqrt(2) * delta) / tau
        records[tuple(row['seam']), row['radius'], tau] = total, err
    outer_upper = outer_charge = outer_written = 0.
    for (seam, radius, tau), (v, eps) in records.items():
        if radius == 1:
            continue
        inside, inner_error = records[seam, 1, tau]
        gain = (inside[1] - inside[2]) / 2
        if gain <= inner_error:
            raise ValueError('unresolved flux calibration')
        outer_charge = max(outer_charge, np.max(abs(v[1:3])) / gain)
        outer_written = max(outer_written, np.max(abs(v[3:5])) / gain)
        outer_upper = max(outer_upper, (np.max(abs(v[1:3])) + eps) / (gain - inner_error))
    n = mobile['shots']
    comparisons = mobile['banks']
    index = max(range(len(comparisons)), key=lambda j: abs(comparisons[j]['receiver_records']['write02']['endpoint_hits']
                                                           - comparisons[j]['receiver_records']['write11']['endpoint_hits']))
    selected = comparisons[index]
    # Two-sided simultaneous exact binomial coverage for both compared rates in every appointment.
    tail = .01 / (4 * len(comparisons))
    endpoint_intervals = {}
    for source in ('write02', 'write11'):
        k = selected['receiver_records'][source]['endpoint_hits']
        endpoint_intervals[source] = [float(betaincinv(k, n - k + 1, tail)) if k else 0.,
                                      float(betaincinv(k + 1, n - k, 1 - tail)) if k < n else 1.]
    return dict(real_quadrature_banks=bank_count, simultaneous_flux_coverage=.99,
                maximum_abs_charged_outer_ratio=float(outer_charge),
                maximum_abs_written_extreme_outer_ratio=float(outer_written),
                simultaneous_upper_on_abs_charged_outer_ratio=float(outer_upper),
                endpoint_selection={k: selected[k] for k in ('seam', 'sign', 'time', 'launch')},
                simultaneous_selected_endpoint_intervals=endpoint_intervals)


if __name__ == '__main__':
    data = json.loads(Path(__file__).with_name('results.json').read_text())
    print(json.dumps(intervals(data['boundary_work'], data['mobile_probe']), indent=2))
