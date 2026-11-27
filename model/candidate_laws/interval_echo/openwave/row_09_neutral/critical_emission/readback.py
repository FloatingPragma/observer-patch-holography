"""Count-based report of parent depletion, distal acceptance and its bright-only alternative."""

import json
from pathlib import Path

import numpy as np
from scipy.special import betaincinv

HERE = Path(__file__).resolve().parent


def summary(data):
    nshot = data['shots']
    rows = []
    for panel in data['panels']:
        if panel['source_link_fraction'] != .5:
            continue
        n = panel['dimers']
        parent, sham, packet = panel['sources']
        closed = next(r for r in data['panels'] if r['dimers'] == n and r['source_link_fraction'] == 0)
        seen = next(q for q in parent['local_records'] if q['time'] == 4*n)['joint_end_counts']
        if sum(seen) != nshot or sum(parent['joint_counts']) != nshot:
            raise ValueError('incomplete terminal records')
        k = parent['joint_counts'][0]
        lower = float(betaincinv(k, nshot-k+1, .01/24))
        if not lower > panel['bright_only_tag_upper']:
            raise ValueError('low-energy population unresolved')
        rows.append(dict(dimers=n, parent_excited_at_four_n=sum(seen[2:])/nshot,
                         depleted_and_absorbed=k, all_launches=nshot, lower_99=lower,
                         uninspected_absorption=parent['uninspected_counts'][0],
                         ground_sham=sham['joint_counts'][0], incident_packet=packet['joint_counts'][0],
                         closed_parent=closed['sources'][0]['joint_counts'][0],
                         closed_parent_incident_packet=closed['sources'][2]['joint_counts'][0]))
    return rows


def release_summary(data):
    levels = np.asarray(data['clock_levels'])
    starts = data['original_starts_per_contact_and_arm']
    rows = []
    for panel in data['panels']:
        for mount in panel['mounts']:
            contrast = variance_bound = 0.
            for term in mount['contacts']:
                u, v = term['path_positions']
                width = 1 if u//2 == v//2 else 33
                for key, sign in [('daughter',1),('host',-1)]:
                    counts = np.array(term[key])
                    if sum(counts) != starts or min(counts) < 0:
                        raise ValueError('incomplete static contact counts')
                    accepted = counts[2:].sum()
                    contrast += sign * (levels @ counts[2:])/accepted
                    variance_bound += width**2/accepted
            radius = np.sqrt(np.log(800)*variance_bound/2)
            np.testing.assert_allclose([contrast,radius],
                                       [mount['measured_increment'],mount['joint_99_radius']],atol=5e-12)
            rows.append(dict(dimers=panel['dimers'],mount=mount['mount'],
                             increment=contrast,joint_99_radius=radius,
                             low_weight_upper=mount['static_ground_spectral_bound']['low_weight_upper']))
    return rows


if __name__ == '__main__':
    print(json.dumps(dict(decay=summary(json.loads((HERE/'emitter.json').read_text())),
                          release=release_summary(json.loads((HERE/'release.json').read_text()))),indent=2))
