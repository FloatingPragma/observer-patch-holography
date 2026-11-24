"""Signed fixed-area diagnostic on the protected native annular range instrument."""
import argparse
import json
from math import pi

import numpy as np
from scipy.sparse.csgraph import shortest_path

from model.executable_model import federation
from model.candidate_laws.interval_echo import candidate
from ..collar_pairing.measure import compare
from . import measure

HERE = measure.HERE


def masked_load():
    # On separate copies the observer measures the energy of a blank contact
    # with the ordinary sharp load parked, then performs the protected route survey.
    law = candidate.native_law
    h = law.seam_term(candidate.KAPPA)
    u, v = [np.eye(h.shape[0])[law.SRC.index(x)] for x in [(1,0),(0,1)]]
    empty = np.eye(h.shape[0])[law.SRC.index((0,0))]
    for vec, energy in [(empty,0),(u+v,8),(u-v,9)]:
        np.testing.assert_allclose(h@vec,energy*vec,rtol=0,atol=2e-13)
    board = federation(1)
    sites = np.arange(board.n)
    allowed = sites[sites%12!=4]
    if np.any(np.array(board.external)%12 == 4):
        raise ValueError("a masked port is an external gateway")
    # An independent shortest-path survey on the uncompressed port graph.
    reduced = board.adjacency()[allowed][:,allowed]
    destinations = np.searchsorted(allowed,sites[sites%12==0])
    ranges = shortest_path(reduced, directed=False, unweighted=True, indices=destinations[0])
    lengths = ranges[destinations]
    return dict(contact_gap_mean=float(u@h@u-empty@h@empty),local_energy_levels=[8,9],
                masked_ports=board.cells, maximum_range=int(lengths.max()),
                home_annulus_stations=int(np.sum((lengths>3)&(lengths<=6))))


def residual(charge, width, radius):
    # Every function of the protected range tape is identical in both source arms.
    # Dividing the signed area residual by positive kappa leaves charge/(2*pi).
    return dict(charge=charge, simultaneous_radius=width,
                per_coupling_residual_over_radius_fourth=charge/(2*pi*radius**4),
                per_coupling_lower=(charge-width)/(2*pi*radius**4))


def run():
    global_bank = json.loads((HERE/'results.json').read_text())
    ancestor_bank = json.loads((HERE/'ancestor.json').read_text())
    local = masked_load()
    rows = []
    for src, adjacent in zip(global_bank['panels'],ancestor_bank['panels'],strict=True):
        geometry = src['geometry']
        primary = src['settings'][0]
        arm = primary['arms']['separated']['contrasts']['home']
        ix = next(i for i, s in enumerate(adjacent['geometry']['settings'])
                  if (s['scale'], s['partner']) == (1,'far'))
        region = adjacent['arms'][ix]['contrasts']['home']
        radius = primary['relative_radius']
        if geometry['cells'] == local['masked_ports']:
            assert geometry['calibration_range'] == local['maximum_range']
            assert geometry['annulus_counts'][0][0] == local['home_annulus_stations']
        rows.append(dict(cells=geometry['cells'], relative_routed_radius=radius,
                         global_tag=residual(arm['measured_shared_change'],arm['simultaneous_radius'],radius),
                         ancestor_tag=residual(region['shared'],region['simultaneous_radius'],radius)))
    return dict(masked_ordinary_contact=local, boards=rows)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    data=run()
    dest=HERE/'signed_area.json'
    if parser.parse_args().write:
        dest.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    else:
        compare(data,json.loads(dest.read_text()))
    print(data['boards'][-1])
