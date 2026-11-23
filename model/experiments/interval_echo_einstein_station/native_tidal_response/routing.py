"""Classify every sharp stationary seam and replay a shared native routing response."""

import argparse
import hashlib
import json

import numpy as np
import sympy as sp

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_null_balance.regional import verify
from .measure import HERE, SHOTS, TAGS, TIMES, build_station, evolve, local_calibration


SEED = 26112363


def exact_seam():
    law = candidate.native_law
    W = sp.Matrix(law._write_pattern())
    reversal, squared = map(sp.Matrix, law._label_ops())
    # The signed record reversal and squared contrast preserve the columns of W.
    record = sp.kronecker_product((sp.eye(law.N_LAB)-reversal)/2 + int(candidate.KAPPA)*squared,
                                  sp.eye(law.N_SRC))
    complement = sp.eye(W.rows) - W*W.T/2
    assert complement*record*W == sp.zeros(W.rows, W.cols)
    native = sp.diag(W.T*record*W/2, complement*record*complement)
    exchange, cost = law.local_term_parts()
    np.testing.assert_allclose(np.asarray(native, dtype=float), exchange+candidate.KAPPA*cost,
                               atol=2e-13, rtol=0)
    blank = sp.diag(sp.eye(law.N_SRC), sp.zeros(law.N_LAB*law.N_SRC))
    assert native*blank == blank*native
    return native, record


def stationary_census():
    """Full memory-ray intersections, with exact arithmetic over the rational field."""
    law = candidate.native_law
    native, written_record = exact_seam()
    energies = (0, 8, 9, 32, 33)  # spectrum of native record generator
    assert set(native.eigenvals()) == set(energies)
    fibres = {}
    identity = sp.eye(law.DIM)
    for word in law.SRC:
        index = law.SRC.index(word)
        fibre = identity[:, index::law.N_SRC]
        fibres[word] = {e: native*fibre - e*fibre for e in energies}
    census = []
    for orientation in ("left", "right"):
        for resident in range(3):
            spaces = [((), sp.eye(1+law.N_LAB))]
            for q in (0, *TAGS):
                word = (resident, q) if orientation == "left" else (q, resident)
                following = []
                for history, basis in spaces:
                    for e, constraints in fibres[word].items():
                        directions = (constraints*basis).nullspace()
                        if directions:
                            following.append(((*history, e), basis*sp.Matrix.hstack(*directions)))
                spaces = following
            census.append(dict(orientation=orientation, resident=resident,
                               spaces=[dict(energy=list(h), shifts=[h[1]-h[0], h[2]-h[0]],
                                            dimension=b.cols) for h, b in spaces]))

    # Double-tag flight at an empty endpoint can use only the blank private memory.
    i, j = [law.SRC.index(w) for w in ((2, 0), (0, 2))]
    hop = native[j::law.N_SRC, i::law.N_SRC]
    expected = sp.zeros(1+law.N_LAB)
    expected[0, 0] = -sp.Rational(1, 2)
    assert hop == expected
    # With one double tag on an empty material background, an arbitrary written
    # memory evolves by its record generator independently of the tag position.
    for word in ((0, 0), (2, 0), (0, 2)):
        word_index = law.SRC.index(word)
        columns = list(range(law.N_SRC + word_index, law.DIM, law.N_SRC))
        assert native[law.N_SRC:, columns] == written_record[:, word_index::law.N_SRC]
        assert native[:law.N_SRC, columns] == sp.zeros(law.N_SRC, law.N_LAB)
    # The blank link has exactly the same exchange for a unit tag.
    i, j = [law.SRC.index(w) for w in ((1, 0), (0, 1))]
    assert native[j::law.N_SRC, i] == hop[:, 0]
    shared = sorted({s["shifts"][0] for row in census for s in row["spaces"]
                     if s["shifts"][0] == s["shifts"][1] and s["shifts"][0] != 0})
    return dict(resident_census=census, shared_nonzero_shifts=shared, blank_flag_commutes=True,
                double_tag_hop=dict(blank_eigenvalue=-.5, written_eigenvalue=0, written_dimension=10,
                                    all_written_memories_factor_on_empty_material=True))


def route_counts():
    rng = np.random.default_rng(SEED)
    blanks, meter = local_calibration()
    odd = meter["coin_odd"]
    assert odd["source_energy"] == 32 and odd["native_potential"] == [0, 0]
    hosts = []
    for level in (0, 1, 2):
        board = federation(level)
        original, changed = (build_station(board, blanks, alternate) for alternate in (False, True))
        assert original[1] == changed[1]
        np.testing.assert_array_equal(original[2], changed[2])
        old_edge = min(original[3] - changed[3])
        new_edge = min(changed[3] - original[3])
        assert all(port < 12 for edge in (old_edge, new_edge) for port in edge)
        enabled = original[2]
        a, b = map(lambda p: int(np.flatnonzero(enabled == p)[0]), old_edge)
        prepared = np.zeros((len(enabled), 1))
        prepared[a] = 1
        records = []
        for route, (free, _, _, _, _) in (("original", original), ("changed", changed)):
            # The inserted dark memory has a scalar gap, whether an endpoint is loaded or empty.
            distribution = evolve(free, a, prepared, 0)[:, :, 0]
            for q in TAGS:
                for t, law in zip(TIMES, distribution, strict=True):
                    at_receiver = float(law[b])
                    at_exterior = float(law[enabled//12 != 0].sum())
                    records.append(dict(route=route, tag=q, time=t,
                                        measured_receiver=int(rng.binomial(SHOTS, at_receiver)),
                                        exact_receiver=at_receiver,
                                        measured_exterior=int(rng.binomial(SHOTS, at_exterior)),
                                        exact_exterior=at_exterior))
        hosts.append(dict(cells=board.cells, original_energy_mark=list(old_edge),
                          replacement_energy_mark=list(new_edge), records=records))
    n = sum(2*len(h["records"]) for h in hosts)
    radius = float(np.sqrt(np.log(2*n/.01)/(2*SHOTS)))
    for h in hosts:
        for r in h["records"]:
            for key in ("exterior", "receiver"):
                assert abs(r["measured_"+key]/SHOTS - r["exact_"+key]) <= radius
    return dict(seed=SEED, shots=SHOTS, primary_frequencies=n, simultaneous_confidence=.99,
                record_radius=radius, hosts=hosts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="save the exact census and routing records")
    args = parser.parse_args()
    result = dict(protocol_sha256=hashlib.sha256((HERE/"ROUTING_PROTOCOL.md").read_bytes()).hexdigest(),
                  law=candidate.NAME, stationary_classification=stationary_census(), routing=route_counts())
    path = HERE/"routing.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    else:
        verify(json.loads(path.read_text()), result)
    census = result["stationary_classification"]
    print("Shared nonzero stationary shifts", census["shared_nonzero_shifts"])
    for row in census["resident_census"][:3]:
        print(row["resident"], row["spaces"])
    for h in result["routing"]["hosts"]:
        print(h["cells"], [(r["route"], r["measured_receiver"]/SHOTS, r["measured_exterior"]/SHOTS)
                           for r in h["records"] if r["tag"] == 1 and r["time"] == 2])


if __name__ == "__main__":
    main()
