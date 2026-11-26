"""Material hand, conditional launch and stationary belt transfer."""

import numpy as np
from scipy.linalg import eigh
from .. import experiment as belt
from model.candidate_laws.interval_echo import candidate


def contact_certificate():
    """Check the sharp public position effect, its unwrite and source stationarity."""
    law = candidate.native_law
    echo = law.echo()
    for q in (1, 2):
        for word in ((q, 0), (0, q), (0, 0)):
            column = echo[:, law.SRC.index(word)]
            for wanted in (-q, 0, q):
                mask = np.r_[
                    np.zeros(9), np.repeat([int(d == wanted) for _, d in law.LABELS], 9)
                ]
                actual = echo.T @ (column * mask)
                expected = np.eye(law.DIM)[:, law.SRC.index(word)] * int(
                    word[0] - word[1] == wanted
                )
                if np.linalg.norm(actual - expected) > 1e-12:
                    raise ValueError(
                        "Signed sharp-position read or inverse echo failed"
                    )
    blocks, inserts = belt.calibration()
    g = belt.geometry()
    ports = g["perimeter"]
    if len(set(ports)) != 10 or set(g["source_contact_sites"]) != {ports[2], ports[8]}:
        raise ValueError("Unexpected source or transport geometry")
    return blocks, inserts, g


def histories(k, q, arm, blocks, inserts):
    h = belt.generator(k, q, arm, blocks, inserts)
    clocks = [belt.generator(k, q, arm, blocks, inserts, bond=b) for b in (2, 4)]
    return h, clocks


def hand(h, t):
    e, u = eigh(h)
    return abs(np.exp(-1j * np.asarray(t)[..., None] * e) @ (u[0] * u[1])) ** 2


def impulse(h, source, delays, parked=False):
    if parked:
        h = np.diag(np.diag(h))
    e, u = eigh(h)
    return abs((np.exp(-1j * np.outer(delays, e)) * u[source]) @ u.T) ** 2


def output(hclock, hguide, emit, t, d, parked=False):
    """Each draw retains a silent source or one of ten sharp terminal positions."""
    emission = hand(hclock, t[:, None] - d)
    flight = impulse(hguide, emit, d, parked)
    ports = emission @ flight / len(d)
    return np.column_stack((1 - emission.mean(axis=1), ports))


def formed_source():
    """Construct each collar from the native echo, with the stated odd-coin control."""
    from scipy.linalg import expm
    from model.experiments.material_inertia.dark_collar import label

    law = candidate.native_law
    gate = law.echo()
    pulse = expm(-1j * np.pi / 2 * law.seam_term(8))
    sharp = law.SRC.index((2, 0))
    seed = pulse[:, sharp]
    # An isolated private 8 n^2 phase gives a relative shift of 3 pi/2.
    seed *= np.exp(-1j * (3 * np.pi / 64) * 8 * np.array([u * u for u, v in law.SRC]))
    blank = np.eye(11)[0]
    even = (gate @ np.kron(blank, seed)).reshape(11, 9)
    material = law.SRC.index((1, 1))
    if abs(abs(np.vdot(label(2, 1), even[:, material])) - 1) > 1e-12:
        raise ValueError("Native even source formation failed")
    odd = even.copy()
    for i, (c, d) in enumerate(law.LABELS):
        if c:
            odd[1 + i] *= -1
    if abs(abs(np.vdot(label(2, 1, -1), odd[:, material])) - 1) > 1e-12:
        raise ValueError("Supplied relative coin flip failed")
    sham = (gate @ np.kron(blank, np.eye(9)[0])).reshape(11, 9)
    if abs(np.vdot(label(0, 1), sham[:, 0]) - 1) > 1e-12:
        raise ValueError("Native zero-direction formation failed")
    return {
        "even_echo_material_word": [1, 1],
        "odd_control": "supplied pi between stationary coin branches",
        "separated_source_energy_per_contact": [32, 32, 0],
        "material_staging": "separate shuttered routing exchanges park the load created at each belt site",
    }


def terminal_tree(port_chart):
    # The nine perimeter links form a spanning path. The pole memories are
    # outside that path. The least reconstructed load on ten sites is zero.
    edges = list(zip(port_chart[:-1], port_chart[1:], strict=True))
    for address in range(10):
        for tag in (1, 2):
            word = np.zeros(10, dtype=int)
            word[address] = tag
            differences = word[:-1] - word[1:]
            relative = np.r_[0, -np.cumsum(differences)]
            if not np.array_equal(relative - relative.min(), word):
                raise ValueError("Terminal public position decoder lost a word")
    return edges


def preparation_routes(chart):
    """Find fresh storage for the known loads created while forming the two writes."""
    import networkx as nx
    from model.executable_model import federation

    board = federation(0)
    pole = chart["parked_unit"]
    if board.degree[pole] != 5:
        raise ValueError("Missing private phase contact")
    arms = set(chart["perimeter"]) | {pole}
    used = set()
    graph = nx.Graph(board.seams)
    for point in chart["source_contact_sites"]:
        graph.remove_edge(pole, point)
    paths = []
    # After the first echo park its two singles, bring a separate sharp
    # double from storage to the pole, form the second write, evacuate its
    # outgoing belt single and retain the final sharp pole unit.
    first, second = chart["source_contact_sites"]
    for start in (first, pole, pole, second):
        # An endpoint transposition restores every intermediate load.
        options = nx.single_source_shortest_path(graph, start)
        path = min(
            (v for i, v in options.items() if i not in used | arms),
            key=lambda x: (len(x), x),
        )
        used.add(path[-1])
        paths.append(path)
    return paths
