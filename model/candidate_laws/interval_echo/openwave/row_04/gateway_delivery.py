"""Commission a passive two-gateway delivery circuit using one-load transfer records."""

import argparse
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations
import json

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .corridor_pair import certify_native_collar, neutral_column, seal_path


def link(a, b):
    return (a, b) if a < b else (b, a)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def station(board, cell):
    """The tidal station's port rule, applied before any station is installed."""
    local = set(range(12 * cell, 12 * (cell + 1)))
    gateways = {p for p in local if board.degree[p] == 6}
    source = min(
        local - gateways, key=lambda p: (-len(gateways.intersection(board.neighbours[p])), p)
    )
    relays = tuple(sorted(gateways.intersection(board.neighbours[source])))
    sinks = []
    if len(relays) == 2:
        a, c = relays
        sinks = sorted(
            local.intersection(board.neighbours[a], board.neighbours[c]) - {source, a, c}
        )
    return source, relays, sinks


def wiring(board, source, relays, sink):
    a, c = relays
    live = {link(p, q) for p in relays for q in (source, sink)}
    incident = {link(p, q) for p in (source, a, c) for q in board.neighbours[p]}
    require(live <= incident and len(live) == 4, "missing delivery path")
    return live, incident - live


def destinations(board, port, written):
    return {q for q in board.neighbours[port] if link(port, q) not in written}


def delivery_sink(board, source, relays, written, *, allow_direct=False):
    """Challenge all first contacts and both single surviving relays. Return the common receiver."""
    left, right = (destinations(board, p, written) - {source} for p in relays)
    if left != right or len(left) != 1:
        return None
    sink = next(iter(left))
    if sink // 12 != source // 12 or sink in relays:
        return None
    first = destinations(board, source, written)
    if allow_direct:
        first.discard(sink)
    if first != set(relays):
        return None
    return sink


def exact_contact_checks():
    """Bind the sharp moves, writes and one-load endpoint readings to native columns."""
    certify_native_collar()
    law = candidate.native_law
    echo = law.echo()
    swap, contrast = law.local_term_parts()
    pulse = expm(-1j * np.pi * (swap + candidate.KAPPA * contrast))
    blank = np.eye(law.DIM)
    np.testing.assert_allclose(
        echo @ blank[:, law.SRC.index((0, 0))], neutral_column((0, 0)), atol=1e-13, rtol=0
    )
    for word in ((0, 1), (1, 0)):
        initial = blank[:, law.SRC.index(word)]
        opposite = blank[:, law.SRC.index(word[::-1])]
        np.testing.assert_allclose(swap @ initial, (initial - opposite) / 2, atol=1e-13, rtol=0)
        np.testing.assert_allclose(contrast @ initial, initial, atol=1e-13, rtol=0)
        np.testing.assert_allclose(
            pulse @ initial, np.exp(-1j * np.pi * candidate.KAPPA) * opposite, atol=1e-13, rtol=0
        )
        np.testing.assert_allclose(
            pulse @ neutral_column(word), neutral_column(word), atol=1e-13, rtol=0
        )
    # Read an endpoint on a blank incident seam. In the one-load sector a
    # positive oriented contrast means that the reader's endpoint is occupied.
    labels = [None] + [d for _, d in law.LABELS]
    for endpoint, sign in ((0, 1), (1, -1)):
        flag = np.repeat([d == sign for d in labels], len(law.SRC))
        for word in ((0, 0), (0, 1), (1, 0)):
            ket = blank[:, law.SRC.index(word)]
            np.testing.assert_allclose(
                echo.T @ (flag * (echo @ ket)), word[endpoint] * ket, atol=1e-13, rtol=0
            )


def challenge_contacts(board, ports, written, pulse):
    """Challenge each addressed link on a fresh one-load copy and read its far endpoint."""
    law = candidate.native_law
    records = []
    for p in ports:
        for q in board.neighbours[p]:
            edge = link(p, q)
            # A separate blank incident link reads q against a vacant reference.
            witnesses = [r for r in board.neighbours[q] if r != p and link(q, r) not in written]
            require(witnesses, "no blank endpoint meter")
            word = (1, 0) if p < q else (0, 1)
            initial = (
                neutral_column(word) if edge in written else np.eye(law.DIM)[:, law.SRC.index(word)]
            )
            state = (pulse @ initial).reshape(-1, law.N_SRC)
            # The signed echo read above realizes this endpoint projector exactly.
            effect = [loads[int(p < q)] == 1 for loads in law.SRC]
            yes = float(np.linalg.norm(state[:, effect]) ** 2)
            click = round(yes)
            require(
                click in (0, 1)
                and abs(yes - click) < 1e-12
                and abs(np.linalg.norm(state[:, np.logical_not(effect)]) ** 2 - (1 - click))
                < 1e-12,
                "non-sharp calibration or missing outcome",
            )
            records.append((p, q, min(witnesses), click))
    return records


def contact_receipts(board, ports, written, records):
    """Bind each record to its dispatched contact and to a distinct blank meter."""
    expected = {(p, q) for p in ports for q in board.neighbours[p]}
    seen = {}
    for p, q, witness, click in records:
        require(
            all(type(x) is int for x in (p, q, witness, click))
            and (p, q) in expected
            and click in (0, 1),
            "invalid contact record",
        )
        require((p, q) not in seen, "duplicate contact record")
        meters = [r for r in board.neighbours[q] if r != p and link(q, r) not in written]
        require(meters and witness == min(meters), "wrong endpoint meter")
        seen[p, q] = click
    require(seen.keys() == expected, "incomplete contact panel")
    require(
        all(seen.get((q, p), click) == click for (p, q), click in seen.items()),
        "inconsistent reciprocal contact",
    )
    return seen


def switches_from_records(board, source, records):
    """Infer the relays, common receiver and unwanted contacts from successful blank trials."""
    seen = contact_receipts(board, (source, *board.neighbours[source]), set(), records)
    require(all(seen.values()), "blank contact did not transfer")
    arrivals = {}
    for p, q in seen:
        arrivals.setdefault(p, set()).add(q)
    require(source in arrivals and arrivals[source] <= arrivals.keys(), "missing neighbour trials")
    relays = sorted(q for q in arrivals[source] if any(r // 12 != q // 12 for r in arrivals[q]))
    require(len(relays) == 2, "calibration needs two gateway neighbours")
    a, c = relays
    sinks = (arrivals[a] & arrivals[c]) - {a, source, c}
    sinks = {d for d in sinks if d // 12 == source // 12}
    require(len(sinks) == 1, "calibration has no unique common receiver")
    sink = sinks.pop()
    allowed = {source: {a, c}, a: {source, sink}, c: {source, sink}}
    return (
        tuple(relays),
        sink,
        {link(p, q) for p, keep in allowed.items() for q in arrivals[p] - keep},
    )


def commission(source, relays, live, cut):
    """Clear both ends for each write using only the two original loads and the four live links."""
    initial = link(source, min(relays))
    state = initial
    transcript = []

    def move_until(accept):
        nonlocal state
        queue = deque([(state, ())])
        seen = {state}
        while queue:
            word, moves = queue.popleft()
            if accept(word):
                transcript.extend(("swap", e) for e in moves)
                state = word
                return
            for edge in sorted(live):
                if len(set(word).intersection(edge)) == 1:
                    after = tuple(sorted(set(word).symmetric_difference(edge)))
                    if after not in seen:
                        seen.add(after)
                        queue.append((after, (*moves, edge)))
        raise ArithmeticError("no vacant write route")

    for e in sorted(cut):
        move_until(lambda word: set(word).isdisjoint(e))
        transcript.append(("echo", e))
    move_until(lambda word: word == initial)
    return transcript


def replay_installation(board, initial, transcript):
    require(
        len(initial) in (1, 2)
        and len(set(initial)) == len(initial)
        and all(type(p) is int and 0 <= p < board.n for p in initial),
        "invalid sharp loads",
    )
    loads = set(initial)
    memories = set()
    for op, edge in transcript:
        a, b = edge
        edge = (a, b)
        require(
            type(a) is int
            and type(b) is int
            and 0 <= a < b < board.n
            and b in board.neighbours[a]
            and edge not in memories,
            "missing or previously written contact",
        )
        if op == "echo":
            require(loads.isdisjoint(edge), "nonvacant memory write")
            memories.add(edge)
        elif op == "swap":
            require(len(loads.intersection(edge)) == 1, "invalid sharp swap")
            loads.symmetric_difference_update(edge)
        else:
            raise ValueError("unknown local operation")
    return tuple(sorted(loads)), memories


def two_attempts(board, source, relays, sink, written, *, stop_on_receipt=False):
    """Choose one incident contact uniformly twice, retaining every refusal and return."""
    population = {source: Fraction(1)}
    for _ in range(2):
        after = Counter()
        for p, mass in population.items():
            if stop_on_receipt and p == sink:
                after[p] += mass
                continue
            for q in board.neighbours[p]:
                after[p if link(p, q) in written else q] += mass / len(board.neighbours[p])
        population = after
    bins = (
        population.get(source, 0),
        sum(population.get(p, 0) for p in relays),
        population.get(sink, 0),
    )
    return [str(x) for x in (*bins, 1 - sum(bins))]


def free_packet(board, source, relays, sink, cut):
    """Native one-load flight. All terminal ports are read in independent preparations."""
    u, v = np.array([e for e in board.seams if e not in cut]).T
    adjacency = sparse.csr_matrix(
        (np.ones(2 * len(u)), (np.r_[u, v], np.r_[v, u])), shape=(board.n,) * 2
    )
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    require(np.all(degree >= 1), "endpoint missing a blank reading seam")
    h = (sparse.diags(degree) - adjacency) / 2 + sparse.diags(
        candidate.KAPPA * (degree + 6 - board.degree)
    )
    start = np.zeros(board.n)
    start[source] = 1
    first = h @ start
    second = h @ first
    others = sorted(set(range(board.n)) - {source, sink, *relays})
    require(
        np.all(first[list(relays)] == -0.5) and first[sink] == 0 and np.all(first[others] == 0),
        "first-hop crosstalk",
    )
    require(second[sink] == 0.5 and np.all(second[others] == 0), "second-hop crosstalk")
    prob = abs(expm_multiply(-0.5j * h, start)) ** 2
    require(np.isfinite(prob).all() and abs(prob.sum() - 1) < 1e-11, "one-load flight lost weight")
    return prob.tolist()


def all_masks(board, source, relays, incident, *, allow_direct=False):
    """Evaluate every local switching pattern using only the one-load delivery requirement."""
    solutions = []
    for count in range(len(incident) + 1):
        for mask in combinations(sorted(incident), count):
            sink = delivery_sink(board, source, relays, set(mask), allow_direct=allow_direct)
            if sink is not None:
                solutions.append((sink, set(mask)))
    return solutions


def compare_nomination(board):
    # This is the consumer. The station, compiler and switch tests take no
    # pair addresses, strength, persistence threshold or nominated mask as inputs.
    from . import open_pair

    source, relays, sinks = station(board, 13)
    require(len(sinks) == 1, "ambiguous filed receiver")
    sink = sinks[0]
    law = candidate.native_law
    swap, contrast = law.local_term_parts()
    pulse = expm(-1j * np.pi * (swap + candidate.KAPPA * contrast))
    before = challenge_contacts(board, (source, *board.neighbours[source]), set(), pulse)
    measured_relays, measured_sink, cut = switches_from_records(board, source, before)
    require(
        (measured_relays, measured_sink) == (relays, sink), "record and carrier selectors disagree"
    )
    live, predicted = wiring(board, source, relays, sink)
    require(cut == predicted, "calibration selected different switches")
    solutions = all_masks(board, source, relays, live | cut)
    require(solutions == [(sink, cut)], "delivery admits a different mask")
    initial = link(source, relays[0])
    actions = commission(source, relays, live, cut)
    pair, written = replay_installation(board, initial, actions)
    after = challenge_contacts(board, (source, *relays), written, pulse)
    installed = contact_receipts(board, (source, *relays), written, after)
    require(
        all(click == (link(p, q) in live) for (p, q), click in installed.items()),
        "installed circuit misroutes",
    )
    a, c = relays
    # The filed vacancy schedule can also consume the synthesized switches.
    # Comparing its ordered record tests stricter equality than equal support.
    setup = seal_path((a, source, c), sorted(written))
    require(
        pair == setup["final_pair"] and set(setup["neutral_edges"]) == written,
        "different prepared state",
    )
    require(
        (a, source, c) == open_pair.PATH
        and sink == open_pair.BRIDGE
        and setup == open_pair.seal_from_sharp(),
        "different filed installation",
    )
    require(
        all(delivery_sink(board, source, relays, cut - {edge}) is None for edge in cut),
        "an omitted write preserves passive delivery",
    )
    # An ordinary receiver can accept the direct source-to-destination hop.
    # Test that option explicitly; the filed circuit requires gateway mediation.
    shortcut = link(source, sink)
    require(shortcut in written, "no direct delivery control in the filed cell")
    direct = written - {shortcut}
    with_bypass = all_masks(board, source, relays, live | cut, allow_direct=True)
    require(
        with_bypass == [(sink, direct), (sink, written)], "unexpected direct-receipt alternatives"
    )
    matched = {
        name: {
            rule: two_attempts(board, source, relays, sink, mask, stop_on_receipt=stop)
            for rule, stop in (("two_attempts", False), ("stop_on_receipt", True))
        }
        for name, mask in (("blank", set()), ("ten_writes", direct), ("eleven_writes", written))
    }
    return {
        "source": source,
        "relays": relays,
        "receiver": sink,
        "tested_masks": 2 ** len(live | cut),
        "accepted_masks": len(solutions),
        "neutral_writes": sorted(written),
        "actions": actions,
        "filed_setup": setup,
        "challenges_before_and_after": {"blank": before, "installed": after},
        "two_attempt_bins_source_relays_receiver_other": {
            "passive": matched["eleven_writes"]["two_attempts"],
            "blank": matched["blank"]["two_attempts"],
            "first_relay_unavailable": two_attempts(
                board, source, relays, sink, written | {link(source, relays[0])}
            ),
            "second_relay_unavailable": two_attempts(
                board, source, relays, sink, written | {link(source, relays[1])}
            ),
        },
        "stop_on_receipt_with_direct_link": matched["ten_writes"]["stop_on_receipt"],
        "matched_delivery_controls": matched,
        "addressed_blank_delivery": [
            replay_installation(
                board, (source,), [("swap", link(source, p)), ("swap", link(p, sink))]
            )[0]
            for p in relays
        ],
        "port_probabilities_at_half_exposure": free_packet(board, source, relays, sink, cut),
    }


def survey(board):
    counts = Counter()
    homes = []
    for cell in range(board.cells):
        source, relays, sinks = station(board, cell)
        if len(relays) != 2:
            status = "wrong relay count"
        elif not sinks:
            status = "no common receiver"
        elif len(sinks) > 1:
            status = "ambiguous receiver"
        else:
            status = "unique"
        counts[status] += 1
        layouts = []
        for sink in sinks:
            live, cut = wiring(board, source, relays, sink)
            require(delivery_sink(board, source, relays, cut) == sink, "wrong passive destination")
            initial = link(source, relays[0])
            actions = commission(source, relays, live, cut)
            require(
                replay_installation(board, initial, actions) == (initial, cut),
                "incomplete installation",
            )
            layouts.append({"sink": sink, "writes": len(cut), "swaps": len(actions) - len(cut)})
        homes.append({"source": source, "relays": relays, "status": status, "layouts": layouts})
    return {"cells": board.cells, "counts": dict(counts), "homes": homes}


def verified_installation(board, comparison):
    """Recover the state actually prepared before passing it to either pair flight."""
    b, station_relays, sinks = station(board, 13)
    require(type(comparison["source"]) is int and comparison["source"] == b, "wrong launch station")
    panels = comparison["challenges_before_and_after"]
    relays, d, cut = switches_from_records(board, b, panels["blank"])
    require(
        relays == station_relays
        and sinks == [d]
        and all(type(p) is int for p in comparison["relays"])
        and tuple(comparison["relays"]) == relays
        and type(comparison["receiver"]) is int
        and comparison["receiver"] == d,
        "wrong relay layout",
    )
    initial = link(b, relays[0])
    pair, written = replay_installation(board, initial, comparison["actions"])
    require(
        pair == initial
        and written == cut
        and json.dumps(sorted(cut)) == json.dumps(comparison["neutral_writes"]),
        "wrong commissioned state",
    )
    installed = contact_receipts(board, (b, *relays), written, panels["installed"])
    require(
        all(click == (link(p, q) not in cut) for (p, q), click in installed.items()),
        "wrong installed observation",
    )
    setup = seal_path((relays[0], b, relays[1]), sorted(cut))
    require(
        json.dumps(setup, sort_keys=True) == json.dumps(comparison["filed_setup"], sort_keys=True),
        "wrong filed preparation transcript",
    )
    return b, relays, d, cut, setup


def replay_filing(board, comparison):
    """Hand the generated installation to the registered acquisition and full independent replay."""
    from tools.strict_json import loads
    from . import open_pair, open_pair_check, open_pair_release
    from .corridor_release import monitored_source

    b, (a, c), _, _, setup = verified_installation(board, comparison)
    require((a, b, c) == open_pair.PATH and setup == open_pair.seal_from_sharp(), "source map")
    exchange, contrast = open_pair.source_parts()
    require(np.array_equal(contrast, 2 * np.eye(3)), "nonuniform acquisition phase")
    pulse = expm((-2j * np.pi / 3) * exchange)
    rho, exhausted, reach, first = monitored_source(
        pulse, open_pair_release.ROUNDS, open_pair_release.ATTEMPTS
    )
    floor, _ = open_pair_release.source_bound()
    home = [e for e in board.seams if e[0] // 12 == e[1] // 12 == b // 12]
    strengths, times = open_pair_check.STRENGTHS, open_pair_check.TIMES
    report = {
        "setup": setup,
        "home_seams_in_order": home,
        "guaranteed_joint_home": str(floor),
        "first_vacancy_in_one_attempt": first,
        "attempt_reach": reach,
        "exhaustion": exhausted,
        "appointments": [
            {"kappa": k, "exposure": t, "bins": open_pair_release.terminal(rho, exhausted, k, t)}
            for k in strengths
            for t in times
        ],
    }
    return open_pair_check.verify(
        report, strengths, times, loads(open_pair_check.ARCHIVE.read_text(encoding="utf-8"))
    )


def direct_link_pair_control(board, comparison):
    """Retain the acquired density and detector; restore the ordinary receiver shortcut."""
    from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
    from model.executable_model.carrier import Carrier
    from . import open_pair_check as filed

    b, (a, c), d, cut, _ = verified_installation(board, comparison)
    shortcut = link(b, d)
    initial = link(b, a)
    # Leave the shortcut blank from the start. The same sharp swaps restore
    # the source loads without erasing a memory before the free flight.
    actions = [(op, e) for op, e in comparison["actions"] if (op, tuple(e)) != ("echo", shortcut)]
    pair, written = replay_installation(board, initial, actions)
    require(pair == initial and written == cut - {shortcut}, "wrong direct-link preparation")
    live = Carrier("blank delivery contacts", board.n, [e for e in board.seams if e not in written])
    pairs = PairSector(live)
    exchange = pairs.exchange()
    weights = live.degree + 6 - board.degree
    contrast = weights[pairs.aa] + weights[pairs.bb] - 2 * pairs.adjacent
    launch = np.zeros((pairs.P, 2))
    source = [pairs.index(b, a), pairs.index(b, c)]
    launch[source, [0, 1]] = 1
    # Squared norm off the two supported words for (X-Z)/sqrt(2).
    escape = exchange @ (launch[:, 0] - launch[:, 1])
    escape[source] = 0
    require(escape @ escape / 2 == 0.25, "wrong shortcut escape")
    home = [pairs.index(*e) for e in board.seams if e[0] // 12 == e[1] // 12 == b // 12]
    quiet = np.ones(pairs.P, dtype=bool)
    quiet[home] = False
    output = []
    for k in filed.STRENGTHS:
        amplitudes = filed.flight_amplitudes(exchange, contrast, k, filed.TIMES, launch, home)
        state = np.sqrt(float(filed.RETRY)) * (launch @ filed.VECTOR)
        exhaustion = 1 - float(state @ state)
        h = exchange + sparse.diags(k * (contrast - contrast[source[0]]), dtype=np.float64)
        previous = 0
        for t, z in zip(filed.TIMES, amplitudes, strict=True):
            state = expm_multiply(-1j * (t - previous) * h, state)
            previous = t
            bins = float(filed.RETRY) * abs(z @ filed.VECTOR) ** 2
            expected = np.r_[bins, 1 - bins.sum()]
            direct = np.r_[
                abs(state[home]) ** 2, np.vdot(state[quiet], state[quiet]).real + exhaustion
            ]
            require(
                np.isfinite(direct).all() and min(direct) >= 0 and abs(direct.sum() - 1) < 2e-10,
                "invalid direct-link probabilities",
            )
            np.testing.assert_allclose(direct, expected, rtol=0, atol=2e-10, equal_nan=False)
            output.append({"kappa": k, "exposure": t, "bins": direct.tolist()})
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--filed-replay",
        action="store_true",
        help="repeat all four filed arms at all nine appointments",
    )
    args = parser.parse_args()
    exact_contact_checks()
    output = {
        "carriers": [survey(federation(i)) for i in range(3)],
        "comparison": compare_nomination(federation(0)),
    }
    if args.filed_replay:
        output["filed_replay"] = replay_filing(federation(0), output["comparison"])
        output["direct_link_pair_control"] = direct_link_pair_control(
            federation(0), output["comparison"]
        )
    print(json.dumps(output, allow_nan=False))
