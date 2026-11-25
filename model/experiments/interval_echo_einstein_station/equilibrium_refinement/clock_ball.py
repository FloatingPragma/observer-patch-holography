"""Compare finite boundary counts, radial shell readings and joint probe flags on one timed ball."""
from collections import Counter, deque
from itertools import combinations
import json
import math

import numpy as np

from model.candidate_laws.interval_echo.clock import Registers
from model.executable_model import federation
from .measure import apparatus, generator
from ..proper_ball_ruler.measure import pulse_calibration

SHOTS = 2**19
CENTER = 0
RADIUS = 2
# This permutation preserves the entire monitored Hamiltonian. The receiver stays at CENTER.
REFLECTION = (9, 1, 2, 10, 11, 5, 6, 8, 7, 0, 3, 4)
SPREAD = ((0, 1), (2, 9))
SPLIT = ((0, 2), (1, 9))


def shortest_chains(ports, contacts, center):
    neighbours = {p: [] for p in ports}
    for a, b in contacts:
        neighbours[a].append(b)
        neighbours[b].append(a)
    chains = {center: ()}
    queue = deque([center])
    while queue:
        p = queue.popleft()
        for q in sorted(neighbours[p]):
            if q not in chains:
                chains[q] = (*chains[p], (p, q))
                queue.append(q)
    return chains


def read_word_menu(reg, border, ports, surveys, word):
    """Use only native binary comparisons. Each destructive transport is reversed."""
    reg.loads.fill(0)
    reg.loads[list(word)] = 2
    owner = {0}
    gateway = min(border)
    external = border[gateway]
    witness = {0, external // 12}

    def read(site):
        swaps = reg.route(site, gateway, owner)
        for a, b in swaps:
            reg.swap(a, b, owner)
        bit = int(reg.flag(gateway, external, witness))
        for a, b in reversed(swaps):
            reg.swap(a, b, owner)
        return bit

    flags = [read(p) for p in ports]
    # Transport the binary reading at the center. A local terminal comparison is
    # scored separately from the signal leg. This is a scheduled relay clock.
    for center, chains in surveys.items():
        for p, chain in chains.items():
            for a, b in chain:
                reg.swap(a, b, owner)
            if read(p) != flags[center]:
                raise ValueError("Routed message differs from the origin's registered bit")
            for a, b in reversed(chain):
                reg.swap(a, b, owner)
    if tuple(p for p, flag in zip(ports, flags, strict=True) if flag) != word:
        raise ValueError("Checkpoints failed to report the prepared word")
    if tuple(np.flatnonzero(reg.loads)) != word:
        raise ValueError("Reversible survey failed to restore the preparation")
    return flags


def two_flag_capacity(m, occupied_probability, both_occupied_probability):
    # A uniformly chosen single and an unordered pair of distinct shell receivers.
    # With at most two loads, log binom(m,k) is exactly quadratic in k.
    return (m * occupied_probability * math.log(m)
            + math.comb(m, 2) * both_occupied_probability
            * math.log(math.comb(m, 2) / m**2))


def measure():
    car = federation(0)
    ports, contacts, border, words = apparatus(car, [0])
    chains = shortest_chains(ports, contacts, CENTER)
    radii = [len(chains[p]) for p in ports]
    inside = [p for p in ports if radii[p] < RADIUS]
    rim = [p for p in ports if radii[p] == RADIUS]
    ball = set(inside + rim)
    h = generator(car, words, contacts)
    perm = [words.index(tuple(sorted(REFLECTION[p] for p in w))) for w in words]
    # The entries here are exact representable integers and half-integers.
    np.testing.assert_array_equal(2*h, np.rint(2*h))
    np.testing.assert_array_equal(h, h[np.ix_(perm, perm)])
    if {tuple(sorted(REFLECTION[p] for p in w)) for w in SPREAD} != set(SPLIT):
        raise ValueError("The symmetry misses one of the held preparations")
    reg = Registers(car, seed=31709)
    surveys = {p: shortest_chains(ports, contacts, p) for p in [CENTER, *rim]}
    flags = {w: np.asarray(read_word_menu(reg, border, ports, surveys, w)) for w in words}
    pulse = pulse_calibration()
    allowed = [w for w in words if set(w) <= ball]
    completions = Counter(tuple(p for p in inside if flags[w][p]) for w in allowed)
    m = len(rim)
    pairs = list(combinations(rim, 2))
    rng = np.random.default_rng(194055)
    # Independent baseline distribution of interior records, used in the first-law score.
    labels = sorted(completions)
    baseline = rng.multinomial(SHOTS, [completions[label]/len(allowed) for label in labels])
    entropy_score = dict(zip(labels, -np.log(baseline / SHOTS), strict=True))
    samples = {}
    for name, stock in (("split", SPLIT), ("spread", SPREAD)):
        aggregate = np.mean([flags[w] for w in stock], axis=0)
        word_tables = []
        for w in stock:
            bits = flags[w]
            b = [tuple(int(bits[p]) for p in pair) for pair in pairs]
            word_tables.append(Counter(sum(pair) for pair in b))
        pair_law = [sum(row[k] for row in word_tables)/(len(stock)*len(pairs)) for k in range(3)]
        probabilities = [float(aggregate[p]) for p in ports]
        local_counts = [int(rng.binomial(SHOTS, p)) for p in probabilities]
        pair_counts = rng.multinomial(SHOTS, pair_law)
        single_counts = int(rng.binomial(SHOTS, aggregate[rim].mean()))
        counts = rng.multinomial(SHOTS, [1/len(stock)] * len(stock))
        labels_arm = [tuple(p for p in inside if flags[w][p]) for w in stock]
        capacity = sum(math.log(completions[key]) for key in labels_arm) / len(stock)
        samples[name] = dict(words=stock, one_site_occupied_probability=probabilities,
            one_site_counts=local_counts, shell_pair_occupied_count_probability=pair_law,
            shell_pair_counts=pair_counts.tolist(), independent_single_shell_occupied=single_counts,
            volume_vacancy_units=float(sum(1-aggregate[p] for p in ball)),
            shell_vacancy_units=float(sum(1-aggregate[p] for p in rim)),
            signed_leg_count=counts.tolist(),
            entropy_score_mean=float(sum(n*entropy_score[key] for n,key in zip(counts, labels_arm, strict=True))/SHOTS),
            capacity_nats=capacity,
            capacity_from_pair_flags=two_flag_capacity(m, float(aggregate[rim].mean()), pair_law[2]),
            sampled_capacity_from_pair_flags=two_flag_capacity(m, single_counts/SHOTS, pair_counts[2]/SHOTS))
    s, t = samples["spread"], samples["split"]
    if not np.array_equal(s["one_site_occupied_probability"], t["one_site_occupied_probability"]):
        raise ValueError("Held words failed to preserve a one-station reading")
    # Each rim observer accepts a resident report from its own site or one
    # exchange away. This independently defined aperture can respond to pairs.
    panels = {p: [q for q in rim if len(surveys[p][q]) <= 1] for p in rim}
    def lit_panels(w):
        return sum(any(flags[w][p] for p in panel) for panel in panels.values())
    zero = tuple(sorted(inside[:2]))
    one = tuple(sorted([inside[0], rim[0]]))
    def capacity(w):
        return math.log(completions[tuple(p for p in inside if flags[w][p])])
    gain = (capacity(one) - capacity(zero)) / (lit_panels(one) - lit_panels(zero))
    panel_change = (sum(map(lit_panels, SPREAD)) - sum(map(lit_panels, SPLIT))) / 2
    shadows = dict(panels=panels, calibration_words=[zero, one], nats_per_lit_panel=gain,
        held_rim_pairs=[dict(word=w, lit_panels=lit_panels(w), capacity_nats=capacity(w))
                       for w in pairs], spread_minus_split_lit_panels=panel_change,
        predicted_capacity_contrast=gain*panel_change)
    return dict(center=CENTER, radius=RADIUS, unit="one calibrated outgoing exchange",
        single_exchange_hand=round(pulse["decoded_single_step"], 12),
        hand_hops_to_stations=radii, interior=inside, shell=rim, total_word_records=len(words),
        retained_ball_word_records=len(allowed), interior_reference_records=labels,
        independent_interior_reference_counts=baseline.tolist(), shots_per_meter=SHOTS,
        exact_energy_symmetry=list(REFLECTION), native_gates=dict(reg.gates),
        arms=samples, timed_rim_apertures=shadows, spread_minus_split=dict(
            capacity_nats=s["capacity_nats"]-t["capacity_nats"],
            flag_reconstruction_nats=s["sampled_capacity_from_pair_flags"]-t["sampled_capacity_from_pair_flags"],
            interior_entropy_first_variation_nats=s["entropy_score_mean"]-t["entropy_score_mean"],
            two_rim_flags_probability=s["shell_pair_occupied_count_probability"][2]-t["shell_pair_occupied_count_probability"][2]))


if __name__ == "__main__":
    print(json.dumps(measure(), indent=2, allow_nan=False))
