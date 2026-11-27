"""Release two kinds of load on observer zero's settled repair pattern."""

from collections import Counter
from functools import cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from matched_einstein_limit import observer_field
from model.candidate_laws.interval_echo import candidate
from model.experiments.material_inertia.dark_collar import label
from observer_zero_capacity import component
from strict_history.wide_replay import SCHEME, WideReplay


HERE = Path(__file__).resolve().parent
SEED = 20261126
SHOTS = 2**30
WIDTHS = (1 / 16, 1 / 32, 1 / 64)
TIMES = sorted(
    {0.0, 0.5, 1.0, 1.5, 2.0, *(1 + sign * h for h in WIDTHS for sign in (-1, 1))}
)


def settled_lab(model):
    links = cache(model.links)
    home = component(links, 0)
    ports = sorted(home)
    index = {p: i for i, p in enumerate(ports)}
    last = max(t for p in ports for _, t in links(p))
    if last >= 96:
        raise ValueError(
            "A material flight may cross an active write at the calibration time"
        )
    # Each message uses the settled material relay, with a separately routed reply.
    outbound = np.array(
        [[row[q] for q in ports] for row in (component(links, p) for p in ports)],
        dtype=int,
    )
    np.testing.assert_array_equal(outbound, outbound.T)
    roundtrips = outbound + outbound.T
    archived = json.loads((HERE / "level202_echo_reach.json").read_bytes())
    observed = {
        (12 * int(owner) + port): wait
        for (owner, port), wait in zip(
            archived["source_labels"],
            archived["surveys"][-1]["echo_wait_by_label"],
            strict=True,
        )
        if wait is not None
    }
    assert observed == {
        p: int(roundtrips[index[0], i]) for i, p in enumerate(ports) if p // 12
    }
    boundary = [(p, q) for p in ports for q, _ in links(p) if q not in index]
    inner, outer = min(boundary, key=lambda edge: (home[edge[0]], *edge))
    open_matrix = np.zeros_like(roundtrips)
    for p in ports:
        for q, written in links(p):
            if not written:
                open_matrix[index[p], index[q]] = 1
    degree = np.array([len(links(p)) for p in ports])
    for _, endpoint in boundary:
        owner, slot = divmod(endpoint, 12)
        assert 0 < model.history(owner).ends[slot] < 96
    return (
        ports,
        index[inner],
        outer,
        open_matrix,
        degree,
        roundtrips,
        last,
        sorted(boundary),
    )


def strict_releases(model, ports, home_ranges):
    def load(station):
        owner, slot = divmod(station, 12)
        h = model.history(owner)
        return 1 if h.ends[slot] else 2 * ((h.mask >> slot) & 1)

    launches = Counter()
    terminal_reads = []
    for i, p in enumerate(ports):
        initial = load(p)
        if initial == 1:
            continue
        for q, written in model.links(p):
            other = load(q)
            assert not (not written and abs(initial - other) == 2)
            assert not (not written and abs(1 - other) == 2)
        launches["add_one" if initial == 0 else "remove_one"] += 1
        distance = int(home_ranges[i])
        terminal_reads.append(
            dict(
                terminal=i,
                sham_load=initial,
                released_load=1,
                query_and_return_ticks=[[s, s + 2 * distance] for s in (112, 128)],
            )
        )
    return dict(
        launches=dict(launches),
        launches_with_an_enabled_write=0,
        sites_read=len(ports),
        unchanged_range_count=sum(launches.values()),
        terminal_reads=terminal_reads,
    )


def source_and_clock():
    law = candidate.native_law
    exchange, cost = law.local_term_parts()
    h = exchange + candidate.KAPPA * cost
    reversal, squared = law.generator_parts()
    echo = law.echo()
    loads = np.eye(law.N_SRC)
    blank = np.eye(law.DIM, law.N_SRC)
    np.testing.assert_allclose(
        h @ blank, blank @ law.seam_term(candidate.KAPPA), atol=1e-12
    )
    changed = (
        blank @ law.seam_term(candidate.KAPPA) @ blank.T
        + reversal
        + candidate.KAPPA * squared
    )
    # The actual signed strict write has the symmetric coin at both input orientations.
    writes = [law.write()[:, law.SRC.index(word)] for word in ((0, 2), (2, 0))]
    balanced = loads[:, law.SRC.index((1, 1))]
    modes = {}
    for name, sign in (("even", 1), ("odd", -1)):
        memory = label(2, sign)
        lift = np.kron(memory[:, None], loads)
        pointer = echo @ lift
        fields = []
        for meter in (squared, reversal):
            weights = np.sum(pointer * (meter @ pointer), axis=0)
            np.testing.assert_allclose(meter @ pointer, pointer * weights, atol=2e-12)
            fields.append(weights)
        levels = fields[0] * candidate.KAPPA + fields[1]
        np.testing.assert_allclose(h @ lift, lift * levels, atol=2e-12)
        changed_levels = np.sum(lift * (changed @ lift), axis=0)
        np.testing.assert_allclose(changed @ lift, lift * changed_levels, atol=2e-12)
        signed_projection = np.kron(memory[1:], balanced)
        pre = law.SRC.index((0, 1))
        modes[name] = dict(
            enrollment_probability_by_written_sign=[
                round(float(abs(signed_projection @ w) ** 2), 12) for w in writes
            ],
            source_squared_contrast=round(float(fields[0][pre]), 12),
            source_odd_reversal=round(float(fields[1][pre]), 12),
            energies_by_mobile_load=[
                round(float(levels[law.SRC.index((q, 1))]), 12) for q in (0, 1, 2)
            ],
            energies_without_resident=[
                round(float(levels[law.SRC.index((q, 0))]), 12) for q in (0, 1, 2)
            ],
            changed_law_with_resident=[
                round(float(changed_levels[law.SRC.index((q, 1))]), 12)
                for q in (0, 1, 2)
            ],
        )
    a, b = [law.SRC.index(w) for w in ((1, 0), (0, 1))]
    energies, vectors = np.linalg.eigh(h[np.ix_([a, b], [a, b])])
    clock = []
    for t in TIMES:
        p = abs(np.sum(vectors[0] * vectors[1] * np.exp(-1j * energies * t))) ** 2
        clock.append(float(2 * np.arcsin(np.sqrt(p))))
    np.testing.assert_allclose(clock, TIMES, atol=2e-13)
    return modes, clock


def receiver_source_meter():
    """Read sign parity and transition energy through the reachable end of the collar."""
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    levels, vectors = np.linalg.eigh(h)
    holds = (0.0, np.pi / 128, np.pi / 33, np.pi / 32, np.pi)
    gates = [(vectors * np.exp(-1j * t * levels)) @ vectors.T for t in holds]
    load_basis = np.eye(law.N_SRC)
    refs = [law.SRC.index(w) for w in ((0, 1), (1, 0))]
    values, axes = np.linalg.eigh(h[np.ix_(refs, refs)])
    splitter = (axes * np.exp(-0.5j * np.pi * values)) @ axes.T
    modes = {"even": label(2, 1), "odd": label(2, -1)}

    def terminal(memory, gate):
        paths = []
        for inside in (0, 1):
            local = law.SRC.index((inside, 1))
            output = (gate @ np.kron(memory, load_basis[:, local])).reshape(
                -1, law.N_SRC
            )
            np.testing.assert_allclose(
                np.linalg.norm(output[:, local]), 1.0, atol=2e-12
            )
            paths.append(output[:, local])
        return splitter.conj().T @ (splitter[:, :1] * paths)

    admission = []
    for word in ((0, 2), (2, 0)):
        written = law.write()[:, law.SRC.index(word)].reshape(law.N_LAB, law.N_SRC)
        signed = np.r_[0.0, written[:, law.SRC.index((1, 1))]]
        outcomes = terminal(signed, gates[-1])
        probs = np.sum(abs(outcomes) ** 2, axis=1)
        fidelity = [
            float(abs(memory @ state) ** 2 / p)
            for memory, state, p in zip(modes.values(), outcomes, probs, strict=True)
        ]
        np.testing.assert_allclose(probs, 0.5, atol=2e-12)
        np.testing.assert_allclose(fidelity, 1.0, atol=2e-12)
        admission.append(
            dict(
                signed_input=word,
                probabilities=probs.tolist(),
                post_read_fidelities=fidelity,
            )
        )
    rng = np.random.default_rng(20261129)
    shots = 2**26
    observed = {}
    for name, memory in modes.items():
        returns = [float(np.linalg.norm(terminal(memory, g)[0]) ** 2) for g in gates]
        clicks = int(rng.binomial(shots, returns[1]))
        phase_energy = float(np.arccos(2 * clicks / shots - 1) / holds[1])
        observed[name] = dict(
            exact_reference_returns=returns,
            short_hold_return_count=clicks,
            estimated_energy=phase_energy,
        )
    assert max(levels) * holds[1] < np.pi
    return dict(
        holds=holds, enrollment=admission, short_hold_shots=shots, readings=observed
    )


def wall_and_volume(adj, ranges, source):
    law = candidate.native_law
    h = law.seam_term(candidate.KAPPA)
    vacant = float(h[law.SRC.index((0, 0)), law.SRC.index((0, 0))])
    positive = float(h[law.SRC.index((1, 0)), law.SRC.index((1, 0))])
    walls = []
    for radius in (2, 3, 4):
        inside = ranges[source] <= radius
        cuts = int(adj[np.ix_(inside, ~inside)].sum())
        assert cuts > 0
        walls.append(
            dict(
                radius=radius,
                responding_ports=int(inside.sum()),
                switched_blank_contacts=cuts,
                switch_off_work_each_source=-cuts * vacant,
                single_loaded_contact_control_work=-positive,
            )
        )
    triples = []
    for neighbours in combinations(np.flatnonzero(adj[source]), 3):
        lengths = ranges[source, list(neighbours)] ** 2
        gram = (
            lengths[:, None]
            + lengths[None, :]
            - ranges[np.ix_(neighbours, neighbours)] ** 2
        ) / 2
        eigen = np.linalg.eigvalsh(gram)
        triples.append(
            dict(
                terminals=list(map(int, neighbours)),
                twice_gram=(2 * gram).astype(int).tolist(),
                classification="incompatible"
                if eigen[0] < -1e-10
                else ("degenerate" if eigen[0] < 1e-10 else "positive"),
            )
        )
    return dict(
        chambers=walls,
        volume_checks=triples,
        proper_volume=None,
        matter_pressure_in_physical_units=None,
    )


def transverse_echoes(ranges, home, source):
    regions = []
    for radius in (2, 3, 4, 8, 16, None):
        keep = (
            np.arange(len(ranges))
            if radius is None
            else np.flatnonzero(ranges[source] <= radius)
        )
        distances = ranges[np.ix_(keep, keep)]
        squared = distances**2
        flat = -0.5 * (
            squared
            - squared.mean(axis=0)
            - squared.mean(axis=1)[:, None]
            + squared.mean()
        )
        spherical = np.cos(distances / 100.0)
        regions.append(
            dict(
                radius=radius,
                sites=len(keep),
                flat_min_gram_eigenvalue=float(np.linalg.eigvalsh(flat)[0]),
                sphere_100_min_gram_eigenvalue=float(np.linalg.eigvalsh(spherical)[0]),
            )
        )
    opposite = []
    for triple in combinations(np.flatnonzero(ranges[home] == 1), 3):
        if all(ranges[a, b] == 2 for a, b in combinations(triple, 2)):
            opposite.append(list(map(int, triple)))
    return dict(regions=regions, home_unit_radius_pairwise_diameter_two=opposite)


def flight(h, launches, ranges):
    energies, vectors = np.linalg.eigh(h)
    # Global phase drops out of all terminal counts.
    energies -= energies[0]
    starts = vectors.T[:, launches]
    means, histograms = [], {}
    for time in TIMES:
        coeff = np.exp(-1j * time * energies)[:, None] * starts
        psi = vectors @ coeff
        prob = abs(psi) ** 2
        np.testing.assert_allclose(prob.sum(axis=0), 1, atol=1e-12)
        first, companion = prob[:, 0], prob[:, 1:].mean(axis=1)
        means.append(float(first @ ranges @ companion))
        if time in (1 - 1 / 16, 1.0, 1 + 1 / 16):
            joint = np.outer(first, companion)
            hist = np.bincount(ranges.ravel(), weights=joint.ravel())
            histograms[time] = hist / hist.sum()
        if time == 1:
            vel = vectors @ ((-1j * energies)[:, None] * coeff)
            acc = vectors @ ((-(energies**2))[:, None] * coeff)
            dp = 2 * (psi.conj() * vel).real
            ddp = 2 * ((psi.conj() * acc).real + abs(vel) ** 2)
            d2 = (
                ddp[:, 0] @ ranges @ companion
                + 2 * dp[:, 0] @ ranges @ dp[:, 1:].mean(axis=1)
                + first @ ranges @ ddp[:, 1:].mean(axis=1)
            )
    second = [
        (
            means[TIMES.index(1 + h)]
            - 2 * means[TIMES.index(1.0)]
            + means[TIMES.index(1 - h)]
        )
        / h**2
        for h in WIDTHS
    ]
    return means, histograms, float(d2), second


def count_ranges(rng, shots, law):
    counts = rng.multinomial(shots, law)
    bins = np.arange(len(law))
    mean = float(counts @ bins / shots)
    se = float(np.sqrt((counts @ ((bins - mean) ** 2)) / (shots * (shots - 1.0))))
    return counts.tolist(), mean, se


def probe_banks(adj, degree, ranges, source, modes, spectators):
    opens = adj.sum(axis=1)
    launches = [source, *map(int, np.flatnonzero(adj[source]))]
    rng = np.random.default_rng(SEED)
    retained_rng = np.random.default_rng(20261130)
    repeat_rng = np.random.default_rng(20261128)
    repeat_shots = 2**38
    banks = []
    for tag in (1, 2):
        base = (
            np.diag(opens) / 2
            - adj / 2
            + np.diag(candidate.KAPPA * tag**2 * (opens - degree))
        )
        cache_by_shift = {}
        for arm in (
            "resident",
            "resident_removed",
            "changed_written_law",
            "settled_exterior",
        ):
            records = {}
            for name, mode in modes.items():
                energy_reading = {
                    "resident": "energies_by_mobile_load",
                    "resident_removed": "energies_without_resident",
                    "changed_written_law": "changed_law_with_resident",
                    "settled_exterior": "energies_by_mobile_load",
                }
                values = mode[energy_reading[arm]]
                shift = values[tag] - values[0]
                key = (shift, arm == "settled_exterior" and tag == 1)
                if key not in cache_by_shift:
                    h = base.copy()
                    h[source, source] += shift
                    if key[1]:
                        for port, parity in spectators:
                            levels = modes[parity]["energies_by_mobile_load"]
                            h[port, port] += levels[tag] - levels[0]
                    cache_by_shift[key] = flight(h, launches, ranges)
                means, histograms, d2, second = cache_by_shift[key]
                main_rng = retained_rng if arm == "settled_exterior" else rng
                counts, mean, se = count_ranges(main_rng, SHOTS, histograms[1.0])
                records[name] = dict(
                    local_energy_shift=shift,
                    mean_ranges=means,
                    time_one_second_derivative=d2,
                    centered_second_differences=second,
                    time_one_range_histogram=counts,
                    measured_mean=mean,
                    mean_se=se,
                )
                if arm in ("resident", "settled_exterior"):
                    observations = [
                        count_ranges(repeat_rng, repeat_shots, histograms[t])
                        for t in (1 - 1 / 16, 1.0, 1 + 1 / 16)
                    ]
                    weights = np.array([1.0, -2.0, 1.0]) * 16**2
                    records[name]["precision_repeat"] = dict(
                        shots_per_time=repeat_shots,
                        histograms=[o[0] for o in observations],
                        acceleration=float(weights @ [o[1] for o in observations]),
                        acceleration_se=float(
                            np.linalg.norm(weights * [o[2] for o in observations])
                        ),
                    )
            even, odd = records["even"], records["odd"]
            banks.append(
                dict(
                    tag=tag,
                    source_mount=arm,
                    source_readings=records,
                    measured_odd_minus_even=odd["measured_mean"]
                    - even["measured_mean"],
                    difference_se=float(np.hypot(odd["mean_se"], even["mean_se"])),
                    acceleration_difference=odd["time_one_second_derivative"]
                    - even["time_one_second_derivative"],
                )
            )
            if "precision_repeat" in even:
                e, o = even["precision_repeat"], odd["precision_repeat"]
                banks[-1]["measured_acceleration_difference"] = (
                    o["acceleration"] - e["acceleration"]
                )
                banks[-1]["acceleration_difference_se"] = float(
                    np.hypot(e["acceleration_se"], o["acceleration_se"])
                )
    return launches, banks


def gravity_control(vacuum):
    radar = observer_field
    coupling, mass = radar.COUPLING, radar.MASS_SQUARED

    def forward(launch):
        hx, hy, field, velocity = launch
        density = (velocity**2 + mass * field**2) / 2
        hz = (vacuum + coupling * density - hx * hy) / (hx + hy)
        return radar.propagate(
            np.r_[1.0, 1.0, 1.0, hx, hy, hz, field, velocity, 0.0, 0.0, 0.0],
            0.0,
            1.8,
            mass,
            coupling,
            vacuum,
        )

    sources = [forward(s) for s in radar.LAUNCHES]
    results = []
    for spacing in radar.SPACINGS:
        calibration_tape = radar.records(sources[0], 0.32, spacing)
        fit = radar.calibrate(calibration_tape)
        held = []
        for source in sources[1:]:
            for t in (0.4, 0.72, 1.15):
                tape = radar.records(source, t, spacing)
                _, _, _, tidal, field, vel, _ = radar.local_state(tape)
                potential = fit["mass_squared"] * field**2 / 2
                rho, pressure = vel**2 / 2 + potential, vel**2 / 2 - potential
                inferred = float(
                    sum(tidal) + fit["coupling"] * (rho + 3 * pressure) / 2
                )
                row = dict(
                    tidal_trace=float(sum(tidal)),
                    matter_density=rho,
                    pressure=pressure,
                    inferred_lambda=inferred,
                )
                if spacing == radar.SPACINGS[-1]:
                    row["clock_and_scalar_records"] = tape
                held.append(row)
        results.append(
            dict(
                spacing=spacing,
                calibration=fit,
                held=held,
                max_held_error=max(abs(r["inferred_lambda"] - vacuum) for r in held),
            )
        )
    return dict(
        input_lambda=vacuum,
        refinements=results,
        fine_calibration_records=calibration_tape,
    )


def measure():
    model = WideReplay(202, 110, SEED)
    ports, src, outer, adj, degree, trips, last, boundary = settled_lab(model)
    ranges = trips // 2
    modes, clock = source_and_clock()
    outcomes = np.random.default_rng(20261127)
    index = {p: i for i, p in enumerate(ports)}
    spectators = [
        (index[p], ("even", "odd")[int(outcomes.integers(2))])
        for p, q in boundary
        if (p, q) != (ports[src], outer)
    ]
    launches, banks = probe_banks(adj, degree, ranges, src, modes, spectators)

    def address(p):
        owner, port = divmod(p, 12)
        return [str(owner), port]

    def original_load(p):
        owner, slot = divmod(p, 12)
        return 2 * ((model.history(owner).mask >> slot) & 1)

    inner = ports[src]
    write_time = dict(model.links(inner))[outer]
    assert write_time == dict(model.links(outer))[inner] > 0
    return dict(
        level=202,
        seed=SEED,
        receiver=["0", 0],
        archive_scheme=SCHEME,
        protocol_sha256=sha256(
            (HERE / "LEVEL202_GRAVIMETER.md").read_bytes()
        ).hexdigest(),
        clock_archive_sha256=sha256(
            (HERE / "level202_echo_reach.json").read_bytes()
        ).hexdigest(),
        settled_component=dict(
            ports=len(ports),
            owners=len({p // 12 for p in ports}),
            last_incident_write=last,
            written_boundary_contacts=len(boundary),
        ),
        source_contact=[address(inner), address(outer)],
        archived_write_tick=write_time,
        inside_minus_outside_at_write=original_load(inner) - original_load(outer),
        terminal_labels=[address(p) for p in ports],
        roundtrip_ticks=trips.tolist(),
        launch_indices=launches,
        strict_one_unit_releases=strict_releases(model, ports, ranges[ports.index(0)]),
        prepared_source_energy=modes,
        receiver_side_energy_meter=receiver_source_meter(),
        exterior_spectator_readouts=spectators,
        spectator_word_probability=2.0 ** -len(spectators),
        exposure_times=TIMES,
        measured_exchange_clock=clock,
        second_difference_widths=WIDTHS,
        simulated_pairs_per_histogram=SHOTS,
        banks=banks,
        wall_and_volume=wall_and_volume(adj, ranges, src),
        transverse_echoes=transverse_echoes(ranges, ports.index(0), src),
        controlled_einstein_scalar=[gravity_control(v) for v in (0.0, 0.12)],
        native_gravitational_coupling=None,
        observer_vacuum_energy_density_J_m3=None,
        scope="The strict release preserves the settled stock apart from one test load. The finite free-flight apparatus retains the written pattern, prepares the ordinary loads and uses exchange-time and material radar units. The Einstein-scalar controls supply their own metric, proper clock and matter action.",
    )


if __name__ == "__main__":
    print(json.dumps(measure(), separators=(",", ":"), allow_nan=False))
