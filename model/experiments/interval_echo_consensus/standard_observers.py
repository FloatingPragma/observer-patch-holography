"""Evaluate supplied contrast gaps on the registered patch clock and gauge menus.

The concentration of signed amplitudes at a parked rail provides a separate
interference reading. All routines use existing carrier, gauge and clock
implementations. The ordinary candidate and its module-level caches are untouched.
"""

import importlib.util
import json
from math import pi, remainder, sin
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import eye

from model.candidate_laws.contrast_memory_spherical_v1 import candidate as contrast
from model.candidate_laws.interval_echo import clock, currents
from model.executable_model import cell, federation
from model.experiments.interval_echo_clock import publication


def local_copy(module):
    """Keep candidate substitutions and cached clock gates private to this replay."""
    spec = importlib.util.spec_from_file_location(
        f"{module.__package__}._standard_{module.__name__.split('.')[-1]}",
        Path(module.__file__),
    )
    copy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(copy)
    return copy


def seam_gates(e1, e2):
    law = clock.law
    echo = law.echo()
    exchange, _ = law.local_term_parts()
    record_energies = [0.0, *((0.0, e1, e2)[abs(d)] for _, d in law.LABELS)]
    h = exchange + (echo * np.repeat(record_energies, law.N_SRC)) @ echo.T
    basis = np.eye(law.DIM)
    rails = [law.SRC.index(word) for word in ((2, 0), (0, 2))]
    if np.max(abs(h[:, rails] - basis[:, rails] @ h[np.ix_(rails, rails)])) > 1e-12:
        raise ValueError("The clock rails do not form a reducing subspace.")
    pulses = {t: expm(-1j * t * h) for t in (pi, pi / 2, pi / 16)}
    dimer = {t: matrix[np.ix_(rails, rails)] for t, matrix in pulses.items()}
    phase = np.exp(-1j * pi * remainder(e2, 2.0))
    rates = {}
    column_error = 0.0
    for a in (0, 2):
        for b in (0, 2):
            source, target = (law.SRC.index(word) for word in ((a, b), (b, a)))
            swap = pulses[pi][:, source]
            if np.max(abs(swap - (phase if a != b else 1) * basis[:, target])) > 1e-11:
                raise ValueError(
                    "The native contact did not exchange these sharp rails."
                )
            column_error = max(
                column_error, float(np.max(abs(swap - basis[:, target])))
            )
            if a != b:
                rates[a] = float(abs(pulses[pi / 2][target, source]) ** 2)
    if max(abs(rate - 0.5) for rate in rates.values()) > 1e-12:
        raise ValueError("The half-contact stopped giving fair transfers.")
    return echo, rates, column_error, dimer, phase


def carried_hands(gates):
    echo, rates, error, dimer, phase = gates
    engine = local_copy(clock)
    engine._primitives = lambda: (echo, rates, error)
    engine._dimer_pulse = dimer.__getitem__
    swap_counts = []

    class TransportRegisters(engine.Registers):
        def __init__(self, board, seed):
            super().__init__(board, seed)
            self.coherent_swaps = [0, 0]
            swap_counts.append(self.coherent_swaps)

        def swap(self, a, b, owner):
            tagged = [(p, *self.coherent[p]) for p in (a, b) if p in self.coherent]
            if len(tagged) == 2:
                if tagged[0][1] is not tagged[1][1]:
                    raise ValueError(
                        "A contact joining two coherent hands requires a joint-state simulation."
                    )
                tagged[0][1].state *= phase
                self.coherent_swaps[1] += 1
            elif tagged:
                port, hand, rail = tagged[0]
                other = self.loads[b if port == a else a]
                if other not in (0, 2):
                    raise ValueError(
                        "The companion rail is outside the prepared sharp menu."
                    )
                # Exactly one endpoint holds a double in this coherent branch.
                hand.state[rail if other == 0 else 1 - rail] *= phase
                self.coherent_swaps[0] += 1
            super().swap(a, b, owner)

        def fair(self, site, owner):
            if any(12 * site + i in self.coherent for i in (2, 3, 4)):
                raise ValueError(
                    "The fair-source contact or detector touches a coherent hand."
                )
            return super().fair(site, owner)

    engine.Registers = TransportRegisters
    archived = json.loads((publication.HERE / "publication.json").read_text())
    event_differences = bit_differences = 0
    probability_error = 0.0
    for specimen in archived["specimens"]:
        board = federation({20: 0, 80: 1, 320: 2}[specimen["cells"]])
        for analyzer in (0, 1):
            run = engine.material_replay(
                board,
                publication.SCHEDULE,
                age_bits=specimen["width"],
                seed=specimen["seed"],
                analyzer=analyzer,
            )
            for field, stored in (
                ("primary", "primary"),
                ("tick", "moving_public"),
                ("home_tick", "home_public"),
            ):
                event_differences += sum(
                    event[field] != expected
                    for event, expected in zip(
                        run["local_events"], specimen[stored], strict=True
                    )
                )
            for name, stored in specimen["endpoint"].items():
                terminal = run["oscillators"][name]
                delta = (
                    np.asarray(terminal["first_rail_probabilities"])
                    - stored["probabilities"]
                )
                probability_error = max(probability_error, float(np.max(abs(delta))))
                bit_differences += int(
                    terminal["bit"] != stored["published_bits_at_settings"][analyzer]
                )
                if (
                    terminal["owners"] != [terminal["bit"]] * 2
                    or run["reunion"][name] != stored["public_count"]
                ):
                    raise ValueError(
                        "The two patch owners disagree at fixed-home reunion."
                    )
            if run["reunion"]["traveller"] != run["home_retained_traveller"]:
                raise ValueError("The home patch did not retain the visiting word.")
    return dict(
        material_trips=len(swap_counts),
        public_event_differences=event_differences,
        seeded_pointer_differences=bit_differences,
        max_endpoint_probability_error=round(probability_error, 12),
        published_probability_tolerance_pass=probability_error <= 1e-10,
        max_single_rail_swaps=max(x[0] for x in swap_counts),
        max_two_rail_swaps=max(x[1] for x in swap_counts),
    )


def parked_arm(dimer):
    """Split one double, expose only one arm, recombine and retain all three positions."""
    ports = (6, 7, 2)
    board = cell()
    if any(b not in board.neighbours[ports[0]] for b in ports[1:]):
        raise ValueError("The parked arm leaves its native cell.")
    state = np.array([1.0, 0.0, 0.0], complex)
    for indices, t in (([0, 2], pi / 2), ([0, 1], pi / 16), ([0, 2], pi / 2)):
        state[indices] = dimer[t] @ state[indices]
    probabilities = abs(state) ** 2
    if abs(probabilities.sum() - 1) > 1e-12:
        raise ValueError("A parked-arm outcome was lost.")
    return dict(ports=list(ports), probabilities=np.round(probabilities, 12).tolist())


def observer_tests(e1, e2):
    """Test magnitude-one and magnitude-two gaps, in native exchange splitting units."""
    e1, e2 = float(e1), float(e2)
    if not np.all(np.isfinite([e1, e2])):
        raise ValueError("Both contrast energies must be finite.")
    gates = seam_gates(e1, e2)
    gauge = local_copy(currents)
    gauge.KAPPA = e1  # Binary pairs access only contrasts zero and one.
    gauge_result = gauge._response(cell(), True, 1.0)
    board = cell()
    levels = np.linalg.eigvalsh(board.laplacian().toarray() / 2)
    for screen in (board, federation(0), federation(1), federation(2)):
        difference = (
            contrast.generator(screen, e1)
            - screen.laplacian() / 2
            - 6 * e1 * eye(screen.n)
        )
        if np.max(abs(difference.data), initial=0) > 1e-12:
            raise ValueError(
                "The centered one-load law differs from the spatial generator."
            )
    return dict(
        slow_rank=int(np.isclose(levels, (5 - np.sqrt(5)) / 2).sum()),
        gauge=gauge_result["guarded"],
        sharp_swap_phase_error=float(
            f"{2 * abs(sin(pi * remainder(e2, 2.0) / 2)):.5g}"
        ),
        published_sharp_primitive_tolerance_pass=gates[2] <= 1e-11,
        carried_hands=carried_hands(gates),
        parked_arm=parked_arm(gates[3]),
    )
