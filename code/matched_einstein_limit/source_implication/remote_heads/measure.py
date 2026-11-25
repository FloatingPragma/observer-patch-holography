"""Identify both retained parity blocks using one unmarked material return port."""

import json
from math import comb
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_einstein_station.collar_pairing.phase_geometry import (
    remote,
)
from ..full_algebra import native_precision, station


def derivative_rows(generators, port, degree, prime):
    """Native return-effect derivatives over a prime field, in both parity blocks."""
    symmetric = np.triu_indices(4)
    skew = np.triu_indices(4, 1)
    address = slice(4 * port, 4 * port + 4)
    columns = [np.zeros((h.shape[0], 4), dtype=np.int64) for h in generators]
    powers = [[] for _ in generators]
    for column in columns:
        column[address] = np.eye(4, dtype=np.int64)
    for n in range(degree + 1):
        coords = np.zeros((2, 16), dtype=np.int64)
        for b, h in enumerate(generators):
            blocks = powers[b]
            blocks.append(columns[b][address].copy())
            effect = (
                sum(
                    ((-1) ** j * comb(n, j) % prime) * (blocks[j].T @ blocks[n - j])
                    for j in range(n + 1)
                )
                % prime
            )
            if n % 2:
                coords[b, 10:] = effect[skew]
            else:
                coords[b, :10] = effect[symmetric]
            columns[b] = (h @ columns[b]) % prime
        yield n, coords.ravel()


def modular_minor(rows, prime):
    """Select independent rows and return their signed determinant modulo prime."""
    echelon = {}
    orders = []
    pivots = []
    diagonal = []
    for n, row in rows:
        for j in sorted(echelon):
            row = (row - row[j] * echelon[j]) % prime
        nz = np.flatnonzero(row)
        if len(nz):
            j = int(nz[0])
            a = int(row[j])
            echelon[j] = row * pow(a, -1, prime) % prime
            diagonal.append(a)
            pivots.append(j)
            orders.append(n)
    determinant = None
    if len(pivots) == 32:
        inversions = sum(a > b for i, a in enumerate(pivots) for b in pivots[i + 1 :])
        determinant = pow(-1, inversions, prime)
        for a in diagonal:
            determinant = determinant * a % prime
    return orders, determinant


def generators():
    board = remote.formation.federation(1)
    route = remote.formation.region(board)
    return [native_precision.generator(board, route, 1, odd) for odd in (False, True)]


def exchange_control():
    """Rebuild twice the flight generator with contrast turned off on the same edges."""
    board = remote.formation.federation(1)
    route = remote.formation.region(board)
    marked = [tuple(sorted(e)) for e in (route[:2], route[2:])]
    a, b = np.array([e for e in board.seams if e not in marked]).T
    ends = np.r_[a, b]
    adjacency = sparse.csr_matrix(
        (np.ones(len(ends), dtype=np.int64), (ends, np.r_[b, a])),
        shape=(board.n, board.n),
    )
    laplacian = (
        sparse.diags(np.bincount(ends, minlength=board.n), dtype=np.int64) - adjacency
    )
    eye = np.eye(2, dtype=np.int64)
    hs = []
    for odd in (False, True):
        h = sparse.kron(laplacian, np.eye(4, dtype=np.int64), format="csr")
        for slot, edge in enumerate(marked):
            parity = -1 if slot == 0 and odd else 1
            local = 2 * station.iso.ranging.parity_blocks(1, 0, parity)
            np.testing.assert_array_equal(local, local.astype(np.int64))
            local = local.astype(np.int64)
            for i in range(2):
                for j in range(2):
                    part = local[1 + i, 1 + j] - int(i == j) * local[0, 0]
                    coin = np.kron(part, eye) if slot == 0 else np.kron(eye, part)
                    site = sparse.csr_matrix(
                        ([1], ([edge[i]], [edge[j]])), shape=adjacency.shape
                    )
                    h += sparse.kron(site, coin, format="csr")
        stride = 2 if odd else 4
        spectator = np.eye(stride, dtype=np.int64)
        assert (
            h - sparse.kron(h[::stride, ::stride], spectator, format="csr")
        ).nnz == 0
        hs.append(h)
    orders, _ = modular_minor(derivative_rows(hs, 5, 40, 1009), 1009)
    assert len(orders) == 5
    return dict(rank=len(orders), independent_orders=orders)


def detector_effects(hs, ports, quarters):
    """Separate sharp launch and terminal position reads at quarter-tick appointments."""
    addresses = np.array([4 * p + c for p in ports for c in range(4)])
    blocks = []
    for h in hs:
        source = np.zeros((h.shape[0], len(addresses)), complex)
        source[addresses] = np.eye(len(addresses))
        columns = expm_multiply(
            -0.5j * h, source, start=0, stop=quarters / 4, num=quarters + 1
        )
        restricted = columns[1:, addresses].reshape(
            quarters, len(ports), 4, len(ports), 4
        )
        transfers = restricted.transpose(0, 1, 3, 2, 4).reshape(-1, 4, 4)
        blocks.append(transfers.conj().transpose(0, 2, 1) @ transfers)
    return np.stack(blocks, axis=1)


def unread_flight(hs, state, port, time):
    """Leave the source heads in place and remove the single probe with its position unread."""
    after = []
    for generator, block in zip(hs, state):
        start = np.zeros((generator.shape[0], 4), dtype=complex)
        start[4 * port : 4 * port + 4] = np.eye(4)
        wave = expm_multiply(-0.5j * time * generator, start).reshape(-1, 4, 4)
        after.append(np.einsum("axi,ij,ayj->xy", wave, block, wave.conj()))
    return np.array(after)


def counted_survey(hs):
    ports, quarters = (5, 40, 0, 38), 32
    effects = detector_effects(hs, ports, quarters)
    archived = json.loads(Path(station.__file__).with_name("results.json").read_text())
    source = np.array(archived["x_real"]) + 1j * np.array(archived["x_imag"])
    probabilities = station.pair(effects, source)
    shots, seed = 2**36, 815027
    hits = np.random.default_rng(seed).binomial(shots, probabilities)

    center = np.stack([np.eye(4) / 8] * 2)
    basis = station.directions()
    offset = station.pair(effects, center)
    design = np.einsum("ibjk,mbkj->im", effects, basis).real
    inverse = np.linalg.pinv(design)
    estimated = center + np.einsum(
        "m,mbij->bij", inverse @ (hits / shots - offset), basis
    )
    odd_coordinates = basis[:, 1].trace(axis1=1, axis2=2).real
    energy_weights = odd_coordinates @ inverse
    # For a sum of independent binomial frequencies, Hoeffding uses the
    # squared l2 norm of its fixed weights, divided by shots per setting.
    radius = np.linalg.norm(energy_weights) * np.sqrt(np.log(200) / (2 * shots))
    low_dose = dict(
        shots_per_setting=shots,
        seed=seed,
        yes_counts=hits.tolist(),
        hilbert_schmidt_error=float(np.linalg.norm(source - estimated)),
        estimated_min_eigenvalue=float(np.linalg.eigvalsh(estimated).min()),
        estimated_parked_energy=16 + float(estimated[1].trace().real),
        energy_99pct_halfwidth=float(radius),
    )

    # Precondition by a variance envelope from the detector alone.
    upper = np.linalg.eigvalsh(effects).max(axis=(1, 2))
    scale = np.sqrt(upper)
    weighted = np.linalg.pinv(design / scale[:, None]) / scale[None, :]
    high_shots = 2**52
    high_seed = 815028
    rng = np.random.default_rng(high_seed)
    counts = {}

    def acquire(label, state):
        counts[label] = rng.binomial(high_shots, station.pair(effects, state)).tolist()
        return np.array(counts[label]) / high_shots

    def bernstein(weights, log_level, independent_arms):
        variance = independent_arms * ((weights**2) @ upper) / high_shots
        bound = np.max(abs(weights), axis=-1) / high_shots
        return np.sqrt(2 * log_level * variance) + 2 * log_level * bound / 3

    reference = acquire("reference", source)
    state_fit = center + np.einsum(
        "m,mbij->bij", weighted @ (reference - offset), basis
    )
    state_error = float(np.linalg.norm(bernstein(weighted, np.log(2 * 31 / 0.005), 1)))
    smallest = float(np.linalg.eigvalsh(state_fit).min())
    if state_error >= smallest:
        raise ArithmeticError(
            "Reference precision does not resolve the modular logarithm."
        )
    log_error = -np.log1p(-state_error / smallest)
    modular_weights = -station.pair(basis, station.log_blocks(state_fit)) @ weighted
    energy_weights = odd_coordinates @ weighted
    log_level = np.log(2 * 4 / 0.005)
    right = np.trace(source.sum(axis=0).reshape(2, 2, 2, 2), axis1=0, axis2=2)
    left = np.trace(source.reshape(2, 2, 2, 2, 2), axis1=2, axis2=4)
    arms = {
        "ordinary": np.stack([source.sum(axis=0) / 2] * 2),
        "pairing": np.stack([np.kron(block, right) for block in left]),
    }
    product = arms["pairing"]
    collar = station.log_blocks(product) - station.log_blocks(source)
    forward = -float(station.pair(collar, source))
    backward = float(station.pair(collar, product))
    # First port to the same port at time two, with the clock schedule above.
    return_effect = effects[(4 * 2 - 1) * len(ports) ** 2]
    commutator = collar @ return_effect - return_effect @ collar
    after_flight = unread_flight(hs, source, port=5, time=2)
    collar_charge = dict(
        state_to_product_relative_entropy=forward,
        product_to_state_relative_entropy=backward,
        unpair_slope=forward + backward,
        unread_flight_charge_change=float(station.pair(collar, after_flight - source)),
        unread_flight_parked_energy_change=float(
            (after_flight - source)[1].trace().real
        ),
        at_time_two_port_five_return_commutator_hs=float(np.linalg.norm(commutator)),
    )
    held = acquire("held_baseline", source)
    signed = {}
    for label, arm in arms.items():
        contrast = acquire(label, arm) - held
        signed[label] = dict(
            modular_nats=float(modular_weights @ contrast),
            modular_halfwidth=float(
                bernstein(modular_weights, log_level, 2) + 2 * log_error
            ),
            parked_energy=float(energy_weights @ contrast),
            energy_halfwidth=float(bernstein(energy_weights, log_level, 2)),
        )
    return dict(
        ports=list(ports),
        quarter_ticks=quarters,
        order="appointment, terminal port, launch port (last index fastest)",
        collar_charge=collar_charge,
        low_dose=low_dose,
        singular_min=float(np.linalg.svd(design, compute_uv=False)[-1]),
        high_dose=dict(
            shots_per_setting=high_shots,
            seed=high_seed,
            counts=counts,
            simultaneous_confidence=0.99,
            smallest_estimated_eigenvalue=smallest,
            hilbert_schmidt_error=float(np.linalg.norm(source - state_fit)),
            hilbert_schmidt_radius=state_error,
            held_contrasts=signed,
        ),
    )


def run():
    hs = generators()
    prime = 1009
    orders, determinant = modular_minor(derivative_rows(hs, 5, 40, prime), prime)
    # A second prime checks the very same ordered minor.
    selected = set(orders)
    check_prime = 1013
    check = (
        (n, row) for n, row in derivative_rows(hs, 5, 40, check_prime) if n in selected
    )
    check_orders, check_determinant = modular_minor(check, check_prime)
    assert len(orders) == 32 and check_orders == orders
    return dict(
        port=5,
        probe_charge=1,
        orders=orders,
        determinant_residues={
            str(prime): determinant,
            str(check_prime): check_determinant,
        },
        exchange_only=exchange_control(),
        finite_read=counted_survey(hs),
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
