"""Compare two vacuum-calibrated written hands exposed to a freely moving material source."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_null_balance.regional import verify
from model.experiments.material_inertia.dark_collar import label

HERE = Path(__file__).resolve().parent
SHOTS = 2**36
SEED = 26112327
TIMES = (.5, 1., 2.)
LEVELS = (0, 1, 2)
SOURCES = ("vacuum", "contact", "one_bridge", "two_bridges")
ARMS = ("empty_meter", "occupied_meter", "restored_contact")
GAPS = (32., 33.)


def hand_columns():
    """Compile the written channels and both Ramsey analyses from the 99-state law."""
    law = candidate.native_law
    exchange, cost = law.local_term_parts()
    native = exchange + candidate.KAPPA * cost
    pair = law.SRC.index((1, 1))
    identity = np.eye(law.N_SRC)
    write = law.echo()
    polarizations = [np.kron(label(2, sign)[:, None], identity) for sign in (1, -1)]
    repair = np.zeros_like(native)
    for gap, basis in zip(GAPS, polarizations, strict=True):
        energies = np.full(law.N_SRC, gap)
        energies[pair] = 0
        np.testing.assert_allclose(native @ basis, basis * energies, rtol=0, atol=1e-13)
        record = basis[:, pair]
        repair += gap * np.outer(record, record)
    for gap, basis in zip(GAPS, polarizations, strict=True):
        np.testing.assert_allclose((native + repair) @ basis, gap * basis, rtol=0, atol=1e-13)
    guard = np.kron(label(0, 1)[:, None], identity)
    np.testing.assert_allclose(native @ guard, 0, rtol=0, atol=1e-13)
    np.testing.assert_array_equal(repair[:, :law.N_SRC], 0)
    # Positive sine is checked at the physical analyzer, before use on arbitrary mixtures.
    sharp = law.SRC.index((2, 0))
    levels = np.column_stack([basis[:, pair] for basis in polarizations])
    initial = levels.T @ write[:, sharp]
    analyzed = write @ levels
    sine = expm(-1.5j * np.pi * native) @ analyzed
    for time in TIMES:
        wave = initial * np.exp(-1j * time * np.asarray(GAPS))
        yes = abs(np.array([analyzed[sharp], sine[sharp]]) @ wave)**2
        np.testing.assert_allclose(2*yes-1, [np.cos(time), np.sin(time)], atol=2e-13, rtol=0)
    singles = [law.SRC.index(word) for word in ((1, 0), (0, 1))]
    return native[np.ix_(singles, singles)]


def instruments(board, blank):
    b, near = min(edge for edge in board.external if edge[0] < 12)
    p = min(i for i in board.neighbours[b] if i < 12)
    first_ring = {v//12 for u, v in board.external if u < 12}
    far = min(v for a, c in board.external for u, v in ((a, c), (c, a))
              if u//12 == near//12 and v//12 not in first_ring | {0})
    base = candidate.generator(board).tolil()
    for i in board.neighbours[p]:
        edge = [p, i]
        base[np.ix_(edge, edge)] -= blank
    base = base.tocsr()
    active = np.arange(board.n) != p
    np.testing.assert_allclose(base[active][:, p].toarray(), 0, rtol=0, atol=1e-13)
    # A common scalar phase is invisible in a number or hand reading.
    base -= 6*candidate.KAPPA * sparse.eye(board.n, format="csr")
    return base, p, b, near, far


def propagate(base, gateway, wave, removed):
    collision = sparse.csr_matrix(([removed], ([gateway], [gateway])), shape=base.shape)
    h = -1j*(base-collision)
    return expm_multiply(h, wave, start=0., stop=2., num=5, traceA=h.diagonal().sum())


def readout(low, high, t, gateway, vacuum):
    if vacuum:
        pop = np.array([0., 0., 1.])
        fringe = pop * np.exp(1j*t)
    else:
        entry = (abs(low)**2 + abs(high)**2)/2
        overlap = high.conj()*low*np.exp(1j*t)
        def bins(w):
            patch = w[:12].sum()
            return np.array([w[gateway], patch-w[gateway], w.sum()-patch])
        pop, fringe = bins(entry), bins(overlap)
    yes = (pop + np.stack([fringe.real, fringe.imag]))/2
    outcomes = np.stack([yes, pop-yes], axis=-1)
    np.testing.assert_allclose(outcomes.sum(axis=(1,2)), 1, rtol=0, atol=2e-12)
    if outcomes.min() < -1e-13:
        raise ArithmeticError("negative joint reading")
    return outcomes.clip(0., 1.)


def phase(xy):
    return float(np.arctan2(xy[1], xy[0]))


def compare(left, right):
    shift = np.angle(np.exp(1j*(left["phase_rad"]-right["phase_rad"])))
    radius = left["angle_radius"] + right["angle_radius"]
    return dict(shift_rad=float(shift), simultaneous_interval=[float(shift-radius),float(shift+radius)])


def acquire(outcomes, rng, coordinate_radius):
    counts = np.stack([rng.multinomial(SHOTS, row.ravel()/row.sum()).reshape(3,2)
                       for row in outcomes])
    xy = 2*counts[:,:,0].sum(axis=1)/SHOTS-1
    angle_radius = np.arcsin(np.sqrt(2)*coordinate_radius / np.linalg.norm(xy))
    if not np.isfinite(angle_radius):
        raise ArithmeticError("the measured fringe does not resolve an angle")
    exact_xy = 2*outcomes[:,:,0].sum(axis=1)-1
    return dict(counts=counts.tolist(), phase_rad=phase(xy), angle_radius=float(angle_radius),
                coherence=float(np.linalg.norm(xy)), exact_phase_rad=phase(exact_xy),
                source_at_gateway=float(outcomes[0,0].sum()),
                source_in_cell=float(outcomes[0,:2].sum()))


def run():
    blank = hand_columns()
    rng = np.random.default_rng(SEED)
    axes = 2*len(LEVELS)*len(ARMS)*len(SOURCES)*len(TIMES)
    coordinate_radius = np.sqrt(2*np.log(2*axes/.01)/SHOTS)
    boards = []
    for level in LEVELS:
        board = federation(level)
        base, parked, gate, near, far = instruments(board, blank)
        prepared = np.zeros((board.n, 3))
        prepared[[gate, near, far], range(3)] = 1
        mobile = {g: propagate(base, gate, prepared, g) for g in (0., *GAPS)}
        source_energy = {}
        for name, port in (("one_bridge", near), ("two_bridges", far)):
            degree = int(board.degree[port])
            odd = rng.binomial(SHOTS, .5, degree)
            source_energy[name] = dict(port=int(port), shots_per_seam=SHOTS, odd_counts=odd.tolist(),
                                       native_energy=candidate.KAPPA*6 + float(odd.sum()/SHOTS))
        arms = {}
        for arm in ARMS:
            gaps = GAPS if arm == "occupied_meter" else (0., 0.)
            rows = []
            for j, source in enumerate(SOURCES):
                for time in TIMES:
                    tick = int(2*time)
                    low, high = (mobile[g][tick,:,max(j-1,0)] for g in gaps)
                    result = acquire(readout(low, high, time, gate, source=="vacuum"), rng, coordinate_radius)
                    rows.append(dict(source=source, appointment=tick, **result))
            arms[arm] = rows
        differences = []
        for empty, occupied, restored in zip(*arms.values(), strict=True):
            vacuum = arms["empty_meter"][TIMES.index(empty["appointment"]/2)]
            differences.append(dict(source=empty["source"], appointment=empty["appointment"],
                                    resident_minus_empty=compare(occupied, empty),
                                    empty_minus_vacuum=compare(empty, vacuum),
                                    restored_minus_empty=compare(restored, empty)))
        boards.append(dict(cells=board.cells, parked_port=int(parked), gateway=int(gate),
                           source_energy=source_energy, arms=arms, differences=differences))
    return dict(candidate=candidate.NAME, design_sha256=hashlib.sha256((HERE/"DESIGN.md").read_bytes()).hexdigest(),
                seed=SEED, shots_per_setting=SHOTS, simultaneous_coverage=.99, screens=boards)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    else:
        verify(json.loads((HERE/"results.json").read_text()), result)
    for board in result["screens"]:
        for row in board["differences"]:
            if row["appointment"]==4:
                print(board["cells"], row["source"], row["resident_minus_empty"])


if __name__ == "__main__":
    main()
