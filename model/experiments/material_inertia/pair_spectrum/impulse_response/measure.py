"""A signed centre-velocity test on the same free bound pair, with every filter loss."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .. import dispersion, measure as flight

HERE = Path(__file__).resolve().parent
THETA = 0.25
SHOTS = 2**51


def response(h, energies, basis, x, impulse, cluster_error, rest_error, residual):
    """Restrict the full and projected double commutators to the initial doublet."""
    edges = h.tocoo(copy=True)
    edges.data *= -(impulse[edges.row] - impulse[edges.col]) * (
        x[edges.row] - x[edges.col]
    )
    kernel = edges.tocsr()
    full = basis[:, :2].T @ (kernel @ basis[:, :2])
    full_norm = np.asarray(abs(kernel).sum(axis=1)).max()
    xb = basis.T @ (x[:, None] * basis)
    gb = basis.T @ (impulse[:, None] * basis)
    commutator = (energies[:, None] - energies) * xb
    projected = gb @ commutator - commutator @ gb

    def interval(matrix, error):
        low, high = np.linalg.eigvalsh(matrix)[[0, -1]]
        return [float(low - error), float(high + error)]

    # Rotate the numerical 24-frame onto the exact spectral subspace. The polar
    # overlap differs from identity by <= eps². H has no cross-band block.
    width = np.ptp(energies)
    rotation = 2 * cluster_error
    energy_error = residual + 3 * 80 * cluster_error**2
    # X has norm <=1, the signed impulse has norm <=2. Both commutators are bilinear.
    band_error = 8 * energy_error + 32 * rotation * (width + energy_error)
    band_norm = np.linalg.norm(projected, 2)
    return dict(
        all_ports=interval(full, 4 * rest_error * full_norm + 1e-12),
        within_cluster=interval(
            projected[:2, :2],
            band_error + 4 * (rest_error + rotation) * (band_norm + band_error) + 1e-14,
        ),
    )


def count_outcomes(state, pointer, admission, seed, population_bank):
    """Record filter rejection and every coordinate bin after the two preparations."""
    values, inverse = np.unique(pointer.round(12), return_inverse=True)
    amplitude = np.abs(state) ** 2
    probabilities = np.bincount(inverse, weights=amplitude, minlength=len(values))
    live = float(probabilities.sum())
    if live > 1 + 2e-12:
        raise ValueError("Computed detector gain exceeds one")
    rng = np.random.default_rng(seed)
    selected = int(rng.binomial(SHOTS, admission))
    # Condition only on first (rest) filter admission. Failed second reads score zero.
    probs = np.r_[probabilities, max(0.0, 1 - live)]
    normalized = probs / probs.sum()
    if population_bank is not None:
        population_bank.append(
            dict(
                seed=list(seed), admission=admission, probabilities=normalized.tolist()
            )
        )
    counts = rng.multinomial(selected, normalized)
    return dict(
        enrolled=SHOTS,
        failed_preparation=SHOTS - selected,
        failed_band_filter=int(counts[-1]),
        values=values.tolist(),
        counts=counts[:-1].tolist(),
        mean=float(counts[:-1] @ values / selected),
        exact_mean=float(amplitude @ pointer),
        halfwidth_99_simultaneous=float(
            np.sqrt(2 * np.log(2 * 24 / 0.01) / selected) + 1e-12
        ),
    )


def clocks(h, energies, basis, state, x, impulse, admission, seed, population_bank):
    output = {}
    # Matrices are shifted by a scalar energy to avoid irrelevant fast phases.
    shifted = h - flight.sparse.eye(h.shape[0], format="csr") * flight.CENTER
    columns = state[:, None] * np.exp(1j * THETA * impulse[:, None] * np.array([-1, 1]))
    for name, time in (("free", 1 / 16), ("band_selected", 32768)):
        if name == "free":
            waves = expm_multiply(-1j * time * shifted, columns)
        else:
            waves = basis @ (
                np.exp(-1j * time * (energies - energies[0]))[:, None]
                * (basis.T @ columns)
            )
        trials = [
            count_outcomes(
                waves[:, i], x, admission, [*seed, int(time), i], population_bank
            )
            for i in range(2)
        ]
        scale = 2 * THETA * time
        output[name] = dict(
            time=time,
            records=trials,
            signed_position_response=(trials[1]["mean"] - trials[0]["mean"]) / scale,
            halfwidth_99_simultaneous=sum(
                r["halfwidth_99_simultaneous"] for r in trials
            )
            / scale,
            nominal_response=(trials[1]["exact_mean"] - trials[0]["exact_mean"])
            / scale,
        )
    return output


def run(population_bank=None):
    board = flight.federation(0)
    excess, vacancy, adjacent, kinetic, h, cost = flight.defect.generators(board)
    original = json.loads((flight.HERE / "results.json").read_text())
    reduced = flight.PairReduction(h, cost)
    energies, basis, residual, _ = reduced.lowest(h, 24)
    cut = (energies[1] + energies[2]) / 2
    if reduced.inertia(cut)["upper_lines_below"] != 2:
        raise ValueError("The enrolled rest doublet cannot be separated")
    cluster_error = max(
        original["cluster_projector_error"], residual / (flight.CUT - energies[-1])
    )
    rest_error = residual / (
        min(cut - energies[1], flight.CUT - energies[-1]) - residual
    )
    f0 = dispersion.geometry(board)[0][1]
    panels = []
    for launch in original["flights"]:
        if launch["admission"] < 0.001:
            continue
        owner = launch["owner"]
        axis = min(
            v // 12
            for edge in board.external
            for u, v in (edge, edge[::-1])
            if u // 12 == owner
        )
        f = f0[:, axis] / max(abs(f0[:, axis]))
        x = (f[excess // 12] + f[vacancy // 12]) / 2
        relative = f[excess // 12] - f[vacancy // 12]
        s = flight.defect.ordered_index(*launch["ports"], board.n)
        g = basis[s, :2]
        admission = float(g @ g)
        state = basis[:, :2] @ (g / np.sqrt(admission))
        controls = {}
        for label, generator in (("native", h), ("contrast_removed", kinetic)):
            evolved = expm_multiply(
                -1j
                * (
                    generator
                    - flight.sparse.eye(h.shape[0], format="csr")
                    * generator.diagonal().mean()
                ),
                state,
            )
            controls[label] = float(np.linalg.norm(evolved[adjacent]) ** 2)
        # The same coordinate instrument on a real uniform launch has positive
        # instantaneous response under the exchange hopping, by its Dirichlet form.
        edges = kinetic.tocoo()
        controls["uniform_common_response"] = float(
            -np.mean(
                np.bincount(
                    edges.row,
                    weights=edges.data * (x[edges.row] - x[edges.col]) ** 2,
                    minlength=kinetic.shape[0],
                )
            )
        )
        cases = {}
        for k, (label, impulse) in enumerate(
            (("common", x), ("signed_load", relative))
        ):
            cases[label] = dict(
                derivative_intervals=response(
                    h, energies, basis, x, impulse, cluster_error, rest_error, residual
                ),
                finite=clocks(
                    h,
                    energies,
                    basis,
                    state,
                    x,
                    impulse,
                    admission,
                    [174920, owner, k],
                    population_bank,
                ),
            )
        panels.append(
            dict(
                owner=owner,
                axis=axis,
                rest_admission=admission,
                cases=cases,
                controls=controls,
            )
        )
    return dict(
        theta=THETA,
        preparations=panels,
        residual=residual,
        rest_projector_error=rest_error,
        band_projector_error=cluster_error,
        native_energy_units="exchange",
        spatial_units="normalized owner graph harmonic",
    )


if __name__ == "__main__":
    cmd = argparse.ArgumentParser(description=__doc__)
    cmd.add_argument("--write", action="store_true")
    args = cmd.parse_args()
    from . import replay

    if args.write:
        # Preserve the published realization. Regeneration is possible on a
        # runtime that exactly reproduces the original full acquisition law.
        from hashlib import sha256

        bank = []
        result = run(bank)
        outputs = [
            ("results.json", result, replay.RESULT_SHA256),
            ("acquisition.json", bank, replay.POPULATIONS_SHA256),
        ]
        serialized = [
            (name, json.dumps(data, indent=2, allow_nan=False) + "\n", digest)
            for name, data, digest in outputs
        ]
        if any(
            sha256(data.encode()).hexdigest() != digest
            for _, data, digest in serialized
        ):
            raise ValueError(
                "Fresh proposal differs from the original acquisition; use the default physical replay to verify it"
            )
        for name, data, _ in serialized:
            (HERE / name).write_text(data)
    else:
        replay.verify()
    print(
        "Verified all original impulse counts and independently recomputed the recoil gates."
    )
