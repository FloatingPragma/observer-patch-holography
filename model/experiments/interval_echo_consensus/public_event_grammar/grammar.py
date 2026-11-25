"""Finite complete native/adaptive program grammar and effect signatures."""

from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import argparse
import inspect
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import clock, cohorts, observations, replay
from ..signed_effects.sequential import compare

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RECEIPT = HERE / "results.json"


def require(test, message):
    if not test:
        raise ValueError(message)


@dataclass(frozen=True)
class Halt:
    label: str = "complete"


@dataclass(frozen=True)
class Step:
    name: str
    branches: tuple  # (public label, Kraus matrix, next program)


def compile_program(program, dimension):
    """Flatten the entire tree; hidden records can be mapped away afterwards."""
    if isinstance(program, Halt):
        return {(program.label,): np.eye(dimension, dtype=complex)}
    require(isinstance(program, Step), "unknown grammar constructor")
    require(
        len({label for label, _, _ in program.branches}) == len(program.branches),
        "duplicate primitive labels",
    )
    gram = np.zeros((dimension, dimension), complex)
    result = {}
    for label, k, next_program in program.branches:
        require(
            k.shape == (dimension, dimension) and np.isfinite(k).all(),
            "invalid native map",
        )
        gram += k.conj().T @ k
        for history, following in compile_program(next_program, dimension).items():
            result[(label, *history)] = following @ k
    require(
        np.max(np.abs(gram - np.eye(dimension))) < 2e-12,
        "incomplete primitive: retain refusal outcomes",
    )
    return result


def event_effect(branches, event):
    dim = next(iter(branches.values())).shape[0]
    return sum(
        (k.conj().T @ k for label, k in branches.items() if event(label)),
        np.zeros((dim, dim), complex),
    )


def runtime_law(branches, rho):
    require(
        np.max(np.abs(rho - rho.conj().T)) < 1e-12
        and np.linalg.eigvalsh(rho).min() > -1e-12,
        "positive initial density required",
    )
    require(abs(np.trace(rho) - 1) < 1e-12, "normalized initial law required")
    weights = {
        label: float(np.trace(k @ rho @ k.conj().T).real)
        for label, k in branches.items()
    }
    require(
        min(weights.values()) > -1e-12 and abs(sum(weights.values()) - 1) < 2e-11,
        "runtime mass lost",
    )
    return weights


def unitary(name, u, following):
    return Step(name, ((name, u, following),))


def material_read(following):
    return Step(
        "sharp material read",
        tuple(
            (
                f"material:{a}{b}",
                np.diag(np.tile(np.arange(9) == j, 11)).astype(complex),
                following(a, b),
            )
            for j, (a, b) in enumerate(law.SRC)
        ),
    )


def coin(probabilities, following, dimension=99):
    require(
        all(p >= 0 for p in probabilities) and sum(probabilities) == 1,
        "normalized classical law required",
    )
    return Step(
        "classical random stream",
        tuple(
            (f"coin:{j}", np.sqrt(float(p)) * np.eye(dimension), following(j))
            for j, p in enumerate(probabilities)
        ),
    )


def native_program(refined=False, post=False):
    """Full99 dirty-collar execution, free exposure and record-based stopping."""
    h0, h1 = law.local_term_parts()
    free = expm(-0.173j * (h0 + 8 * h1))

    def after_read(a, b):
        if a + b == 0:
            return Halt("refused:empty")
        end = Halt("early:unequal" if a != b else "complete:balanced")
        if post:
            end = unitary("postrecord echo", law.echo(), end)
        if a < b:
            end = unitary("adaptive free contact", free, end)
        if refined:
            end = coin((Fraction(1, 4), Fraction(3, 4)), lambda _: end)
        return end

    return unitary(
        "native free contact",
        free,
        unitary("native echo", law.echo(), material_read(after_read)),
    )


def probes(effect):
    """State-independent Hermitian signature: d diagonals and complex pairs."""
    d = len(effect)
    output = list(np.diag(effect).real)
    for i in range(d):
        for j in range(i + 1, d):
            output.append(
                float((effect[i, i] + effect[j, j]).real / 2 + effect[i, j].real)
            )
            output.append(
                float((effect[i, i] + effect[j, j]).real / 2 - effect[i, j].imag)
            )
    return np.asarray(output)


def reconstruct(signature, dimension):
    signature = np.asarray(signature)
    require(
        signature.shape == (dimension**2,),
        "complete state-independent probe signature required",
    )
    result = np.diag(signature[:dimension]).astype(complex)
    pos = dimension
    for i in range(dimension):
        for j in range(i + 1, dimension):
            half = (result[i, i] + result[j, j]).real / 2
            result[i, j] = signature[pos] - half + 1j * (half - signature[pos + 1])
            result[j, i] = result[i, j].conj()
            pos += 2
    return result


def exact_counterexample():
    # All effects are scalar: this obstruction precedes tomography or physics.
    fine = [Fraction(1, 4), Fraction(1, 4), Fraction(1, 2)]
    coarse = [Fraction(1, 2), Fraction(1, 2)]

    def changed(p):
        norm = sum(x * x for x in p)
        return [x * x / norm for x in p]

    a, b = changed(fine), changed(coarse)
    return dict(
        equal_event_effect="I/2",
        fine_weights=list(map(str, a)),
        coarse_weights=list(map(str, b)),
        normalized_each_context=sum(a) == sum(b) == 1,
        fine_union=str(a[0] + a[1]),
        coarse_same_effect=str(b[0]),
        discrepancy=str(b[0] - a[0] - a[1]),
    )


def bounded_address_reader(q, tries=12):
    """Exact compressed record law, including the address-timeout leaf."""
    require(
        type(q) is int and q >= 1 and type(tries) is int and tries >= 1,
        "bounded address request required",
    )
    addresses = 1 << (q - 1).bit_length()
    reject = Fraction(addresses - q, addresses)
    result = {
        f"accepted:{j}": sum(
            (reject**k / Fraction(addresses) for k in range(tries)), Fraction()
        )
        for j in range(q)
    }
    result["timeout"] = reject**tries
    require(sum(result.values()) == 1, "address refusals lost")
    return {k: str(v) for k, v in result.items()}


def source_initial_law():
    """Literal bulk-entry independent tag/matter preparation, before optics."""
    import itertools

    ps = (Fraction(1, 25), Fraction(1, 25), Fraction(4, 25))
    weights = {
        "".join(map(str, bits)): np.prod(
            [p if bit else 1 - p for bit, p in zip(bits, ps, strict=True)]
        )
        for bits in itertools.product((0, 1), repeat=3)
    }
    require(
        sum(weights.values()) == 1 and min(weights.values()) >= 0,
        "source preparation law incomplete",
    )
    return dict(
        source="observations.measure: independent left/right tag offers and matter offer",
        exact_complete_weights={key: str(value) for key, value in weights.items()},
        initial_sharp_material_words="all other specified occupied rails and blank/source reference registers are deterministic",
        quantum_density="diagonal preparation law followed by the deposited native tag write; positivity/normalization follow before any event-effect quotient",
    )


def source_audit():
    files = [law, clock, cohorts, observations, replay]
    hashes = {
        str(Path(module.__file__).relative_to(ROOT)): sha256(
            Path(module.__file__).read_bytes()
        ).hexdigest()
        for module in files
    }
    # Semantic guards detect changes to actual entrypoint draw contracts.
    checks = {
        "candidate_entry_uses_material_replay": "clock.material_replay"
        in inspect.getsource(replay.generate),
        "coherent_terminal_born_draw": "rng.random() < probabilities[analyzer]"
        in inspect.getsource(clock.SeamOscillator.terminal),
        "native_fair_draw": "self.rng.random() < self.rates"
        in inspect.getsource(clock.Registers.fair),
        "finite_rejection_bound": clock.ADDRESS_TRIES == 12 and clock.STEPS == 12,
        "cohort_independent_primary_draws": "source.integers(2"
        in inspect.getsource(cohorts.sample),
        "cohort_fair_bits": "getrandbits(1)" in inspect.getsource(cohorts.sample),
        "bulk_explicit_preparation_probabilities": "< 0.04"
        in inspect.getsource(observations.measure)
        and "< 0.16" in inspect.getsource(observations.measure),
        "overflow_is_exception_not_native_record": "raise OverflowError"
        in inspect.getsource(clock.Registers.post),
    }
    require(all(checks.values()), "normal entrypoint contract changed")
    return dict(
        source_sha256=hashes,
        checks=checks,
        quantum="native quantum state, unitary/Kraus update and squared-amplitude draw are primitive architecture semantics",
        classical="independent uniform random streams and declared preparation probabilities are supplied; fixed PRNG seeds replay draws",
        bounded_rejection="all failed address draws contribute to explicit timeout or release; no success-only renormalization",
        abort="invalid allocation, ownership, unprepared read, numeric invariant and exhausted elapsed register raise host exceptions; no physical observer refusal record is claimed",
        finite_machine="ideal real probability semantics; a 53-bit uniform threshold differs by at most2^-53 per scalar draw; a finite tree with N draws has coupling error at most N*2^-53 before propagation error",
    )


def report():
    base = compile_program(native_program(), 99)
    refined = compile_program(native_program(refined=True), 99)
    post = compile_program(native_program(post=True), 99)
    rng = np.random.default_rng(8303)
    # Arbitrary complex full-support initial law, including inter-charge blocks.
    vectors = rng.normal(size=(99, 7)) + 1j * rng.normal(size=(99, 7))
    rho = vectors @ vectors.conj().T
    rho /= np.trace(rho)
    runtime = runtime_law(base, rho)
    material_labels = [f"material:{a}{b}" for a, b in law.SRC]
    signature_error = 0.0
    refinement_error = 0.0
    post_error = 0.0
    pairing_error = 0.0
    for label in material_labels:

        def event(history, label=label):
            return label in history

        effect = event_effect(base, event)
        other = event_effect(refined, event)
        twisted = event_effect(post, event)
        signature_error = max(
            signature_error,
            float(np.max(np.abs(reconstruct(probes(effect), 99) - effect))),
        )
        refinement_error = max(refinement_error, float(np.max(np.abs(effect - other))))
        post_error = max(post_error, float(np.max(np.abs(effect - twisted))))
        value = sum(p for history, p in runtime.items() if event(history))
        pairing_error = max(
            pairing_error, abs(value - float(np.trace(rho @ effect).real))
        )
    gram = event_effect(base, lambda _: True)
    require(
        max(signature_error, refinement_error, post_error, pairing_error) < 1e-11,
        "event bridge failed",
    )
    return dict(
        source_audit=source_audit(),
        source_initial_law=source_initial_law(),
        producer_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        domain="every finite complete native/adaptive branch program on the finite registered carrier; host exceptions outside valid execution domains are separately reported",
        dirty_collar_dimension=99,
        initial_density_rank=int(np.linalg.matrix_rank(rho)),
        initial_density_min_eigenvalue=float(np.linalg.eigvalsh(rho).min()),
        initial_density_trace=float(np.trace(rho).real),
        initial_charge_coherence=float(
            np.max(
                np.abs(
                    rho
                    * (
                        np.tile([sum(x) for x in law.SRC], 11)[:, None]
                        - np.tile([sum(x) for x in law.SRC], 11)[None, :]
                    )
                )
            )
        ),
        complete_raw_records={str(k): v for k, v in runtime.items()},
        base_histories=len(base),
        refined_histories=len(refined),
        complete_effect_error=float(np.max(np.abs(gram - np.eye(99)))),
        state_independent_probe_coordinates=99**2,
        probe_reconstruction_error=signature_error,
        refined_event_operator_error=refinement_error,
        postrecord_effect_error=post_error,
        runtime_effect_pairing_error=pairing_error,
        classical_rejection_q3=bounded_address_reader(3),
        counterexample=exact_counterexample(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    actual = report()
    if args.write:
        with RECEIPT.open("x") as f:
            json.dump(actual, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        compare(json.loads(RECEIPT.read_text()), actual)
        print(
            "public event grammar: complete native histories, refinements, independent signatures and sampler boundary replay"
        )


if __name__ == "__main__":
    main()
