"""Exact Gaussian-rational replay of the modal Maxwell-shaped factorization.

The Lean theorem inserts the OPH FZ-12 frequency.  This independent companion
replays the vector algebra with exact pairs of ``Fraction`` values: the
Fourier multiplier ``i * scale * (k x)`` with opposite E/B signs squares to
the oscillator wave law, whereas same signs or a wrong normalization do not.
It is not a simulation and makes no physical electromagnetic identification.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Final, TypeAlias

Real: TypeAlias = Fraction
Gaussian: TypeAlias = tuple[Real, Real]
Momentum: TypeAlias = tuple[Real, Real, Real]
ComplexVec3: TypeAlias = tuple[Gaussian, Gaussian, Gaussian]
ComplexPairedState: TypeAlias = tuple[ComplexVec3, ComplexVec3]

GZERO: Final[Gaussian] = (Fraction(0), Fraction(0))
VZERO: Final[ComplexVec3] = (GZERO, GZERO, GZERO)


def gadd(x: Gaussian, y: Gaussian) -> Gaussian:
    return x[0] + y[0], x[1] + y[1]


def gneg(x: Gaussian) -> Gaussian:
    return -x[0], -x[1]


def gmul(x: Gaussian, y: Gaussian) -> Gaussian:
    return x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def vadd(u: ComplexVec3, v: ComplexVec3) -> ComplexVec3:
    return tuple(gadd(x, y) for x, y in zip(u, v, strict=True))  # type: ignore[return-value]


def vneg(v: ComplexVec3) -> ComplexVec3:
    return tuple(gneg(x) for x in v)  # type: ignore[return-value]


def vscale(c: Gaussian, v: ComplexVec3) -> ComplexVec3:
    return tuple(gmul(c, x) for x in v)  # type: ignore[return-value]


def real_dot(u: Momentum, v: Momentum) -> Real:
    return sum((x * y for x, y in zip(u, v, strict=True)), start=Fraction(0))


def momentum_dot(momentum: Momentum, v: ComplexVec3) -> Gaussian:
    result = GZERO
    for k, component in zip(momentum, v, strict=True):
        result = gadd(result, gmul((k, Fraction(0)), component))
    return result


def momentum_cross(momentum: Momentum, v: ComplexVec3) -> ComplexVec3:
    k0, k1, k2 = ((x, Fraction(0)) for x in momentum)
    return (
        gadd(gmul(k1, v[2]), gneg(gmul(k2, v[1]))),
        gadd(gmul(k2, v[0]), gneg(gmul(k0, v[2]))),
        gadd(gmul(k0, v[1]), gneg(gmul(k1, v[0]))),
    )


def fourier_curl(
    curl_scale: Real, momentum: Momentum, v: ComplexVec3
) -> ComplexVec3:
    """Exact multiplier ``i * curl_scale * (momentum x v)``."""

    return vscale((Fraction(0), curl_scale), momentum_cross(momentum, v))


def spatial_action(frequency: Real, v: ComplexVec3) -> ComplexVec3:
    return vscale((frequency * frequency, Fraction(0)), v)


def maxwell_shaped_generator(
    curl_scale: Real,
    momentum: Momentum,
    state: ComplexPairedState,
) -> ComplexPairedState:
    """Opposite-sign pair ``(D B, -D E)`` for ``D=i*scale*(k x)``."""

    electric, magnetic = state
    return (
        fourier_curl(curl_scale, momentum, magnetic),
        vneg(fourier_curl(curl_scale, momentum, electric)),
    )


def same_sign_mutation(
    curl_scale: Real,
    momentum: Momentum,
    state: ComplexPairedState,
) -> ComplexPairedState:
    """Adversarial pair ``(D B, D E)`` with the wrong relative sign."""

    electric, magnetic = state
    return (
        fourier_curl(curl_scale, momentum, magnetic),
        fourier_curl(curl_scale, momentum, electric),
    )


def second_order_residual(
    generator,
    curl_scale: Real,
    momentum: Momentum,
    frequency: Real,
    state: ComplexPairedState,
) -> ComplexPairedState:
    """Return ``G(G(state)) + omega^2 state`` componentwise."""

    second = generator(curl_scale, momentum, generator(curl_scale, momentum, state))
    return (
        vadd(second[0], spatial_action(frequency, state[0])),
        vadd(second[1], spatial_action(frequency, state[1])),
    )


def _serialize_vec(v: ComplexVec3) -> list[list[str]]:
    return [[str(x[0]), str(x[1])] for x in v]


def exact_replay() -> dict[str, object]:
    """Run an off-axis witness and two independent mutation controls."""

    momentum: Momentum = (Fraction(3), Fraction(4), Fraction(0))
    frequency = Fraction(7)
    curl_scale = Fraction(7, 5)
    electric: ComplexVec3 = (
        (Fraction(4), Fraction(0)),
        (Fraction(-3), Fraction(0)),
        (Fraction(2), Fraction(1)),
    )
    magnetic: ComplexVec3 = (
        (Fraction(0), Fraction(4)),
        (Fraction(0), Fraction(-3)),
        (Fraction(5), Fraction(0)),
    )
    state = (electric, magnetic)

    assert curl_scale * curl_scale * real_dot(momentum, momentum) == frequency**2
    assert momentum_dot(momentum, electric) == GZERO
    assert momentum_dot(momentum, magnetic) == GZERO
    assert momentum_dot(momentum, fourier_curl(curl_scale, momentum, electric)) == GZERO
    assert momentum_dot(momentum, fourier_curl(curl_scale, momentum, magnetic)) == GZERO

    wave_residual = second_order_residual(
        maxwell_shaped_generator, curl_scale, momentum, frequency, state
    )
    assert wave_residual == (VZERO, VZERO)

    same_sign_residual = second_order_residual(
        same_sign_mutation, curl_scale, momentum, frequency, state
    )
    expected_same_sign_residual = (
        vscale((2 * frequency * frequency, Fraction(0)), state[0]),
        vscale((2 * frequency * frequency, Fraction(0)), state[1]),
    )
    assert same_sign_residual == expected_same_sign_residual

    wrong_scale_residual = second_order_residual(
        maxwell_shaped_generator, Fraction(2), momentum, frequency, state
    )
    assert wrong_scale_residual != (VZERO, VZERO)

    return {
        "arithmetic": "exact_gaussian_fraction",
        "momentum": [str(x) for x in momentum],
        "frequency": str(frequency),
        "curl_scale": str(curl_scale),
        "maxwell_shaped_wave_residual": [
            _serialize_vec(wave_residual[0]),
            _serialize_vec(wave_residual[1]),
        ],
        "same_sign_mutation_rejected": True,
        "same_sign_residual_exactly_twice_spatial_action": True,
        "wrong_scale_mutation_rejected": True,
        "computes_fz12_symbol": False,
        "physical_claim": False,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(exact_replay(), indent=2, sort_keys=True))
