"""Exact native column checks and the invisible ordered correlation on two private seams."""
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_einstein_station import backreaction as seam

HERE = Path(__file__).resolve().parent
ORIENTED = np.array([[0, 1, -1], [-1, 0, 1], [1, -1, 0]], dtype=int)


def exact_columns():
    law = seam.law
    w = law._write_pattern()  # sqrt(2) W
    reflection, d2 = law._label_ops()
    dim = law.N_SRC*law.N_LAB
    twice_projection = 2*np.eye(dim, dtype=int) - w@w.T
    twice_generator = np.kron(np.eye(law.N_LAB, dtype=int) - reflection + int(2*seam.KAPPA)*d2,
                              np.eye(law.N_SRC, dtype=int))
    legs = np.eye(law.N_SRC, dtype=int)
    vacuum = legs[law.SRC.index((0, 0))]
    sym = legs[law.SRC.index((0, 1))] + legs[law.SRC.index((1, 0))]
    for vector, energy in [(vacuum, 0), (sym, 8)]:
        intermediate = twice_generator @ w @ vector
        np.testing.assert_array_equal(w.T @ intermediate, 4*energy*vector)
        np.testing.assert_array_equal(twice_projection @ intermediate, np.zeros(dim, dtype=int))

    def eigenword(magnitude, parity, phase, load, energy):
        word = seam.coin_basis(magnitude, parity) @ np.array([1, phase])
        vector = np.kron(word, load)
        intermediate = twice_generator @ twice_projection @ vector
        blank = w.T @ intermediate
        occupied_residual = twice_projection @ intermediate - 8*energy*vector
        np.testing.assert_array_equal(blank, np.zeros(law.N_SRC, dtype=int))
        np.testing.assert_array_equal(occupied_residual, np.zeros(dim, dtype=int))
        return max(int(np.max(np.abs(blank))), int(np.max(np.abs(occupied_residual))))

    protected_modes = [(0, 1, -1, 0), (0, 1, 1, 0), (2, -1, -1, 33), (2, 1, -1, 32)]
    words = np.column_stack([seam.coin_basis(m, parity) @ np.array([1, phase])
                             for m, parity, phase, _ in protected_modes])
    gram = words.T @ words
    np.testing.assert_array_equal(gram, np.diag([2, 2, 4, 4]))
    residual = max(eigenword(m, parity, phase, load, e)
                   for m, parity, phase, e in protected_modes for load in legs)
    for magnitude, energies in [(0, (0, 0)), (2, (32, 0))]:
        for ports, energy in zip([(2, 0), (1, 1)], energies, strict=True):
            eigenword(magnitude, 1, 1, legs[law.SRC.index(ports)], energy)
    return residual


def distribution(parameter):
    return (1 + parameter*ORIENTED)/9


def audit():
    scalar_residual = exact_columns()
    energies = (0, 32, 33)
    # Two logical registers in one protected memory. The second bit distinguishes
    # coin phase at magnitude zero and reversal parity at magnitude two.
    table = np.array([[0, 0], [33, 32]])
    difference = np.array([[-1, 1], [1, -1]], dtype=int)
    energy_contrast = Fraction(int((table*difference).sum()), 4)
    assert energy_contrast == Fraction(1, 4)
    sigma = np.exp(-table/16)
    sigma /= sigma.sum()
    modular = np.log(np.outer(sigma.sum(axis=1), sigma.sum(axis=0))/sigma)
    measured_energy = float(4*np.sum(modular*difference))
    np.testing.assert_allclose(measured_energy, float(energy_contrast), atol=2e-14)
    blind = dict(logical_energy_table=table.tolist(), energy_contrast=str(energy_contrast),
                 gibbs_energy_contrast=measured_energy, max_integer_scalar_action_residual=scalar_residual,
                 protected_input_columns=4*seam.law.N_SRC)
    samples = []
    for t in [Fraction(-1, 2), Fraction(0), Fraction(1, 2)]:
        probs = distribution(t)
        assert list(probs.sum(axis=0)) == list(probs.sum(axis=1)) == [Fraction(1,3)]*3
        histogram = {}
        for i, u in enumerate(energies):
            for j, v in enumerate(energies):
                histogram[u+v] = histogram.get(u+v, 0) + probs[i,j]
        samples.append(histogram)
    assert samples[0] == samples[1] == samples[2]
    rows = []
    t0 = 0.5
    for radius in (0.5, 0.25, 0.125, 0.0625):
        f = radius**2
        reference = distribution(f*t0)
        tangent = f*ORIENTED/9
        modular = float(np.sum(-np.log(reference*9)*tangent))
        analytic = -2*f*np.arctanh(f*t0)/3
        np.testing.assert_allclose(modular, analytic, atol=1e-16, rtol=0)
        rows.append(dict(radius=radius, inclusion=f, modular_slope=modular,
                         rescaled_stress_slope=15*modular/(8*np.pi**2*radius**4)))
    return dict(protected_logical_registers=blind, arm_parameters=[-0.5,0,0.5], exact_source_energies=list(energies),
                common_energy_histogram={str(e):str(w) for e,w in sorted(samples[0].items())},
                full_inclusion_modular_slope=-np.log(3)/3,
                rescaled_limit=-5/(8*np.pi**2), regional=rows)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    rows = audit()
    dest = HERE/'protected.json'
    if parser.parse_args().write:
        dest.write_text(json.dumps(rows, indent=2, allow_nan=False)+'\n')
    else:
        from cosmology.rar_deep_regime.lensing_profile import check
        check(json.loads(dest.read_text()), rows)
    print(rows['full_inclusion_modular_slope'], rows['rescaled_limit'])
