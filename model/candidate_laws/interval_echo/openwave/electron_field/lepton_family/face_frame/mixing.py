"""Mixing from one face per mass eigenbasis and an identity family current."""

import json
from collections import Counter
from itertools import combinations, product

import numpy as np
import sympy as sp

from .measure import character_frames, faces, theorem_projectors


def exact_moduli():
    """Enumerate overlaps from the face normals, vertices of the dual dodecahedron."""
    phi = (1 + sp.sqrt(5)) / 2
    normals = list(product((-1, 1), repeat=3))
    for a, b in product((-1, 1), repeat=2):
        normal = (0, a / phi, b * phi)
        normals.extend(normal[i:] + normal[:i] for i in range(3))

    # All normals have length sqrt(3). Transitivity allows (1,1,1) as first axis.
    cosines = Counter(sp.simplify(sp.Abs(sp.Add(*n)) / 3) for n in normals)
    assert cosines == {sp.Integer(1): 2, sp.Rational(1, 3): 12, sp.sqrt(5) / 3: 6}
    moduli = set()
    for c in cosines:
        entries = (c, sp.sqrt((1 - c**2) / 2), (1 + c) / 2, (1 - c) / 2)
        moduli.update(sp.simplify(v) for v in entries)
    return sorted(moduli, key=float)


def quartets(mixing):
    for i, j in combinations(range(3), 2):
        for k, l in combinations(range(3), 2):
            direct = mixing[i, k] * mixing[j, l]
            crossed = mixing[i, l] * mixing[j, k]
            yield direct * crossed.conjugate()


def run():
    bands, _ = theorem_projectors()
    bases = []
    for face in faces():
        columns = []
        for projector in character_frames(face, bands[1][0]):
            vector = projector[:, np.argmax(projector.diagonal().real)]
            columns.append(vector / np.linalg.norm(vector))
        basis = np.column_stack(columns)
        np.testing.assert_allclose(
            basis.conj().T @ basis, np.eye(3), atol=1e-13, rtol=0
        )
        bases.append(basis)

    moduli = exact_moduli()
    expected = np.array(moduli, dtype=float)
    hits = np.zeros(len(moduli), dtype=int)
    imaginary = []
    for up, down in product(bases, repeat=2):
        mixing = up.conj().T @ down
        magnitudes = abs(mixing)
        # In character order 0,+,- the circular diagonal and off-diagonal
        # entries respectively have magnitudes (1+c)/2 and (1-c)/2.
        cosine = magnitudes[1, 1] - magnitudes[1, 2]
        diagonal, opposite = (1 + cosine) / 2, (1 - cosine) / 2
        mixed = np.sqrt(max(0.0, (1 - cosine**2) / 2))
        # Square before comparing so that roundoff at parallel axes stays O(eps).
        predicted = np.array([
            [abs(cosine), mixed, mixed],
            [mixed, diagonal, opposite],
            [mixed, opposite, diagonal],
        ])
        np.testing.assert_allclose(magnitudes**2, predicted**2, atol=1e-13, rtol=0)
        for value in magnitudes.flat:
            nearest = np.argmin(abs(expected - value))
            assert abs(expected[nearest] - value) < 1e-13
            hits[nearest] += 1
        imaginary.extend(abs(q.imag) for q in quartets(mixing))
    assert np.all(hits > 0)
    assert max(imaginary) < 1e-13
    return dict(
        face_pairs=len(bases)**2,
        quartets=len(imaginary),
        maximum_absolute_imaginary_quartet=float(max(imaginary)),
        entry_moduli=[str(v) for v in moduli],
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
