"""A coherent three-face phase with unit coefficients, using the registered band.

Run this module from the repository with python -m and its package name.
The first three outward faces in the apparatus order select one fixed example.
The positive helicity and the two operators below are design choices. No
quark mass or mixing measurement enters. A physical CKM interpretation needs
independent production, mass and charged-current evidence on the same sector.

For V=[a b c], G=V†V and rank-one P_i=v_i v_i†, let A=I+P_a+P_b and
B=I+P_b+P_c. Their commutator is V K V† where K_ij=(u_i d_j-u_j d_i)G_ij,
u=(1,1,0), d=(0,1,1). K has zero diagonal and is skew-Hermitian. Expanding
its determinant gives -2i Im(G_ab G_bc G_ca). The traceless 3x3 commutator
has tr([A,B]^3)=3 det([A,B]). Thus a nonzero directed three-projector phase
and linear independence give a CP-odd invariant of this particular A,B.
"""

import json

import numpy as np

from .measure import (
    character_frames,
    faces,
    theorem_projectors,
)


def witness(vectors):
    projectors = [np.outer(vector, vector.conj()) for vector in vectors.T]
    a, b, c = projectors
    loop = np.trace(a @ b @ c)
    up = np.eye(3) + a + b
    down = np.eye(3) + b + c
    commutator = up @ down - down @ up
    invariant = np.trace(np.linalg.matrix_power(commutator, 3)).imag
    gram = vectors.conj().T @ vectors
    expected = -6 * np.linalg.det(gram).real * loop.imag
    np.testing.assert_allclose(invariant, expected, atol=2e-14, rtol=0)

    mass_up, frame_up = np.linalg.eigh(up)
    mass_down, frame_down = np.linalg.eigh(down)
    mixing = frame_up.conj().T @ frame_down
    j = (mixing[0, 0] * mixing[1, 1] * (mixing[0, 1] * mixing[1, 0]).conj()).imag

    # A projected path through c, b and a has amplitude loop for input a.
    # The idle path returns to a with unit amplitude. After the two paths
    # interfere, report the two Y outcomes and every declined projection.
    yes = abs(1 - 1j * loop)**2 / 4
    no = abs(1 + 1j * loop)**2 / 4
    declined = (1 - abs(loop)**2) / 2
    np.testing.assert_allclose([yes + no + declined, yes - no], [1, loop.imag], atol=2e-14)
    return dict(
        loop=[loop.real, loop.imag], invariant=invariant, quartet=j,
        squared_mass_spectra=[mass_up.tolist(), mass_down.tolist()],
        Y_plus_minus_refused=[yes, no, declined],
    )


def main():
    bands, _ = theorem_projectors()
    band = bands[1][0]
    coordinates = np.linalg.eigh(band)[1][:, -3:]
    enrolled_faces = list(faces())[:3]
    vectors = []
    for face in enrolled_faces:
        projector = character_frames(face, band)[1]
        vector = projector[:, np.argmax(projector.diagonal().real)]
        vectors.append(coordinates.conj().T @ vector / np.linalg.norm(vector))
    vectors = np.column_stack(vectors)
    actual = witness(vectors)
    conjugate = witness(vectors.conj())
    repeated = witness(vectors[:, [0, 1, 0]])
    for spectrum in actual["squared_mass_spectra"]:
        assert min(spectrum) > 0 and min(np.diff(spectrum)) > 0.1
    assert actual["invariant"] > 0.008 and actual["quartet"] > 0.004
    np.testing.assert_allclose(conjugate["invariant"], -actual["invariant"], atol=2e-14)
    np.testing.assert_allclose(repeated["invariant"], 0, atol=2e-14)
    print(json.dumps(dict(faces=enrolled_faces, comparison=actual,
                          conjugate_invariant=conjugate["invariant"],
                          repeated_face_invariant=repeated["invariant"]), indent=2))


if __name__ == "__main__":
    main()
