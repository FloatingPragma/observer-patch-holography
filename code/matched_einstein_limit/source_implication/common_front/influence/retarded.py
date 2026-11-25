"""Exact spatial support of the head self-energy, with a nominal finite-frequency replay."""

import argparse
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from . import measure
from ...directional_stress.bridge import tangents


def sites(rows):
    return (4 * np.asarray(rows)[:, None] + np.arange(4)).ravel()


def expect(matrix, state):
    return np.einsum("bijkm,bmk->ij", matrix, state)


def tensor(p):
    n = p.shape[0] // 4
    return p.reshape(n, 4, n, 4).transpose(0, 2, 1, 3)


def compute(q, board, route, states):
    native = measure.action.native
    # Subtract the energy of a localized pilot at blank port 5, in each charge
    # experiment. This is unchanged by either source or any frequency setting.
    hs = [native.generator(board, route, q, o) / 2 for o in (False, True)]
    blank = sites([5])
    energy = hs[0][blank, :][:, blank].toarray()[0, 0]
    assert all(
        np.array_equal(h[blank, :][:, blank].toarray(), energy * np.eye(4)) for h in hs
    )
    hs = [h - energy * sparse.eye(h.shape[0], format="csr") for h in hs]
    outer = np.array([p for p in range(board.n) if p not in route])
    ext = sites(outer)
    mark = sites(route)
    full = sparse.eye(hs[0].shape[0], format="csc")
    h0 = hs[0][ext, :][:, ext][::4, ::4].tocsr()
    contact = hs[0][mark, :][:, ext][::4, ::4].tocsr()
    eye = sparse.eye(4, format="csr")
    for h in hs:
        assert (h[ext, :][:, ext] - sparse.kron(h0, eye)).count_nonzero() == 0
        assert (h[mark, :][:, ext] - sparse.kron(contact, eye)).count_nonzero() == 0
    hit = np.array([h[mark, :][:, mark].toarray() for h in hs])
    hbar = {k: expect(np.array([tensor(v) for v in hit]), s) for k, s in states.items()}
    exterior_read_port = 5
    remote = np.flatnonzero(outer == exterior_read_port)[0]
    pole = {}
    for omega, eta in ((0, 0.25), (0, 1), (0, 4), (-1, 1), (1, 1)):
        z = complex(omega, eta)
        operator = z * sparse.eye(len(outer), format="csc") - h0
        exterior = splu(operator.tocsc()).solve(contact.T.toarray())
        surface = contact @ exterior
        rn = np.array(
            [
                np.linalg.inv(z * np.eye(16) - h - np.kron(surface, np.eye(4)))
                for h in hit
            ]
        )
        kernels = np.array([tensor(x) for x in rn])
        rr = {k: expect(kernels, s) for k, s in states.items()}
        inv = np.linalg.inv(rr["reference"])
        sigma = z * np.eye(4) - surface - inv - hbar["reference"]
        # Scalar exterior propagation gives delta G_BB through the same
        # four-by-four reduced block, even for arbitrarily distant terminals.
        v = exterior[remote]
        panel = {}
        for k in ("ordinary", "pairing"):
            delta = inv @ rr[k] @ inv - hbar[k]
            amplitude = sigma[3, 0]
            slope = (amplitude.conj() * delta[3, 0]).real / abs(amplitude)
            panel[k] = dict(
                connected_39_from_0=[delta[3, 0].real, delta[3, 0].imag],
                modulus_slope=slope,
                blank_return_slope=[(v @ rr[k] @ v).real, (v @ rr[k] @ v).imag],
                full_marked_matrix_frobenius_slope=np.linalg.norm(delta),
            )
        rhs = np.zeros((hs[0].shape[0], 20))
        rhs[np.r_[mark, sites([exterior_read_port])]] = np.eye(20)
        direct = []
        err = 0
        for h, k in zip(hs, rn):
            a = (z * full - h).tocsc()
            y = splu(a).solve(rhs)
            err = max(
                err, np.linalg.norm(a @ y - rhs), np.linalg.norm(y[mark, :16] - k)
            )
            direct.append(y[sites([exterior_read_port]), -4:])
        for label in ("ordinary", "pairing"):
            exact = sum(np.trace(u @ s) for u, s in zip(direct, states[label]))
            go = v @ rr[label] @ v
            assert abs(go - exact) < 1e-11
        assert err < 1e-10
        denom = panel["ordinary"]["modulus_slope"]
        pole[f"{omega},{eta}"] = dict(
            omega=omega,
            width=eta,
            reference_memory_39_from_0=[sigma[3, 0].real, sigma[3, 0].imag],
            ordinary=panel["ordinary"],
            pairing=panel["pairing"],
            modulus_collar_over_ordinary=(
                panel["pairing"]["modulus_slope"] / denom if denom else None
            ),
            maximum_sparse_residual_and_schur_difference=err,
        )
    return dict(
        blank_reference_port=5, exterior_receiver=exterior_read_port, settings=pole
    )


def run():
    n = measure.action.native
    b = n.iso.pairing.formation.federation(1)
    route = n.iso.pairing.formation.region(b)
    return dict(
        exterior_ports=b.n - len(route),
        marked_ports=list(route),
        charges={str(q): compute(q, b, route, tangents()) for q in (-1, 1)},
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "retarded.json").write_text(out)
    else:
        print(out, end="")
