"""Source maps from full seam and collar-fuel instruments to the used columns."""

from functools import lru_cache
from math import pi, sin, hypot

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import collar
from model.candidate_laws.interval_echo import candidate
from ..photon_identity.native_emission.measure import four_h


def require(check, reason):
    if not check:
        raise ArithmeticError(reason)


def fuel_probability(k, *, matched):
    # On the two vacuum words the diagonal difference is zero when the
    # occupied fuel rail has its native reference. Bare fuel has level zero.
    if matched:
        return .5
    ratio = hypot(1, k/(2*collar.G))
    return (sin(pi*ratio/4)/ratio)**2


@lru_cache(maxsize=1)
def check_local():
    law = candidate.native_law
    n, m = law.N_SRC, law.DIM
    eye = np.eye(m)
    words = ((0,0), (1,0), (0,1))

    def ix(c, d, word):
        return (1+law.LABELS.index((c,d)))*n + law.SRC.index(word)

    # Integer columns for sqrt(2) times each normalized arrow state.
    v = np.column_stack([eye[:, ix(c, -1, w)]+eye[:, ix(c, 1, w)] for c in (0,1) for w in words]).astype(int)
    active = np.zeros(6, dtype=int)
    active[[2,4]] = 1
    twice_a = np.outer(active, active)
    h0, h1 = four_h(0), four_h(1)-four_h(0)
    np.testing.assert_array_equal(h0 @ v, 0*v)
    np.testing.assert_array_equal(h1 @ v, v @ (4*np.eye(6, dtype=int)-2*twice_a))
    # On written inputs 2 U P_blank U is (sqrt(2) W)(sqrt(2) W)^T.
    w = law._write_pattern()
    np.testing.assert_array_equal((w @ w.T) @ v[n:], v[n:] @ twice_a)
    flag = np.arange(m)<n
    u = law.echo()
    for seen in (False,True):
        actual = u.T @ ((flag == seen)[:,None]*(u @ v))
        np.testing.assert_allclose(actual, v @ (twice_a/2 if seen else np.eye(6)-twice_a/2),
                                   atol=2e-14, rtol=0)
    # Copy just the coin bit, never the unresolved signed contrast. After
    # either blank-query branch the two head outcomes on a sharp active input
    # are equiprobable sharp words, with the even sign superposition intact.
    coin_one = np.repeat([False]+[c==1 for c,_ in law.LABELS],n)
    for index in (2,4):
        for seen in (False,True):
            after = u.T @ ((flag==seen)*(u @ v[:,index]))/np.sqrt(2)
            for coin,j in ((0,2),(1,4)):
                part = (coin_one==coin)*after
                np.testing.assert_allclose(abs(np.vdot(v[:,j]/np.sqrt(2),part))**2,.25,atol=2e-14,rtol=0)
                np.testing.assert_allclose(np.vdot(part,part).real,.25,atol=2e-14,rtol=0)
    # Physical endpoint read on a blank private reference with vacant far end.
    nonzero = np.repeat([False]+[d != 0 for _,d in law.LABELS], n)
    ref = eye[:, [law.SRC.index((0,0)), law.SRC.index((1,0))]]
    fuel_rail = ref.T @ h1 @ ref / 4
    np.testing.assert_array_equal(fuel_rail, np.diag([0,1]))
    np.testing.assert_array_equal(h1 @ ref, 4*ref @ fuel_rail)
    for seen in (0,1):
        np.testing.assert_allclose(u.T @ ((nonzero==seen)[:,None]*(u @ ref)),
                                   ref @ np.diag([int(seen==0),int(seen==1)]), atol=2e-14, rtol=0)
    # A sharp move along a blank contact is a SWAP at pi, times a common
    # phase for ONE sharp unit. No such phase is used in the five coherent queries.
    blank = eye[:, [law.SRC.index(w) for w in words[1:]]]
    np.testing.assert_array_equal(h0 @ blank, 2*blank @ np.array([[1,-1],[-1,1]]))
    np.testing.assert_array_equal(h1 @ blank, 4*blank)
    # Reflection of the complete native seam: exchange endpoint slots AND
    # reverse sign and coin. Mounting the donor converter at the other end
    # puts its even-sign word at the opposite physical head.
    perm = [law.SRC.index(w[::-1]) for w in law.SRC]
    perm += [ix(1-c,-d,w[::-1]) for c,d in law.LABELS for w in law.SRC]
    np.testing.assert_array_equal(np.sort(perm), np.arange(m))
    np.testing.assert_allclose(u[np.ix_(perm,perm)], u, atol=2e-14, rtol=0)
    for part in (h0,h1):
        np.testing.assert_array_equal(part[np.ix_(perm,perm)], part)
    trans = np.zeros((m,m))
    for w0 in law.SRC:
        for d in (-1,1):
            trans[ix(0,d,w0),law.SRC.index(w0)] = 1/np.sqrt(2)
    vacuum = eye[:,law.SRC.index((0,0))]
    for c, actor in enumerate((trans, trans[np.ix_(perm,perm)])):
        target = v[:,3*c]/np.sqrt(2)
        np.testing.assert_allclose(actor @ vacuum, target, atol=2e-14, rtol=0)
        np.testing.assert_allclose(actor.T @ target, vacuum, atol=2e-14, rtol=0)
        np.testing.assert_array_equal(h0 @ v[:,3*c], 0*v[:,0])
        np.testing.assert_array_equal(h1 @ v[:,3*c], 4*v[:,3*c])
        embedding = np.zeros((2*m,2))
        embedding[m:,0], embedding[:m,1] = vacuum, target
        for k in (0., .5, 8., -128., 128.):
            h = (h0+k*h1)/4
            for matched in (False,True):
                ef = k if matched else 0
                free_drift = np.kron(np.eye(2),h)
                if matched:
                    free_drift += k*np.kron(fuel_rail,eye)
                zero = np.zeros_like(h)
                full = free_drift+collar.G*np.block([[zero,actor],[actor.T,zero]])
                small = np.array([[ef,collar.G],[collar.G,k]])
                np.testing.assert_allclose(full @ embedding, embedding @ small, atol=1e-12, rtol=0)
                state = embedding @ expm((-1j*pi/(4*collar.G))*small)[:,0]
                np.testing.assert_allclose([np.vdot(state[:m],state[:m]).real,
                                            np.vdot(state[m:],state[m:]).real],
                                           [fuel_probability(k, matched=matched),1-fuel_probability(k, matched=matched)],
                                           atol=1e-14, rtol=0)
    return True
