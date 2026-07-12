"""Hadron two-point contractions from point propagators.

Channels follow the reference-spec export granularity:

    pi_iso          pseudoscalar pi correlator
    N_iso_direct    nucleon direct Wick contraction
    N_iso_exchange  nucleon exchange Wick contraction
    N_iso = N_iso_direct - N_iso_exchange   (repo-side combination rule)

The nucleon interpolator is eps_abc (u_a^T C gamma5 d_b) u_c with degenerate
u = d propagators, projected with P_plus = (1 + gamma_4)/2. The contraction
is the direct Wick enumeration of the two u-quark pairings; the exchange
pairing carries the odd-permutation sign, realized by the combination rule.
The source-side barred spin matrix equals -C gamma5, a global phase common
to both terms, so it is dropped. The free-field test anchors the result:
the projected N_iso effective mass approaches three times the free quark
pole mass.
"""

from __future__ import annotations

import numpy as np

from .core import CCONJ, GAMMA, GAMMA5

CG5 = CCONJ @ GAMMA5
P_PLUS = 0.5 * (np.eye(4, dtype=complex) + GAMMA[3])

_EPS = np.zeros((3, 3, 3))
for _i, _j, _k, _s in (
    (0, 1, 2, 1.0), (1, 2, 0, 1.0), (2, 0, 1, 1.0),
    (0, 2, 1, -1.0), (2, 1, 0, -1.0), (1, 0, 2, -1.0),
):
    _EPS[_i, _j, _k] = _s


def pion_correlator(prop: np.ndarray) -> np.ndarray:
    """C_pi(t) = sum_x Tr[S(x) S(x)^dag], the pi^+ pseudoscalar correlator
    via gamma5-hermiticity. Positive by construction."""
    dens = np.einsum("txyzscud,txyzscud->t", prop, np.conj(prop))
    return np.real(dens)


def _quark_lines(prop_t: np.ndarray) -> np.ndarray:
    """Reorder one timeslice to Q[x, c_sink, c_src, s_sink, s_src]."""
    v3 = prop_t.shape[0] * prop_t.shape[1] * prop_t.shape[2]
    q = prop_t.reshape(v3, 4, 3, 4, 3)
    return q.transpose(0, 2, 4, 1, 3)


def nucleon_correlators(prop: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(direct, exchange) nucleon two-point terms, P_plus projected.

    With A = C gamma5 and Q the quark line, the two u-quark pairings give

      direct   = eps_abc eps_ijk A_{ef} A_{gh} (P_plus)_{st}
                 Q_{ai,eg} Q_{bj,fh} Q_{ck,ts}
      exchange = eps_abc eps_ijk A_{ef} A_{gh} (P_plus)_{st}
                 Q_{ak,es} Q_{bj,fh} Q_{ci,tg}
    """
    t_extent = prop.shape[0]
    direct = np.zeros(t_extent)
    exchange = np.zeros(t_extent)
    for t in range(t_extent):
        q = _quark_lines(prop[t])
        d_t = np.einsum(
            "abc,ijk,ef,gh,st,xaieg,xbjfh,xckts->",
            _EPS, _EPS, CG5, CG5, P_PLUS, q, q, q, optimize=True)
        e_t = np.einsum(
            "abc,ijk,ef,gh,st,xakes,xbjfh,xcitg->",
            _EPS, _EPS, CG5, CG5, P_PLUS, q, q, q, optimize=True)
        direct[t] = np.real(d_t)
        exchange[t] = np.real(e_t)
    return direct, exchange


def effective_mass(corr: np.ndarray) -> np.ndarray:
    """am_eff(t) = log(|C(t)| / |C(t+1)|) where both entries are nonzero."""
    out = np.full(len(corr) - 1, np.nan)
    for t in range(len(corr) - 1):
        if corr[t] != 0.0 and corr[t + 1] != 0.0 and abs(corr[t]) > abs(corr[t + 1]):
            out[t] = np.log(abs(corr[t]) / abs(corr[t + 1]))
    return out
