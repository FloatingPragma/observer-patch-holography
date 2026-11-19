"""Local charge hopping and the radiative pair converter on the OPH carrier."""

import math

import numpy as np
from scipy import sparse as sp

from ..radiative_repair import source


MASS = 1
EXTERNAL_WEIGHT = 1 / 8


def charged_operator(car, external_weight=EXTERNAL_WEIGHT, *, parked=False):
    """Use difference squares and one unit of onsite charged energy.

    The parked control keeps the diagonal costs of the external terms.
    """
    if not math.isfinite(external_weight) or external_weight < 0:
        raise ValueError("external seam weight must be finite and nonnegative")
    ends = np.asarray(car.seams, dtype=int).reshape(-1, 2)
    left, right = ends.T
    external = left // 12 != right // 12
    weights = np.where(external, external_weight, 1.0)
    degree = np.bincount(ends.ravel(), weights=np.repeat(weights, 2), minlength=car.n)
    hops = weights * (~external if parked else 1)
    off_diagonal = sp.coo_matrix(
        (-np.tile(hops, 2), (np.r_[left, right], np.r_[right, left])),
        shape=(car.n, car.n),
    )
    return (sp.diags(MASS + degree) + off_diagonal).tocsr()


def pair_generator(car, external_weight=EXTERNAL_WEIGHT, *, parked=False, blocked_conversion=False):
    """Use the same converter and neutral field as radiative_repair at K=2."""
    original_charged, _ = source.one_body(car)
    charged_delta = charged_operator(car, external_weight, parked=parked) - original_charged
    pair_delta = sp.kronsum(charged_delta, charged_delta)
    old_law = "dephased_contact" if blocked_conversion else "exchange"
    base = source.generator(car, old_law)
    return base + sp.block_diag((pair_delta, sp.csr_matrix(pair_delta.shape)), format="csr")
