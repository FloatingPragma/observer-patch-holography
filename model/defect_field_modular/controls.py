"""Frozen emitted source words; every replay checks their realized target.

These are compiler outputs for declared reference data, not new pulse kinds.
Keeping their actual angles fixes resources and avoids platform-dependent
phase choices during QR compilation of an equivalent isometry.
"""
from pathlib import Path
from functools import lru_cache

import numpy as np

from ..observer.source_scalar import compiler

PATH=Path(__file__).with_name('source_controls.npz')
KINDS=('I','X','Y','Z','H','J')


def encode(word):
    return np.array([[KINDS.index(e['kind']),e['angle'],e['targets'][0],
                      e['targets'][1] if len(e['targets'])==2 else -1] for e in word],float)


def decode(array):
    a=np.asarray(array)
    if a.ndim!=2 or a.shape[1]!=4 or not np.isfinite(a).all():
        raise ValueError('complete finite frozen source-word rows required')
    result=[]
    for kind,angle,first,second in a:
        if kind!=int(kind) or not 0<=kind<len(KINDS) or first!=int(first) or second!=int(second):
            raise ValueError('exact pulse kinds and target indices required')
        targets=[int(first)] if second==-1 else [int(first),int(second)]
        result.append(compiler.pulse(KINDS[int(kind)],float(angle),targets))
    return result


@lru_cache(maxsize=1)
def load():
    with np.load(PATH,allow_pickle=False) as data:
        if set(data.files)!={'preparation','recovery'}:
            raise ValueError('both complete frozen source words are required')
        return {key:decode(data[key]) for key in data.files}


def freeze():
    from .modular import reference
    from ..common_source.mass_response.history import preparation_word
    from ..field_thermal.thermal import car_restriction
    from ..reference_recovery.recovery import RecoveryService
    data=reference()
    preparation=preparation_word(data['root'].T.reshape(-1))
    recovery=RecoveryService(data['rho'],car_restriction(4,(0,1)))
    np.savez_compressed(PATH,preparation=encode(preparation),recovery=encode(recovery.word))
    load.cache_clear()
