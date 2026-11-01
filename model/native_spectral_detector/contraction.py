"""Auditable fused complex128 contraction; all rows and columns are retained."""
import ctypes
import atexit
from functools import lru_cache
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

_BUILD_DIRECTORY=None


@lru_cache(maxsize=1)
def library():
    global _BUILD_DIRECTORY
    source=Path(__file__).with_name('contract.cpp')
    compiler=shutil.which(os.environ.get('CXX','c++'))
    if not compiler:raise RuntimeError('a C++ compiler is required for the full-row native contraction')
    key=hashlib.sha256(source.read_bytes()+sys.platform.encode()+compiler.encode()).hexdigest()
    # A cold process compiles the currently bound source afresh; a stale
    # or externally replaced shared-library cache is never evidence.
    _BUILD_DIRECTORY=tempfile.TemporaryDirectory(prefix='oph-native-spectral-contraction-')
    atexit.register(_BUILD_DIRECTORY.cleanup)
    root=Path(_BUILD_DIRECTORY.name)
    target=root/(key+('.dylib' if sys.platform=='darwin' else '.so'))
    command=[compiler,'-std=c++17','-O3','-fPIC','-shared',str(source),'-o',str(target)]
    result=subprocess.run(command,capture_output=True,text=True)
    if result.returncode:raise RuntimeError('native contraction compilation failed: '+result.stderr)
    lib=ctypes.CDLL(str(target));function=lib.oph_apply
    array=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    integer=np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS')
    function.argtypes=[array,integer,ctypes.c_int64,integer,ctypes.c_int,ctypes.c_int,array,array]
    function.restype=ctypes.c_int
    return lib


def layout(bits,targets):
    if (type(bits) is not int or not 1<=bits<=20 or type(targets) is not tuple or
            not 1<=len(targets)<=3 or len(set(targets))!=len(targets) or
            any(type(j) is not int or not 0<=j<bits for j in targets)):
        raise ValueError('bounded complete tensor-factor layout required')
    return _layout(bits,targets)


@lru_cache(maxsize=256)
def _layout(bits,targets):
    mask=sum(1<<j for j in targets)
    rows=np.arange(1<<bits,dtype=np.int64)
    bases=np.ascontiguousarray(rows[(rows&mask)==0])
    offsets=np.array([sum(((index>>j)&1)<<target for j,target in enumerate(targets))
                      for index in range(1<<len(targets))],dtype=np.int64)
    bases.flags.writeable=False;offsets.flags.writeable=False
    return bases,offsets


def apply_inplace(state,matrix,targets,bits):
    if (state.dtype!=np.complex128 or not state.flags.c_contiguous or not state.flags.writeable or
            state.ndim!=2 or state.shape[0]!=1<<bits or not 1<=state.shape[1]<=32):
        raise ValueError('complete writable contiguous complex128 ensemble required')
    dimension=1<<len(targets)
    if matrix.shape!=(dimension,dimension):raise ValueError('local matrix shape differs from targets')
    matrix=np.ascontiguousarray(matrix,dtype=np.complex128)
    if not np.isfinite(matrix).all():raise ValueError('finite local matrix required')
    bases,offsets=layout(bits,targets);density=np.zeros((dimension,dimension),complex)
    code=library().oph_apply(state.view(np.float64),bases,len(bases),offsets,state.shape[1],dimension,
                             matrix.view(np.float64),density.view(np.float64))
    if code:raise ArithmeticError('full-row native contraction rejected its dimensions')
    if not np.isfinite(density).all():raise ArithmeticError('nonfinite native contraction density')
    return density
