"""Exact logical algebra of the native written-contrast memory and one controlled extension."""
import numpy as np
import sympy as sy

from .. import law
from .source import require


def transfer_algebra():
    # X = sqrt(2) W, split by the two repaired output loads y. A Hermitian memory
    # generator preserves the written code exactly iff it commutes with W W*.
    X = law._write_pattern().reshape(10,9,9)
    blocks = np.einsum('ayk,bzk->yzab', X,X).reshape(-1,10,10)
    I = np.eye(10,dtype=np.int64)
    equations = [np.kron(I,B.T)-np.kron(B,I) for B in blocks if np.any(B)]
    constraints = sy.Matrix(np.unique(np.vstack(equations),axis=0))
    kernel = constraints.nullspace()
    images = []
    odd = [law.SRC.index(v) for v in [(0,1),(1,0)]]
    second = [law.SRC.index(v) for v in [(1,2),(2,1)]]
    doubled = [law.SRC.index(v) for v in [(0,2),(2,0)]]
    for vec in kernel:
        G = np.asarray(vec,dtype=object).reshape(10,10)
        M = sum(B.T@G@B for B in X.transpose(1,0,2))  # twice the induced matrix
        form = np.zeros((9,9),dtype=object)
        for i in (0,4,8):
            form[i,i] = M[0,0]
        for rows,src in ((odd,odd),(second,odd),(doubled,doubled)):
            form[np.ix_(rows,rows)] = M[np.ix_(src,src)]
        require(np.array_equal(form,M), 'unexpected memory-only logical action')
        images.append(list(M.flat))
    dimension = sy.Matrix(images).rank()
    require(dimension == 9, 'memory-only algebra missing a logical matrix unit')
    return dict(label_commutant_dimension=len(kernel), logical_dimension=dimension,
                blocks=['equal loads: common scalar', 'unit contrast: repeated 2x2', 'double contrast: 2x2'])


def local_conversion():
    """Prove balanced-output memory control leaves no written record on any logical input."""
    X = sy.Matrix(law._write_pattern().tolist())
    W = X/sy.sqrt(2)
    def even_label(d):
        return sy.Matrix([int(t==d) for _,t in law.LABELS])/sy.sqrt(2)
    v = even_label(0)-(even_label(2)+even_label(-2))/sy.sqrt(2)
    F = sy.zeros(9)
    F[4,4] = 1
    G = sy.kronecker_product(v*v.T,F)
    u = sy.zeros(9,1)
    u[4],u[2],u[6] = 1,-1/sy.sqrt(2),-1/sy.sqrt(2)
    V = u*u.T
    require((G*W-W*V).applyfunc(sy.simplify)==sy.zeros(90,9),'controlled memory leaks')
    # Removing the output predicate lets the vacuum leak out of the native code.
    bare = sy.kronecker_product(v*v.T,sy.eye(9))
    vac = bare*W[:,0]
    leakage = vac-W*(W.T*vac)
    require(sy.simplify(leakage.dot(leakage))==1, 'missing vacuum counterexample')
    return np.asarray(V,dtype=float)
