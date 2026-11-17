"""Complete-band positive flattening, in dimensionless proper R units."""
import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
PAULI = np.array([X, Y, Z])
ALPHA = np.array([np.kron(s, X) for s in PAULI])
BETA = np.kron(I2, Z)
SPIN = np.array([np.kron(s, I2) / 2 for s in PAULI])
I4 = np.eye(4, dtype=complex)


def symbol(p, R=1., mass=0.):
    """All input momenta, full four-component source, no band deletion."""
    p = np.asarray(p, float)
    r = np.linalg.norm(p, axis=-1)
    D = np.einsum('...j,jab->...ab', p, ALPHA)
    b = mass + R**3 * r**4
    E_uv = np.sqrt(r*r + b*b)
    # At the massless isolated zero mode choose the beta reference filling.
    den = np.where(E_uv > 0, E_uv, 1.)
    Hhat = (D + b[..., None, None] * BETA) / den[..., None, None]
    Hhat = np.where((E_uv == 0)[..., None, None], BETA, Hhat)
    E = np.sqrt(r*r + mass*mass)
    H = E[..., None, None] * Hhat
    Pminus = (I4 - Hhat) / 2
    return H, Pminus, E


def rotation(p, R=1., mass=0.):
    p = np.asarray(p, float)
    r = np.linalg.norm(p, axis=-1)
    den = np.where(r > 0, r, 1.)
    A = np.einsum('...j,jab->...ab', p / den[..., None], ALPHA)
    # Difference of the two polar angles; mass must be nonnegative here.
    phi = np.arctan2(mass + R**3*r**4, r) - np.arctan2(mass, r)
    phi = np.where(r == 0, 0., phi)
    return np.cos(phi/2)[..., None, None]*I4 + np.sin(phi/2)[..., None, None]*(BETA @ A)


def rotation_derivatives(p, R=1.):
    """Analytic first and second momentum derivatives, for nonzero p."""
    p = np.asarray(p, float)
    r = np.linalg.norm(p)
    n = p/r
    A = np.einsum('j,jab->ab', n, ALPHA)
    Ai = np.array([(ALPHA[i] - A*n[i])/r for i in range(3)])
    Aij = np.empty((3, 3, 4, 4), complex)
    for i in range(3):
        for j in range(3):
            Aij[i,j] = (-ALPHA[i]*n[j]-ALPHA[j]*n[i]-A*(i==j)+3*A*n[i]*n[j])/r**2
    phi = np.arctan(R**3*r**3)
    F = 3*R**3*r/(1+R**6*r**6)
    Fp = 3*R**3*(1-5*R**6*r**6)/(1+R**6*r**6)**2
    phii = F*p
    phiij = F*np.eye(3)+(Fp/r)*np.outer(p,p)
    c, s = np.cos(phi/2), np.sin(phi/2)
    ci, si = -s*phii/2, c*phii/2
    cij = -c*np.outer(phii,phii)/4-s*phiij/2
    sij = -s*np.outer(phii,phii)/4+c*phiij/2
    V = c*I4+s*(BETA@A)
    Vi = np.array([ci[i]*I4+si[i]*(BETA@A)+s*(BETA@Ai[i]) for i in range(3)])
    Vij = np.empty((3,3,4,4),complex)
    for i in range(3):
        for j in range(3):
            Vij[i,j] = cij[i,j]*I4+sij[i,j]*(BETA@A)+si[i]*(BETA@Ai[j])+si[j]*(BETA@Ai[i])+s*(BETA@Aij[i,j])
    return V, Vi, Vij


def tilted_covariance(q, speed, R=1.):
    """Proper-coordinate plane t=u x, boost along x; exact real bands.

    The source fields are evaluated on the plane; no transported-net
    assumption is made. Canonical whitening is performed on the full Gram.
    """
    q = np.asarray(q, float)
    shape = q.shape[:-1]
    Eobs = np.linalg.norm(q, axis=-1)
    gamma = 1/np.sqrt(1-speed**2)
    G = np.zeros(shape+(4,4),complex)
    N = np.zeros_like(G)
    records = []
    for sign in [-1,1]:
        k = q.copy()
        k[...,0] = gamma*(q[...,0]+sign*speed*Eobs)
        Elab = gamma*(Eobs+sign*speed*q[...,0])
        _, Pm, _ = symbol(k,R)
        Ps = Pm if sign == -1 else I4-Pm
        velocity = np.divide(k[...,0],Elab,out=np.zeros_like(Elab),where=Elab>0)
        jac = gamma*(1-speed*sign*velocity)
        contribution = Ps/jac[...,None,None]
        G += contribution
        if sign == -1:
            N = contribution
        records.append((k,jac))
    e, U = np.linalg.eigh(G)
    if np.min(e) <= 0:
        raise ValueError('Nonpositive full CAR Gram')
    invsqrt = (U*e[...,None,:]**(-.5))@U.conj().swapaxes(-1,-2)
    C = invsqrt@N@invsqrt
    return C,G,records


def entropy(C):
    e = np.clip(np.linalg.eigvalsh((C+C.conj().T)/2),0,1)
    from scipy.special import xlogy
    return float(np.sum(-xlogy(e,e)-xlogy(1-e,1-e)))
