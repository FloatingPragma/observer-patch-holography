"""Complete finite control of the refinement theorem's exact identities."""
from pathlib import Path
import json
import numpy as np
from scipy.special import xlogy

HERE = Path(__file__).resolve().parent


def entropy(e):
    e = np.clip(e, 0, 1)
    return float(np.sum(-xlogy(e, e)-xlogy(1-e, 1-e)))


def main():
    N, length, region_sites = 64, 32., 32
    x = np.arange(N)*length/N
    momentum = 2*np.pi*(np.arange(-N//2, N//2)+.5)/length
    F = np.exp(1j*x[:, None]*momentum[None, :])/np.sqrt(N)
    X = np.array([[0., 1.], [1., 0.]])
    Z = np.diag([1., -1.])
    blocks = []
    for p in momentum:
        H = p*X+p**4*Z
        blocks.append((np.eye(2)-H/np.sqrt(p*p+p**8))/2)
    P = np.einsum('xk,kab,yk->xayb', F, np.array(blocks), F.conj()).reshape(2*N, 2*N)
    beta_reference = np.kron(np.eye(N), np.diag([0., 1.]))
    chi = np.diag(np.repeat(np.arange(N)<region_sites, 2).astype(float))
    D = P-beta_reference
    assert np.linalg.norm(P@P-P, 2) < 1e-12
    A = chi@P@(np.eye(2*N)-chi)
    C = P[:2*region_sites, :2*region_sites]
    exact_entropy = entropy(np.linalg.eigvalsh(C))
    records = []
    previous = None
    for width in (32, 16, 8, 4, 2, 1):
        ncell = region_sites//width
        W = np.zeros((N, ncell))
        for k in range(ncell):
            W[k*width:(k+1)*width, k] = 1/np.sqrt(width)
        U = np.kron(W, np.eye(2))
        Q = U@U.T
        assert np.linalg.norm(U.T@U-np.eye(2*ncell), 2) < 1e-13
        assert np.linalg.norm(Q@beta_reference-beta_reference@Q, 2) < 1e-13
        assert np.linalg.norm(chi@Q-Q, 2) < 1e-13
        if previous is not None:
            assert np.linalg.norm(Q@previous-previous, 2) < 1e-13
        previous = Q
        compressed = U.T@P@U
        An = Q@P@(np.eye(2*N)-Q)
        assert np.linalg.norm(An-Q@D@(np.eye(2*N)-Q), 2) < 1e-13
        variance = U@(compressed-compressed@compressed)@U.T
        identity_error = float(np.linalg.norm(variance-An@An.conj().T, 2))
        assert identity_error < 1e-12
        singular = np.linalg.svd(An, compute_uv=False)
        roots = (1-np.sqrt(np.clip(1-4*singular**2, 0, 1)))/2
        direct = entropy(np.linalg.eigvalsh(compressed))
        by_offdiagonal = entropy(roots)
        assert abs(direct-by_offdiagonal) < 2e-11
        sdelta = np.linalg.svd(An-A, compute_uv=False)
        records.append({'cell_width_in_grid_sites': width, 'observed_one_particle_dimension': 2*ncell,
            'entropy': direct, 'full_finite_region_entropy': exact_entropy,
            'S_3_over_2_offdiagonal_error': float(np.sum(sdelta**1.5)**(2/3)),
            'variance_identity_error': identity_error,
            'entropy_identity_error': abs(direct-by_offdiagonal)})
    assert abs(records[-1]['entropy']-exact_entropy) < 1e-13
    out = {'schema': 'oph.cell_entropy_limit.finite_control.v1', 'records': records,
        'all_modes_retained': True, 'continuum_error_rate_certified': False,
        'entropy_monotonicity_assumed': False, 'physical_gravity_solved': False}
    (HERE/'RESULT.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2), flush=True)


if __name__ == '__main__':
    main()
