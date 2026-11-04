"""Finite compact packets acted on by the literal bounded radial source.

Arrays cover the packet support and its one-hop boundary. A distant anchor
is an integer coordinate, not an omitted volume of occupied rows.
"""
from itertools import product
from math import fsum

import numpy as np

from model.native_mass_channels import algebra as source
from .bounds import exterior_packet_bound


def packet(cells, displacement, momentum, spacing=1, velocity=1, mass=1,
           energy_sign=1, seed=None):
    bounds = exterior_packet_bound(spacing, velocity, mass, cells, displacement, momentum)
    if type(energy_sign) is not int or energy_sign not in (-1, 1):
        raise ValueError('energy sign must be exactly -1 or1')
    if cells > 64:
        raise ValueError('numerical diagnostic limited to64 cells per axis')
    if displacement > 2**52:
        raise ValueError('numerical displacement must be at most2^52; exact bounds have no such cap')
    a, v, m = map(float, (spacing, velocity, mass))
    p = np.asarray(tuple(map(float, momentum)))
    gamma, beta = source.folded_generators()
    h = sum(v/a*np.sin(a*p[j])*np.asarray(gamma[j], complex) for j in range(3))
    h += m*np.asarray(beta[3], complex)
    energy = energy_sign*np.sqrt(m*m+(v/a)**2*np.sum(np.sin(a*p)**2))
    if seed is None:
        seed = np.array([complex(j+1, (-1)**j*(17-j)) for j in range(16)])
    else:
        seed = np.asarray(seed)
        if seed.shape != (16,) or not np.isfinite(seed).all():
            raise ValueError('finite sixteen-row seed required')
    spinor = (seed+h@seed/energy)/2
    norm = np.linalg.norm(spinor)
    if norm < 1e-12:
        raise ValueError('seed has zero or unresolved projection onto the requested energy branch')
    spinor /= norm
    length = 2*cells+2
    psi = np.zeros((length, length, length, 2), complex)
    line = np.sqrt(2/(cells+1))*np.sin(np.pi*(np.arange(cells)+1)/(cells+1))
    envelope = line[:, None, None]*line[None, :, None]*line[None, None, :]
    cell = np.indices((cells, cells, cells)).transpose(1, 2, 3, 0)
    for bits in product((0, 1), repeat=3):
        coefficient = sum((-1)**sum(b*x for b, x in zip(valley, bits, strict=True))*spinor[2*j:2*j+2]
                          for j, valley in enumerate(source.SITES))/np.sqrt(8)
        # Remove the anchor's common phase before numerical range reduction.
        # Formula(9) differs only by exp(i*a*p.anchor), an irrelevant phase.
        points = 2*cell+np.array(bits)
        amplitude = envelope*np.exp(1j*a*(points@p))
        target = tuple(slice(1+b, 1+2*cells+b, 2) for b in bits)
        psi[target] = amplitude[..., None]*coefficient
    return {'state': psi, 'energy': float(energy), 'fiber_spinor': spinor,
            'source_fiber': h, 'bounds': bounds,
            'anchor': (0, 0, displacement), 'spacing': a, 'velocity': v, 'mass': m}


def act(value, *, freeze_direction=False):
    """Literal i*v*sigma_j/(2*a) hopping plus original radial mass."""
    if type(freeze_direction) is not bool:
        raise ValueError('Boolean diagnostic profile switch required')
    psi = value['state']
    length = psi.shape[0]
    n = np.indices((length, length, length)).transpose(1, 2, 3, 0)-1
    n += np.array(value['anchor'])
    pauli = [np.asarray(sigma, complex) for sigma in source.SIGMA]
    result = np.zeros_like(psi)
    for j in range(3):
        head = [slice(None)]*3; tail = [slice(None)]*3
        head[j] = slice(None, -1); tail[j] = slice(1, None)
        head, tail = tuple(head), tuple(tail)
        result[head] += 1j*value['velocity']/(2*value['spacing'])*np.einsum('ij,...j->...i', pauli[j], psi[tail])
        result[tail] -= 1j*value['velocity']/(2*value['spacing'])*np.einsum('ij,...j->...i', pauli[j], psi[head])
    if freeze_direction:
        directions = np.zeros_like(n, dtype=float)
        directions[..., 2] = 1
    else:
        norms = np.linalg.norm(n, axis=-1)
        directions = np.divide(n, norms[..., None], out=np.zeros_like(n, dtype=float), where=norms[..., None] != 0)
    for j in range(3):
        coefficient = value['mass']*directions[..., j]*np.where(n[..., j] % 2, -1, 1)
        result += coefficient[..., None]*np.einsum('ij,...j->...i', pauli[j], psi)
    return result


def diagnostic(cells=8, displacement=10000, momentum=(1, 0, 1), spacing=1, velocity=1, mass=1):
    value = packet(cells, displacement, momentum, spacing, velocity, mass)
    psi, energy = value['state'], value['energy']
    actual, frozen = act(value), act(value, freeze_direction=True)
    return {'cells': cells, 'displacement': displacement,
            'norm_error': abs(fsum(float(z.real)**2+float(z.imag)**2 for z in psi.flat)-1),
            'fiber_eigenvector_error': float(np.linalg.norm(value['source_fiber']@value['fiber_spinor']-energy*value['fiber_spinor'])),
            'kinetic_envelope_residual': float(np.linalg.norm((frozen-energy*psi).ravel())),
            'texture_residual': float(np.linalg.norm((actual-frozen).ravel())),
            'actual_radial_residual': float(np.linalg.norm((actual-energy*psi).ravel())),
            'bounds': value['bounds'], 'floating_diagnostic': True,
            'all_packet_and_boundary_rows_retained': True}
