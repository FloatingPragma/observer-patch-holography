"""Faithful Gibbs response for the same nonquadratic finite Hamiltonian."""
import numpy as np
from scipy.linalg import expm, expm_frechet

from .finite import numerical_block


def gibbs(hamiltonian, beta):
    h = np.asarray(hamiltonian,dtype=complex)
    if h.ndim != 2 or h.shape[0] != h.shape[1] or not np.isfinite(h).all() or not np.isfinite(beta) or beta <= 0:
        raise ValueError('finite square Hamiltonian and positive beta required')
    if np.linalg.norm(h-h.conj().T)>1e-12:
        raise ValueError('Hermitian Hamiltonian required')
    energy, basis = np.linalg.eigh(h)
    weights = np.exp(-beta*(energy-energy.min()))
    probability = weights/weights.sum()
    if probability.min() <= 0:
        raise ValueError('floating arithmetic lost faithfulness; rescale the numerical witness')
    log_z = -beta*energy.min()+np.log(weights.sum())
    rho = (basis*probability)@basis.conj().T
    return {'rho':rho,'probability':probability,'basis':basis,'energy':energy,'log_z':log_z,
            'entropy':float(-probability@np.log(probability)),
            'mean_energy':float(probability@energy),'free_energy':float(-log_z/beta)}


def response(hamiltonian, perturbation, beta):
    value = gibbs(hamiltonian,beta)
    x = np.asarray(perturbation,dtype=complex)
    if x.shape != np.asarray(hamiltonian).shape or not np.isfinite(x).all() or np.linalg.norm(x-x.conj().T)>1e-12:
        raise ValueError('finite Hermitian perturbation in the same algebra required')
    p,u = value['probability'],value['basis']
    mean = float(np.trace(value['rho']@x).real)
    centered = u.conj().T@x@u-mean*np.eye(len(p))
    logarithmic_mean = np.empty((len(p),len(p)))
    for i in range(len(p)):
        for j in range(len(p)):
            if abs(np.log(p[i])-np.log(p[j])) < 1e-12:
                logarithmic_mean[i,j] = (p[i]+p[j])/2
            else:
                logarithmic_mean[i,j] = (p[i]-p[j])/(np.log(p[i])-np.log(p[j]))
    drho = -beta*u@(logarithmic_mean*centered)@u.conj().T
    covariance = float(np.sum(logarithmic_mean*np.abs(centered)**2))
    denergy = float(np.trace(drho@hamiltonian).real+mean)
    log_rho = (u*np.log(p))@u.conj().T
    dentropy = float(-np.trace(drho@log_rho).real)
    return {**value,'drho':drho,'mean_perturbation':mean,'BKM_variance':covariance,
            'energy_derivative':denergy,'entropy_derivative':dentropy,
            'free_energy_second_derivative':-beta*covariance}


def control():
    h,x,_ = numerical_block(1.,.5)
    beta = .7
    value = response(h,x,beta)
    exponential, derivative = expm_frechet(-beta*h,-beta*x)
    z = np.trace(exponential)
    independent = derivative/z-exponential*np.trace(derivative)/(z*z)
    eps = 1e-4
    minus,plus = gibbs(h-eps*x,beta),gibbs(h+eps*x,beta)
    finite_second = (plus['free_energy']-2*value['free_energy']+minus['free_energy'])/(eps*eps)
    errors = {'Frechet_derivative':np.linalg.norm(independent-value['drho']),
              'entropy_energy_work_identity':abs(value['entropy_derivative']-beta*(value['energy_derivative']-value['mean_perturbation'])),
              'free_energy_curvature_difference':abs(finite_second-value['free_energy_second_derivative'])}
    if errors['Frechet_derivative'] > 1e-12 or errors['entropy_energy_work_identity'] > 1e-12 or errors['free_energy_curvature_difference'] > 1e-6:
        raise ArithmeticError('finite Gibbs response control failed')
    if np.linalg.norm(h@x-x@h) < 1 or value['BKM_variance'] <= 0:
        raise ArithmeticError('response witness must be noncommuting and nonzero')
    return {'Hamiltonian': 'm=1, lambda=1/2 in the exact16-state block',
            'beta': '7/10', 'faithful_state': True,
            'BKM_variance':format(value['BKM_variance'],'.12g'),
            'free_energy_second_derivative':format(value['free_energy_second_derivative'],'.12g'),
            'noncommuting_derivative_checked_against_matrix_exponential_Frechet_derivative':True,
            'free_energy_curvature_checked_against_centered_difference':True,
            'state_preparation': 'declared spectral Gibbs reset with all input outcomes retained',
            'native_CE_produces_this_Gibbs_state':False,
            'Gibbs_stationarity_implies_geometric_entropy_stationarity':False}
