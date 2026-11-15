"""Finite instance of particle_primary_mass_readout/PROOF.md coupling (1).

This evaluates source-marginal spectral laws, not a physical mass or apparatus.
"""
import numpy as np

def spectral_coupling(hamiltonian,density,probe,tolerance=1e-10):
    h,rho,b=map(lambda x:np.asarray(x,dtype=complex),(hamiltonian,density,probe))
    if h.ndim!=2 or h.shape[0]!=h.shape[1] or rho.shape!=h.shape or b.shape!=h.shape:raise ValueError('equal square source operators required')
    if not all(np.isfinite(x).all() for x in (h,rho,b)):raise ValueError('finite operators required')
    if np.max(abs(h-h.conj().T))>tolerance or np.max(abs(rho-rho.conj().T))>tolerance:raise ValueError('Hermitian H and rho required')
    if abs(np.trace(rho)-1)>tolerance or np.linalg.eigvalsh(rho).min()<-tolerance:raise ValueError('density required')
    if np.max(abs(h@rho-rho@h))>tolerance:raise ValueError('same-clock stationary density required')
    energies,u=np.linalg.eigh(h)
    # Diagonalize rho within each degenerate H block, retaining all initial weights.
    groups=[]
    for i,e in enumerate(energies):
        if not groups or abs(e-energies[groups[-1][0]])>tolerance:groups.append([i])
        else:groups[-1].append(i)
    for ids in groups:
        _,v=np.linalg.eigh(u[:,ids].conj().T@rho@u[:,ids]);u[:,ids]=u[:,ids]@v
    probabilities=np.diag(u.conj().T@rho@u).real
    matrix=u.conj().T@b@u
    weights=abs(matrix)**2*probabilities[None,:]
    transfer=energies[:,None]-energies[None,:]
    excitation=np.broadcast_to((energies-energies[0])[:,None],weights.shape)
    excess=float(np.trace(rho@(h-energies[0]*np.eye(len(h)))).real)
    cap=float(np.linalg.norm(b,2)**2)
    cost=float(np.sum(weights*abs(excitation-transfer)))
    if cost>cap*excess+100*tolerance:raise ValueError('same-source transport inequality failed')
    return {'energies':energies.tolist(),'initial_weights':probabilities.tolist(),'joint_original_weights':weights.tolist(),
            'energy_transfers':transfer.tolist(),'final_excitations':excitation.tolist(),'original_mass':float(weights.sum()),
            'initial_ground_excess':excess,'transport_cost':cost,'transport_cost_upper':cap*excess,
            'success_probability_division':False,'reference_counterevolution_used':False,'physical_mass_identified':False}
