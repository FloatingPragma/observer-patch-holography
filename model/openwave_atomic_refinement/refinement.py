"""Lowest parity sectors of the same finite Dirichlet Gauss Hamiltonian."""
from functools import lru_cache
import numpy as np
from scipy.sparse.linalg import eigsh
from model.openwave_field_observables.orbitals import reduced_hamiltonian,dirichlet_source


@lru_cache(maxsize=None)
def spectrum(n,length):
    rows=[]
    for odd in ((False,False,False),(True,False,False)):
        h,a=reduced_hamiltonian(n,length,odd)
        values,vectors=eigsh(h,k=2 if not any(odd) else 1,which='SA',tol=3e-13,
            v0=np.random.default_rng(18067).normal(size=h.shape[0]))
        residuals=[float(np.linalg.norm(h@vectors[:,i]-values[i]*vectors[:,i])) for i in range(len(values))]
        if max(residuals)>2e-9:raise ArithmeticError('eigenpair residual exceeds declared diagnostic budget')
        row={'parity':list(odd),'dimension':h.shape[0],
             'energies':[round(float(x),12) for x in values],
             'residual_norms':[float(x) for x in residuals]}
        if not any(odd):
            # Exact positive-vector theorem; these float ratios are diagnostics,
            # not outward-rounded bounds for the exact Green source.
            u=np.abs(vectors[:,0]);u=np.maximum(u,np.finfo(float).tiny)
            ratios=(h@u)/u
            row['positive_vector_ratio_diagnostic']=[float(np.min(ratios)),float(np.max(ratios))]
            row['minimum_trial_component']=float(np.min(u))
        rows.append(row)
    return {'grid':n,'length':length,'spacing':a,'full_spatial_sites':n**3,
            'sectors':rows,'ground_to_odd':round(rows[1]['energies'][0]-rows[0]['energies'][0],12),
            'second_shell_splitting':round(abs(rows[1]['energies'][0]-rows[0]['energies'][1]),12)}


def volume_control():
    small,pot_s,a=dirichlet_source(7,8.)
    large,pot_l,b=dirichlet_source(11,12.)
    assert a==b==1
    nested=pot_l[2:9,2:9,2:9]
    if np.max(nested-pot_s)>1e-12:raise ArithmeticError('Dirichlet Green domain monotonicity failed')
    return {'fixed_spacing':a,'small_grid':7,'large_grid':11,
            'maximum_nested_potential_change':float(np.max(nested-pot_s)),
            'strictly_more_attractive_in_larger_box':bool(np.max(nested-pot_s)<0)}


def run():
    # Shared a=1/3 tests volume independently of spacing. Shared L=24 tests spacing.
    volume=[spectrum(71,24.),spectrum(95,32.),spectrum(119,40.)]
    spacing=[volume[0],spectrum(95,24.),spectrum(127,24.),spectrum(159,24.)]
    for j in range(2):
        for sector in range(2):
            if volume[j+1]['sectors'][sector]['energies'][0]>volume[j]['sectors'][sector]['energies'][0]+2e-9:
                raise ArithmeticError('nested-box variational order failed')
    if not all(spacing[i+1]['second_shell_splitting']<spacing[i]['second_shell_splitting'] for i in range(3)):
        raise ArithmeticError('second-shell refinement diagnostic failed')
    return {'schema':'oph.openwave.atomic_refinement.v1','volume_rows':volume,'spacing_rows':spacing,
            'independent_nested_Green_control':volume_control(),
            'hydrogen_reference_gap':.375,'mass_input':1,'coupling_squared_input':'4*pi',
            'continuum_spacing_error_certified':False,'physical_hydrogen_dictionary_derived':False,
            'residual_meaning':'Every residual bounds distance to some eigenvalue of its numerical finite matrix. It does not alone certify eigenvalue index, floating error, volume error or continuum error.',
            'positive_vector_ratio_meaning':'Collatz bounds are exact for strictly positive trial vectors and exact matrices; displayed floating ratios are diagnostics, with no outward-rounded claim.',
            'scope':'Independent volume and spacing studies of the declared nonrelativistic Dirichlet Green atom; all parity/source conventions unchanged. No physical precision spectrum or spontaneous atom formation is claimed.',
            'verdict':'PASS_SCOPED_ATOMIC_VOLUME_SPACING_CONTROL'}
