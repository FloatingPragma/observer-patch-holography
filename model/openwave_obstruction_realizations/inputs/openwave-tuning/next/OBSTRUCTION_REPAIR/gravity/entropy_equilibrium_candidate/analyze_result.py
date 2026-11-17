"""Post-solve controls; no geometry target enters the preparation solve."""
from pathlib import Path
import json
import numpy as np
from run import region,source,vacuum_and_pair,mixed_entropy,bernoulli_entropy,gaussian_entropy_modular,L,NC,coords,neighbor

HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'RESULT.json').read_text())
a=np.load(HERE/'SOURCE_ARRAYS.npz')
K=a['K_all'];J=a['J_all'];sigma=a['sigma'];H=a['H_reference'];P=a['P_reference'];Q=a['P_pair']
train=K[:NC];jt=J[:NC]
selector_K=np.zeros_like(K)
for axis in range(3):
    for c in range(NC):
        d=neighbor(c,axis);ids=region([c]);ids2=region([c,d])
        S0=gaussian_entropy_modular(P[np.ix_(ids,ids)])[0]
        S1=gaussian_entropy_modular(P[np.ix_(ids2,ids2)])[0]
        selector_K[axis*NC+c,c]=-3*(S1-S0)
        selector_K[axis*NC+c,d]=-1.5*(S1-S0)
field_K=K-selector_K
fourier=[]
for mode in coords:
    wave=np.exp(2j*np.pi*(coords@mode)/L)
    eig=[np.vdot(wave,m[:NC]@wave)/NC for m in (K,field_K,selector_K)]
    fourier.append(dict(mode=mode.tolist(),total=[eig[0].real,eig[0].imag],field=[eig[1].real,eig[1].imag],
                        selector_probability=[eig[2].real,eig[2].imag],selector_Shannon_linear=0,
                        Fourier_residual=float(np.linalg.norm(train@wave-eig[0]*wave))))
project=np.eye(NC)-np.ones((NC,NC))/NC
x,_,rank,_=np.linalg.lstsq(train@project,-jt,rcond=1e-11)
sz=project@x
lap=np.zeros((NC,NC))
for c in range(NC):
    for axis in range(3):
        for sign in (-1,1):
            d=neighbor(c,axis,sign);lap[c,c]+=1;lap[c,d]-=1
dC=Q-P
energy_density=np.real(np.diag((H@dC+dC@H)/2)).reshape(NC,4).sum(axis=1)
rho=energy_density-energy_density.mean()
curvature=lap@sigma
slope=float(rho@curvature/(rho@rho))
remainder=float(np.linalg.norm(curvature-slope*rho)/max(1e-30,np.linalg.norm(curvature)))
eta=json.loads((HERE.parent/'curved_ball/PLANAR_NORMALIZATION.json').read_text())['entropy_per_coordinate_face']
eps=r['validations'][0]['epsilon']
controls=[]
for name,sg in [('no_geometry',np.zeros(NC)),('wrong_sign_geometry',-eps*sigma)]:
    PP,QQ,_,_,_=vacuum_and_pair(source(sg));vol=np.exp(3*sg)
    for c in (0,13):
        d=neighbor(c,0);ids=region([c]);ids2=region([c,d])
        s0,_=mixed_entropy(PP[np.ix_(ids,ids)],QQ[np.ix_(ids,ids)],eps)
        s1,_=mixed_entropy(PP[np.ix_(ids2,ids2)],QQ[np.ix_(ids2,ids2)],eps)
        b0,_=mixed_entropy(P[np.ix_(ids,ids)],Q[np.ix_(ids,ids)],0)
        b1,_=mixed_entropy(P[np.ix_(ids2,ids2)],Q[np.ix_(ids2,ids2)],0)
        w=(1.5-vol[c])/vol[d]
        change=bernoulli_entropy(w)+(1-w)*s0+w*s1-(np.log(2)+.5*(b0+b1))
        controls.append(dict(name=name,cell=c,entropy_change=float(change),linear_expected=float(eps*jt[c]*(2 if name=='wrong_sign_geometry' else 1))))
validation=r['validations']
ratios={key:[validation[i][key]/validation[i+1][key] for i in (0,1)] for key in ['trained_max_abs','heldout_y_max_abs','heldout_z_max_abs','reversed_x_max_abs']}
unused_radius=[]
entropy_decomposition=[]
for case in validation:
    epsilon=case['epsilon'];sg=epsilon*sigma
    pp,qq,_,_,_=vacuum_and_pair(source(sg));vol=np.exp(3*sg)
    changes=[]
    for c in range(NC):
        ids=region([c])
        actual,_=mixed_entropy(pp[np.ix_(ids,ids)],qq[np.ix_(ids,ids)],epsilon)
        ref,_=mixed_entropy(P[np.ix_(ids,ids)],Q[np.ix_(ids,ids)],0)
        w=.5/vol[c]
        changes.append(bernoulli_entropy(w)+w*actual-(np.log(2)+.5*ref))
    unused_radius.append(dict(epsilon=epsilon,protocol='concentric radii .1 and .3 at x+e_x/4: empty versus one cell; expected volume .5',
                              maximum_abs_entropy_change=max(abs(x) for x in changes),changes=changes))
    for c in (0,13):
        d=neighbor(c,0);ids=region([c]);ids2=region([c,d])
        S0=mixed_entropy(pp[np.ix_(ids,ids)],qq[np.ix_(ids,ids)],epsilon)[0]
        S1=mixed_entropy(pp[np.ix_(ids2,ids2)],qq[np.ix_(ids2,ids2)],epsilon)[0]
        B0=gaussian_entropy_modular(P[np.ix_(ids,ids)])[0]
        B1=gaussian_entropy_modular(P[np.ix_(ids2,ids2)])[0]
        w=(1.5-vol[c])/vol[d]
        field_change=.5*((S0-B0)+(S1-B1))
        selector_reference=(w-.5)*(B1-B0)
        selector_cross=(w-.5)*((S1-S0)-(B1-B0))
        shannon=bernoulli_entropy(w)-np.log(2)
        # If the preparation flag were also in the observer algebra, its
        # cq entropy would contain H(epsilon), not the mixture entropy.
        s0cq=(1-epsilon)*gaussian_entropy_modular(pp[np.ix_(ids,ids)])[0]+epsilon*gaussian_entropy_modular(qq[np.ix_(ids,ids)])[0]
        s1cq=(1-epsilon)*gaussian_entropy_modular(pp[np.ix_(ids2,ids2)])[0]+epsilon*gaussian_entropy_modular(qq[np.ix_(ids2,ids2)])[0]
        all_flags=bernoulli_entropy(epsilon)+bernoulli_entropy(w)+(1-w)*s0cq+w*s1cq-(np.log(2)+.5*(B0+B1))
        entropy_decomposition.append(dict(epsilon=epsilon,cell=c,field_change_at_fixed_selector=field_change,
                selector_probability_on_reference=selector_reference,selector_cross=selector_cross,
                selector_Shannon_change=shannon,total=field_change+selector_reference+selector_cross+shannon,
                including_external_preparation_flag_in_regional_algebra=all_flags))
out=dict(zero_mean_constrained_training_rank=int(rank),zero_mean_constrained_residual_max=float(np.max(abs(train@sz+jt))),
         zero_mode_retained_in_original=True,finite_epsilon_halving_ratios=ratios,negative_controls=controls,
         unused_radius=unused_radius,
         response_decomposition=dict(field_norm=float(np.linalg.norm(field_K)),selector_probability_norm=float(np.linalg.norm(selector_K)),
                selector_Shannon_linear_norm=0,uniform_mode_fourier=fourier[0]),
         all_periodic_Fourier_modes=fourier,finite_entropy_decomposition=entropy_decomposition,
         postsolve_geometry_comparison=dict(periodic_energy_background_subtracted_for_comparison_only=True,
             neutralized_energy_density=rho.tolist(),laplacian_sigma=curvature.tolist(),
             fitted_curvature_over_energy_coefficient=slope,relative_non_Poisson_remainder=remainder,
             conditional_tower_coefficient_pi_over_eta=float(np.pi/eta),
             ratio_to_conditional_tower_coefficient=slope/(np.pi/eta),
             fit_used_to_prepare_sigma=False,lapse_rule_supplied=False),
         expected_volume_only=True,matched_gravity_solved=False)
(HERE/'INSPECTION.json').write_text(json.dumps(out,indent=2)+'\n')
np.savez_compressed(HERE/'RESPONSE_DECOMPOSITION.npz',K_total=K,K_field=field_K,K_selector_probability=selector_K,
                    K_selector_Shannon_linear=np.zeros_like(K))
print(json.dumps(out,indent=2))
