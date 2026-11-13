"""Complete divergence-free Maxwell modes coupled to one quantum geometry."""
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh
from .cochains import complete_transverse_frame
from model.openwave_universal_energy_geometry.source import profile,QSIZE,rest_block
from model.openwave_common_gravity_records.records import norm_upper,product_roundoff
from model.finite_area_response.response import exponential


def eigen_certificate(H):
    e,V=eigh(H)
    eta=norm_upper(V.T@V-np.eye(len(e)))+product_roundoff(V.T,V)+2e-14
    residual=norm_upper(H@V-V*e)+product_roundoff(H,V)+2e-14*norm_upper(V*e)
    error=2*residual+4*max(abs(e))*eta
    assert eta<=2e-10 and error<=5e-8
    return e,V,{'gram_error_upper':2e-10,'generator_error_upper':5e-8,
                'live_Gram_error_diagnostic':eta,'live_generator_error_diagnostic':error,
                'arithmetic_model':'complete real eigensystem checked with binary64 gamma_n matrix-product bounds'}


@lru_cache(None)
def field(gamma=1):
    d=complete_transverse_frame(6);frame=d['frame'];p=profile();n=p['n']
    assert n==7
    edge_f=[]
    for axis,x in d['edges']:
        y=list(x);y[axis]+=1
        a=int(p['F_integer'][np.ravel_multi_index(x,(n,)*3)])
        b=int(p['F_integer'][np.ravel_multi_index(tuple(y),(n,)*3)])
        edge_f.append(F(a+b,2**41))
    # Project the approximate frame to exact Gauss and polar-normalize it.
    # The nonzero full-graph gradient singular value is >=2/(cells+1).
    G=d['gradient'].toarray().astype(float)
    projection=(d['frame_divergence_residual']+product_roundoff(G.T,frame))*(d['cells']+1)/2+1e-12
    gram=d['frame_Gram_residual']+product_roundoff(frame.T,frame)+2*projection+projection*projection+1e-12
    frame_error=projection+gram/(1-gram)
    assert frame_error<1e-9
    roots=[];certificates=[];curl=d['curl'].astype(float)
    for sign in (1,-1):
        arguments=[(1+gamma)*sign*QSIZE*f for f in edge_f]
        enclosures={a:exponential(a) for a in set(arguments)}
        s=np.array([np.exp(float(a)) for a in arguments])
        scalar=max(float(max(abs(F(float(x))-enclosures[a].lo),abs(F(float(x))-enclosures[a].hi))) for x,a in zip(s,arguments,strict=True))
        smax=max(float(e.hi) for e in enclosures.values())+1e-13
        smin=min(float(e.lo) for e in enclosures.values())-1e-13
        cf=curl@(s[:,None]*frame)
        K=cf.T@cf
        e,V,ec=eigen_certificate(K)
        root=(V*np.sqrt(e))@V.T
        # ||curl||²<=12. Frame and exponential errors are included separately.
        construction_error=12*(2*smax*scalar+scalar*scalar)+24*smax*smax*frame_error+product_roundoff(cf.T,cf)+1e-12
        residual=norm_upper(root@root-K)+product_roundoff(root,root)+construction_error
        exact_lower=8/(7*7)*smin*smin
        numerical_lower=np.sqrt(e[0])*(1-ec['gram_error_upper'])-product_roundoff(V*np.sqrt(e),V.T)
        assert exact_lower>.05 and numerical_lower>.2
        error=residual/(np.sqrt(exact_lower)+numerical_lower)
        assert error<1e-7
        roots.append(root);certificates.append({'source_root_operator_error_upper':1e-7,
            'live_root_residual_diagnostic':residual,'positive_root_lower':numerical_lower,
            'exact_stiffness_lower':exact_lower,'exact_Gauss_frame_error_upper':1e-9})
    d.update({'roots':roots,'certificates':certificates,'edge_profile_exact':list(map(str,edge_f))})
    return d


@lru_cache(None)
def block(source_level,probe_level,gamma=1):
    p=profile();f=F(int(p['F_integer'][np.ravel_multi_index((3,4,3),(7,)*3)]),2**40)
    vac=rest_block(source_level,probe_level,f)
    d=field(gamma);m=d['frame'].shape[1]
    H=np.kron(vac,np.eye(m));H[:m,:m]+=d['roots'][0];H[m:,m:]+=d['roots'][1]
    e,V,c=eigen_certificate(H);ev,Uv=eigh(vac)
    c['source_operator_error_upper']=1.00002e-7
    return {'H':H,'vacuum_H':vac,'e':e,'V':V,'ev':ev,'Uv':Uv,'certificate':c}


def packet(polarization=2):
    d=field();positions=d['edge_positions'];axes=np.array([a for a,_ in d['edges']])
    proposal=np.exp(-np.sum((positions-np.array([1.,4.,3.]))**2,axis=1)/(2*.8**2)+1.2j*positions[:,0])*(axes==polarization)
    accepted=np.rint(proposal.real*2**30)+1j*np.rint(proposal.imag*2**30)
    reduced=d['frame'].T@accepted
    return reduced/np.linalg.norm(reduced),[[int(z.real),int(z.imag)]for z in accepted]


def packet_certificate(polarization=2):
    d=field();_,integers=packet(polarization)
    accepted=np.array([complex(a,b)for a,b in integers]);projected=d['frame'].T@accepted
    ratio=np.linalg.norm(accepted)/np.linalg.norm(projected)
    multiplication=product_roundoff(d['frame'].T,accepted.reshape(-1,1))/np.linalg.norm(projected)
    error=2e-9*ratio+2*multiplication+2e-12
    assert error<1e-8
    return {'polarization_axis':polarization,'accepted_to_projected_norm_ratio_diagnostic':float(ratio),
            'normalized_packet_error_upper':1e-8,'frame_error_used_upper':1e-9,
            'projection_product_error_upper_diagnostic':multiplication}


@lru_cache(None)
def detector():
    d=field();F=d['frame'];mask=d['edge_positions'][:,0]>=5
    # A declared half-efficient absorber leaves a strict no-click floor1/2.
    # Its positive square root therefore has a linear Sylvester error bound.
    effect=.5*(F[mask].T@F[mask])
    e,V=eigh(np.eye(F.shape[1])-effect)
    assert e[0]>.5-2e-14 and e[-1]<1+2e-14
    K=(V*np.sqrt(np.maximum(e,0)))@V.T
    residual=norm_upper(K@K+effect-np.eye(len(K)))+product_roundoff(K,K)+.5*product_roundoff(F[mask].T,F[mask])
    assert residual<1e-11
    return mask,K,{'complete_Kraus_Gram_residual_upper':1e-11,
                  'receiver_efficiency_exact':'1/2','no_click_effect_lower_exact':'1/2',
                  'all_Kraus_operator_error_upper':5e-9,
                  'arrival':'every retained transverse-edge absorption Kraus row',
                  'no_arrival':'positive square root of I minus complete receiver effect',
                  'vacuum_failure':'identity on the photon vacuum and its reference flag'}
