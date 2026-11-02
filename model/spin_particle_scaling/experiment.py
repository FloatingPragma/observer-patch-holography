"""Complete 256-state finite controls, not an infinite-volume approximation claim."""
import numpy as np
from model.spin_goldstone_sector.experiment import matrices


def permutation_matrix(site_permutation):
    size=256;out=np.zeros((size,size),complex)
    for mask in range(size):
        target=sum(((mask>>x)&1)<<site_permutation[x] for x in range(8))
        out[target,mask]=1
    return out


def rotation(axis,angle):
    sigma=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]))[axis]
    local=np.cos(angle/2)*np.eye(2)-1j*np.sin(angle/2)*sigma
    out=np.ones((1,1),complex)
    for _ in range(8):out=np.kron(out,local)
    return out


def check():
    h,ops,order=matrices(.2);values,vectors=np.linalg.eigh(h)
    if (np.linalg.norm(vectors.conj().T@vectors-np.eye(256))>1e-10
            or np.linalg.norm(h@vectors-vectors*values)>1e-10):
        raise ArithmeticError('finite eigensolver controls failed')
    g=vectors[:,0];gaps=values-values[0];Q=sum((o[0] for o in ops),np.zeros_like(h))
    if gaps[1]<1e-6:raise ArithmeticError('finite field ground is not isolated')
    if np.linalg.norm(Q@g)>1e-10:raise ArithmeticError('wrong residual charge of ground')
    A=[sum((o[a] for o in ops),np.zeros_like(h))/np.sqrt(8) for a in (1,2)]
    v=[a@g for a in A];z=[vectors.conj().T@x for x in v]
    errors=[];significant_clusters=0;start=0
    for end in range(1,257):
        if end<256 and gaps[end]-gaps[end-1]<1e-9:continue
        block=np.array([[np.vdot(z[a][start:end],z[b][start:end]) for b in range(2)] for a in range(2)])
        errors.extend([abs(block[0,0]-block[1,1]),abs(block[0,1])])
        if block.trace().real>1e-10:significant_clusters+=1
        start=end
    if significant_clusters<3 or np.linalg.norm(v[0])**2<1e-4:
        raise ArithmeticError('finite positive transverse response absent')
    vp=(v[0]+1j*v[1])/np.sqrt(2);vm=(v[0]-1j*v[1])/np.sqrt(2)
    errors.extend([np.linalg.norm(Q@vp-vp),np.linalg.norm(Q@vm+vm),abs(np.vdot(vp,vm))])
    # Origin-centered C4 about axis1 swaps the y,z coordinates of a mod2 cube.
    p=[(x&1)|((x&2)<<1)|((x&4)>>1) for x in range(8)]
    c4=permutation_matrix(p)@rotation(0,np.pi/2)
    phase=np.vdot(g,c4@g)
    errors.extend([np.linalg.norm(c4@g-phase*g),np.linalg.norm(c4@vp+1j*phase*vp),
                   np.linalg.norm(c4@vm-1j*phase*vm)])
    for j in (1,2):
        D=permutation_matrix([x^(1<<j) for x in range(8)])@rotation(j,np.pi)
        phase=np.vdot(g,D@g)
        errors.append(np.linalg.norm(D@g-phase*g))
        for a in range(2):
            sign=1 if a+1==j else -1
            errors.append(np.linalg.norm(D@v[a]-sign*phase*v[a]))
    # The two original-spin components have different physical Fourier signs.
    for a in (1,2):
        carrier=7^(1<<a)
        source_ops=[(-1)**((x&carrier).bit_count())*ops[x][a] for x in range(8)]
        reconstructed=sum(((-1)**((x&carrier).bit_count())*source_ops[x] for x in range(8)),np.zeros_like(h))/np.sqrt(8)
        errors.append(np.linalg.norm((reconstructed-A[a-1])@g))
    if max(errors)>2e-11:raise ArithmeticError('finite symmetry/spectral control failed')
    # Global spin symmetry alone does not imply equal ± norms in a polarized state.
    up=np.array([1,0],complex);sx=np.array([[0,1],[1,0]],complex)/2;sy=np.array([[0,-1j],[1j,0]],complex)/2
    plus=(sx+1j*sy)@up;minus=(sx-1j*sy)@up
    if np.linalg.norm(plus)!=0 or abs(np.linalg.norm(minus)-1)>1e-15:raise ArithmeticError('polarized countercontrol failed')
    return {'full_spin_dimension':256,'periodic_bonds_with_parallel_edges':24,
            'field':'.2','residual_charge_of_ground':0,
            'spectral_matrix_scalar_at_each_cluster':True,
            'circular_charges':[1,-1],'circular_channels_orthogonal':True,
            'actual_translation_and_C4_actions_checked':True,
            'every_finite_energy_eigenspace_checked':True,
            'significant_positive_transverse_response_checked':True,
            'one_spin_U1_alone_does_not_force_equal_weights':True,
            'finite_control_proves_infinite_volume_limit':False,
            'native_pulses_executed':False}
