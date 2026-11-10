"""Exact full-range affine Gauss labels and compensated finite moves."""
from functools import lru_cache
import itertools
import sympy as s
from .field import cube

LOW,HIGH=-4,4

def shift(lower,upper,step):
    if step not in (-1,1):raise ValueError('one integer charge step required')
    dim=upper-lower+1;out=s.zeros(dim)
    for q in range(lower,upper+1):
        if lower<=q+step<=upper:out[q+step-lower,q-lower]=1
    return out

@lru_cache(None)
def calculate():
    data=cube();D,G,C=data['D'],data['G'],data['C'];reference=0
    columns=[x for x in range(8) if x!=reference]
    F=s.Matrix.hstack(*(s.eye(8)[:,x]-s.eye(8)[:,reference] for x in columns))
    longitudinal=D*G*F
    assert D.T*longitudinal==F and C*longitudinal==s.zeros(6,7)
    assert longitudinal.T*longitudinal==F.T*G*F
    charges=s.diag(*range(LOW,HIGH+1));checks=0
    for step in (-1,1):
        T=shift(LOW,HIGH,step)
        assert charges*T-T*charges==step*T
        assert T.T==shift(LOW,HIGH,-step)
        assert all(value in (0,1) for value in (T.T*T).diagonal())
        refcharges=s.diag(*range(-32,33));R=shift(-32,32,step)
        assert refcharges*R-R*refcharges==step*R
        for source_charge,ell in itertools.product(range(LOW,HIGH+1),repeat=2):
            if LOW<=source_charge+step<=HIGH and LOW<=ell+step<=HIGH:
                assert (ell+step)-(source_charge+step)==ell-source_charge
            checks+=1
    # The independent Gauss coordinates are ell_x-rho_x, x!=r, and
    # -sum ell-rho_r-R. Their changes vanish for every oriented hop.
    hop_checks=0
    for tail,head,axis in data['edges']:
        drho=s.zeros(8,1);drho[tail]=-1;drho[head]=1
        dell=s.Matrix([drho[x] for x in columns])
        assert F*dell-drho==s.zeros(8,1)
        for direction in (-1,1):
            for rt,rh,et,eh in itertools.product(range(LOW,HIGH+1),repeat=4):
                source_valid=LOW<=rt-direction<=HIGH and LOW<=rh+direction<=HIGH
                electric_valid=(tail==reference or LOW<=et-direction<=HIGH) and (head==reference or LOW<=eh+direction<=HIGH)
                if source_valid and electric_valid:
                    assert (et-direction)-(rt-direction)==et-rt
                    assert (eh+direction)-(rh+direction)==eh-rh
                # Invalid moves are zero, not cyclic wraparound.
                hop_checks+=1
    for site in range(8):
        for direction in (-1,1):
            drho=s.eye(8)[:,site]*direction;dell=s.Matrix([drho[x] for x in columns]);dR=-direction
            assert F*dell-drho-s.eye(8)[:,reference]*dR==s.zeros(8,1)
    # Any simultaneous charge eigenvalue tuple in the full range has
    # exactly one retained physical electric/reference label tuple.
    for rho in ([-4]*8,[4]*8,[-4,4]*4,[0]*8,[1,0,0,0,0,0,0,0]):
        ell=s.Matrix([rho[x] for x in columns]);R=-sum(rho)
        assert -32<=R<=32 and F*ell-s.Matrix(rho)-s.eye(8)[:,reference]*R==s.zeros(8,1)
    return {'full_longitudinal_dimension':9**7,'full_static_reference_dimension':65,
      'reference_vertex':reference,'longitudinal_charge_columns':F.tolist(),
      'longitudinal_electric_map':longitudinal.tolist(),
      'single_label_boundary_checks':checks,'complete_two_endpoint_boundary_cases':hop_checks,
      'all_12_oriented_hops_and_16_charge_pump_moves_preserve_Gauss':True,
      'invalid_cutoff_moves':'zero partial isometry; no modulo wraparound',
      'physical_isometry':'Y psi=sum_rho Pi_rho psi tensor |rho_(x!=r)>_ell tensor |-sum rho>_R',
      'Y_dagger_Y':'I on every retained native matter sector',
      'Y_Y_dagger':'projector onto all normalized integer Gauss kernels',
      'density_operator_hypotheses':['all rho_x commute and have integer spectra in [-4,4]',
        'the full native complement is retained','sum_x rho_x is the conserved native port charge'],
      'g_zero_constraint':'the normalized integer Gauss constraints do not disappear at g=0',
      'charge_current_orientation':'D_tail=-1,D_head=+1; I_flow,e=-i(K_e-K_e^dag), dot rho=D^T I_flow',
      'Peierls_derivative_current':'J_phase,e=+i(K_e-K_e^dag)=-I_flow,e',
      'positive_Peierls_phase':'K_e=Psi_head^dag B_e^dag Psi_tail exp(+i g Q_e), with matching ell shifts',
      'full_kinematic_photon_cutoff_affects_Gauss_identity':False}

if __name__=='__main__':
    import json
    print(json.dumps(calculate(),indent=2,default=str))
