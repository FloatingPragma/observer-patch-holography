"""Complete spin-zero meson relative fibres of the finite-range strong member."""
from fractions import Fraction as F
import sympy as sy
from model.openwave_hadron_refinement_flow.analytic import constants
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_hadron_refinement_flow.certified_rest import eigen_intervals
from model.openwave_capacity_stress.source import enclosure
from model.openwave_strange_composites.current import mass_parameters
from model.openwave_baryon_dispersion.source import repair_matrix


def data(strange=False):
    # constants applies dn to mbar and dU to U before returning p. Thus
    # every pair has +2dn and split/outer states also have +8dU/3.
    # The native two-cell control below deliberately uses unshifted inputs.
    _,_,_,p,k,u=constants(*degree_coefficients());dm,gs=mass_parameters()
    c=F(3,4)*(p['Jw']+p['eta']);light=F(3,4)*(p['Jw']+p['eta'])*(1 if strange else 2)
    # Charged pi-/K- have one up-type and one down-type constituent: D=0.
    base=2*p['mbar']+2*p['t']+light-c+(dm+gs if strange else 0)
    edge=base+2*p['A']/9+8*p['U']/3+F(3,2)*p['Js']+3*p['spin_gradient']
    outside=2*p['mbar']-2*p['t']+light+(dm if strange else 0)+2*p['A']/9+8*p['U']/3+4*k/3+F(3,2)*(p['Js']+p['spin_gradient'])
    return base,edge,outside,p,u


def degree_ancestry():
    from model.openwave_hadron_refinement_flow.positivity import inputs
    p,_,_=inputs();dU,dn=degree_coefficients()
    corrected=data(False)[3]
    assert corrected['mbar']-p['mbar']==dn and corrected['U']-p['U']==dU
    return {'bulk_colocated_pair_shift_source_exact':str(2*dn),
            'bulk_split_and_tail_pair_shift_source_exact':str(2*dn+8*dU/3),
            'two_cell_source_shift_exact':'0',
            'implementation':'constants(*degree_coefficients()) applies the common bulk corrections before each meson diagonal is formed.'}


def matrices(cosine,sine,strange=False,repair=True):
    a,b,C,p,u=data(strange);z=sy.Rational(cosine)-sy.I*sy.Rational(sine)
    t=p['t'];A=sy.Matrix([[a,-t*(1+z),-t*(1+z)],[-t*(1+sy.conjugate(z)),b,0],[-t*(1+sy.conjugate(z)),0,b]])
    if repair:
        M=sy.Matrix([[a,-2*t,-2*t],[-2*t,b,0],[-2*t,0,b]])
        zeta=repair_matrix()[2]
        A+=sy.Rational(2*(1-F(cosine))*zeta)*M.inv()
    assert A==A.conjugate().T
    if sine:
        R,I=A.applyfunc(sy.re),A.applyfunc(sy.im);A=R.row_join(-I).col_join(I.row_join(R))
    return A,C,p,u


def ground(cosine,sine,strange=False,repair=True):
    A,C,p,u=matrices(cosine,sine,strange,repair)
    bounds=eigen_intervals(A,sy.eye(A.rows))
    lo,hi=bounds[0];assert hi<C
    # Complete two tails r<=-2 and r>=2. Each has endpoint coupling <=2t.
    error=4*p['t']**2/(C-hi)
    return {'cosine':str(cosine),'sine':str(sine),'strange':strange,'correlated_kinetic':repair,
       'full_ground_GeV_interval':[str((lo-error)*u),str(hi*u)],
       'full_ground_GeV_display':enclosure(((lo-error)*u,hi*u)),
       'complete_relative_tail_floor_source':str(C),'compact_ground_source_interval':[str(lo),str(hi)],
       'full_tail_Schur_error_source':str(error),'all_relative_positions_retained':True}


def native_two_cell_check(strange=False):
    import numpy as np
    from model.openwave_strange_composites.source import source,INDEX
    from model.openwave_hadron_refinement_flow.positivity import inputs
    p,_,u=inputs();dm,gs=mass_parameters();d=source('meson');lookup={tuple(z):i for i,z in enumerate(d['labels'])}
    # pi-: d anti-u; K-: s anti-u. Complete four position states at S=0.
    P=np.zeros((256,4));positions=((0,0),(0,1),(1,0),(1,1))
    for j,(x,y)in enumerate(positions):
        for sq,sa,sign in((0,1,1),(1,0,-1)):
            q=INDEX[x,int(strange),sq,1];a=INDEX[y,0,sa,0]
            P[lookup[q,a],j]=sign/np.sqrt(2)
    H=d['H0'].copy()
    import scipy.sparse as sp
    H+=sp.diags(float(dm)*d['Nstrange']+float(gs)*d['strange_pair'])
    c=F(3,4)*(p['Jw']+p['eta']);base=2*p['mbar']+2*p['t']+c*(1-int(strange))+(dm+gs if strange else 0)
    edge=base+2*p['A']/9+8*p['U']/3+F(3,2)*p['Js']+3*p['spin_gradient']
    expected=np.diag([float(base),float(edge),float(edge),float(base)])
    for i,j in((0,1),(0,2),(3,1),(3,2)):expected[i,j]=expected[j,i]=-float(p['t'])
    residual=np.max(abs(H@P-P@expected));assert residual<1e-9
    return {'strange':strange,'complete_native_spin_zero_position_rank':4,'intertwiner_residual_diagnostic':float(residual),'raw_meson_singlets_retained':256}
