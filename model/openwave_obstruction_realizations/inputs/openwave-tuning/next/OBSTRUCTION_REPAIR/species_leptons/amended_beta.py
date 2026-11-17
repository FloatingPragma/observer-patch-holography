"""Three-dimensional beta records of the current-square colour amendment.

No empirical hadron, axial, or neutrino parameters enter this source. One
massless neutral and an electron of unit rest mass are declared on its clock.
The original 432-state lepton/common-field transfer is not silently imported.
"""
import ast
from functools import lru_cache
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sy

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def main():
    current_path=HERE/'AMENDED_CURRENT.json'
    d=json.loads(current_path.read_text())
    vbox=list(map(F,d['vector_interval']));abox=list(map(F,d['axial_interval']))
    gbox=list(map(F,d['axial_over_vector_interval']));qbox=list(map(F,d['source_gap_exact']))
    assert F(999,1000)<vbox[0]<=vbox[1]<1
    assert F(166,100)<gbox[0]<=gbox[1]<F(167,100)
    assert 1<qbox[0]<qbox[1]<F(14,10)
    fv=d['vector_float'];g=d['ratio_float'];Q=float(sum(qbox,F())/2)
    # Reuse only the source-independent exact 64-state current theorem. No
    # calibrated source builder, current table or empirical input is imported.
    norm_path=ROOT/'reverse-engineering-reality/model/openwave_common_mass_source/source.py'
    tree=ast.parse(norm_path.read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='local_norm_certificate')
    namespace={'sy':sy,'product':product,'lru_cache':lru_cache}
    exec(compile(ast.Module(body=[fn],type_ignores=[]),str(norm_path),'exec'),namespace)
    norm=namespace['local_norm_certificate']()
    assert norm['exact_operator_norm']=='6'
    # Only one lepton family: six per quark, three quarks. Contact positions
    # are orthogonal blocks for each quark, so there is no extra factor two.
    full_norm=18
    G=F(1,64);t=F(1,100000);gt=G*t;vflo=F(999,1000);glo=F(166,100)
    rows=[]
    for L,w,chiral in product((4,8),(0,1),(1,-1,0)):
        th=np.array(list(product(2*np.pi*np.arange(L)/L,repeat=3)));N=len(th)
        p=2*np.sin(th);p[abs(p)<1e-14]=0
        wmass=w*np.sum(1-np.cos(th),axis=1)
        me=1+wmass;mn=wmass
        Ee=np.sqrt((p*p).sum(axis=1)+me*me);En=np.sqrt((p*p).sum(axis=1)+mn*mn)
        ae=p/Ee[:,None];bn=-np.divide(p,En[:,None],out=np.zeros_like(p),where=En[:,None]>0)
        ve=(p*2*np.cos(th)+me[:,None]*w*np.sin(th))/Ee[:,None]
        vh=np.divide(p*2*np.cos(th)+mn[:,None]*w*np.sin(th),En[:,None],out=np.zeros_like(p),where=En[:,None]>0)
        ve[abs(ve)<1e-13]=0;vh[abs(vh)<1e-13]=0
        form=np.cos((th[:,0,None]-th[None,:,0])/2)**2
        spin=fv*fv*(1+3*g*g+(1-g*g)*(ae@bn.T)+chiral*(2*g*(1-g)*ae[:,2,None]+2*g*(1+g)*bn[None,:,2]))/2
        assert spin.min()>-1e-12 and spin.max()<16
        mismatch=Ee[:,None]+En[None,:]-Q;assert np.max(abs(mismatch))<14
        pr=(float(gt)/N)**2*form*spin*np.sinc(mismatch*float(t)/(2*np.pi))**2
        pe=pr.sum(axis=1);pn=pr.sum(axis=0);total=float(pr.sum())
        odd_e=float(pe@np.sign(ve[:,2]));odd_n=float(pn@np.sign(-vh[:,2]))
        delta=(full_norm*gt)**2/2;err=(8*gt+delta)*delta
        finite_time_error=16*gt**2*(14*t)**2/12
        emin=gt**2*vflo**2*glo*(glo-1)*F(2,7*L)-finite_time_error-err
        nmin=gt**2*vflo**2*glo*(glo+1)*F(2,7*L)-finite_time_error-err
        assert emin>0 and nmin>0
        if w==0 or chiral==0:
            assert abs(odd_e)<1e-27 and abs(odd_n)<1e-27
        else:
            assert -chiral*odd_e>float(emin) and chiral*odd_n>float(nmin)
        rows.append(dict(cube_side=L,wilson=w,chirality_control=chiral,
            chirality_control_meaning='left' if chiral==1 else 'right' if chiral==-1 else 'equal classical mixture of left/right; fair coin retained',
            first_order_selected_probability=total,electron_velocity_odd=odd_e,
            antineutrino_velocity_odd=odd_n,electron_conditional_velocity_asymmetry=odd_e/total,
            antineutrino_conditional_velocity_asymmetry=odd_n/total,
            original_preparation_failure='119/120',complete_quantum_complement_retained=True,
            exact_original_input_electron_negative_margin=str(emin/120) if w and chiral==1 else None,
            exact_original_input_antineutrino_positive_margin=str(nmin/120) if w and chiral==1 else None,
            full_Fock_selected_probability_error_upper=str(err),
            finite_time_weight_error_upper=str(finite_time_error)))
    # A definite nonzero on-shell current point in the continuous momentum
    # family: electron p=0; neutral along y, q_x=q_z=0. E_nu(q)=Q-1.
    lo,hi=0.,np.pi/2
    for _ in range(70):
        mid=(lo+hi)/2;energy=np.hypot(2*np.sin(mid),1-np.cos(mid))
        if energy<Q-1:lo=mid
        else:hi=mid
    q=(lo+hi)/2
    surface=dict(electron_momentum=[0,0,0],neutral_momentum=[0,q,0],
        energy_error=float(abs(1+np.hypot(2*np.sin(q),1-np.cos(q))-Q)),
        spin_weight_at_point=(fv*fv+3*(g*fv)**2)/2,
        continuous_momentum_family_only=True,
        statement='Nonzero current on a regular five-dimensional free energy shell; positive leading-order spectral density. This is not a full-interacting lifetime theorem or a finite-cube on-shell decay claim.')
    assert surface['energy_error']<1e-14 and surface['spin_weight_at_point']>0
    out=dict(schema='oph.current_square_3d_beta.v1',current_sha256=hashlib.sha256(current_path.read_bytes()).hexdigest(),
        norm_source_sha256=hashlib.sha256(norm_path.read_bytes()).hexdigest(),
        declared_source=d['source']+' + one electron/neutral family, 3D Wilson coefficient1, local native weak current G=1/64',
        empirical_hadron_or_axial_or_neutrino_inputs=[],electron_mass=1,neutral_mass=0,
        neutral_zero_mode_reference='Chosen negative beta eigenvectors occupied at zero energy; both active spin weights1/2, all zero-mode complements retained.',
        source_time=str(t),weak_coefficient=str(G),full_Fock_current_norm_upper=full_norm,
        reused_exact_one_quark_one_family_norm=norm,source_Q_interval=d['source_gap_exact'],
        vector_interval=d['vector_interval'],axial_interval=d['axial_interval'],
        records=rows,free_onshell_spectral_density=surface,
        complete_weak_native_execution_claimed=False,physical_beta_row_solved=False,
        missing=['native preparation of the amended ground states','complete common-source compiler','shared physical clock and mass mapping','physical mass splitting and axial comparison','controlled long-time interacting decay'])
    (HERE/'AMENDED_BETA.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'source_Q':Q,'fV':fv,'gA':g,'records':len(rows),
        'positive_density_point':surface,'L8_left':[r for r in rows if r['cube_side']==8 and r['wilson']==1 and r['chirality_control']==1][0]},indent=2))


if __name__=='__main__':main()
