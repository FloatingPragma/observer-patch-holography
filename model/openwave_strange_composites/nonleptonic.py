"""Native nonleptonic vertex and a complete same-source closed-channel bound."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import comb
import numpy as np
from .source import source,select,INDEX,native,parent
from .current import mass_parameters


def word(mask,operators):
    sign=1
    for create,i in reversed(operators):
        if bool(mask>>i&1)==create:return None,0
        sign*=(-1)**((mask&((1<<i)-1)).bit_count())
        mask=mask|(1<<i) if create else mask^(1<<i)
    return mask,sign


def raw_state(kind,state):
    d=source(kind);v=d['B']@state
    return {mask:complex(z)for mask,z in zip(d['full'],v)if z!=0}


def wedge(left,right):
    out=defaultdict(complex)
    for a,u in left.items():
        ops=[(True,i)for i in range(96)if a>>i&1]
        for b,v in right.items():
            mask,sign=word(b,ops)
            if sign:out[mask]+=sign*u*v
    return dict(out)


def certificate():
    dm,pair=mass_parameters();unit=F(parent.INPUTS['electron_unit_GeV']);p=parent.parameters()
    Lam=select('baryon',(1,1,0,1),(0,0,0,0),.5,0,float(dm),float(pair))
    pro=select('baryon',(2,1,0,0),(0,0,0,0),.5,.5,float(dm),float(pair))
    pion=select('meson',(0,1,0,0),(1,0,0,0),0,.5,float(dm),float(pair))
    initial=raw_state('baryon',Lam['state']);target=wedge(raw_state('baryon',pro['state']),raw_state('meson',pion['state']))
    norm=sum(abs(v)**2 for v in target.values());assert norm>0
    target={k:v/np.sqrt(norm)for k,v in target.items()}
    pauli=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.])]
    bar=[np.eye(2),*(-x for x in pauli)];eta=(1,-1,-1,-1);epsilon=np.array([[0,1],[-1,0]])
    s=F('0.22431');cab=float(s)*np.sqrt(1-float(s*s))/2
    joins=native()[2];out=defaultdict(complex)
    for x,su,ss,sd,sa,c,d,dd in product(range(2),range(2),range(2),range(2),range(2),range(3),range(3),range(3)):
        spin=sum(eta[mu]*bar[mu][su,ss]*(bar[mu]@epsilon)[sd,sa]for mu in range(4))
        assert spin in(0,2,-2)
        coefficient=cab*spin*joins[0][d,dd]
        if coefficient==0:continue
        up=3*INDEX[x,0,su,0]+c;strange=3*INDEX[x,1,ss,1]+c
        down=3*INDEX[x,0,sd,1]+d;antiup=48+3*INDEX[x,0,sa,0]+dd
        ops=[(True,up),(False,strange),(True,down),(True,antiup)]
        for mask,value in initial.items():
            final,sign=word(mask,ops)
            if sign:out[final]+=coefficient*sign*value
    phi={k:v for k,v in out.items()if v!=0}
    amplitude=sum(target.get(k,0).conjugate()*v for k,v in phi.items())
    phi2=sum(abs(v)**2 for v in phi.values());assert phi2>0 and abs(amplitude)>1e-6
    # The complete two-cell Nq4/Na1 source retains its positive colour,
    # kinetic and spin terms. K>=nq*na-min(nq,na)=3, net N=3 density>=5.
    c=F(3,4)*(p['Jw']+p['eta'])
    floor=5*p['mbar']-p['delta']/2+5*p['A']/9+3*p['pair_triplet']-c
    target_lambda=F('1.115683')/unit
    assert floor>target_lambda
    # An exact native co-located Lambda trial: the light ud spin/isospin
    # singlet and a spin-up strange quark. Its two orbital coordinates have
    # coefficients +1 and -1; the symmetric-cube colour join is unchanged.
    baryon=source('baryon');lookup={tuple(z):i for i,z in enumerate(baryon['labels'])}
    trial=np.zeros(len(baryon['labels']),complex)
    trial[lookup[tuple(sorted((INDEX[0,0,0,0],INDEX[0,0,1,1],INDEX[0,1,0,1])))]]=1
    trial[lookup[tuple(sorted((INDEX[0,0,1,0],INDEX[0,0,0,1],INDEX[0,1,0,1])))]]=-1
    assert np.linalg.norm(baryon['S']@trial-float(F(3,4))*trial)<1e-12
    assert np.linalg.norm(baryon['I']@trial)<1e-12
    upper=3*p['mbar']+dm+p['delta']/2+p['A']+F(3,4)*(p['Js']+p['spin_gradient'])+3*p['t']
    rayleigh=np.vdot(trial,baryon['H0']@trial+float(dm)*baryon['Nstrange']*trial).real/2
    assert abs(rayleigh-float(upper))<2e-9
    assert floor>upper>target_lambda
    # Full unfiltered native weak-current square on48 left quark modes:
    # each local current norm<=12, two cells, four Lorentz components and
    # three weak axes, coefficient2. The canonical Nambu transformation is
    # unitary on the complete96-mode CAR carrier.
    weak_norm=F(6912);G=F(1,4096)
    assert floor-upper>2*G*weak_norm
    dimension=sum(comb(6,u)*comb(6,2-u)*comb(6,d)*comb(6,2-d)*6
                  for u,d,a in product(range(3),range(3),range(2))if u+d+a==3)
    return {'native_N5_vertex_norm_squared':float(phi2),
            'normalized_proton_pion_packet_vertex_amplitude':[float(amplitude.real),float(amplitude.imag)],
            'complete_native_N5_flavour_spin_occupation_dimension':dimension,
            'antisymmetrized_proton_pion_packet_norm_squared':float(norm),
            'weak_coefficient_source_units':str(G),'full_unfiltered_weak_norm_upper':str(weak_norm),
            'strong_N4q_1antiquark_channel_floor_GeV':str(floor*unit),
            'strong_Lambda_to_channel_gap_GeV':str((floor-target_lambda)*unit),
            'strong_Lambda_ground_variational_upper_GeV':str(upper*unit),
            'rigorous_channel_gap_above_Lambda_ground_GeV':str((floor-upper)*unit),
            'gap_after_two_full_weak_norm_shifts_GeV':str((floor-upper-2*G*weak_norm)*unit),
            'fitted_Lambda_gap_role':'The Lambda-to-channel comparison uses the mass calibration target; the separate exact trial upper bound proves the closed-channel conclusion independently of fit rounding.',
            'strong_channel_is_kinematically_open':False,
            'outgoing_channel_counts':{'u_quarks':2,'d_quarks':2,'anti_u':1,'total_charge':0,'total_baryon_number':1},
            'source':'The full native left-quark current square2 sum_x,mu,a eta_mu J_x,mu,a² under psi_L=(b+epsilon a_dagger)/sqrt2; its N5 projection is computed from the literal CAR word, while every other-number channel is retained by the full operator norm.',
            'interpretation':'The nonzero vertex is a virtual finite-channel response. The same-source five-body floor excludes an on-shell two-cell Lambda→proton+pion decay for the declared weak coefficient. No lifetime or physical strange-decay validation is claimed.',
            'energy':'The autonomous full Hamiltonian conserves its exact total energy. A proton/pion packet filter generally performs work; it cannot convert the closed channel into an energy-conserving outgoing decay.'}


if __name__=='__main__':
    import json
    print(json.dumps(certificate(),indent=2),flush=True)
