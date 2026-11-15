"""Supplied spectral record energies on the implemented finite integer sector.

J03's static identity is audited separately from its dynamical use. On transit
states E=(F_a-F_b)/2 is a declared half-record extension. Adding this diagonal
energy changes the generator; its massive-square identity is not retained.
"""
import math
import numpy as np
from sympy import Matrix, Rational
from .construction import Construction
from .quantum_records import JointRecords


class RecordEnergy:
    def __init__(self,construction,kind='first_order',mass=0.):
        if kind not in ('first_order','quadratic','massive'):raise ValueError('unknown record energy')
        if not math.isfinite(mass) or mass<0:raise ValueError('invalid energy mass')
        self.construction=construction; self.kind=kind; self.mass=mass
        self.b=construction.incidence
        u,s,vt=np.linalg.svd(self.b,full_matrices=True)
        tolerance=max(self.b.shape)*np.finfo(float).eps*(s[0] if len(s) else 1)
        rank=int(np.count_nonzero(s>tolerance))
        self.gradient=vt[:rank].T@vt[:rank]
        self.cycle=np.eye(self.b.shape[1])-self.gradient
        nonzero=s[:rank]**2
        phi=np.ones(rank) if kind=='quadratic' else np.sqrt(nonzero+(mass**2 if kind=='massive' else 0))
        self.phi_zero=1. if kind=='quadratic' else mass if kind=='massive' else 0.
        self.kernel=(vt[:rank].T*phi)@vt[:rank]+self.phi_zero*self.cycle
        self.charge_kernel=(u[:,:rank]*(phi/nonzero))@u[:,:rank].T

    def components(self,records):
        e=np.asarray(records,float)
        if e.shape!=(self.b.shape[1],) or not np.all(np.isfinite(e)):raise ValueError('invalid records')
        rho=self.b@e
        total=float(e@self.kernel@e/2)
        charge=float(rho@self.charge_kernel@rho/2)
        loop=float(self.phi_zero*np.linalg.norm(self.cycle@e)**2/2)
        return {'total':total,'charge':charge,'loop':loop,'identity_error':abs(total-charge-loop)}

    def joint_diagonal(self,joint):
        if joint.construction!=self.construction:raise ValueError('construction mismatch')
        return np.array([self.components((np.asarray(f[::2])-f[1::2])/2)['total']
                         for _,f in joint.states])

    def augment(self,joint,coupling):
        """Explicit supplied diagonal interaction, with unchanged Gauss sector."""
        if not math.isfinite(coupling) or coupling<0:raise ValueError('invalid coupling')
        return joint.hamiltonian+coupling*np.diag(self.joint_diagonal(joint))


def run():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    first=RecordEnergy(c);quadratic=RecordEnergy(c,'quadratic');massive=RecordEnergy(c,'massive',.23)
    strings=[np.array([1,0,0])+k*np.ones(3,dtype=int) for k in range(-2,3)]
    values={name:[op.components(e) for e in strings] for name,op in
            (('first_order',first),('quadratic',quadratic),('massive',massive))}
    # Independent exact rational solution of the grounded graph system.
    b=Matrix(c.incidence.astype(int).tolist());lap=b*b.T
    rho=Matrix([1,-1,0]);potential=Matrix([*(lap[:2,:2].inv()*rho[:2,0]),0])
    real_min=(rho.T*potential)[0]/2
    integer_min=min(sum(int(x)**2 for x in e)/2 for e in strings)
    j=JointRecords(c,cutoff=1);op=massive.augment(j,.29)
    diagonal=massive.joint_diagonal(j)
    base=j.hamiltonian
    square_change=float(np.max(np.abs(op@op-base@base)))
    from .record_preparation import finite_evolve
    initial=np.eye(len(j.states))[0]
    free_state=finite_evolve(base,initial,.37)
    interacting_state=finite_evolve(op,initial,.37)
    def populations(state):
        return [float(sum(abs(state[i])**2 for i,s in enumerate(j.states) if s[0]==pos))
                for pos in range(c.ports+len(c.seams))]
    free_pop=populations(free_state);interacting_pop=populations(interacting_state)
    observable_change=float(np.max(np.abs(np.array(free_pop)-interacting_pop)))
    # Check transit extension against exact integer Gauss charges.
    g=np.array(j.gauss(j.states[0])[:c.ports])
    half_record_error=0.
    for pos,f in j.states:
        e=(np.asarray(f[::2])-f[1::2])/2
        effective=g-np.array([int(pos==x) for x in range(c.ports)],float)
        if pos>=c.ports:
            a,d=c.seams[pos-c.ports];effective[a]-=.5;effective[d]-=.5
        half_record_error=max(half_record_error,float(np.max(np.abs(c.incidence@e-effective))))
    return {'scope':'Supplied spectral energy on finite integer current records; triangle static identity and declared transit extension. No dynamical Coulomb or confinement claim.',
            'checks':{'spectral_charge_cycle_identity':max(v['identity_error'] for rows in values.values() for v in rows)<1e-12,
                      'first_order_loop_independent':bool(np.ptp([v['total'] for v in values['first_order']])<1e-12),
                      'exact_real_minimum':real_min==Rational(1,3),
                      'integer_minimum':integer_min==.5,
                      'transit_half_record_gauss':half_record_error==0,
                      'augmented_hermitian':bool(np.allclose(op,op.T)),
                      'augmented_probability_normalized':abs(sum(interacting_pop)-1)<1e-12},
            'mutations':{'quadratic_distinguishes_loops':bool(np.ptp([v['total'] for v in values['quadratic']])>1),
                         'massive_adds_loop_cost':bool(np.ptp([v['loop'] for v in values['massive']])>.1),
                         'bare_massive_square_transfer_rejected':square_change>.01},
            'full_criterion_met':False,
            'supplied':{'graph':'oriented triangle','mass':.23,'cutoff':1,'energy_coupling':.29,
                        'energy_choice':'massive spectral record energy',
                        'transit_extension':'E_e=(F_ea-F_eb)/2; effective charge includes half a unit at each endpoint of occupied seam'},
            'measurements':{'unit_dipole_real_quadratic_minimum':str(real_min),
                            'unit_dipole_integer_quadratic_minimum':integer_min,
                            'integer_cycle_energy_polynomial':'(3*k^2+2*k+1)/2',
                            'real_cycle_minimizer':'-1/3',
                            'same_charge_loop_shifts':[-2,-1,0,1,2],
                            'energy_components':values,'augmented_square_change':square_change,
                            'joint_diagonal_range':[float(diagonal.min()),float(diagonal.max())],
                            'evolution_time':.37,'free_position_probabilities':free_pop,
                            'augmented_position_probabilities':interacting_pop,
                            'maximum_position_probability_change':observable_change},
            'physical_rows':{'20':'static energy choice explicit; no force-law/strength validation',
                             '22':'integer-vs-real cost diagnostic; no physical string tension or string breaking',
                             '30':'interaction changes H; free massive-square cannot be reused'}}
