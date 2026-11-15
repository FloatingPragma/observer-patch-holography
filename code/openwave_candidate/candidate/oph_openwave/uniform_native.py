"""Exact local-coarse reduction for a separately declared uniform-record preparation.

No definite-record substitution: V|z> uniformly superposes the complete modular
Gauss fiber at z. Kraus intertwiners prove the reduction, including references.
The canonical reduced execution stores 87 amplitudes rather than 87*N^40.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

PACKAGE=Path(__file__).resolve().parent
from .construction import Construction
from .geometry import flat_ring
from .modular_native import LocalCoarseModular
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator,expm_multiply


class UniformFiberReduction:
    """Bare position channel exact on the uniform modular fiber isometry image."""
    def __init__(self,construction,modulus=6):
        if type(modulus) is not int or modulus<2:raise ValueError('integer modulus>=2 required')
        if not construction.seams:raise ValueError('native event instrument requires at least one seam')
        self.construction=construction;self.modulus=modulus
        self.n=construction.ports+len(construction.seams)
        adjacency=[[] for _ in range(self.n)]
        self.half_edges=[]
        for edge,(a,b) in enumerate(construction.seams):
            for side,port in enumerate((a,b)):
                midpoint=construction.ports+edge;idx=len(self.half_edges)
                self.half_edges.append((port,midpoint))
                adjacency[port].append((midpoint,idx));adjacency[midpoint].append((port,idx))
        seen={0};queue=[0]
        for vertex in queue:
            for neighbor,_ in adjacency[vertex]:
                if neighbor not in seen:seen.add(neighbor);queue.append(neighbor)
        if len(seen)!=self.n:raise ValueError('connected graph required for uniform fibers at every position')
        self.adjacency=adjacency
        self.cycle_rank=len(construction.seams)-construction.ports+1
        self.fiber_cardinality=modulus**self.cycle_rank
        self.hamiltonian=construction.hamiltonian
        self.maps=[]
        for seam,(a,b) in enumerate(construction.seams):
            grouped={}
            for law in ('C','E'):
                for coin in (0,1):
                    for pos in range(self.n):
                        contrast=int(pos==a)-int(pos==b)
                        current=(contrast-(2*coin-1)*(abs(contrast)%2))//2
                        destination=b if current==1 else a if current==-1 else pos
                        value=current if law=='C' else contrast
                        grouped.setdefault((law,coin,value),[]).append((pos,destination))
            operations=[]
            for key,pairs in sorted(grouped.items()):
                src,dst=map(np.array,zip(*pairs))
                if len(set(dst))!=len(dst):raise ValueError('noninjective local outcome')
                operations.append(({'seam':seam,'law':key[0],'coin':key[1],'record_value':key[2]},src,dst))
            self.maps.append(operations)
        self.eigenvalues,self.eigenvectors=np.linalg.eigh(self.hamiltonian)

    def initial_density(self):
        """Bare port0 represents the UNIFORM record fiber at port0, not F=0."""
        result=np.zeros((self.n,self.n),complex);result[0,0]=1;return result

    def density(self,rho):
        rho=np.asarray(rho,complex)
        if rho.shape!=(self.n,self.n) or not np.all(np.isfinite(rho)) or not np.allclose(rho,rho.conj().T) or not np.isclose(np.trace(rho),1) or np.linalg.eigvalsh(rho).min()<-1e-10:
            raise ValueError('finite normalized positive density required')
        return rho

    def kraus_for(self,seam):
        result=[]
        for record,src,dst in self.maps[seam]:
            k=np.zeros((self.n,self.n));k[dst,src]=.5;result.append((record,k))
        return result

    def _seam_output(self,rho,seam):
        output=np.zeros_like(rho);ledger=[]
        for record,src,dst in self.maps[seam]:
            branch=rho[np.ix_(src,src)]/4
            output[np.ix_(dst,dst)]+=branch
            ledger.append({**record,'probability':float(np.trace(branch).real)})
        return output,ledger

    def commit(self,rho,seam=None):
        rho=self.density(rho)
        seams=range(len(self.maps)) if seam is None else (seam,)
        if not len(seams) or any(type(e) is not int or not 0<=e<len(self.maps) for e in seams):raise ValueError('invalid seam')
        output=np.zeros_like(rho);ledger=[]
        for e in seams:
            branch,records=self._seam_output(rho,e);output+=branch/len(seams)
            for record in records:record['probability']/=len(seams)
            ledger.extend(records)
        return output,ledger

    def step(self,rho,time,probability):
        if not math.isfinite(time) or not math.isfinite(probability) or not 0<=probability<=1:raise ValueError('invalid split step')
        rho=self.density(rho)
        v=self.eigenvectors;u=(v*np.exp(-1j*time*self.eigenvalues))@v.T
        evolved=u@rho@u.conj().T
        committed,ledger=self.commit(evolved)
        for record in ledger:record['probability']*=probability
        ledger.append({'status':'no_commit','probability':1-probability})
        return (1-probability)*evolved+probability*committed,ledger

    def generator_action(self,rho,rate):
        """Exact per-seam GKSL action; no large Liouville matrix assembled."""
        if not math.isfinite(rate) or rate<0:raise ValueError('invalid per-seam rate')
        rho=self.density(rho)
        result=-1j*(self.hamiltonian@rho-rho@self.hamiltonian)
        for seam in range(len(self.maps)):result+=rate*(self._seam_output(rho,seam)[0]-rho)
        return result

    def half_record_marginals(self,rho):
        """Readout probabilities of F_h moduloN on V rho V†, without measuring it.

        A cycle-edge coordinate is uniform; a bridge coordinate is fixed by
        charge on its cut. Executing record measurement generally leaves V's image.
        """
        probabilities=np.diag(self.density(rho)).real;result=[]
        for edge,(tail,head) in enumerate(self.half_edges):
            component={tail};queue=[tail]
            for vertex in queue:
                for neighbor,index in self.adjacency[vertex]:
                    if index!=edge and neighbor not in component:component.add(neighbor);queue.append(neighbor)
            if head in component:
                marginal=np.ones(self.modulus)/self.modulus;kind='cycle_uniform'
            else:
                marginal=np.zeros(self.modulus)
                initial=int(0 in component)
                for pos,p in enumerate(probabilities):marginal[(initial-int(pos in component))%self.modulus]+=p
                kind='bridge_cut_charge'
            result.append({'half_seam':edge,'type':kind,'residue_probabilities':marginal.tolist()})
        return result

    def boundary_charge_readout(self,rho,vertices):
        """Joint net-boundary residue and whole-load regional charge, fixed G."""
        vertices=frozenset(vertices)
        if any(type(v) is not int or not 0<=v<self.n for v in vertices):raise ValueError('invalid region')
        probabilities=np.diag(self.density(rho)).real;g=int(0 in vertices)
        inside=float(sum(p for pos,p in enumerate(probabilities) if pos in vertices))
        return [{'regional_charge':q,'net_boundary_record_residue':(g-q)%self.modulus,
                 'regional_gauss_residue':g,'probability':inside if q else 1-inside} for q in (0,1)]


    def per_seam_operator(self,rate):
        """Matrix-free GKSL forward/adjoint and analytic Liouville supertrace."""
        if not math.isfinite(rate) or rate<0:raise ValueError('finite nonnegative per-seam rate required')
        n=self.n;h=csr_matrix(self.hamiltonian)
        def action(vector,adjoint=False):
            matrix=np.asarray(vector).reshape((n,n),order='F')
            commutator=h@matrix-(h.T@matrix.T).T
            result=(1j if adjoint else -1j)*commutator
            for operations in self.maps:
                channel=np.zeros((n,n),complex)
                for _,source,target in operations:
                    if adjoint:channel[np.ix_(source,source)]+=matrix[np.ix_(target,target)]/4
                    else:channel[np.ix_(target,target)]+=matrix[np.ix_(source,source)]/4
                result+=rate*(channel-matrix)
            return result.reshape(-1,order='F')
        trace=rate*sum(sum((np.count_nonzero(source==target)/2)**2 for _,source,target in operations)-n*n
                       for operations in self.maps)
        return LinearOperator((n*n,n*n),matvec=action,rmatvec=lambda vector:action(vector,True),dtype=complex),float(trace)

    def evolve_per_seam(self,rho,time,rate):
        """Numerical exponential action of the exact reduced finite generator.

        No dense Liouville matrix or certified numerical-error bound is used.
        """
        if not math.isfinite(time) or time<0:raise ValueError('finite nonnegative duration required')
        rho=self.density(rho);generator,trace=self.per_seam_operator(rate)
        return expm_multiply(time*generator,rho.reshape(-1,order='F'),traceA=time*trace).reshape((self.n,self.n),order='F')


def tiny_intertwiner_audit(modulus):
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    full=LocalCoarseModular(c,modulus);reduced=UniformFiberReduction(c,modulus)
    v=np.zeros((full.n,reduced.n))
    counts=[]
    for pos in range(reduced.n):
        indices=[i for i,s in enumerate(full.states) if s[0]==pos];counts.append(len(indices))
        v[indices,pos]=1/np.sqrt(len(indices))
    errors=[]
    for seam in range(3):
        bare={tuple(sorted(record.items())):k for record,k in reduced.kraus_for(seam)}
        for record,k in full.operators[seam]:errors.append(float(np.max(np.abs(k@v-v@bare[tuple(sorted(record.items()))]))))
    rho=reduced.initial_density();lifted=v@rho@v.T;definite=full.initial_density()
    for _ in range(8):
        rho,_=reduced.step(rho,.37,.29);lifted,_=full.step(lifted,.37,.29)
        definite,_=full.step(definite,.37,.29)
    fixed_gauss=full.gauss_residue(full.states[0])
    # Fully source-monitored outcome must resolve an individual joint basis state.
    monitored,_=full.commit(v@reduced.initial_density()@v.T,fully_monitored_control=True)
    code_projection=v@v.T
    leakage=float(np.trace(monitored@(np.eye(full.n)-code_projection)).real)
    diagonal=np.diag(lifted).real
    actual_marginal=[float(sum(diagonal[i] for i,s in enumerate(full.states) if s[1][0]==residue)) for residue in range(modulus)]
    expected_marginal=reduced.half_record_marginals(rho)[0]['residue_probabilities']
    return {'modulus':modulus,'full_states':full.n,'fiber_counts':counts,
            'predicted_fiber_cardinality':reduced.fiber_cardinality,
            'isometry_error':float(np.max(np.abs(v.T@v-np.eye(reduced.n)))),
            'H_intertwining_error':float(np.max(np.abs(full.hamiltonian@v-v@reduced.hamiltonian))),
            'every_coarse_Kraus_intertwining_error':max(errors),
            'eight_step_lifted_density_error':float(np.max(np.abs(lifted-v@rho@v.T))),
            'definite_record_substitution_density_error':float(np.max(np.abs(definite-v@rho@v.T))),
            'global_source_monitor_code_leakage':leakage,
            'half_record_readout_error':max(abs(a-b) for a,b in zip(actual_marginal,expected_marginal)),
            'fixed_modular_gauss':all(full.gauss_residue(s)==fixed_gauss for s in full.states)}


def run(full_ledger=False):
    tiny=[tiny_intertwiner_audit(n) for n in (3,6)]
    ring=UniformFiberReduction(flat_ring(1,1,.23),6)
    rho=ring.initial_density();initial_mean_energy=float(np.trace(rho@ring.hamiltonian).real)
    step_ledgers=[]
    for step in range(8):
        rho,ledger=ring.step(rho,.37,.29)
        if full_ledger:step_ledgers.append({'step':step+1,'outcomes':ledger})
    boundary=ring.boundary_charge_readout(rho,range(12))
    derivative=ring.generator_action(rho,.29/(63*.37))
    rate=.29/(3*.37)
    continuous=ring.evolve_per_seam(ring.initial_density(),.1,rate)
    result={'schema':'oph.uniform_modular_fiber_reduction.v1','preparation_variant':'uniform_modular_Gauss_fiber',
            'scope':'Exact coarse-instrument reduction under explicitly different uniform-record preparation; no definite-record substitution or physical criterion.',
            'checks':{'equal_composite_modulus_fibers':all(len(set(r['fiber_counts']))==1 and r['fiber_counts'][0]==r['predicted_fiber_cardinality'] for r in tiny),
                      'every_Kraus_intertwines':all(r['every_coarse_Kraus_intertwining_error']<1e-12 for r in tiny),
                      'whole_instrument_repeated_evolution':all(r['eight_step_lifted_density_error']<1e-12 for r in tiny),
                      'record_marginal_readout':all(r['half_record_readout_error']<1e-12 for r in tiny),
                      'ring_trace_positive':bool(abs(np.trace(rho)-1)<1e-12 and np.linalg.eigvalsh(rho).min()>-1e-12),
                      'per_seam_generator_trace':bool(abs(np.trace(derivative))<1e-12),
                      'same_local_rate_ring_trace_positive':bool(abs(np.trace(continuous)-1)<1e-12 and np.linalg.eigvalsh(continuous).min()>-1e-12)},
            'mutations':{'definite_preparation_substitution_rejected':all(r['definite_record_substitution_density_error']>.01 for r in tiny),
                         'global_full_source_monitor_rejected':all(r['global_source_monitor_code_leakage']>.1 for r in tiny)},
            'full_criterion_met':False,'tiny_full_sector_audits':tiny,
            'per_seam_ring':{'geometry':'flat_ring(1,1)','dimension':87,'per_seam_rate':rate,'duration':.1,
                             'dense_liouville_assembled':False,'trace':float(np.trace(continuous).real),
                             'minimum_eigenvalue':float(np.linalg.eigvalsh(continuous).min()),
                             'position_probabilities':np.diag(continuous).real.tolist(),
                             'numerical_error':'expm_multiply numerical action; no certified error claimed'},
            'canonical_ring':{'geometry':'flat_ring(1,1)','ports':24,'seams':63,'reduced_dimension':ring.n,
                              'cycle_rank':ring.cycle_rank,'fiber_cardinality':str(ring.fiber_cardinality),
                              'formal_full_sector_dimension':str(ring.n*ring.fiber_cardinality),
                              'final_position_probabilities':np.diag(rho).real.tolist(),
                              'port_probability':float(np.diag(rho).real[:24].sum()),
                              'transit_probability':float(np.diag(rho).real[24:].sum()),
                              'final_purity':float(np.trace(rho@rho).real),
                              'boundary_readout_first_carrier_ports':boundary,
                              'ledger_probability':sum(row['probability'] for row in ledger),
                              'last_step_outcome_count':len(ledger),
                              'ledger_scope':'step-marginal local records; joint multistep history not reconstructed',
                              'initial_mean_H_energy':initial_mean_energy,
                              'final_mean_H_energy':float(np.trace(rho@ring.hamiltonian).real)},
            'supplied':{'preparation':'uniform pure modular fixed-Gauss fiber at initial port0, DISTINCT from definite F=0',
                        'modulus':6,'mass':.23,'time':.37,'split_probability':.29,'steps':8,
                        'record_energy':'none','instrument':'local coarse C/E equal mixture, transit elsewhere',
                        'ring_split_equivalent_per_seam_rate':.29/(63*.37),
                        'continuous_per_seam_rate':rate,'continuous_duration':.1,
                        'physical_charge_unit_mapping':None,'physical_species_mapping':None},
            'sources':{str(p.relative_to(PACKAGE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       (Path(__file__).resolve(),PACKAGE/'modular_native.py',PACKAGE/'geometry.py',PACKAGE/'construction.py')}}
    if full_ledger:result['complete_step_marginal_ledgers']=step_ledgers
    return result



if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--full-ledger',action='store_true');args=p.parse_args();result=run(args.full_ledger)
    if not all(result['checks'].values()) or not all(result['mutations'].values()):raise SystemExit('uniform reduction failed')
    payload=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(payload)
    print(payload)
