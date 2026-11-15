"""Complete regional charge apparatus on the bounded integer-record instrument."""
from pathlib import Path
import argparse
import json
import operator
import numpy as np

from .construction import Construction
from .quantum_records import JointRecords
from .interacting_records import InteractingRecords
from .observer_protocol import ObserverProtocol


class RegionalCharge:
    """A region is a supplied subset of vertices of the subdivided graph.

    Transit seams are assigned as whole vertices. The half-record effective
    charge used by RecordEnergy is a different observable.
    """
    def __init__(self,joint,ports=(),transit_seams=()):
        self.joint=joint;c=joint.construction
        try:
            ports=tuple(operator.index(v) for v in ports)
            seams=tuple(operator.index(v) for v in transit_seams)
        except TypeError as exc:raise ValueError('integer regional vertices required') from exc
        if len(set(ports))!=len(ports) or len(set(seams))!=len(seams):raise ValueError('duplicate regional vertex')
        if any(not 0<=v<c.ports for v in ports) or any(not 0<=v<len(c.seams) for v in seams):raise ValueError('invalid regional vertex')
        self.vertices=frozenset((*ports,*(c.ports+v for v in seams)))
        self.charge=np.array([int(s[0] in self.vertices) for s in joint.states],dtype=int)
        self.transit=np.array([int(s[0]>=c.ports) for s in joint.states],dtype=int)
        # C is incidence of half-seams oriented from port to seam midpoint.
        incidence=np.zeros((c.ports+len(c.seams),2*len(c.seams)),dtype=int)
        for e,(a,b) in enumerate(c.seams):
            incidence[a,2*e]=1;incidence[b,2*e+1]=1
            incidence[c.ports+e,2*e:2*e+2]=-1
        indicator=np.array([int(v in self.vertices) for v in range(len(incidence))])
        self.boundary_coefficients=indicator@incidence
        self.boundary=np.array([sum(int(a)*int(b) for a,b in zip(self.boundary_coefficients,s[1]))
                                for s in joint.states],dtype=int)
        self.gauss=np.array([sum(joint.gauss(s)[v] for v in self.vertices) for s in joint.states])

    def instrument(self,rho,refine_transit=False):
        """Complete projective instrument, retaining impossible outcomes too."""
        rho=self.joint._density(rho)
        result=np.zeros_like(rho);outcomes=[]
        for charge in (0,1):
            for transit in ((0,1) if refine_transit else (None,)):
                mask=self.charge==charge
                if transit is not None:mask=mask & (self.transit==transit)
                branch=rho*mask[:,None]*mask[None,:]
                probability=float(np.trace(branch).real)
                result+=branch
                outcomes.append({'sharp_charge':charge,'transit':transit,
                                 'probability':probability,
                                 'conditional_density':branch/probability if probability>0 else None,
                                 'unnormalized_density':branch})
        return result,outcomes

    def readout(self,rho):
        rho=self.joint._density(rho);p=np.diag(rho).real
        mean=float(p@self.charge)
        census={}
        for weight,q,b,g in zip(p,self.charge,self.boundary,self.gauss):
            key=(int(q),int(b),int(g))
            census[key]=census.get(key,0.)+float(weight)
        return {'mean_charge':mean,'variance':float(p@(self.charge**2)-mean**2),
                'complement_mean_charge':float(p@(1-self.charge)),
                'joint_charge_boundary_gauss_probabilities':[
                    {'charge':q,'boundary_record':b,'gauss':g,'probability':prob}
                    for (q,b,g),prob in sorted(census.items())],
                'sharp_charge_spectrum':sorted(set(map(int,self.charge))),
                'port_vertices':sorted(v for v in self.vertices if v<self.joint.construction.ports),
                'transit_seams':sorted(v-self.joint.construction.ports for v in self.vertices if v>=self.joint.construction.ports)}

    def continuity_audit(self,rho):
        rho=self.joint._density(rho);j=self.joint
        h=j.hamiltonian;q=np.diag(self.charge);b=np.diag(self.boundary)
        charge_rate=float((1j*np.trace(rho@(h@q-q@h))).real)
        boundary_rate=float((1j*np.trace(rho@(h@b-b@h))).real)
        _,ledger=j.commit(rho)
        max_error=0;rejected=0
        for outcome in ledger:
            source=outcome['source'];target=outcome['target']
            delta_q=int(self.charge[target]-self.charge[source])
            delta_b=int(self.boundary[target]-self.boundary[source])
            e=outcome['seam'];current=outcome['applied_current']
            expected_delta_b=current*(self.boundary_coefficients[2*e]-self.boundary_coefficients[2*e+1])
            max_error=max(max_error,abs(delta_q+delta_b),abs(delta_b-expected_delta_b))
            if outcome['status']=='boundary_rejected':
                rejected+=1
                if current!=0 or delta_q!=0 or delta_b!=0:raise ValueError('rejection changed regional charge')
        return {'basis_gauss_identity_error':int(np.max(np.abs(self.charge+self.boundary-self.gauss))),
                'coherent_charge_rate':charge_rate,'coherent_boundary_record_rate':boundary_rate,
                'coherent_continuity_error':abs(charge_rate+boundary_rate),
                'event_continuity_max_error':int(max_error),'complete_event_outcomes':len(ledger),
                'rejected_outcomes_retained':rejected}


def run():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    variants=[]
    for kind in ('free','first_order','massive'):
        j=JointRecords(c) if kind=='free' else InteractingRecords(c,energy_kind=kind,energy_coupling=.29,energy_mass=.23)
        # Assign occupied seam0 wholly to the region along with port0.
        regional=RegionalCharge(j,ports=(0,),transit_seams=(0,))
        rho,_=j.step(j.initial_density(),.37,.29)
        measured,branches=regional.instrument(rho)
        refined,fine=regional.instrument(rho,True)
        variants.append({'variant':kind,'readout':regional.readout(rho),
                         'outcomes':[{'sharp_charge':r['sharp_charge'],'probability':r['probability']} for r in branches],
                         'refined_outcomes':[{'sharp_charge':r['sharp_charge'],'transit':r['transit'],'probability':r['probability']} for r in fine],
                         'trace_after_charge_measurement':float(np.trace(measured).real),
                         'trace_after_refinement':float(np.trace(refined).real),
                         'continuity':regional.continuity_audit(rho)})
    j=JointRecords(c);region=RegionalCharge(j,ports=(0,),transit_seams=(0,))
    inside=0;outside=next(i for i,s in enumerate(j.states) if s[0]==1)
    psi=np.zeros(len(j.states),complex);psi[inside]=psi[outside]=1/np.sqrt(2)
    rho=np.outer(psi,psi.conj());measured,outcomes=region.instrument(rho)
    fragmentation={'state_basis_indices':[inside,outside],'readout':region.readout(rho),
                   'outcomes':[{'sharp_charge':r['sharp_charge'],'probability':r['probability']} for r in outcomes],
                   'initial_purity':float(np.trace(rho@rho).real),
                   'post_charge_purity':float(np.trace(measured@measured).real),
                   'total_one_load_probability':float(np.trace(rho).real),
                   'half_eigenvalue_present':bool(np.any(region.charge==.5)),
                   'interpretation':'mean one-half has variance one-quarter; sharp outcomes are zero or one. One-load sector does not classify physical species or multi-load fragmentation.'}
    return {'schema':'oph.regional_charge.v1',
            'scope':'Finite operational one-load charge instrument and continuity; no physical row1 pass.',
            'supplied':{'graph':'oriented triangle','cutoff':1,'propagation_mass':.23,'time':.37,
                        'commit_probability':.29,'interaction_coupling':.29,'record_energy_mass':.23,
                        'region_ports':[0],'region_transit_seams':[0],
                        'charge_unit':'one supplied load unit, physical electric-charge mapping absent'},
            'variants':variants,'fragmentation_control':fragmentation,
            'checks':{'complete_complementary_outcomes':all(abs(sum(o['probability'] for o in v['outcomes'])-1)<1e-12 for v in variants),
                      'complete_transit_refinement':all(abs(sum(o['probability'] for o in v['refined_outcomes'])-1)<1e-12 for v in variants),
                      'regional_continuity':all(v['continuity']['basis_gauss_identity_error']==0 and v['continuity']['event_continuity_max_error']==0 and v['continuity']['coherent_continuity_error']<1e-12 for v in variants),
                      'sharp_half_absent':not fragmentation['half_eigenvalue_present']},
            'mutations':{'expectation_is_not_sharp_half':abs(fragmentation['readout']['mean_charge']-.5)<1e-12 and abs(fragmentation['readout']['variance']-.25)<1e-12,
                         'charge_monitor_disturbs_cross_region_coherence':abs(fragmentation['initial_purity']-fragmentation['post_charge_purity'])>.4},
            'physical_species_mapping':None,'physical_charge_unit_mapping':None,
            'full_criterion_met':False,'official_cells_awarded':0}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=run();payload=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(payload)
    print(payload)
