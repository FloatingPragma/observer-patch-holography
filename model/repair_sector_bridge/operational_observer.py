"""Actual bounded native source and two addressed-record observers."""
from dataclasses import dataclass,replace
from pathlib import Path
import hashlib,json
from .record_observer import Observer,tables,encode,memory_bound

def transition_table():
    from model.observer.native import native
    states,seams,_=tables();index={s:i for i,s in enumerate(states)};rows=[]
    for i,state in enumerate(states):
        for address,seam in enumerate(seams):
            for coin in (0,1):
                out,(_,current)=native(state,seam,coin)
                rows.append({'state':i,'branch':2*address+coin,'next':index[out],'letter':encode(address,coin,current)})
    return rows

@dataclass(frozen=True)
class Event:
    sequence:int
    letter:int

@dataclass(frozen=True)
class Source:
    initial_index:int
    state_index:int
    capacity:int
    attempts:int=0
    halted:bool=False
    @classmethod
    def prepare(cls,index,capacity):
        Observer.prepare(index,capacity)
        return cls(index,index,capacity,halted=capacity==0)
    def attempt(self,address,coin):
        if self.halted:return self,None
        if type(address) is not int or not 0<=address<30 or type(coin) is not int or coin not in (0,1):raise ValueError('valid seam/coin required')
        from model.observer.native import native
        states,seams,_=tables();out,(_,j)=native(states[self.state_index],seams[address],coin);n=self.attempts+1
        return replace(self,state_index=states.index(out),attempts=n,halted=n==self.capacity),Event(self.attempts,encode(address,coin,j))

@dataclass(frozen=True)
class Pair:
    source:Source
    left:Observer
    right:Observer
    status:str='active'
    @classmethod
    def prepare(cls,index,capacity):
        return cls(Source.prepare(index,capacity),Observer.prepare(index,capacity),Observer.prepare(index,capacity),status='halted' if capacity==0 else 'active')
    @property
    def overlap(self):
        return tuple(Event(i,letter) for i,letter in enumerate(self.left.tape))
    def step(self,address,coin):
        if self.status!='active':return self
        if self.left.initial_index!=self.right.initial_index or self.left.initial_index!=self.source.initial_index or self.left.state_index!=self.right.state_index or self.left.capacity!=self.right.capacity or self.left.capacity!=self.source.capacity or self.left.tape!=self.right.tape or len(self.left.tape)!=self.source.attempts or self.left.status!='active' or self.right.status!='active':
            return replace(self,source=replace(self.source,halted=True),status='observer_mismatch')
        source,event=self.source.attempt(address,coin)
        if event is None:return replace(self,source=source,status='halted')
        left=self.left.consume(event.sequence,event.letter);right=self.right.consume(event.sequence,event.letter)
        status='halted' if source.halted else 'active'
        if left.status not in ('active','full') or right.status not in ('active','full'):source=replace(source,halted=True);status='record_error'
        return replace(self,source=source,left=left,right=right,status=status)
    def readback(self):
        return {'left':list(self.left.tape),'right':list(self.right.tape),'shared':[(e.sequence,e.letter) for e in self.overlap]}

def resources(capacity=39):
    per=memory_bound(capacity)
    return {'per_observer_registers':per,'two_observer_retained_bits_upper':2*per['retained_register_bits_upper'],'source_state_index_bits':7,'source_preparation_label_bits':7,'pair_coordinator_status_bits':2,'source_working_binary_load_bits':12,'source_counter_bits':capacity.bit_length(),'source_halt_bit':1,'shared_incoming_packet_bits':7+capacity.bit_length(),'excluded_resources':['immutable program/geometry','known initial-state preparation apparatus','seam selector and coin generation','native attempt clock physical calibration','optional external archive beyond capacity','hardware implementation/error correction'],'shared_overlap_storage':'derived from retained left tape and implicit sequence indices; no third persistent tape','readback_output_workspace_bits_upper':capacity*(7+capacity.bit_length()),'Python_object_overhead_included':False,'external_archive_used_by_execution':False,'source_halts_at_capacity':True}

def validate_resource_contract(contract,capacity=39):
    expected=resources(capacity)
    for name in ('per_observer_registers','two_observer_retained_bits_upper','source_state_index_bits','source_preparation_label_bits','pair_coordinator_status_bits','source_working_binary_load_bits','source_counter_bits','source_halt_bit','shared_incoming_packet_bits','readback_output_workspace_bits_upper','excluded_resources'):
        if contract.get(name)!=expected[name]:raise ValueError('missing or altered resource obligation: '+name)
    if contract.get('external_archive_used_by_execution') is not False or contract.get('source_halts_at_capacity') is not True:raise ValueError('unadmitted continuing archive/source resource')
    return True

def report():
    rows=transition_table();noops=0
    for row in rows:
        p=Pair.prepare(row['state'],1).step(row['branch']//2,row['branch']%2)
        if not(p.left.state_index==p.right.state_index==p.source.state_index==row['next'] and p.left.tape==p.right.tape==(row['letter'],) and p.step(0,0)==p):raise ValueError('actual continuation mismatch')
        noops+=int(row['next']==row['state'])
    p=Pair.prepare(0,2).step(0,0);a=p.step(0,0);b=p.step(1,0)
    if a.overlap[:1]!=b.overlap[:1] or a.overlap[1].letter==b.overlap[1].letter:raise ValueError('scheduler witness failed')
    full=Pair.prepare(0,39)
    for t in range(39):full=full.step(t%30,t%2)
    if full.source.attempts!=39 or full.status!='halted' or full.step(0,0)!=full or full.left.state_index!=full.source.state_index or full.right.state_index!=full.source.state_index:raise ValueError('39-attempt source/readback mismatch')
    root=Path(__file__).resolve().parents[2];paths=[Path(__file__),Path(__file__).with_name('record_observer.py'),root/'model/observer/native.py',root/'model/observer/quantum_mixing.py',root/'model/repair_sector_bridge/bridge.py',root/'Lean/Tower/OperationalObserver.lean',root/'Lean/QFT/OperationalOverlapEvidence.lean',root/'Lean/Tower/RepairRecordObserver.lean',Path(__file__).with_name('test_operational_observer.py'),Path(__file__).with_name('OPERATIONAL_OBSERVER.md')]
    return {'schema':'oph.addressed.operational_observer.v1','bounded_progress_status':'partial executable binding; full OperationalObserver tower attachment open','canonical_program_update_boundary':'Both queue/current-renderer implementations are frozen by model/completion/component_interfaces.json source bindings; no campaign rebind or full-completion claim is made','scheduler_contract':'Declared uniform choice of thirty seams and fair coin gives the 1/60 Markov generator; executable step accepts each selector output as explicit input, not an autonomous random device.','actual_transition_table_sha256':hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest(),'all_transition_rows':rows,'checks':{'all3960_actual_two_observer_continuations':len(rows)==3960,'noops_retained':noops>0,'source_halt_capacity':True,'own_record_readback':p.readback()['left']==p.readback()['right']},'bounded39_attempt_execution':full.readback(),'no_op_source_branches':noops,'resource_contract':resources(),'self_description':{'domain':'Own bounded tape readback, not immutable program reconstruction','same_prefix':p.readback(),'indistinguishable_prefix_scheduler_variants':'Both schedules emit first branch0; one continues branch0, the other branch2. Same native law and observed prefix, different unobserved continuation.','next_letters':[a.overlap[1].letter,b.overlap[1].letter]},'Lean_status':'Finite table binding separate; no full OperationalObserver tower attachment claimed','source_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'physical_or_unrestricted_internal_closure':False}

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    result=report();path=Path(__file__).with_name('operational_observer_results.json')
    if args.check:
        if json.loads(path.read_text())!=json.loads(json.dumps(result)):raise SystemExit('stale operational observer source/result receipt')
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['checks']))
