"""Complete native spin-half neutral transfer on the side-four source."""
from collections import defaultdict
from functools import lru_cache
import sympy as s
from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration, KineticTerm, kinetic_terms, yukawa_words, cells, cell_index
from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as native
from model.defect_mass_join.source import car


def clean(v):return {k:s.simplify(z)for k,z in v.items()if s.simplify(z)!=0}
def dot(a,b):return s.simplify(sum(s.conjugate(z)*b.get(k,0)for k,z in a.items()))
def local_action(poly,state):
 out=defaultdict(lambda:s.Integer(0))
 for (mask,h),z in state.items():
  for key,w in car.column(poly,mask,h).items():out[key]+=z*w
 return clean(out)

@lru_cache(maxsize=1)
def spinors():
 masks=[sum(1<<i for i in q)for q in [(10,11,24),(9,10,26),(9,11,25)]]
 up={(m,0):v/s.sqrt(6)for m,v in zip(masks,(2,-1,1))}
 lower=car.add(car.spin_current(0),car.scale(car.spin_current(1),-s.I))
 down=local_action(lower,up);states=[up,down]
 for state in states:
  assert dot(state,state)==1
  for a in range(12):assert not local_action(car.gauge_current(a),state)
  sq=defaultdict(lambda:s.Integer(0))
  for a in range(3):
   for k,v in local_action(car.spin_current(a),local_action(car.spin_current(a),state)).items():sq[k]+=v
  assert clean(sq)==clean({k:s.Rational(3,4)*v for k,v in state.items()})
 assert dot(up,down)==0
 return states


def two_site(state,site):return {(mask<<(30*site),(0,),(0,0)):value for(mask,_),value in state.items()}
def apply_hop(state,term):
 out=defaultdict(lambda:s.Integer(0))
 for key,value in state.items():
  for target,z in native.hopping_column(key,(term,)).items():out[target]+=value*z
 return clean(out)

@lru_cache(maxsize=1)
def edge_moments():
 matrices=[];supports=[]
 for axis in range(3):
  term=KineticTerm(1,0,axis,1,(0,),s.Integer(1));M=s.zeros(2)
  for col,state in enumerate(spinors()):
   phi=two_site(state,0);powers=[phi]
   for k in range(3):powers.append(apply_hop(powers[-1],term))
   for row,target in enumerate(spinors()):
    receiver=two_site(target,1)
    assert dot(receiver,powers[0])==dot(receiver,powers[1])==dot(receiver,powers[2])==0
    M[row,col]=dot(receiver,powers[3])
   supports.append([len(p)for p in powers])
  matrices.append(M)
 sigma=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
 assert all(M==-6*sig for M,sig in zip(matrices,sigma))
 return matrices,supports


def certificate():
 config=Configuration(side=4);terms=kinetic_terms(config)
 x=0;y=cell_index(4,(0,0,1));selected=[t for t in terms if{t.left_cell,t.right_cell}=={x,y}]
 assert len(selected)==2
 amplitudes=[s.conjugate(t.coefficient)if t.left_cell==x else t.coefficient for t in selected]
 assert amplitudes==[-s.I*s.sqrt(2)/8]*2
 matrices,supports=edge_moments()
 coefficient=s.simplify(-6*sum(z**3 for z in amplitudes))
 assert coefficient==-3*s.I*s.sqrt(2)/64
 # A scalar/CAR matrix unit and compressed group coordinate are contractions.
 onsite=sum(abs(word[-1])for word in yukawa_words())
 B=s.simplify(64*(onsite+1)+192+60*sum(abs(t.coefficient)for t in terms))
 upper=s.Integer(s.ceiling(B));time=s.Rational(1,10**20)
 lead=3*s.sqrt(2)*time**3/(64*6);remainder=upper**4*time**4/24
 assert s.simplify(lead-remainder).is_positive
 return {'schema':'oph.side_four_native_neutral_transfer.v1','source':native.descriptor(config),
  'physical_preparation':'normalized local e^c LL gauge singlet with spin1/2, all other matter/scalars empty and links blank; idle opposite-L odd reference retained',
  'local_spin_up_mask_amplitudes':{str(mask):str(v)for(mask,_),v in spinors()[0].items()},
  'exact_local_gauge_commutators_per_spin':12,'local_total_spin_Casimir':'3/4','electric_charge':'0','raw_lepton_number':'1',
  'all_axes_complete_link_third_moments':[[[str(z)for z in row]for row in M.tolist()]for M in matrices],
  'complete_two_site_power_supports':supports,
  'source_paths':[{'left':t.left_cell,'right':t.right_cell,'path':list(t.path),'coefficient':str(t.coefficient)}for t in selected],
  'full_source_H0_H1_H2_receiver_amplitudes':['0','0','0'],'full_source_H3_receiver_amplitude':str(coefficient),
  'source_operator_norm_upper':str(upper),'time':str(time),'Taylor_degree':3,
  'leading_receiver_amplitude_modulus':str(lead),'full_source_state_remainder_upper':str(remainder),
  'original_receiver_probability_interval':[str(s.expand((lead-remainder)**2)),str(s.expand((lead+remainder)**2))],
  'original_receiver_probability_decimal_interval':[float((lead-remainder)**2),float((lead+remainder)**2)],
  'complete_records':'128 mutually orthogonal site/spin packet projectors plus their full-carrier complement; source/reference preparation deterministic by declared parity-even unitary',
  'all_192_PW119_links_and_1920_native_CAR_modes_retained':True,
  'boundary_control':'stripping both seam signs makes the two path amplitudes opposite and cancels this third moment',
  'boundary_control_third_moment':str(s.simplify(-6*((s.I*s.sqrt(2)/8)**3+(-s.I*s.sqrt(2)/8)**3))),
  'same_side_two_low_energy_band_or_ground_transferred':False,'relativistic_mass_or_observed_neutrino_identified':False,
  'verdict':'PASS_NATIVE_SIDE_FOUR_NEUTRAL_TRANSFER'}
