"""Attach the finite collar bound to authenticated native port windows."""
from collections import deque
from fractions import Fraction as F
from model.native_rooted_refinement import routing
from . import collar
from .certificate import require


def row(level,radius):
 window=routing.routed_window(level,radius)
 names=window['canonical_ports'];where={name:i for i,name in enumerate(names)}
 cells={c for c,_ in window['routing']}
 boundary=[]
 for cell,neighbours in window['routing']:
  for other,port in neighbours:
   if other not in cells:boundary.append(where[(cell,port)])
 require(boundary and len(set(boundary))==len(boundary),'the source boundary ports are not distinct')
 adjacency=[set()for _ in names]
 for a,b in window['internal_seams']:adjacency[a].add(b);adjacency[b].add(a)
 require(max(len(x)+int(i in boundary)for i,x in enumerate(adjacency))<=6,'the source degree bound failed')
 distance={i:0 for i in window['root_ports']};queue=deque(distance)
 while queue:
  i=queue.popleft()
  for j in adjacency[i]:
   if j not in distance:distance[j]=distance[i]+1;queue.append(j)
 require(len(distance)==len(names),'the routed window is disconnected')
 collar_radius=min(distance[i]for i in boundary);require(collar_radius>=radius,'the actual port collar is shorter than its cell collar')
 V=240*4**level;M=630*4**level
 bound=collar.combined_bound(V,M,V//2,len(names),len(window['internal_seams']),len(boundary),12,6,collar_radius,F(1,100),F(1,2))
 return {'source_level':level,'authenticated_cell_radius':radius,'ports':len(names),'internal_seams':len(window['internal_seams']),'crossing_seams':len(boundary),'actual_port_collar_distance':collar_radius,'absolute_exchange_time':'1/100','time_orientation':'modular parameter 2M*t compares to Heisenberg exchange time -t','density':'1/2','global_ports':V,'global_seams':M,'bounds':{k:str(v)for k,v in bound.items()},'global_state_preparation_claimed':False,'global_graph_constructor_executed':False,'all_level_numeric_routing_claimed':False}
