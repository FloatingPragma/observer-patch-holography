import json, numpy as np
EDGES=[(0,1),(0,5),(0,7),(0,10),(0,11),(1,5),(1,7),(1,8),(1,9),(2,3),(2,4),(2,6),(2,10),(2,11),(3,4),(3,6),(3,8),(3,9),(4,5),(4,9),(4,11),(5,9),(5,11),(6,7),(6,8),(6,10),(7,8),(7,10),(8,9),(10,11)]
rows=[]
for c in (1,2,4,8,16,32,64,128,256,512):
 for w,v in ((1,0),(24,1),(32,1)):
  L=np.zeros((14,14))
  for a,b,k in [(a,b,c) for a,b in EDGES]+[(0,12,w),(0,13,w),(12,13,v)]:
   z=np.eye(14)[a]-np.eye(14)[b];L+=k*np.outer(z,z)
  d,U=np.linalg.eigh(L);f=U[:,1]
  rows.append(dict(core=c,leaf=w,bridge=v,eigenvalues=d.tolist(),lowest_leaf_mass=float(sum(f[12:]**2)),compact_lambda=w+2*v,compact_is_lowest=bool(abs(d[1]-w-2*v)<1e-9)))
print(json.dumps({'scope':'Exploratory spectral design pilot; no confirmatory gates or native-call claim','rows':rows},indent=2))
