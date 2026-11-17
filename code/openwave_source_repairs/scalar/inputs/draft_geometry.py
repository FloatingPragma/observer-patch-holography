"""Exact source quaternion cover and compact-sector finite rotations."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import importlib.util,sys
ROOT=Path(__file__).resolve().parents[4]
spec=importlib.util.spec_from_file_location('symmetric_source_geometry',ROOT/'openwave/shared/paired_source_moments/run.py');m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
q,add,neg,mul,sub,dot,total=m.q,m.add,m.neg,m.mul,m.sub,m.dot,m.total
ZERO=q();ONE=q(1)
def scale(a,x):return tuple(mul(a,v) for v in x)
def cross(x,y):return (sub(mul(x[1],y[2]),mul(x[2],y[1])),sub(mul(x[2],y[0]),mul(x[0],y[2])),sub(mul(x[0],y[1]),mul(x[1],y[0])))
def qm(a,b):
 av,bv=a[1:],b[1:];c=cross(av,bv)
 return (sub(mul(a[0],b[0]),dot(av,bv)),)+tuple(total((mul(a[0],bv[i]),mul(b[0],av[i]),c[i])) for i in range(3))
def qr(a,x):
 v=a[1:];z=tuple(add(mul(a[0],x[i]),cross(v,x)[i]) for i in range(3));t=cross(v,z)
 return tuple(add(x[i],mul(q(2),t[i])) for i in range(3))
def cover(cell):
 points,anti,half,lean=m.geometry(cell)
 alphabet=[q(),q(1),q(-1),q(F(1,2)),q(F(-1,2)),q(0,F(1,2)),q(0,F(-1,2)),q(F(-1,2),F(1,2)),q(F(1,2),F(-1,2))]
 rows=[]
 for a in product(alphabet,repeat=4):
  if total(mul(x,x) for x in a)!=ONE:continue
  moved=[qr(a,p) for p in points]
  if not all(x in points for x in moved):continue
  rows.append((a,tuple(points.index(x) for x in moved)))
 rows.sort();assert len(rows)==120
 lookup={a:i for i,(a,g) in enumerate(rows)};perms={g for a,g in rows};assert len(perms)==60
 products=[]
 for a,g in rows:
  row=[]
  for b,h in rows:
   j=lookup[qm(a,b)];assert rows[j][1]==tuple(g[h[p]] for p in range(12));row.append(j)
  products.append(row)
 assert all(sum(h==g for b,h in rows)==2 for a,g in rows)
 return points,anti,half,rows,products
