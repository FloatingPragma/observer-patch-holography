"""Native zero-momentum rest-frequency and shared-clock control."""
from fractions import Fraction as F
from pathlib import Path
import json,math
import numpy as np
import native_spectroscopy as ns
from run import readout,contrast,pins
HERE=Path(__file__).resolve().parent

def main():
 before=pins();size=2;_,d=ns.b.member(size);initial=np.zeros(1<<d,complex);initial[0]=1
 word=ns.preparation_template(size,(0,0,0));prep=ns.b.complete_run(word,d,initial,96);ns.b.require_complete(prep)
 psi=prep['success'][:,0];norm=float(np.vdot(psi,psi).real);failed=[float(np.trace(g).real) for g in prep['failure_grams']]
 rows=[]
 for mass in (F(0),F(1,2),F(1)):
  walk=ns.walk.Walk(size,h=F(1,size),mass=mass);state=psi.reshape(walk.geometry.shape).copy()
  for _ in range(3):state=walk.step(state)
  c=np.vdot(psi,state.reshape(-1));target=np.exp(-1j*3*float(mass)/size)
  assert abs(c-target)<2e-11
  rho=np.array([[norm,c.conjugate()],[c,norm]])/2;cases=[]
  for q in ('X','Y'):
   for law in ('C','E'):
    meter=readout(rho,q,law);clock=readout(np.ones((2,2))/2,q,law,3,size)
    expected=target.real if q=='X' else target.imag
    expected_clock=math.cos(3/(8*size)) if q=='X' else math.sin(3/(8*size))
    assert abs(contrast(meter)-expected)<2e-11 and abs(contrast(clock)-expected_clock)<2e-11
    assert abs(sum(r['probability'] for r in meter)+sum(failed)-1)<2e-11
    cases.append({'quadrature':q,'law':law,'matter':meter,'clock':clock,'preparation_export_failures':failed})
  rows.append({'mass_exact':str(mass),'mode':[0,0,0],'steps':3,'source_time_in_tau0_exact':'3/2','return':[float(c.real),float(c.imag)],'cases':cases})
 assert pins()==before
 (HERE/'REST_CONTROL.json').write_text(json.dumps({'rows':rows,'same_native_preparation_clock_readout':True,'mass_is_input_not_empirical_prediction':True,'provider_pins_unchanged':True},indent=2,sort_keys=True)+'\n')
 print(json.dumps({'rest_cases':len(rows),'complete_fine_pointer_laws':12,'source_time_in_tau0_exact':'3/2'}),flush=True)
if __name__=='__main__':main()
