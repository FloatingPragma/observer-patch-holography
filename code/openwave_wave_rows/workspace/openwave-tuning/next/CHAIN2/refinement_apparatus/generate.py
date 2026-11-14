"""Emit a finite native macro tape; source primitive expansion is literal."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import apparatus as b

def build(size=2,steps=2,setting='plus',precision=64):
 q,d=b.member(size);source=b.a.reversible_program(size);clock=source['layout']['memory_bits']
 if type(steps)is not int or steps<0:raise ValueError('nonnegative steps required')
 reset=next(b.native_program(size,steps,setting,precision=precision))
 word={'schema':'oph.chain2.refinement_apparatus.native_macro_tape.v1','model_id':b.MODEL,
  'size':size,'steps':steps,'h_exact':str(F(1,size)),'tau0_exact':'1','setting':setting,'clock_quadrature':'X',
  'reset':reset,'source_layout':source['layout'],'clock_wire':clock,
  'preparation_pulses':b.emit(b.selector_template(size,setting)+b.preparation_template(size)+b.H(clock),precision),
  'repeated_source':{'repetitions':steps,'reflection_pulses':b.emit(b.a.reflection_template(),precision),
    'reversible_instructions':source['gates'],
    'native_primitive_templates':{k:b.emit(w,precision)for k,w in b.a.primitive_templates().items()},
    'clock_tick':b.a.clock_tick_pulse(size,F(1,size),precision,clock)},
  'clock_analysis':b.emit(b.H(clock),precision),'terminal_read_bits':[*range(d),clock],
  'terminal_native_pointer_template':b.emit(b.detector_template(),precision),
  'terminal_native_record_operation':{'seam':[0,1],'law':'configured C or E','codes':[0,1,2,3],
    'coarse_bit':'int(2*current+2*coin-1<0)','retained_fields':['coin','label','current','all output configurations']},
  'execution_rule':'Execute chronological lists. For a reversible instruction (kind,targets), substitute targets into every local primitive pulse; preserve all scalar phases and first-export aborts. Repeated_source executes reflection, reversible instructions, then clock_tick. Each terminal read uses its separate ready auxiliary pair from reset.',
  'failure_rule':'At the first M6 complement export or exhausted native reset, retain the failure prefix, full poststate and reference and halt. Native pointer labels are retained. No success normalization.',
  'budget':b.development_budget(size,steps,precision)}
 return word

def main():
 p=argparse.ArgumentParser();p.add_argument('--size',type=int,default=2);p.add_argument('--steps',type=int,default=2);p.add_argument('--output',type=Path,default=Path(__file__).with_name('PROGRAM.json'));args=p.parse_args()
 word=build(args.size,args.steps);args.output.write_text(json.dumps(word,separators=(',',':'))+'\n')
 print(json.dumps({'program':str(args.output),'sha256':hashlib.sha256(args.output.read_bytes()).hexdigest(),'native_source_gates':len(word['repeated_source']['reversible_instructions'])}))
if __name__=='__main__':main()
