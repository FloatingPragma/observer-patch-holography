"""Native terminal pointer on retained unfiltered source-word outputs."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,time
import numpy as np
import native_catalog as nc
import local_flow as lf
from run import decode
from finalize import finish
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    start=time.monotonic();reg=json.loads((HERE/'PREREG_OBSERVER.json').read_text())
    for n,h in reg['files'].items():assert sha(HERE/n)==h,n
    flow=HERE/reg['flow_attempt']
    r=json.loads((flow/'RECEIPT.json').read_text());arrays=np.load(flow/'SOURCE_OUTPUTS.npz')
    basis=np.array(json.loads((HERE/'inputs/certificate.json').read_text())['integer_modes'])
    a=np.column_stack((basis[0]/2,(2*basis[1]-basis[0])/math.sqrt(12)))
    b=np.column_stack((basis[2]/2,(2*basis[3]-basis[2])/math.sqrt(12)))
    pair=np.stack([np.outer(a[:,i],b[:,j])for i in range(2)for j in range(2)],axis=-1)
    effects=np.load(HERE/'inputs/native_pointer.npz')['effects'];pulses={};bounds={}
    for kind,angle in [('Y',-math.pi/4),('Y',math.pi/4),('Z',math.pi/4)]:
        matrix=nc.c.native_pulse(kind,angle);wanted=nc.target(kind,angle)
        sq=sum((nc.iv.scalar_error(matrix[i,j].real,wanted[i][j][0])**2+
                nc.iv.scalar_error(matrix[i,j].imag,wanted[i][j][1])**2 for i in range(6)for j in range(6)),F())
        key=kind+'/'+angle.hex();pulses[key]=matrix;bounds[key]=str(nc.sqrt_upper(sq))
    assert 2*max(F(v)for v in bounds.values())+F(1,10**15)<=F(r['native_catalog']['complete_per_pulse_allowance'])
    prep=pulses['Y/'+(-math.pi/4).hex()];plus=prep[:2,0];prepfail=float(np.vdot(prep[2:,0],prep[2:,0]).real)
    def observe(z,norm):
        rho=np.array([[abs(plus[0])**2,plus[0]*plus[1].conjugate()*z.conjugate()],
                      [plus[1]*plus[0].conjugate()*z,abs(plus[1])**2*norm]],complex)
        result={};conditional=[]
        for setting in ('X','Y'):
            state=rho.copy();failure=0.
            order=[('Z',math.pi/4),('Y',math.pi/4)]if setting=='Y'else[('Y',math.pi/4)]
            for kind,angle in order:
                matrix=pulses[kind+'/'+angle.hex()]
                failure+=float(np.trace(matrix[2:,:2]@state@matrix[2:,:2].conj().T).real)
                state=matrix[:2,:2]@state@matrix[:2,:2].conj().T
            law=[float(np.trace(e@state).real)for e in effects]
            residual=float(np.trace(rho).real-failure-sum(law))
            result[setting]={'pointer_output_centres':law,'reference_preparation_failure_centre':prepfail,
                'analysis_processor_failure_centre':failure,'arithmetic_mass_residual':residual}
            conditional.append((law[0]-law[1])/sum(law))
        result['fringe_centre']=conditional
        result['unresolved_upstream_semantics']='Full source/preparation/route/load/processor/scratch/reset failure laws remain in the complete bound; they are not set to zero or absorbed by normalization.'
        return result
    records={};worst=0.;allres=0.
    for t in ('0.1','4'):
        records[t]={}
        for name,initial in [('single',a),('pair',pair)]:
            final=arrays[name+'_'+t.replace('.','p')];z=np.vdot(initial,final)/initial.shape[-1]
            norm=float(np.vdot(final,final).real/initial.shape[-1]);record=observe(z,norm)
            records[t][name]=record
            worst=max(worst,abs(complex(*record['fringe_centre'])-z))
            allres=max(allres,*(abs(record[k]['arithmetic_mass_residual'])for k in ('X','Y')))
    changed=r['changed_seam']['record'];records['changed_seam']=observe(complex(*changed['fringe_centre']),changed['final_norm_centre'])
    records['full_HOLD']=observe(1+0j,1.);records['source_vacuum']=observe(1+0j,1.)
    graph=json.loads((HERE/'inputs/graph.json').read_text());edges=[tuple(x)for x in graph['seams']]
    steps=r['words']['single']['source_local_edge_words']//(10*len(edges))
    program,_,program_hash=lf.description(edges,4,steps)
    assert program_hash==r['source_program_hashes']['4']
    coordinate=np.zeros((240,1),complex);coordinate[0]=1
    coordinate_out=lf.propagate(coordinate,program,steps,False)
    zcoordinate=np.vdot(coordinate,coordinate_out)
    records['unprepared_coordinate']=observe(zcoordinate,float(np.vdot(coordinate_out,coordinate_out).real))
    coordinate_difference=abs(complex(*records['unprepared_coordinate']['fringe_centre'])-
                              complex(*records['4']['single']['fringe_centre']))
    decoded={}
    for name in ('single','pair'):
        decoded[name]=decode(records['0.1'][name],records['4'][name],float(F(r['complete_budgets'][name]['complete'])))
    one=decoded['single']['energy_interval'];two=decoded['pair']['energy_interval']
    ratio=[two[0]/one[1],two[1]/one[0]];kev=[511*x for x in ratio]
    gates={'native_pointer_centres_match_source_interference':worst<5e-12,
           'complete_reference_mass_residual':allres<1e-12,
           'pair_interval_contains2_and_excludes_old':ratio[0]<2<ratio[1]and kev[1]<1024.45085,
           'same_source_word_resolves_preparation':coordinate_difference>.05,
           'HOLD_and_vacuum_reference_zero':all(abs(complex(*records[k]['fringe_centre'])-1)<5e-12 for k in ('full_HOLD','source_vacuum'))}
    gates={k:bool(v)for k,v in gates.items()}
    result={'prereg_sha256':sha(HERE/'PREREG_OBSERVER.json'),'flow_receipt_sha256':sha(flow/'RECEIPT.json'),
      'records':records,'decoded':decoded,'pair_to_single_interval':ratio,'pair_keV_at511_anchor':kev,
      'native_analysis_pulse_bounds':bounds,'fresh_native_M6_matrices':3,
      'maximum_native_pointer_fringe_difference':worst,'maximum_mass_residual':allres,
      'unprepared_coordinate_fringe_separation':float(coordinate_difference),
      'unprepared_coordinate_source_program_hash':program_hash,
      'gates':gates,'all_gates_pass':all(gates.values()),'elapsed_seconds':time.monotonic()-start}
    out=HERE/'observer';out.mkdir(exist_ok=False)
    np.savez_compressed(out/'NATIVE_TERMINAL_MATRICES.npz',**{str(i):m for i,m in enumerate(pulses.values())})
    result['matrix_keys']=list(pulses)
    exit_code=finish(out/'RECEIPT.json',result)
    print(json.dumps({'gates':gates,'pair_keV':kev,'seconds':result['elapsed_seconds']}))
    raise SystemExit(exit_code)
if __name__=='__main__':main()
