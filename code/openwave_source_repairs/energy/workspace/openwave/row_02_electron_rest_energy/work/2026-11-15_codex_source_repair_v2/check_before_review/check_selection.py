#!/usr/bin/env python3
"""Independent channel reconstruction and all-input convergence certificate.

No producer module is imported. The raw word columns are retained evidence;
this checker separately reconstructs the source law, compact subspace and
complete pump composition. It is not a fresh raw native-word execution.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent;OW=HERE.parents[2]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(record=False):
    start=time.monotonic()
    freeze=json.loads((HERE/'PREREG_CHECK.json').read_text())
    for name,digest in freeze['files'].items():assert sha(HERE/name)==digest
    receipt=json.loads((HERE/'execution/RECEIPT.json').read_text());science=receipt['scientific']
    for name,digest in receipt['artifacts'].items():assert sha(HERE/'execution'/name)==digest
    app=json.loads((HERE/'apparatus/RECEIPT.json').read_text())
    for name,digest in app['artifacts'].items():assert sha(HERE/'apparatus'/name)==digest
    source=np.load(OW/'shared/compact_logarithmic_flow/source_columns.npz')
    counts=np.zeros((240,240),dtype=np.int64)
    # Two native unit branches per edge, independently using floor/ceiling
    # equalization. The 302400 represented columns are not literal calls.
    labels=np.arange(240)
    for a,b in source['seams']:
        for coin in (0,1):
            left=(labels==a).astype(int);right=(labels==b).astype(int)
            total=left+right;leftnew=(total+(2*coin-1)*(total%2))//2
            dest=labels.copy();touch=total==1;dest[touch]=np.where(leftnew[touch],a,b)
            np.add.at(counts,(dest,labels),1)
    assert np.array_equal(counts,source['counts'])
    L=1260*np.eye(240,dtype=np.int64)-counts
    D=np.zeros((240,2),dtype=np.int64);D[[4,5,6,7],0]=[1,-1,-1,1];D[[5,6,8,11],1]=[-1,-1,1,1]
    assert np.array_equal(L@D,6*D)
    P=D[:12]@np.array([[2,-1],[-1,2]],dtype=float)@D[:12].T/6
    assert np.max(abs(P@P-P))<1e-15 and abs(np.trace(P)-2)<1e-15
    data=np.load(HERE/'execution/maps.npz');flow=data['raw_flow_columns'];Q=flow[:12]
    R=flow[12:240];escape_effect=R.conj().T@R;sigma=np.eye(12)/12
    final=[];outs=[]
    for label in range(12):
        rho=np.zeros((12,12),complex);rho[label,label]=1
        escapes=0.
        for n in range(512):
            escaped=float(np.trace(escape_effect@rho).real);escapes+=escaped
            rho=Q@rho@Q.conj().T+escaped*sigma
        final.append(float(np.trace(P@rho).real));outs.append(escapes)
    assert min(final)>1-1e-7
    rho=sigma.copy()
    for n in range(512):
        escaped=float(np.trace(escape_effect@rho).real)
        rho=Q@rho@Q.conj().T+escaped*sigma
    assert np.max(abs(rho-data['terminal_density']))<1e-12
    # Same controller with a complete coordinate dephasing after each epoch.
    deph=sigma.copy()
    for n in range(512):
        escaped=float(np.trace(escape_effect@deph).real)
        deph=Q@deph@Q.conj().T+escaped*sigma;deph=np.diag(np.diag(deph))
    dephweight=float(np.trace(P@deph).real);assert dephweight<.5
    complement=np.eye(12)-P
    # All-input tail operator, without assuming a chosen initial state.
    tail=complement.copy().astype(complex)
    exact_identity_deficit=np.eye(12)-Q.conj().T@Q
    for n in range(64):
        tail=Q.conj().T@tail@Q+np.trace(tail)/12*exact_identity_deficit
        assert np.max(abs(tail))<2
    tail_frobenius=float(np.linalg.norm(tail,'fro'));assert tail_frobenius<.22
    # Standard binary64 matrix-product model: each d<=4 complex dot product
    # has local operator error at most4 gamma32(1+b). No underflow/overflow
    # is used; all retained success amplitudes have magnitude below2.
    assert np.max(abs(flow))<2
    b=F(science['full_primitive_error_upper']);N=int(science['native_pulse_applications'])
    u=F(1,2**53);gamma32=32*u/(1-32*u);local=4*gamma32*(1+b)
    propagation=4*N*local/(1-N*(b+local))**2
    compiled=F(science['flow_complete_halfdiamond_upper'])
    delta=propagation+2*compiled
    channel_error=4*delta+2*delta*delta
    # An entrywise absolute guard for each12x12 complex triple product:
    # every entry has at most4096 unit-scale real operations under the
    # checked magnitude2 bound. Convert entry error to operator error by12,
    # then retain64 complete iterations and their norm gains. The final
    # Frobenius evaluation receives an additional1e-12 outward allowance.
    gamma128=128*u/(1-128*u)
    round_guard=64*12*4096*gamma128/(1-64*channel_error)+F(1,10**12)
    tail_upper=F(math.ceil(tail_frobenius*10**12),10**12)+round_guard+64*channel_error/(1-64*channel_error)
    assert tail_upper<F(1,4)
    apparatus_error=F(app['scientific']['512_epoch_apparatus_halfdiamond_upper'])
    total=512*compiled+apparatus_error
    assert total<F(1,1000)
    # The exact source H spectrum is in[0,-1260log(104/105)] subset[0,13].
    # A bad-sector probability eps bounds E[(H-E6)^2] by169 eps. This is a
    # source-energy window certificate; it is not a rest-mass dictionary.
    eps=F(1,4**8)+total
    result={'all12_initial_basis_compact_weights':final,'all12_expected_escape_counts':outs,
        'minimum_initial_basis_weight':min(final),'dephasing_control_compact_weight':dephweight,
        'all_input64epoch_tail_operator_frobenius':tail_frobenius,
        'all_input64epoch_tail_norm_upper':str(tail_upper),
        'all_input_ideal_convergence':'bad weight after64k epochs <= (1/4)^k for every initial cell density',
        'raw_state_binary64_operator_error_upper':str(propagation),
        'complete512epoch_implementation_error_upper':str(total),
        'all_input_selected_probability_lower':str(1-eps),
        'unresolved_processor_or_apparatus_probability_upper':str(total),
        'unconditional_valid_source_energy_second_moment_about_E6_upper':str(169*eps),
        'conditional_valid_source_energy_second_moment_about_E6_upper':str(169*eps/(1-total)),
        'row2_physical_energy_identified':False,'other_group_G6':False,
        'producer_modules_imported':[],'execution_receipt_sha256':sha(HERE/'execution/RECEIPT.json'),
        'apparatus_receipt_sha256':sha(HERE/'apparatus/RECEIPT.json')}
    target=HERE/'INDEPENDENT_CHECK.json'
    if record:target.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    else:assert json.loads(target.read_text())==result
    print(json.dumps({'all_basis_pass':True,'minimum_weight':min(final),'dephasing_weight':dephweight,
                      'complete_error':float(total),'seconds':time.monotonic()-start}))

if __name__=='__main__':main('--record' in sys.argv)
