"""Independent two-by-two modular algebra/plot audit; no producer edits."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def require(ok,message):
    if not ok:raise ValueError(message)


def array(data):return np.array(data['real'])+1j*np.array(data['imag'])


def run():
    source=HERE/'MODULAR_TIME.json';receipt=json.loads(source.read_text())
    rho=array(receipt['observer']['conditional_density'])
    a=array(receipt['modular_flow']['observable'])
    # Closed 2x2 projectors, independent of producer's eigh/logm algorithms.
    mid=np.trace(rho).real/2
    radius=np.sqrt(((rho[0,0].real-rho[1,1].real)/2)**2+abs(rho[0,1])**2)
    low,high=mid-radius,mid+radius
    require(low>1e-12,'conditional state fails numerical faithfulness')
    p_low=(high*np.eye(2)-rho)/(high-low);p_high=(rho-low*np.eye(2))/(high-low)
    def power(z):return low**z*p_low+high**z*p_high
    def flow(b,z):return power(1j*z)@b@power(-1j*z)
    inv=power(-1);sqrt=power(.5)
    gns=np.kron(inv.T,rho)
    omega=sqrt.reshape(-1,order='F')
    sample=np.array([[1+.2j,2-1j],[.4+.7j,-.6j]])
    sample_b=np.array([[.2,1+.8j],[-.9+.1j,1.3]])
    vec=sample.reshape(-1,order='F')
    gns_vector_error=float(np.max(abs(gns@vec-(rho@sample@inv).reshape(-1,order='F'))))
    gns_omega_error=float(np.max(abs(gns@omega-omega)))
    kms=[];wrong=[]
    for s in [-.4,0.,.7]:
        rhs=np.trace(rho@flow(sample_b,s)@sample)
        kms.append(abs(np.trace(rho@sample@flow(sample_b,s-1j))-rhs))
        wrong.append(abs(np.trace(rho@sample@flow(sample_b,s+1j))-rhs))
    require(max(kms)<1e-11 and min(wrong)>.1,'general complex-observable KMS sign audit failed')
    rows=list(csv.DictReader((HERE/'correlations.csv').open()))
    recorded=receipt['plot_data'];error=0.
    require(len(rows)==len(recorded)==161,'plot data length mismatch')
    for csv_row,json_row in zip(rows,recorded):
        for key in json_row:error=max(error,abs(float(csv_row[key])-json_row[key]))
        s=float(csv_row['s']);evolved=flow(a,s)
        correlation=np.trace(rho@a@evolved)
        error=max(error,abs(correlation-complex(float(csv_row['C_real']),float(csv_row['C_imag']))),
                  abs(evolved[0,1]-complex(float(csv_row['observable_01_real']),float(csv_row['observable_01_imag']))))
    require(error<1e-11,'plot/source correlation mismatch')
    pins=json.loads((HERE/'PLOT_SOURCES.json').read_text())
    require(all(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest for name,digest in pins.items()),'plot source pins mismatch')
    for path,digest in receipt['source_sha256'].items():
        require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'native source pin mismatch: '+path)
    archive=receipt['native_state_archive'];archive_path=HERE/archive['filename']
    require(hashlib.sha256(archive_path.read_bytes()).hexdigest()==archive['sha256'],'native archive hash mismatch')
    with np.load(archive_path,allow_pickle=False) as z:
        full=z['native_density'];P=z['corner_projector'];Q=z['complement_projector']
        inner=z['unnormalized_inside'];outer=z['unnormalized_complement'];terminal=z['nonselective_terminal_density']
        indices=z['position_indices'];corner=z['conditional_corner_density']
        custody_errors=[np.max(abs(P+Q-np.eye(len(full)))),np.max(abs(P@P-P)),np.max(abs(Q@Q-Q)),
                        np.max(abs(P@Q)),np.max(abs(inner-P@full@P)),np.max(abs(outer-Q@full@Q)),
                        np.max(abs(terminal-inner-outer)),np.max(abs(full-full.conj().T)),
                        np.max(abs(corner-rho)),np.max(abs(inner[np.ix_(indices,indices)]/np.trace(inner)-rho))]
        require(max(custody_errors)<1e-11,'full terminal instrument custody failed')
        require(abs(np.trace(full)-1)<1e-11 and np.linalg.eigvalsh(full).min()>-1e-11,'native density positivity/trace failed')
        require(abs(np.trace(inner)-receipt['observer']['corner_probability'])<1e-11 and abs(np.trace(outer)-receipt['observer']['complement_probability'])<1e-11,'terminal branch probability mismatch')
        terminal_disturbance=float(np.linalg.norm(terminal-full))
    complement=receipt['observer']['complement_probability']
    inside=receipt['observer']['corner_probability']
    trace=receipt['native_source']['full_native_trace']
    require(abs(inside+complement-trace)<1e-12,'conditioning/complement probability mismatch')
    require(np.max(abs(gns-array(receipt['modular_flow']['Delta'])))<1e-10,'recorded GNS mismatch')
    require(gns_vector_error<1e-11 and gns_omega_error<1e-11,'GNS conventions failed')
    require(receipt['observer']['depolarizing_floor_added'] is False and receipt['observer']['thermal_or_Gibbs_state_supplied'] is False,'input state ancestry changed')
    return {'schema':'oph.observer.modular_time.independent_review.v1',
            'outcome':'PASS_WITH_SCOPE_NOTES','terminal_archive_max_identity_error':float(max(custody_errors)),'terminal_nonselective_disturbance_frobenius':terminal_disturbance,'conditional_eigenvalues_analytic':[float(low),float(high)],
            'general_complex_KMS_minus_i_max_error':float(max(kms)),
            'wrong_plus_i_min_error':float(min(wrong)),
            'GNS_direct_vectorization_error':gns_vector_error,'GNS_cyclic_vector_error':gns_omega_error,
            'CSV_JSON_analytic_correlation_max_error':float(error),
            'conditioned_plus_complement_minus_native_trace':float(inside+complement-trace),
            'scope_notes':['Input native propagation is numerically generated, without a certified enclosure; this tiny audit checks the supplied conditional state, not an independent full native propagation theorem.',
                           'Full native density and both terminal P/Q branches are archived and independently checked; this terminal instrument is distinct from sequential native monitoring.',
                           'Eigenframe sigma_x and full corner algebra are declared apparatus choices; complex correlations are algebraic, not raw C/E joint probabilities.',
                           'Modular parameter has no physical clock or subjective/experienced-duration identification.',
                           'METHOD wording corrected to declared native evolution parameter; no physical-time conversion is supplied.'],
            'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,HERE/'METHOD.md',HERE/'correlations.csv',Path(__file__)]}}


if __name__=='__main__':print(json.dumps(run(),indent=2))
