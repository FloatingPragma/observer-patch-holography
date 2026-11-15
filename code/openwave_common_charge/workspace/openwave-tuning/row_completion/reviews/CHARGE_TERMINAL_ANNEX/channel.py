"""Fresh source-block and retained-reference tests; no producer evaluator."""
from audit import *


def trace_record(state,cellset,clock_effect):
    result=np.zeros((state.shape[3],state.shape[3]),complex)
    for c in cellset:
        b=state[c]
        result+=np.einsum('ikrj,lk,ilsj->rs',b,clock_effect,b.conj())
    return result


def internal_marginal(state):
    return np.einsum('cikrj,clkrj->il',state,state.conj())


def main():
    compiler=custody();g=compiler.spatial.graph(2)
    coordinates=list(map(tuple,g.cell_coordinates));lookup={x:i for i,x in enumerate(coordinates)}
    direct=[]
    archived=read(frozen('openwave-tuning/next/CHAIN5/common_charge_topology/DYNAMICS.json'))['routes']
    for row in archived:
        c=row['cell'];a=row['axis'];sign=row['signed_rail'];x=list(coordinates[c]);x[a]=(x[a]+sign)%2
        destination=lookup[tuple(x)]
        assert destination==c^(1<<(2-a))==row['destination']
        assert row['native_path'][0]==12*c and row['native_path'][-1]==12*destination
        for v,w in zip(row['native_path'],row['native_path'][1:]):assert w in g.neighbours[v]
        direct.append({'cell':c,'axis':a,'sign':sign,'destination':destination})
    assert len(direct)==48

    # Exact source spin projector algebra. Controlled Weyl unitaries commute
    # with these projectors because they act on the field factor.
    I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1]).astype(complex)
    for A in (np.kron(X,X),np.kron(Y,X),np.kron(Z,X)):
        plus=(np.eye(4)+A)/2;minus=(np.eye(4)-A)/2
        assert np.array_equal(plus@plus,plus) and np.array_equal(minus@minus,minus)
        assert np.array_equal(plus@minus,np.zeros((4,4)))
        assert np.array_equal(plus+minus,np.eye(4))

    rng=np.random.default_rng(297143)
    # A mixed correlated input: rank5 after discarding the last purification
    # coordinate, retaining a 3D quantum reference and a correlated clock.
    initial=rng.normal(size=(8,3,2,3,5))+1j*rng.normal(size=(8,3,2,3,5))
    initial/=np.linalg.norm(initial)
    clock=np.diag(np.exp([-1j/128,1j/128]))
    actual=initial.copy();target=initial.copy();errors=[]
    observables=(X,Y,Z)
    for n in (1,2):
        for axis in range(3):
            out=np.zeros_like(actual)
            for c in range(8):
                M=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
                Q,R=np.linalg.qr(M)
                Q=Q@np.diag(np.diag(R)/abs(np.diag(R)))
                out[c^(1<<(2-axis))]=np.einsum('ij,jkrl->ikrl',Q,actual[c])
            actual=out
        actual=np.einsum('ab,cibrl->ciarl',clock,actual)
        target=target[np.arange(8)^7]
        target=np.einsum('ab,cibrl->ciarl',clock,target)
        maximum=0.
        # Every cell subset, every binary Pauli clock setting, and both
        # outcomes. Reference matrices are compared, not just trace weights.
        for mask in range(256):
            cells=[c for c in range(8) if mask>>c&1]
            for O in observables:
                for sign in (-1,1):
                    effect=(I+sign*O)/2
                    a=trace_record(actual,cells,effect);b=trace_record(target,cells,effect)
                    maximum=max(maximum,float(np.linalg.norm(a-b,2)))
        assert maximum<2e-13
        errors.append({'macros':n,'retained_reference_matrix_comparisons':256*6,
                       'maximum_reference_channel_error':maximum})
    retained_difference=float(np.linalg.norm(internal_marginal(actual)-internal_marginal(target),2))
    assert retained_difference>1e-4
    # L3 is a rejecting control: the two rails no longer have one destination.
    plus=(0+1)%3;minus=(0-1)%3
    assert plus!=minus
    clock_bad=abs(math.sin(F(1,32))-math.sin(F(1,64)))
    assert clock_bad>F(1,100)
    source_files={}
    for name in ('openwave-tuning/common_dynamics/interacting_gaussian/branches.py',
                 'openwave-tuning/common_dynamics/spin/dirac_proposal.py'):
        p=frozen(name);pin(p);source_files[str(p.relative_to(ROOT))]=sha(p)
    save('CHANNEL_AUDIT.json',{'actual_signed_route_destinations':direct,
        'exact_spin_projector_controls':3,'source_files_read':source_files,
        'fresh_mixed_entangled_clock_reference_checks':errors,
        'retained_internal_state_equality_rejected_operator_gap':retained_difference,
        'L3_equal_destination_claim_rejected':{'plus':plus,'minus':minus},
        'wrong_clock_factor_two_Y_mean_gap':clock_bad,
        'theorem_domain':'L2, oneparticle, full g1 internal unitary source, terminal cell-diagonal charge and independent source clock; trace cell/internal outputs, retain arbitrary quantum reference.',
        'generic_internal_unitaries_are_diagnostics_not_replacements_for_actual_field_source':True,
        'all_checks_pass':True})
    save('CHANNEL_SOURCE_PINS.json',PINS)
    print('CHANNEL PASS',errors,'retained internal difference',retained_difference,flush=True)


if __name__=='__main__':main()
