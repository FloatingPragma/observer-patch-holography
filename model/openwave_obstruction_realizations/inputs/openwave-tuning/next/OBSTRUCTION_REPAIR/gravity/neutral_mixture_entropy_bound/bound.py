import numpy as np

def low_rank_log_Z(C,X,signs):
    e,U=np.linalg.eigh(C)
    if np.min(e)<=0 or np.max(e)>=1:raise ValueError('Faithful regional reference required; no eigenvalue floor substituted.')
    Y=U.conj().T@X
    G0=Y.conj().T@(Y/(e*(1-e))[:,None])
    G1=Y.conj().T@(Y/(1-e)[:,None])
    G2=Y.conj().T@(Y*(e/(1-e))[:,None])
    return log_Z_from_grams(X.conj().T@X,G0,G1,G2,signs)

def log_Z_from_grams(M,G0,G1,G2,signs):
    J=np.diag(signs);r=len(signs)
    small=np.eye(2*r)+np.block([[J@G1,J@G0],[-J@G2+J@M@J@G1,-J@G1+J@M@J@G0]])
    phase,logabs=np.linalg.slogdet(small)
    condition=float(np.linalg.cond(small))
    reliable=bool(phase.real>1-1e-7 and abs(phase.imag)<1e-7 and logabs>=-1e-9 and condition<1e12)
    return dict(log_Z=float(logabs),determinant_phase=[float(phase.real),float(phase.imag)],small_condition_number=condition,rank_bound=r,numerical_gate_pass=reliable,interval_arithmetic_certified=False)

def finite_entropy_interval(J,epsilon,log_Z):
    # log(1+epsilon^2 (exp(log_Z)-1)), evaluated stably.
    log_remainder=float(np.logaddexp(np.log1p(-epsilon*epsilon),2*np.log(epsilon)+log_Z))
    log_remainder=max(0.,log_remainder)
    return dict(epsilon=epsilon,linear_entropy=epsilon*J,relative_entropy_upper_bound=log_remainder,entropy_lower=epsilon*J-log_remainder,entropy_upper=epsilon*J)
