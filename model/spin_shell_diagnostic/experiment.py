"""A finite source control; its momenta are outside the low-momentum regime."""
import numpy as np
from model.spin_goldstone_sector.experiment import matrices


def finite_control():
    # matrices(0) counts both periodic length-two parallel bonds. Dividing by2
    # gives the12-edge open cube, with exactly the source coupling J=1 per edge.
    h,ops,_=matrices(0)
    h=h/2
    energy,basis=np.linalg.eigh(h)
    if (not np.all(np.isfinite(energy)) or not np.all(np.isfinite(basis))
            or np.linalg.norm(basis.conj().T@basis-np.eye(256))>1e-10
            or np.linalg.norm(h@basis-basis*energy)>1e-10):
        raise ArithmeticError('finite source eigensystem failed')
    gaps=energy-energy[0];ground=basis[:,0]
    if gaps[1]<1e-6:raise ArithmeticError('finite ground not isolated')
    rows=[]
    for character,expected in [(1,'0.008684'),(3,'0.006647'),(7,'0.031631')]:
        a=sum(((-1)**((x&character).bit_count())*ops[x][2] for x in range(8)))/np.sqrt(8)
        vector=a@ground
        weights=abs(basis.conj().T@vector)**2
        if weights[0]>1e-20:raise ArithmeticError('unexpected vacuum weight')
        norm=float(weights.sum());m1=float(weights@gaps);m2=float(weights@(gaps*gaps))
        inverse=float(weights[1:]@(1/gaps[1:]))
        current=h@a-a@h
        if (abs(np.linalg.norm(current@ground)**2-m2)>1e-10
                or abs(np.vdot(vector,current@ground).real-m1)>1e-10):
            raise ArithmeticError('current moments do not match spectral moments')
        defect=m1*inverse/norm**2-1
        if defect<.006 or abs(defect-float(expected))>6e-7:
            raise ArithmeticError('finite nonzero width control failed')
        rows.append({'cube_character':character,'ratio_moment_defect_rounded':expected,
                     'strictly_positive_defect_above':'0.006',
                     'all256eigenvectors_used':True,'current_moment_crosscheck':True})
    return {'graph':'open2x2x2cube','sites':8,'bonds':12,'J':1,
            'full_spin_dimension':256,'rows':rows,
            'momenta_in_p_le1_regime':False,
            'proves_failure_of_infinite_scaling_shell':False,
            'native_pulses_executed':False,
            'finite_floating_control_not_interval_eigenvalue_certificate':True}
