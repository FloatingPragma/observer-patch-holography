"""Independent point replay by simultaneous inverse-coupling equations.

This module does not import the producer, its implicit solver or its closed
transport kernel. It evaluates the fermion transport by numerical integration.
The producer's interval proof and this numerical crosscheck have distinct roles.
"""
from mpmath.ctx_mp import MPContext


def recompute(design):
    m = MPContext(); m.dps = 55
    phi = (1+m.sqrt(5))/2
    if design == 'sm_color_beta7':
        beta=(m.mpf(41)/10,-m.mpf(19)/6,-m.mpf(7)); k=7
        seed=('51.2','41.2','.00592')
    elif design == 'mssm_beta4_proxy':
        beta=(m.mpf(33)/5,m.mpf(1),-m.mpf(3)); k=4
        seed=('24.3','33.2','.00730')
    else:
        raise ValueError('unknown design')
    def heat(a,rank):
        if rank == 2:
            terms=[(n+1,m.mpf(n*(n+2))/4) for n in range(19)]
        else:
            terms=[((i+1)*(j+1)*(i+j+2)//2,m.mpf(i*i+j*j+i*j+3*i+3*j)/3)
                   for i in range(19) for j in range(19)]
        weights=[dim*m.exp(-4*m.pi*m.pi*a*c) for dim,c in terms]
        return m.fsum(w*m.log(dim) for w,(dim,c) in zip(weights,terms))/m.fsum(weights)
    def source(alpha,A,L):
        P=phi+alpha*m.sqrt(m.pi)
        gauge=[1/(A+b*L/(2*m.pi)) for b in beta]
        ratio=m.sqrt(m.pi*(gauge[1]+3*gauge[0]/5))
        u=1/A
        kap=(u/50)/(1+u/31+u*u/310)
        chi=(u*u/512)/(1-u/77)
        zeta=(alpha*alpha/21)/(1-u/27-u*u/135)
        shape=sorted((1+m.sqrt(2)*m.exp(-chi)*m.cos(m.mpf(2)/9+zeta+2*m.pi*j/3))**2 for j in range(3))
        determinant=m.exp(3*kap)/(2*6**14)
        scale=m.root(determinant/m.fprod(shape),3)
        masses=[scale*r for r in shape]
        quarks=[m.power(6,-n)/m.sqrt(2) for n in (6,3,0,6,4,2)]
        charges=[m.mpf(4)/3]*3+[m.mpf(1)/3]*3
        def kernel(mass):
            z=(ratio/mass)**2
            return 2/m.pi*m.quad(lambda x:x*(1-x)*m.log1p(z*x*(1-x)),[0,m.mpf('.001'),m.mpf('.1'),m.mpf('.5'),m.mpf('.9'),m.mpf('.999'),1])
        # Root iteration uses an independent algebraic primitive; the final
        # transport is checked against quadrature below.
        def closed(mass):
            a=(ratio/mass)**2/4
            integral=-m.mpf(5)/18+1/(6*a)+(2*a-1)*m.sqrt(1+a)*m.asinh(m.sqrt(a))/(6*a*m.sqrt(a))
            return 2*integral/m.pi
        lep=m.fsum(closed(x) for x in masses)
        quark=m.fsum(q*closed(x) for q,x in zip(charges,quarks))
        readout=1/gauge[1]+5/(3*gauge[0])+lep+(1-3*gauge[2]/m.pi)*quark
        return P,gauge,ratio,masses,quarks,charges,readout,kernel
    def residual(A,L,alpha):
        P,gauge,r,*tail=source(alpha,A,L)
        readout=tail[-2]
        return (L+2*m.pi-m.mpf(2)/3*m.log(P)-2*m.pi*A/k+m.log(r),
                heat(gauge[1],2)+heat(gauge[2],3)-P/4,
                alpha*readout-1)
    A,L,alpha=m.findroot(residual,tuple(map(m.mpf,seed)),tol=m.mpf('1e-45'))
    P,gauge,r,masses,quarks,charges,readout,kernel=source(alpha,A,L)
    quadrature=1/gauge[1]+5/(3*gauge[0])+m.fsum(kernel(x) for x in masses)+(1-3*gauge[2]/m.pi)*m.fsum(q*kernel(x) for q,x in zip(charges,quarks))
    if not min(gauge)>0 or not L>0 or not abs(quadrature-readout)<m.mpf('1e-40'):
        raise ArithmeticError('independent positive-branch quadrature check failed')
    return {'P':m.nstr(P,40),'alpha':m.nstr(alpha,40),'alpha_inv':m.nstr(readout,40),
            'alpha_U':m.nstr(1/A,40),'masses_over_v':[m.nstr(x,40) for x in masses],
            'simultaneous_residual_upper':m.nstr(max(abs(x) for x in residual(A,L,alpha)),8),
            'quadrature_difference_upper':m.nstr(abs(quadrature-readout),8)}


def verify_report(report):
    m=MPContext();m.dps=50
    if len(report['closures']) != 2:
        raise ValueError('two declared closures required')
    output={}
    for row in report['closures']:
        result=recompute(row['design'])
        for name,field in (('P','fixed_point_P'),('alpha','fixed_point_alpha'),('alpha_inv','fixed_point_alpha_inv')):
            if not m.mpf(row[field]['lo']) < m.mpf(result[name]) < m.mpf(row[field]['hi']):
                raise ValueError('independent simultaneous root escapes '+field)
        for value,box in zip(result['masses_over_v'],row['implicit_and_mass_enclosures']['masses_over_v']):
            if not m.mpf(box['lo']) < m.mpf(value) < m.mpf(box['hi']):
                raise ValueError('independent mass escapes interval')
        output[row['design']]=result
    return output
