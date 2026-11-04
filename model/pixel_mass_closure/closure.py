"""Declared feedback maps, with explicit census and scale choices.

The imported interval machinery supplies outward rounding, implicit-root
checks and infinite heat-sum tails. No measured mass or coupling is read.
"""
from functools import lru_cache
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BACKEND = 'code/P_derivation/interval_contraction_certificate.py'


@lru_cache(maxsize=1)
def backend():
    spec = importlib.util.spec_from_file_location('_pixel_mass_interval', ROOT/BACKEND)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C = backend()
DESIGNS = ('mssm_beta4_proxy', 'sm_color_beta7')


def configure(b, design):
    if design not in DESIGNS:
        raise ValueError('unknown declared design')
    if design == 'sm_color_beta7':
        b.b_mssm = (b.num(41)/b.num(10), -b.num(19)/b.num(6), -b.num(7))
        b.beta_ew = b.num(7)
    return b


def mass_coefficients(b, pixel, alpha_u):
    """The same conditional MCPR dimensionless law, on arbitrary coordinates."""
    u = alpha_u
    kappa = (u/b.num(50))/(b.one+u/b.num(31)+u*u/b.num(310))
    chi = (u*u/b.num(512))/(b.one-u/b.num(77))
    alpha_p = (pixel-b.phi)/b.sqrt_pi
    zeta = (alpha_p*alpha_p/b.num(21))/(b.one-u/b.num(27)-u*u/b.num(135))
    rho = b.sqrt(b.two)*b.exp(-chi)
    delta = b.num(2)/b.num(9)+zeta
    squares = [b.pow_int(b.one+rho*b.cos(delta+b.two_pi*b.num(j)/b.num(3)),2)
               for j in (1,2,0)]
    scale = b.exp((b.num(3)*kappa-b.log(b.num(2)*b.num(6**14))
                   -sum((b.log(x) for x in squares),b.zero))/b.num(3))
    return [scale*x for x in squares]


def transport(b, pixel, alpha_u, mz, v, a3, *, include_top=True):
    leptons = [v*x for x in mass_coefficients(b,pixel,alpha_u)]
    lepton = sum((C.transport_kernel_exact(b,mz,m,b.one,b.one) for m in leptons),b.zero)
    quarks = C.quark_masses(b,v,{'n_u':(6,3,0),'n_d':(6,4,2)})
    quark = b.zero
    for name,numerator in (('u',4),('d',1),('s',1),('c',4),('b',1),('t',4)):
        if name == 't' and not include_top:
            continue
        quark += C.transport_kernel_exact(b,mz,quarks[name],b.num(numerator)/b.num(9),b.num(3))
    return lepton+(b.one-b.num(3)*a3/b.pi)*quark


class Point(C.PointSolver):
    def __init__(self, design, dps=50, cutoff=18, include_top=True):
        self.design = design
        self.include_top = include_top
        super().__init__(configure(C.MpBackend(dps),design),cutoff,cutoff)

    def solve_mz(self, alpha_u, p, mu_u):
        if self.design != 'sm_color_beta7':
            return super().solve_mz(alpha_u,p,mu_u)
        b = self.b
        v = C.v_transmutation(b,alpha_u,p)
        # Every running inverse coupling is positive on this bracket.
        Lmax = b.two_pi/(b.num(7)*alpha_u)
        lo,hi = mu_u*b.exp(-Lmax),mu_u
        f = lambda mu:C.h_mz(b,mu,alpha_u,mu_u,v)
        return self._illinois(f,lo,hi,f(lo),f(hi))

    def solve_alpha_u(self,p):
        if self.design != 'sm_color_beta7':
            return super().solve_alpha_u(p)
        b = self.b
        def equations(x,y):
            A = (b.num(42)*x-b.num(19)*y)/b.num(23)
            L = b.num(12)*b.pi*(x-y)/b.num(23)
            r = b.sqrt(b.pi*(b.one/x+b.num(69)/(b.num(333)*x-b.num(218)*y)))
            scale = L+b.two_pi-b.num(2)/b.num(3)*b.log(p)-b.two_pi*A/b.num(7)+b.log(r)
            heat = C.ellbar(b,self.su2_terms,b.four_pi_sq/x)+C.ellbar(b,self.su3_terms,b.four_pi_sq/y)-p/b.num(4)
            return scale,heat
        x,y = b.ctx.findroot(equations,(b.num('30.5'),b.num('5.5')),tol=self.rel_tol**2)
        if not 0 < y < x:
            raise C.CertificateError('positive ordered SM running branch required')
        return b.num(23)/(b.num(42)*x-b.num(19)*y),C.mu_u_of_p(b,p)

    def inner(self,alpha,mode=None):
        b = self.b
        p = b.phi+alpha*b.sqrt_pi
        u,mu = self.solve_alpha_u(p)
        mz = self.solve_mz(u,p,mu)
        a = [C.alpha_run(b,u,coef,mz,mu) for coef in b.b_mssm]
        if min(a) <= 0:
            raise C.CertificateError('all three running couplings must be positive')
        v = C.v_transmutation(b,u,p)
        a0 = C.alpha_em_inv(b,a[0],a[1])
        delta = transport(b,p,u,mz,v,a[2],include_top=self.include_top)
        return {'p':p,'alpha_u':u,'mu_u':mu,'mz':mz,'v':v,
                'a1':a[0],'a2':a[1],'a3':a[2],'a0_inv':a0,'delta':delta,
                'alpha_inv':a0+delta,'g':b.one/(a0+delta)}


class Chain(C.IntervalChain):
    def __init__(self, point, dps=50):
        iv = configure(C.IvBackend(dps),point.design)
        super().__init__(iv,point,point.su2_cutoff,point.su3_cutoff)
        configure(self.dual,point.design)

    def g_dual(self, alpha_iv, mode=None):
        iv,d,mpb = self.ivb,self.dual,self.point.b
        p = d.phi+C.Dual(alpha_iv,iv.one)*d.sqrt_pi
        mu = C.mu_u_of_p(d,p)
        width = C._width_mp(mpb,p.x)
        pad = self.scales['dudp_abs']*width*mpb.num(8)+self.scales['u_floor']
        p_mid = C._mid_mp(mpb,p.x)
        candidate,_ = self.point.solve_alpha_u(p_mid)
        u,u_checks = self.verified_alpha_u(p,mu,candidate,pad)
        mz_candidate = self.point.solve_mz(candidate,p_mid,C.mu_u_of_p(mpb,p_mid))
        mz,mz_checks = self.mz_dual(u,p,mu,mz_candidate)
        a = [C.alpha_run(d,u,coef,mz,mu) for coef in d.b_mssm]
        if any(x.x.a <= 0 for x in a):
            raise C.CertificateError('running coupling enclosure is not positive')
        a0 = C.alpha_em_inv(d,a[0],a[1])
        v = C.v_transmutation(d,u,p)
        delta = transport(d,p,u,mz,v,a[2],include_top=self.point.include_top)
        denominator = a0+delta
        if denominator.x.a <= 0:
            raise C.CertificateError('feedback inverse coupling must be positive')
        return d.one/denominator, {'alpha_u':C._iv_pair(u.x),'mz':C._iv_pair(mz.x),
            'v':C._iv_pair(v.x),'masses_over_v':[C._iv_pair(x.x) for x in mass_coefficients(d,p,u)],
            'a0_inv':C._iv_pair(a0.x),'delta':C._iv_pair(delta.x),
            'alpha_inv':C._iv_pair(denominator.x),
            'alpha_u_verification':u_checks,'mz_verification':mz_checks}


def certificate(design,dps=50,cutoff=18):
    point = Point(design,dps,cutoff)
    b = point.b
    alpha = point.fixed_point(None)
    row = point.inner(alpha)
    chain = Chain(point,dps)
    p,u,mz,mu = (row[k] for k in ('p','alpha_u','mz','mu_u'))
    h = b.num('1e-13')
    du = abs(point.solve_alpha_u(p+h)[0]-point.solve_alpha_u(p-h)[0])/(2*h)
    dm_p = abs(point.solve_mz(u,p+h,C.mu_u_of_p(b,p+h))-point.solve_mz(u,p-h,C.mu_u_of_p(b,p-h)))/(2*h*mz)
    dm_u = abs(point.solve_mz(u+h,p,mu)-point.solve_mz(u-h,p,mu))/(2*h*mz)
    floor = b.num(10)**(-(dps-14))
    chain.scales = {'dudp_abs':du+b.num('1e-6'),'dmz_rel_dp':dm_p+b.num('1e-6'),
                    'dmz_rel_du':dm_u+b.num('1e-6'),'u_floor':floor,'mz_floor':floor}
    iv = chain.ivb
    radius = b.num('1e-14')
    domain = iv.hull(iv.thin(alpha-radius).a,iv.thin(alpha+radius).b)
    whole,checks = chain.g_dual(domain)
    middle,_ = chain.g_dual(iv.thin(alpha))
    image = middle.x+whole.d*(domain-iv.thin(alpha))
    lipschitz = iv.sup_abs(whole.d)
    if not lipschitz < iv.one.a or not domain.a < image.a or not image.b < domain.b:
        raise C.CertificateError('Banach contraction or strict self-map failed')
    display = {key:b.ctx.nstr(value,32) for key,value in row.items()}
    return {'design':design,'backend':'private mpmath binary outward interval arithmetic',
            'decimal_precision':dps,'heat_cutoff':cutoff,'infinite_heat_tails_included':True,
            'uniqueness_domain_alpha':C._iv_pair(domain),'fixed_point_alpha':C._iv_pair(image),
            'fixed_point_P':C._iv_pair(iv.phi+iv.sqrt_pi*image),
            'fixed_point_alpha_inv':C._iv_pair(iv.one/image),
            'derivative':C._iv_pair(whole.d),'lipschitz_upper':C._iv_pair(lipschitz)['hi'],
            'strict_self_map':True,'conditional_unique_root_in_domain':True,
            'implicit_and_mass_enclosures':checks,'display':display}
