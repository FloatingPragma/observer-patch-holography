"""Interval Taylor transport and a Volterra remainder, with no fitted target."""
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import sympy as sy
from model.openwave_radiation_capacity.source import Interval as I, exp_point, sigmoid, taylor_step, root_interval, log_interval
from model.openwave_capacity_normalization_transfer.source import integrand, validated
from model.openwave_capacity_transfer.source import background_inputs, constraint_identity
from model.openwave_capacity_stress.source import enclosure

ROOT=Path(__file__).resolve().parents[2]
ZERO=I.point(0)
ONE=I.point(1)

def absolute(v): return max(abs(v.lo),abs(v.hi))
def mid(v): return (v.lo+v.hi)/2
def box(values): return I(min(x.lo for x in values),max(x.hi for x in values))
def mv(matrix,vector): return [sum((a*b for a,b in zip(row,vector)),ZERO) for row in matrix]
def add(a,b): return [x+y for x,y in zip(a,b)]
def polynomial(poly,x,r):
    powers=lambda a,n: ONE if n==0 else a*powers(a,n-1)
    return sum((F(c)*powers(x,i)*powers(r,j) for (i,j),c in poly.terms()),ZERO)

@lru_cache(None)
def mass_polynomials():
    x,r=sy.symbols('x r')
    data=json.loads((ROOT/'model/openwave_capacity_perturbations/results.json').read_text())
    p=sy.Poly(sy.sympify(data['canonical_scalar_completion']['V_second_derivative_over_H_squared_polynomial'],locals={'x':x,'r':r}),x,r)
    return p,p.diff(x),p.diff(r)

def reduction():
    m,r,s,b,h,A,M,hN,Q,U,R,v,w,j,Dm,Dr=sy.symbols('m r s b h A M hN Q U R v w j Dm Dr')
    Phi=R-U;Up=(h-1)*U+R;P=-Phi+(3*m*U+4*r*(U+v)+s*(w+s*U))/2
    Rp=sy.expand(P+Up).subs(h,-(3*m+4*r+s*s)/2).expand()
    y=w+s*U;z=j+b*U+s*Up
    zp=-(3+h)*z-(Q+M)*y+4*s*P-2*A*Phi
    bn=-hN*s-(3+h)*b-M*s+2*h*A
    jp=sy.expand(zp-bn*U-2*b*Up-s*(hN*U+(h-1)*Up+Rp))
    expected=-(3+h)*j-(M-sy.Rational(3,2)*s*s+Q)*w+6*r*s*v-Q*s*U
    residual=sy.expand((jp-expected).subs(b,-(3+h)*s-A).subs(h,-(3*m+4*r+s*s)/2))
    assert residual==0 and sy.expand(Rp-2*r*v-s*w/2)==0
    assert F(8,27)<F(11,20)**2 and 2-F(3,2)*F(8,27)>0
    assert 1+F(3,2)+F(4,27)+F(11,20)<4
    assert -1+3+F(3,2)*F(11,20)<4
    om=F(background_inputs()['Omega_m0']);assert F(5,3)/om<6
    return {'scalar_comoving_equation':'w_N=j; j_N=-(3+h)j-(M-3s^2/2+Q)w+6rs v-Qs U',
        'curvature_equation':'R_N=2 Omega_r v+s w/2; R=Phi+U_m',
        'symbolic_residuals':['0','0'],'weighted_coordinates':['Delta_m','Delta_r','4v','2w','j','U_m','R'],
        'A0_logarithmic_norm_upper_exact':'4','B_operator_norm_upper_exact':'6',
        'norm_domain':'0<=x<=1,0<=radiation_mix<=1,M in[2,6],s^2<=8/27,-2<=h<=0,a<=1; supplied Omega_m0>5/18.'}

def coefficients(t,w,om,orr):
    a=exp_point(-t) if not isinstance(t,I) else I(exp_point(-t.hi).lo,exp_point(-t.lo).hi)
    omega=I(sigmoid(w.lo).lo,sigmoid(w.hi).hi)
    x=root_interval(omega);mix=I.point(orr)/(om*a+orr)
    m=(ONE-mix)*(ONE-omega);r=mix*(ONE-omega)
    s=root_interval(2*omega*(ONE-x))
    h=-(3*m+4*r+s*s)/2
    xp=x*(ONE-x*x)*(ONE+mix+2*x)/2;rp=-mix*(ONE-mix)
    sp=s*xp*(ONE/x-ONE/(2*(ONE-x)))
    mp=(-3*ONE-2*h)*m;rrp=(-4*ONE-2*h)*r
    hp=-(3*mp+4*rrp+2*s*sp)/2
    p,px,pr=mass_polynomials();mass=polynomial(p,x,mix)
    massp=polynomial(px,x,mix)*xp+polynomial(pr,x,mix)*rp
    f=a*a*(ONE-omega)/(om*a+orr);fp=-2*(ONE+h)*f
    return m,r,s,h,mass,f,mp,rrp,sp,hp,massp,fp

def matrices(values):
    m,r,s,h,M,f,mp,rp,sp,hp,Mp,fp=values
    a=[[ZERO for _ in range(8)]for _ in range(8)]
    d=[[ZERO for _ in range(8)]for _ in range(8)]
    def put(i,j,v,dv):a[i][j]=v;d[i][j]=dv
    # (Delta_m,Delta_r,4v,2w,j,U,R) at order k^2; coordinate8 is10T.
    put(0,2,F(3,2)*r,F(3,2)*rp);put(0,3,F(3,4)*s,F(3,4)*sp)
    put(0,7,-f/10,-fp/10)
    put(1,1,ONE,ZERO);put(1,2,-(6*m+2*s*s)/4,-(6*mp+4*s*sp)/4)
    put(1,3,s,sp);put(1,7,-2*f/15,-2*fp/15)
    put(2,1,ONE,ZERO);put(2,2,h,hp);put(3,4,2*ONE,ZERO)
    put(4,2,F(3,2)*r*s,F(3,2)*(rp*s+r*sp))
    put(4,3,-(M-F(3,2)*s*s)/2,-(Mp-3*s*sp)/2)
    put(4,4,-3*ONE-h,-hp);put(4,7,-f*s/10,-(fp*s+f*sp)/10)
    put(5,5,h-ONE,hp);put(5,6,ONE,ZERO)
    put(6,2,r/2,rp/2);put(6,3,s/4,sp/4)
    put(7,7,h-ONE,hp)
    b=[ZERO]*7+[I.point(15)]
    return a,d,b

def lognorm(matrix):
    return max(row[i].hi+sum(absolute(v)for j,v in enumerate(row)if i!=j) for i,row in enumerate(matrix))

@lru_cache(None)
def start_data():
    data=background_inputs();om=F(data['Omega_m0']);orr=F(data['Omega_r0']);ox=1-om-orr
    # The parent quadrature is reused on the exact log-time endpoint t=1.
    step=F(1,512);n=512;w=I(*log_interval(ox/(1-ox)));previous=F(0);rows=[]
    for j in range(12*n+1):
        t=j*step
        if j:w=taylor_step(previous,w,step,om,orr)
        rows.append((t,w,integrand(t,w,'capacity',om,orr,ox)));previous=t
    tail=exp_point(-12).hi**3/(3*root_interval(I.point(orr)).lo)
    integral=I(0,tail)
    for i in range(12*n-1,n-1,-1):
        error=F(3,4)*step**3*rows[i][2].hi
        integral=integral+step*(rows[i][2]+rows[i+1][2])/2+I(-error,error)
    a=exp_point(-1);w1=rows[n][1];omega=I(sigmoid(w1.lo).lo,sigmoid(w1.hi).hi)
    H=root_interval((I.point(om)/(a*a*a)+I.point(orr)/(a*a*a*a))/(ONE-omega))
    T=F(3,2)*H/a*integral;Phi=I.point(F(3,2))-T
    values=coefficients(F(1),w1,om,orr);m,f=values[0],values[5]
    v=[-2*f*Phi/(3*m)]+[ZERO]*6+[10*T]
    assert absolute(v[0])<8
    return om,orr,ox,v,T,Phi

def step_state(center,radius,m0,mb,db,b,h):
    y=list(map(I.point,center));norm=max(sum(absolute(v)for v in row)for row in mb)
    size=max(map(abs,center));bn=max(map(absolute,b))
    sup=(size+h*bn)*exp_point(norm*h).hi
    excursion=h*(norm*sup+bn)
    tube=[I(v-excursion,v+excursion)for v in center]
    slope=add(mv(m0,y),b)
    second=add(mv(db,tube),mv(mb,add(mv(mb,tube),b)))
    image=[v+h*f+h*h*g/2 for v,f,g in zip(y,slope,second)]
    nxt=list(map(mid,image));local=max((v.hi-v.lo)/2 for v in image)
    mu=lognorm(mb)
    full_tube=[I(v-excursion-exp_point(max(F(0),mu)*h).hi*radius,v+excursion+exp_point(max(F(0),mu)*h).hi*radius)for v in center]
    return nxt,I.point(exp_point(mu*h).hi*radius+local).hi,mu,local,full_tube

@lru_cache(None)
def transport(steps=2048):
    if type(steps)is not int or steps<128 or steps>4096:raise ValueError('128..4096 steps required')
    om,orr,ox,initial,T,Phi=start_data();h=F(1,steps)
    # Build complete endpoint and cell background enclosures from the present.
    w=I(*log_interval(ox/(1-ox)));nodes=[w]
    for i in range(steps):w=taylor_step(i*h,w,h,om,orr);nodes.append(w)
    center=list(map(mid,initial));radius=max((v.hi-v.lo)/2 for v in initial)
    mus=[];local_total=F(0);forcing=[]
    for i in range(steps,0,-1):
        # N increases, while t=-N decreases. w decreases monotonically in t.
        wt=I(nodes[i].lo,nodes[i-1].hi);time=I((i-1)*h,i*h)
        m0,_,b=matrices(coefficients(i*h,nodes[i],om,orr))
        vals=coefficients(time,wt,om,orr);mb,db,_=matrices(vals)
        center,radius,mu,local,tube=step_state(center,radius,m0,mb,db,b,h)
        f,s=vals[5],vals[2]
        forcing.append(max(map(absolute,[-f*tube[5],-f*(tube[2]/3+4*tube[5]/3),-f*(tube[3]/2+s*tube[5])])))
        mus.append(mu);local_total+=local
    # Exact finite-k system: y'=A0 y+epsilon B y. The weighted7-coordinate
    # logarithmic norm of A0 is<=4; ||B||<=6 on the complete past branch.
    # ||y0||<=3/2, ||y1(0)||<=8, hence ||y1||<=e^(4t)(8+9t).
    # The k^4 remainder is at most102 epsilon^2 exp(4+6epsilon) for0<=t<=1.
    rows=[];zero=validated()['records'][0]['capacity']['metric_transfer_interval'];phi0=I(*map(F,zero))
    cp=I(center[6]-center[5]-2*radius,center[6]-center[5]+2*radius)
    dm=I(center[0]-radius,center[0]+radius)
    for k in(F(1,100),F(1,50)):
        eps=k*k;remainder=F(0);growth=exp_point((4+6*eps)*h).hi
        for force in forcing:remainder=I.point(growth*(remainder+eps*eps*h*force)).hi
        assert remainder<=102*eps*eps*exp_point(4+6*eps).hi
        metric=phi0+eps*cp+I(-2*remainder,2*remainder)
        density=eps*dm+I(-remainder,remainder)
        rows.append({'k_over_H0_exact':str(k),'metric_transfer_interval':metric.encoded(),
            'metric_transfer_display':enclosure((metric.lo,metric.hi)),
            'comoving_dust_interval':density.encoded(),'comoving_dust_display':enclosure((density.lo,density.hi)),
            'weighted_state_k4_remainder_upper_exact':str(remainder)})
    assert max(mus)<4
    return {'steps':steps,'initial_log_scale_factor_exact':'-1','final_log_scale_factor_exact':'0',
        'coordinate_order':['Delta_m_1','Delta_r_1','4_relative_velocity_1','2_comoving_scalar_1','comoving_scalar_derivative_1','U_m_1','curvature_1','10_T_0'],
        'first_k_squared_coefficient_center_exact':list(map(str,center)),
        'common_weighted_coefficient_radius_exact':str(radius),
        'metric_k_squared_coefficient_display':enclosure((cp.lo,cp.hi)),
        'dust_k_squared_coefficient_display':enclosure((dm.lo,dm.hi)),
        'maximum_cell_logarithmic_norm_upper_exact':str(max(mus)),
        'maximum_k4_forcing_norm_upper_exact':str(max(forcing)),
        'k4_remainder_method':'Cellwise variation of constants with complete enclosed B*y1 forcing, growth exp((4+6*k^2)*h); compared to the independent global102*k^4*exp(4+6*k^2) bound.',
        'initial_T_interval':T.encoded(),'initial_Phi_interval':Phi.encoded(),'finite_wavelength_records':rows}

def run():
    exact=reduction();fine=transport(2048);coarse=transport(1024)
    assert F(fine['common_weighted_coefficient_radius_exact'])<F(coarse['common_weighted_coefficient_radius_exact'])
    for a,b in zip(fine['finite_wavelength_records'],coarse['finite_wavelength_records']):
        for key in ('metric_transfer_interval','comoving_dust_interval'):
            aa=list(map(F,a[key]));bb=list(map(F,b[key]));assert max(aa[0],bb[0])<=min(aa[1],bb[1])
    return {'schema':'oph.capacity_finite_wavelength.v1','verdict':'PASS_VALIDATED_LATE_FINITE_WAVELENGTH_TRANSFER',
        'inputs':background_inputs(),'constraint_identity':constraint_identity(),'exact_reduction':exact,
        'fine':fine,'coarse':coarse,
        'preparation':'At N=-1 the exact regular k=0 solution is supplemented by Delta_m=-2Q Phi/(3 Omega_m); all other comoving corrections vanish. This supplied late-time density correction satisfies the full Hamiltonian constraint exactly.',
        'finite_theorem_shortcut':'Comoving scalar and relative-velocity coordinates annihilate the exact zero-k mode. An8-coordinate affine interval transport gives the k^2 response; a Volterra norm bound supplies the complete k^4 remainder.',
        'scope':{'late_preparation_is_early_primordial_solution':False,'finite_nonzero_wavelength_error_certified':True,
          'native_capacity_horizon_clock_identified':False,'absolute_amplitude_predicted':False,'new_Nature_validation_rows':0,
          'supplied_Einstein_scalar_fluid_equations':True,'full_Boltzmann_or_CMB_match':False},
        'normalization':'H0->lambda H0, density->lambda^2 density, L->L/lambda and k->lambda k leave every dimensionless result invariant. Primordial amplitude and native physical horizon/clock identification are independent inputs.'}
