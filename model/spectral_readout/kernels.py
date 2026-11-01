"""Positive Fejer energy kernels, interval Fourier sums and explicit resolution."""
from fractions import Fraction as F
from functools import lru_cache

from model.certified_reference import interval as iv


def exact(value,name,positive=False):
    if type(value) not in (int,F) or value<0 or (positive and value==0):
        raise ValueError('exact '+('positive ' if positive else 'nonnegative ')+name+' required')
    return F(value)


def count(value):
    if type(value) is not int or not 2<=value<=2**24:
        raise ValueError('integer outcome count in[2,2**24] required')
    return value


def trig_interval(lower,upper):
    """Integer-only sin/cos bounds for any bounded rational angle interval."""
    if type(lower) not in (int,F) or type(upper) not in (int,F) or lower>upper:
        raise ValueError('ordered exact angle interval required')
    midpoint=(F(lower)+F(upper))/2;radius=(F(upper)-F(lower))/2
    divisor=1
    while abs(midpoint)/divisor>1:divisor*=2
    sine,cosine=iv.sine_cosine(midpoint/divisor)
    while divisor>1:
        sine,cosine=2*sine*cosine,cosine*cosine-sine*sine
        divisor//=2
    error=iv.Interval.point(radius).hi
    sine=iv.Interval(max(-iv.SCALE,sine.lo-error),min(iv.SCALE,sine.hi+error))
    cosine=iv.Interval(max(-iv.SCALE,cosine.lo-error),min(iv.SCALE,cosine.hi+error))
    return sine,cosine


@lru_cache(maxsize=4096)
def phase(count_,residue):
    m=count(count_);residue%=m
    # Keep exact cardinal phases, including their signs.
    if (4*residue)%m==0:
        k=(4*residue//m)%4
        return iv.Interval.point((0,1,0,-1)[k]),iv.Interval.point((1,0,-1,0)[k])
    pi=iv.pi();multiple=F(2*residue,m)
    return trig_interval(F(pi.lo,iv.SCALE)*multiple,F(pi.hi,iv.SCALE)*multiple)


def probabilities(characteristics):
    """phi(n*tau)=E exp(-in*tau*E), n=1,...,M-1; return all M effects."""
    m=count(len(characteristics)+1)
    if m>256:raise ValueError('bounded interval calculator supports at most256 outcomes')
    rows=[]
    for k in range(m):
        probability=iv.Interval.point(F(1,m))
        for n,(real,imag) in enumerate(characteristics,1):
            sine,cosine=phase(m,k*n)
            probability+=F(2*(m-n),m*m)*(cosine*real-sine*imag)
        # Nonnegativity and normalization follow from the positive complete
        # kernel theorem, not from replacing or renormalizing observations.
        if probability.hi<0 or probability.lo>iv.SCALE:
            raise ArithmeticError('Fourier data incompatible with a positive energy kernel')
        rows.append(iv.Interval(max(0,probability.lo),min(iv.SCALE,probability.hi)))
    if sum(row.lo for row in rows)>iv.SCALE or sum(row.hi for row in rows)<iv.SCALE:
        raise ArithmeticError('complete kernel probability intervals lose normalization')
    return rows


def energy_centers(count_,tau):
    m=count(count_);tau=exact(tau,'time step',True);pi=iv.pi()
    return [F(2*(k if 2*k<m else k-m),m)/tau*pi for k in range(m)]


def circular_leakage(count_,tau,resolution):
    """P(circular energy error>=resolution), uniform in true energy."""
    m=count(count_);tau=exact(tau,'time step',True)
    width=exact(resolution,'energy resolution',True);angle=tau*width
    if angle>3:raise ValueError('require tau*resolution<=3<pi')
    # pi^2<10. The second estimate sums the two tails of the grid rather
    # than bounding each of the M outcomes by the same worst case.
    bounds=[F(1),10/(m*angle*angle)]
    if m*angle>F(44,7):bounds.append(F(22)/(7*m*angle-44))
    return min(bounds)


def real_energy_leakage(count_,tau,resolution,second_moment,*,support_upper=None):
    """Alias and circular leakage; pi>3 gives a rational Nyquist margin."""
    tau=exact(tau,'time step',True);width=exact(resolution,'resolution',True)
    moment=exact(second_moment,'second moment')
    if tau*width>=3:raise ValueError('positive certified Nyquist margin required')
    margin=3/tau-width
    alias=min(F(1),moment/(margin*margin))
    if support_upper is not None:
        support=exact(support_upper,'absolute energy support bound')
        if support<=margin:alias=F(0)
    circular=circular_leakage(count_,tau,width)
    return {'alias':alias,'circular':circular,'total':min(F(1),alias+circular),
            'certified_Nyquist_interior':margin}


def bin_fourier_error(count_,selected_count,characteristic_error):
    """Error for true POVM laws or a reconstructed event clipped to [0,1].

    Independently noisy Fourier coefficients need not reconstruct positive
    probabilities. The cap by one applies only after clipping that scalar
    event estimate, never to the unrestricted raw Fourier sum.
    """
    m=count(count_)
    if type(selected_count) is not int or not 0<=selected_count<=m:
        raise ValueError('valid number of selected bins required')
    error=exact(characteristic_error,'characteristic error')
    if selected_count in (0,m):return F(0)
    return min(F(1),min(selected_count,m-selected_count)*F(m-1,m)*error)


def clip_event_estimate(value):
    """Classical clipping of one noisy event estimate; not postselection."""
    if type(value) not in (int,F):raise ValueError('exact rational event estimate required')
    return min(F(1),max(F(0),F(value)))
