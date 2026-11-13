# Conserved radiation in capacity-history comparisons

This calculation extends the parent's flat, Planck-anchored Friedmann
comparison. H0=67.4 and Omega_m0=0.315 are inherited empirical inputs.
Omega_r0=0.000092 is a declared conserved ideal-radiation input; the remaining
current fraction follows flatness. Radiation species, thermalization,
entropy-degree changes and the amplitude of the capacity law are not derived.
No additional expansion data are fitted.

## Common conservation equations

For the gamma-one future-event-radius choice, the defining radius equation
gives `w_x=-1/3-2sqrt(Omega_x)/3`. Conserved matter and radiation have
rho_m proportional to a^-3 and rho_r proportional to a^-4. Set
`r=Omega_r/(Omega_m+Omega_r)=Omega_r0/(Omega_m0*a+Omega_r0)`.
The continuity and Friedmann equations then give, exactly,

```
dOmega_x/dlog(a)=Omega_x(1-Omega_x)(1+r+2sqrt(Omega_x)).
```

The positive radiation term cannot be omitted at early times. In the
radiation asymptotic regime r→1 and Omega_x→0, the logarithmic slope tends
to2; the matter-era slope is1. Fixed capacity instead has constant rho_x
and radiation-era fraction proportional to a^4. These are consequences of
the declared histories, not a source-level selection of a horizon.

## Validated integration

Use t=log(1+z), w=log(Omega_x/(1-Omega_x)). Then

```
w'=f(t,w)=-(1+r(t)+2sqrt(sigmoid(w))).
```

On the entire domain, |f|<=4, |f_w|<=2/5, |f_t|<=1/4,
|f_ww|<=2/5, |f_tt|<=1/4 and f_tw=0. The bound for f_w follows from
`max_{0<=s<=1} s(1-s²)=2/(3sqrt3)<2/5`; the other inequalities follow by
direct differentiation. Thus |w''|<=37/20 and

```
|w'''| <= 1/4+(2/5)*16+(2/5)*(37/20)=739/100<8.
```

The second-order Taylor map T_h has derivative in w at least
`1-(2/5)h-(4/5)h²>=3/5` for h<=1/2. It is increasing. Applying T_h to
the two interval endpoints, then enlarging each side by `(4/3)h³`, therefore
encloses the exact next solution. Every exponential, square root and log
uses the existing rational enclosure tools. State endpoints are rounded
outward. The final log-redshift interval contributes the additional
`4*(t_upper-t_lower)` allowance. The receipt is an interval calculation,
not an unbounded-tolerance ordinary differential equation fit. An independent
adaptive floating solver checks selected points without supplying the bounds.

At each redshift the relative expansion is
`H/H_m+r=1/sqrt(1-Omega_x)`, where the denominator uses the same supplied
matter and radiation amplitudes. The present dark fraction is an input
boundary condition, not an inferred cosmic amplitude.

## A retained early-time control

If capacity is proportional to the squared Hubble radius, rho_x is
proportional to H² and its fraction is constant. Separate conservation
then forces its pressure to track the matter/radiation mixture. In radiation
domination its extra expansion-equivalent neutrino density is

```
Delta_Neff=Omega_x/(1-Omega_x)*(1+c_nu*N_eff)/c_nu.
```

The packet retains the declared neutrino conventions and this large control
effect. It is not a CMB perturbation or nucleosynthesis likelihood: a fluid
with the same background density need not share neutrino perturbations.
The calculation reaches z=10^6 in an ideal fixed-radiation model. It does
not claim a particle abundance prediction, a treatment of electron-positron
annihilation, or a measured early-Universe match.
