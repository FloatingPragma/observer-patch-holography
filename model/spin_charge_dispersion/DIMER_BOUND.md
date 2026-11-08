# A variational source energy and current bound

The Hamiltonian is the unchanged cubic spin-one-half source
`H_J=J sum_edges(S_x dot S_y-1/4)`, J>0. Let e denote the negative
spin-dot ground energy per site divided by J, so the ground H_J density
is `-J(e+3/4)`. The constant3/4 counts three source bonds per site and
must not be dropped when comparing the two conventions.

## Exact periodic trial state

On every even cubic torus of side at least four, pair consecutive sites
along coordinate1, with the first coordinate even at the first member.
There are N/2 dimer bonds and5N/2 other bonds among all3N distinct source
bonds, including periodic wrap bonds. At a pair whose first member has
positive staggered sign use

    psi=a|up down>-b|down up>,
    a=(sqrt(3)+sqrt(2))/sqrt(10),
    b=(sqrt(3)-sqrt(2))/sqrt(10).

Use its spin-reversed version at a pair whose first sign is negative.
Then a²+b²=1, 2ab=1/5, and the single-site spin mean is staggered with
magnitude sqrt(6)/5. The source's exact two-spin algebra gives

    <S_x dot S_y>_intra=-1/4-(2ab)/2=-7/20,
    <S_x dot S_y>_inter=-(a²-b²)²/4=-6/25.

The second equality uses the explicit tensor-product state on distinct
dimers, not a claim that the source ground factorizes. More generally,
with u=2ab between0 and1, its unshifted spin-dot energy density is

    E_dot,trial/(JN)=-3/4-u/4+5u²/8
                =-31/40+(5/8)(u-1/5)².                 (1)

Thus the exact variational inequality is

    e>=31/40,        inf spec(H_J)/N<=-61J/40.           (2)

The state is not the selected quantum ground. Its use is only an explicit
upper bound on ground energy. Averaging it over source translations,
cubic rotations and global spin rotations preserves its energy if an
invariant trial density is desired. The ground projector used in the
infrared theorem retains its own invariant definition; the anisotropic
dimer pattern is not substituted for that projector.

The selected B-field state has the same bound after B decreases to zero.
The selection term has norm at most BN/2. Comparing its finite ground to
any zero-field trial gives an H_J expectation at most the trial H_J
energy plus BN. After dividing by JN and taking the established ordered
limits, this correction vanishes. No new ground-state identification or
interchange of the field and volume limits is needed.

## The exact current coefficient

The selected source is invariant under rotations about its ordered axis1
and under the inherited cubic spatial symmetries. Write its common bond
correlations as c_parallel and c_transverse. The existing source theorem
defines the current coefficient

    kappa=-(c_parallel+c_transverse),
    e=-3(c_parallel+2c_transverse).

Hence `kappa=e/3+c_transverse`. Each spin component has norm1/2, so
`c_transverse>=-1/4`. Combining this inequality with (2) yields

    kappa>=1/120.                                        (3)

This kappa is the microscopic current and energy-kick coefficient. It is
different from a detector's retained probability. In particular the
matched complete inverse-current sum rule becomes

    integral E^(-1) nu_j(k,dE)=J kappa/2>=J/240.           (4)

Every frequency in the actual response remains present. Equations(3)-(4)
do not determine the distribution of that current weight, a pole speed,
or a vanishing relative spectral width.

## The same infrared proof gains a stronger order margin

The inherited invariant-ground bound is
`b>=e/9-sqrt(e/18) I`. Its right side is increasing for e>=3/4 and
I<2/5, as proved in the original source order argument. Use (2), the
unchanged exact256-mesh infrared upper bound, and an outward upper
enclosure of sqrt(31/720). This gives a positive rational lower bound

    b_d=(31/40)/9-sqrt(31/720)_upper I_upper.

The fixed t=9 biased-ground construction in [the main proof](PROOF.md)
then gives `m>=18 b_d/(1+81 b_d)>11/100`. The executable retains the
full rational values and their strict comparisons. This strengthens the
source theorem; the readout's default m0=1/10 is a conservative common
constant. No measured magnetization, ground energy or particle datum
enters either inequality.
