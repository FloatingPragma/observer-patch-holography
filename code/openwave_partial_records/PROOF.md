# Partial OpenWave results on native records

Each section below states one partial result for an OpenWave row, the
construction that carries it, the derivation, one scope sentence and the way
it is checked. Rows 3, 4, 6 and 8 use the sealed 240-port signed source of
[`code/openwave_static_signed_source/`](../openwave_static_signed_source/);
row 20 uses residence records on a declared cubic routing; row 31 uses a
native first-hit feedback controller. Every archived receipt in `evidence/`
is a byte copy of a leaf of its original sealed manifest, and
[`EVIDENCE.json`](EVIDENCE.json) records each binding. The official OpenWave
status of every row is untested.

## Native rule and carrier

On a seam (a,b) native repair records

    d = q_a - q_b,
    j = (d - (2c-1)(|d| mod 2))/2,
    q' = q - j e_a + j e_b,

with coin c. The carrier is the sealed graph of 240 ports and 630 seams,
made of twelve-port icosahedral cells. For one positive unit each seam at
its port has exactly one moving coin, so the one-particle transfer is
T = I - L/1260. The compact charged ray f = (|4> - |5> - |6> + |7>)/2
satisfies Lf = 6f, hence Tf = (209/210) f, with energy
E6 = -1260 log(209/210) = 6.014331228348673 in transfer units.

## Row 3: de Broglie clock

Let A embed a position into its 1260 uniform native event arcs,
(Af)(x,e) = f_x/sqrt(1260), and let S reverse each actual arc: a moving event
(seam, coin) from a to b goes to the event (seam, complementary coin) from b
to a, and a held event goes to itself. On the complete padded domain of
2^19 = 524,288 position/event arcs, S^2 = I and A*A = I, so
W = S(2AA* - I) is unitary. Summing the native records gives A*SA = T.
The native currents of the stationary single-unit source give
D/Var(F) = 1 - lambda = 1/210 directly, with lambda = 209/210.

Put u = Af, s = sqrt(1 - lambda^2) and w = (Su - lambda u)/s. Then u and w
are orthonormal and

    W u = lambda u + s w,      W w = -s u + lambda w,

so W^n u = cos(n theta) u + sin(n theta) w with

    theta = acos(209/210) = 0.0976287750 per reflection word.

The native triangle (0,1,5) gives a fixed internal reference in the same
register: put +1 on its three forward moving arcs and -1 on their reversals.
Every port sum vanishes, so A*v = 0; native reversal gives Sv = -v; hence
Wv = v, and v is orthogonal to u and Su. For chi = (u+v)/sqrt(2) every
ordinary position/event coordinate z has

    p_z(n) = C0 + C1 cos(n theta) + S1 sin(n theta) + C2 cos(2n theta) + S2 sin(2n theta),
    C0 = v_z^2/2 + (u_z^2 + w_z^2)/4,   C1 = v_z u_z,   S1 = v_z w_z,
    C2 = (u_z^2 - w_z^2)/4,             S2 = u_z w_z/2.

At the witness arc (port 5, event 3) C1 = -1/(2 sqrt(7560)), so the theta
harmonic is visible in an ordinary coordinate record. The return probability
of u alone is cos^2(n theta), with frequency 2 theta. All 524,288 outcomes
are retained and their probabilities sum to one at every word count. The
sealed run makes 303,660 literal native calls; its harmonic residual over
129 words is 1.2e-15, and a replay regenerates both evidence arrays byte for
byte.

Scope: the reflection walk and the coherent preparation of chi are declared
memory services on the sealed source, and theta is measured in reflection
words.

## Row 4: particle stability

The relations of row 3 make span{u, Su, v} invariant under W, so by
induction W^n chi lies in that span for every integer n. Its position
support is supp(f), the one-step native neighbours of supp(f) and the
triangle, which together are ports 0 to 11: one twelve-port cell. The
confinement is an exact integer identity at every integer time; no sampled
time enters it.

Scope: the statement covers the invariant span of the declared reflection
walk on the sealed source.

## Row 6: spin

For one positive unit put J|p,s> = c^dagger_(p,s)|vacuum>, with two internal
states s. Under the current-only (C) record law the raw record is the native
signed-current code, and a single particle has no occupied intervening mode,
so every actual CAR hop has crossing sign +1. Hence

    V_h J = J (K_h tensor I_2)

for every seam, coin and current h, held and zero-current columns included,
where K_h is the signed-coordinate branch. Fixed-current injectivity makes
the branch effects sum to I. The field spin
S^a = sum_p c^dagger_(p,t) (sigma^a/2)_(t,s) c_(p,s) acts as I tensor
sigma^a/2 on this sector, with Casimir 3/4, and commutes with every V_h. The
compact state f tensor chi_s therefore carries the conserved internal field
spin together with the transfer eigenvalue 209/210 and the energy
certificate E6 = 6.014331228348673 on the same state. The sealed audit makes
302,400 literal native calls and checks 604,800 columns per record law.

Scope: the embedding acts on the positive one-particle sector of the C
record law, and the spin is the internal CAR field spin of that sector.

## Row 8: annihilation

On a seam with loads (j, -j), j = +1 or -1, the load difference 2j is even,
so both coins give current j and the output (0, 0): the signed pair meets
and is removed by native contact. This holds on all 2,520 contact branches
(630 seams, two signs, two coins).

The contact state has a complete energy pointer. Let L_D be the compression
of L tensor I + I tensor L to distinct ordered positions, with the collision
diagonal as an absorbing boundary, and H_c = -1260 log(I - L_D/1260). For the
adjacent basis pair b = |0,1> the moments mu_n = <b, L_D^n b> are exact
integers (1, 11, 130, 1613, ...). Since 0 <= L_D <= 24 I, the log series
through n = 16 and its tail bound give an exact rational enclosure of width
below 10^-27 for the mean contact energy

    <b, H_c b> = 11.051928582384567  (transfer units).

The Bernstein pointer with effects binomial(16,k) X^k (I - X)^(16-k),
X = L_D/24, has 17 nonnegative exact outcomes that sum to one, so every
spectral component is accounted for.

The contact also releases a propagating record field. On the
reciprocal-work ring each site holds two record slots u_i and v_i. The two
swap layers of one cycle give u_i(t+1) = u_(i-1)(t) and v_i(t+1) = v_(i+1)(t),
so E = u + v and B = u - v obey the discrete wave equation

    F_i(t+1) + F_i(t-1) = F_(i-1)(t) + F_(i+1)(t),

and the field norm eps_i = (E_i^2 + B_i^2)/2 obeys local continuity with
bond flux J_(i+1/2) = u_i^2 - v_(i+1)^2. Contact writes the current j into
both slots, so E = 2j, B = 0 and eps = 2, equal to the released load square.
The two characteristics meet at the antipodal observer after L/2 cycles,
whose detector consumes both records and recreates the pair. The sealed run
launches and absorbs 360 of 360 contact flights and satisfies all 147,840
three-time wave equalities, with zero field, budget and native-current
defects.

Scope: the contact energy and the ring field norm are each measured on their
own register and in their own units.

## Row 20: electric force

One unit moves on the declared degree-six cubic routing (twelve ports per
cell, side-48 torus); the observer reconstructs its position from the signed
native current records alone and stops at squared radius 324. Each
trajectory records its visits to the radial shells. With
A_k(p) = mean over the cell centres of shell k of |cell|^-p, the measured
per-port residence fits

    a + b A_k(p),   shells k = 3 to 14,

by generalized least squares with the stratified covariance of the mean:

    p = 1.035144,   approximate 95% profile interval [0.983376, 1.086373],
    b = 0.328320.

The data are 5,040 trajectories with 13,213,995 nonzero currents and no
censoring; 3,168 native branches and 24 complete tapes check the records.
Exponent 1 is inside the interval, and 0.9 and 1.1 are outside it.

Scope: the 1/r form is measured on the declared routing over shells 3 to 14,
and the strength alpha is a separate question.

## Row 31: orbital levels

The inner probe is one native unit on a box of twelve-port cells. Holds and
internal moves keep the cell, and the six external moves are equally likely,
so the target-before-escape probability h_B is the simple cubic harmonic
function: on the 5x5x5 box h_B = 5/22, 1/11 and 1/22 at the face, edge and
corner cells around the target. The outer source accepts an internal seam in
cell x exactly when the probe escapes. Its marginal transfer has
L_eff = w L_external + sum_x (1 - h_B(x)) L_internal,x, and on an internal
band lambda of the cell (5 - sqrt5, 6, 5 + sqrt5 with multiplicities 3, 5, 3)
the spatial operator is

    S_lambda = -Laplacian + lambda (1 - h_B),

an attractive well of depth lambda h_B read from native records alone. On the
27-cell calibration the lowest squared energies are 2.2720093211,
4.2188978370 and 4.6123309041, each below its band threshold, and disabling
the feedback returns each to its threshold. The sealed run makes 305,157
native calls and a fresh replay reproduces it byte for byte.

With a = 1/m, delta = a^2 and probe box B = m^3, the operator
H_delta = delta^-2 L - lambda h_B/delta in the Bohr coordinate y = delta z
converges, shell by shell and counting multiplicity, to -Laplacian - beta/|y|
with beta = lambda c_* and c_* = 1/(4 pi G(0)) = 0.3148702, where G(0) is the
diagonal of the inverse cubic Laplacian. The radial Coulomb equation has
square-integrable solutions exactly when n = q + l + 1 with q >= 0, so

    e_n = -(lambda c_*)^2/(4 n^2),   l = 0, ..., n-1,   sum (2l+1) = n^2,

and shell n carries n^2 spatial orbital states per internal component.
Reflection-positive reconstruction of the same record field gives standing
modes of energy E_a(s) = 2 asinh(a sqrt(s)/2)/a, and the binding measured from
the same free threshold is E_a(lambda) - E_a(s_(n,k)) =
a^2 lambda^(3/2) c_*^2/(8 n^2) + o(a^2). The written controlled-limit proof is
archived in `evidence/row31_coulomb_limit/PROOF.md`.

Scope: routing, target and feedback are declared, and the shell spectrum is
a controlled-limit theorem with a finite native calibration of the well.

## Evidence and reproduction

```sh
python3 code/openwave_partial_records/verify.py
python3 code/openwave_partial_records/recompute.py
python3 tools/build_openwave_partial_records_claim.py --check
```

`verify.py` checks every archived file against its original sealed manifest,
binds the sealed receipt fields of each row and runs `recompute.py`, which
imports none of the sealed experiment scripts. The recompute rebuilds the
native arcs and checks the walk relations, the triangle reference and the
cell confinement in integer arithmetic, and the theta harmonic on all
outcomes for 129 words; it compares the published C-law CAR field code with
an independent Jordan-Wigner count on every endpoint column of the full graph
and every column of the twelve-port cell, including spin conservation and
completeness; it recomputes the contact branches, the Dirichlet moments, the
energy enclosure and all 17 pointer outcomes exactly, and checks the ring
characteristics, wave equation and flux continuity; it refits the residence
exponent and profile interval from the archived trajectory records; and it
recomputes the first-hit probabilities, the cell bands, the three wells and
the n^2 radial shell multiplicities. `verify.py --workspace DIR` also checks
every original sealed leaf in a copy of the OpenWave workspace.
