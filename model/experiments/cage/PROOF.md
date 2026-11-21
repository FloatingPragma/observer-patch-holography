# A cage in the flat branch of the scalar channel: statements and proofs

The scalar law of the [construction manual](../../construction_manual/CONSTRUCTION_MANUAL.md) (law 2) holds an
object on a ring of four seams and couples it to nothing else. This page states and proves that, describes the
three experiments of an observer on such a cage as laws on records, and gives the two operators of the observer
as words on registers. The candidate is `A-native-caged-exciton-v1`: the free law of the scalar channel with
two more operators in the observer's menu. Every statement is checked by the code of this folder, and the
numbers are in [results.json](results.json).

Every matrix, state and evolution parameter below belongs to the simulator and serves to derive what observers
record. An observer holds settings, counts and its own reference.

## 1. Setting

One cell has 12 ports, 30 seams and 120 wedges; a wedge $(a, b, c)$ is the two seams $ab$ and $bc$ at the port
$b$. Every seam and every wedge carries one mode, and a mode holds one object or none. $C$ is the signed
seam-by-wedge incidence: the wedge $(a, b, c)$ has $+1$ on the seam $ab$ and $-1$ on the seam $bc$. $Q$ is the
port-by-seam incidence without signs, and $d_a$ is the number of seams at the port $a$. The law of one object
and the law of many are

$$D = \begin{pmatrix} m & c\,C \\ c\,C^{\top} & -m \end{pmatrix},\qquad H_0 = \sum_{i,j} D_{ij}\, c_i^\dagger c_j,\qquad m = 1,\quad c = \tfrac{1}{\sqrt2} .$$

The hop matrix of the manual is $B = c\,C$ and its seam operator is $K = B B^{\top}$. The reference is the state
with every level of $D$ below zero filled.

## 2. The cage

**Lemma 1 (incidence identity).** On every graph of ports and seams,

$$C\,C^{\top} = \operatorname{diag}(d_a + d_b) - Q^{\top} Q .$$

*Proof.* The entry of $C C^{\top}$ at two seams sums, over the wedges that contain both, the product of their
signs. A seam $ab$ lies in $(d_a - 1) + (d_b - 1)$ wedges, each with the square of a sign. Two different seams
share at most one port, and two seams that share a port lie in exactly one wedge, with opposite signs: their
entry is $-1$. $Q^{\top} Q$ has 2 on the diagonal and 1 for two seams that share a port. ∎

$C C^{\top}$ is therefore the Laplacian of the graph whose points are the seams and whose links join seams that
share a port.

**Lemma 2 (the flat level).** If $Q e = 0$ and every port in the support of $e$ has $k$ seams, then
$C C^{\top} e = 2k\,e$.

*Proof.* Lemma 1 with $Q e = 0$: the diagonal term gives $2k\,e$ on the support. ∎

For a ring of four ports $a, b, c, d$, that is, four ports with the seams $ab$, $bc$, $cd$ and $da$, the profile is

$$e = \tfrac12\,\bigl(\,|ab\rangle - |bc\rangle + |cd\rangle - |da\rangle\,\bigr) .$$

At every port of the ring one seam with $+$ meets one with $-$, so $Q e = 0$, and on the cell
$C C^{\top} e = 10\,e$: the profile lies at the level $5$ of $K$.

**Lemma 3 (the profiles of one cell).** The profiles with $Q e = 0$ form a space of dimension 18. The cell has
30 rings, and their profiles have rank 18. No profile with $Q e = 0$ lives on fewer than four seams, and every
one on four seams is the profile of a ring.

*Proof.* $Q$ has rank 12 on a connected graph that contains a triangle, so its kernel has dimension $30 - 12$.
A profile in the kernel needs two seams at every port it touches. A set of at most three seams with two seams
at every port it touches is a face, and round a face the conditions $e_1 + e_2 = e_2 + e_3 = e_3 + e_1 = 0$
force zero. Four seams with two at every port touch four ports, each with exactly two of them: a ring, with
alternating signs. The counts 30 and 18 are finite checks in integers
([construction.py](construction.py), `flat_branch` and `small_supports`). ∎

**Proposition 4 (the pair).** Let $Q e = 0$ with $|e| = 1$ on ports with $k$ seams, and
$e' = C^{\top} e / \sqrt{2k}$. Then $|e'| = 1$ and

$$D\,(e,\ e') = (e,\ e')\,h,\qquad h = \begin{pmatrix} m & \sigma \\ \sigma & -m \end{pmatrix},\qquad \sigma = \sqrt{2k}\,c,\qquad E = \sqrt{m^2 + \sigma^2} .$$

*Proof.* $|C^{\top} e|^2 = e^{\top} C C^{\top} e = 2k$. $D\,(e, 0) = (m\,e,\ c\,C^{\top} e) = m\,(e, 0) + \sigma\,(0, e')$, and
$D\,(0, e') = (c\,C e',\ -m\,e') = \sigma\,(e, 0) - m\,(0, e')$ by Lemma 2. ∎

On the cell $k = 5$, so $\sigma = \sqrt5$ and $E = \sqrt6$. The pair sits on 4 seams and 28 wedges: at each of
the four ports the wedge of the two seams of the ring and the six wedges of one seam of the ring with one other
seam. The statement holds for every $m$ and $c$ and on every carrier on which the four ports have $k$ seams.
With

$$\cos\vartheta = \sqrt{\tfrac{E + m}{2E}},\quad \sin\vartheta = \sqrt{\tfrac{E - m}{2E}},\qquad \phi_+ = \cos\vartheta\,e + \sin\vartheta\,e',\quad \phi_- = -\sin\vartheta\,e + \cos\vartheta\,e'$$

the two modes of the pair are $D \phi_\pm = \pm E\,\phi_\pm$. In the reference $\phi_-$ holds an object and
$\phi_+$ is empty.

## 3. The excitation and its line

Let $a_+$ and $a_-$ remove an object from $\phi_+$ and $\phi_-$.

**Proposition 5.** With $X^\dagger = a_+^\dagger a_-$ and $P_X = n_+ (1 - n_-)$,

$$[H_0, X^\dagger] = 2E\,X^\dagger,\qquad [H_0, P_X] = 0 .$$

*Proof.* For a law of this form $[H_0, a^\dagger(\phi)] = a^\dagger(D \phi)$. So $[H_0, a_+^\dagger] = E\,a_+^\dagger$
and $[H_0, a_-] = E\,a_-$, which gives the first identity, and $n_+$ and $n_-$ commute with $H_0$. ∎

$X^\dagger$ moves the object of the lower mode into the upper one and leaves the number of objects as it is.
The state it makes lies $2E$ above the reference. The pair of modes has four sectors,

$$P_0 = (1 - n_+)(1 - n_-),\qquad P_{\rm ref} = (1 - n_+)\,n_-,\qquad P_X = n_+\,(1 - n_-),\qquad P_2 = n_+\,n_- ,$$

and each of them commutes with $H_0$.

**Proposition 6 (one line).** For every state and every operator built from $a_+$, $a_-$ and their adjoints
with as many creators as annihilators, the expectation at the parameter value $t$ is a constant plus one term
at the frequency $2E$.

*Proof.* Such an operator is a combination of $1$, $n_+$, $n_-$, $n_+ n_-$, $X^\dagger$ and $X$ and of products
of $X^\dagger$ or $X$ with the commuting numbers. Under the law the first four are constant, $X^\dagger$ carries
$e^{2iEt}$ and $X$ carries $e^{-2iEt}$. ∎

The number of the seam profile is such an operator:

$$n(e) = \cos^2\!\vartheta\; n_+ + \sin^2\!\vartheta\; n_- - \sin\vartheta\cos\vartheta\,\bigl(X + X^\dagger\bigr),\qquad \sin\vartheta\cos\vartheta = \frac{\sigma}{2E} .$$

**Corollary 7 (the swing from rest).** For the object of the pair placed on $e$,

$$P_e(t) = \langle n(e) \rangle_t = 1 - \frac{\sigma^2}{E^2}\,\sin^2(E\,t),\qquad P_e'(0) = 0,\qquad P_e''(0) = -2\sigma^2 .$$

That state has the weight $(1 + m/E)/2$ on the excitation, and its mean level above the lower mode is $E + m$.

*Proof.* $|\langle e|e^{-iht}|e\rangle|^2$ with $h^2 = E^2$ and $\langle e|h|e\rangle = m$. ∎

For the cell the weight is 0.7041, the mean level 3.4495, and the line sits at $2E = 4.898979$.

## 4. The two pulses

The observer writes with one operator: the number $n(e)$ as a potential of strength $v$ for a time $T$, with
the law running. On the pair this is $h + v\,n(e)$.

**Proposition 8.**

| Pulse | $v$ | $T$ | Effect on the lower mode |
| --- | --- | --- | --- |
| coherent | $-2E$ | $\pi / \bigl(2\sqrt{2E(E - m)}\bigr)$ | $\phi_-$ goes to $e$ |
| sharp | $-2E^2/m$ | $\pi m / (2\sigma E)$ | $\phi_-$ goes to $\phi_+$ |

each up to a phase.

*Proof.* Write $h + v\,n(e) = {\rm const} + (m + v/2)\,Z + \sigma\,X$ with the Pauli matrices of the pair
$(e, e')$. A pulse of duration $T$ turns the pair about the axis $(\sigma, 0, m + v/2)$ by the angle
$2T\sqrt{\sigma^2 + (m + v/2)^2}$. The lower mode has the direction $-(\sigma, 0, m)/E$, the profile $e$ has
$(0, 0, 1)$, and the upper mode has $+(\sigma, 0, m)/E$. For $v = -2E$ the axis is $(\sigma, 0, m - E)$, which
lies halfway between the lower mode and $e$ since $\sigma^2 = (E - m)(E + m)$; a half turn exchanges the two.
Its length is $\sqrt{2E(E - m)}$, which gives $T$. For $v = -2E^2/m$ the axis is $(\sigma, 0, -\sigma^2/m)$, at a
right angle to $(\sigma, 0, m)$; a half turn about it reverses the direction of the lower mode. Its length is
$\sigma E/m$. ∎

The coherent pulse does work: $E + m = 3.4495$ on the ideal reference, the mean level of Corollary 7. It also
sets the phase of the swing. The reference by itself has a constant seam record, and so has a sharp excitation.

## 5. Kept at every reading

**Theorem 9.** For every state of the cell or of a larger carrier, correlated between the cage and the rest or
not, and for every value $t$ of the simulator's parameter,

$$\operatorname{Tr}\bigl(\rho_t\,P_X\bigr) = \operatorname{Tr}\bigl(\rho_0\,P_X\bigr),\qquad \rho_t = e^{-itH_0}\,\rho_0\,e^{itH_0},$$

and the same for $P_0$, $P_{\rm ref}$ and $P_2$.

*Proof.* Proposition 5: the four projectors are polynomials in $n_+$ and $n_-$, which commute with $H_0$. ∎

A statement at every value of the parameter is a statement at every reading in any observer's own history. The
theorem uses no clock.

$[H_0, n] = 0$ holds for the number of every mode of every law of this form. What Theorem 9 adds is the place of
the two modes: 4 seams and 28 wedges at four ports, with no tail. That is Proposition 4, and it rests on the
signs of the hop and on five seams at each port of the ring. The same four seams with equal signs, three seams
of a face and a single seam lose their object to the rest of the cell
([test_cage.py](test_cage.py), `test_profiles_without_the_alternating_signs_are_not_kept`).

## 6. The experiments as laws on records

The prepared reference and every state made from it by pulses are states of independent orbitals. Such a state
is given by its occupation matrix $S_{ij} = \langle c_j^\dagger c_i \rangle$, an evolution $U$ of one object acts
as $U S U^\dagger$, and a mode $f$ is found full with the probability $S_{ff}$. After a reading of $f$ with the
probability $p = S_{ff}$ of "full", the other modes have the occupation matrices

$$S^{(0)} = A + \frac{b\,b^\dagger}{1 - p},\qquad S^{(1)} = A - \frac{b\,b^\dagger}{p},\qquad A = S_{\rm rest,rest},\quad b = S_{{\rm rest},f} ,$$

for "empty" and "full". The probability of a later reading is linear in these matrices, so the joint law of two
readings is exact in this description ([instruments.py](instruments.py), `joint_law`). The same law is
computed on all sixteen states of four modes by projection (`fock_joint_law`) and on all 150 modes of the cell
with the fuel; the three routes agree to rounding.

### 6.1 The swing of a cage

One observer holds the profile $e$, a fuel at the level zero and a reference $r$: three transits round one port,
with $D\,r = -m\,r$. Every such circulation is at a right angle to the pair, since $C\,r = 0$ and $e'$ lies in
the image of $C^{\top}$. One trial:

1. the coherent pulse on $e$;
2. a half pulse between $r$ and the fuel, of strength 64 and area $\pi/4$;
3. a wait of $q$ steps of the observer, of equal or of unequal length;
4. one of four phase settings $\theta$ on the fuel, one step, and a second half pulse;
5. the fuel is copied and recorded: the first reading $c$;
6. after one step the fuel, used as found, captures $e$ with a pulse of area $\pi/2$;
7. the fuel is copied and recorded: the second reading $z$.

All four pairs $(c, z)$ are kept, and the law runs through every step. Two readings come from the counts:

$$\Phi_q = \arg \sum_{\theta} e^{-i\theta}\,p(c = 1 \mid q, \theta) \ \text{ up to a constant},\qquad S_q = \tfrac14 \sum_{\theta} p(z = 1 \mid q, \theta) .$$

$\Phi_q$ is the phase of the observer's own reference at its step $q$. The observer takes as given that it
advances by between $1/16$ and $1/4$ per step, since counts cannot show a missed full turn, and the reader
refuses counts whose steps leave these bounds. $S_q$ is the seam record. The fuel's own fringe is a constant
plus one term in $\theta$, and the mean over four quarter turns removes that term. By Proposition 6,

$$S_q = A_0 + A_1 \cos\bigl(\omega\,\Phi_q + \varphi\bigr),\qquad \omega = \frac{2E}{m} ,$$

because the reference advances by $m$ per unit of the simulator's parameter and the seam record by $2E$. The
observer states the line in units of its own reference. No step of the reading uses the simulator's parameter.

### 6.2 The fuel dial

The fuel's level is set by a dial $\epsilon$, the fuel is coupled to $e$ with the strength $g_a = 0.1$ for 128 or
256 steps, and the run is made with the fuel prepared full and prepared empty. On the pair and the fuel the law
is

$$H_\epsilon = \begin{pmatrix} m & \sigma & g_a \\ \sigma & -m & 0 \\ g_a & 0 & \epsilon \end{pmatrix},\qquad p(\text{full} \mid \text{full}) - p(\text{full} \mid \text{empty}) = \bigl|\bigl(e^{-i t H_\epsilon}\bigr)_{ff}\bigr|^2 ,$$

whatever the state of the pair is: the two starts differ by $|f\rangle\langle f|$ in the occupation matrix, and
the probability is linear in it. The fuel empties where $\epsilon$ meets a level of the pair. The reader fits
$\sigma$, with $m$ and $g_a$ known in units of the reference and with the waits read from the observer's own
reference in separate trials (two pulses of area $\pi/4$ between $r$ and the fuel, four phase settings, both
preparations of the fuel). The level difference is $2\sqrt{m^2 + \sigma^2}$.

### 6.3 A cage kept

One trial: one step; one of three preparations (the coherent pulse, the sharp pulse, or no potential for the
time of the coherent pulse); one step; a wait of 0, 8, 32, 128 or 512 steps; one of four analyzers; one step;
the converter of strength 64 and area $\pi/2$ moves the content of $e$ into an empty fuel; the fuel is copied
into two pointers and recorded.

| Analyzer | Potential | Time | Effect on the pair |
| --- | --- | --- | --- |
| none | | 0 | |
| swap | $-2m\,n(e)$ | $\pi / (2\sigma)$ | $e$ and $e'$ change places |
| half swap | $-2m\,n(e)$ | $\pi / (4\sigma)$ | |
| free quarter | 0 | $\pi / (4E)$ | a quarter period of the law |

Let $G$ be the occupation matrix of the pair before the analyzer, in the basis $(e, e')$. Analyzer $j$ with the
step and the capture is a row vector $r_j$ on the pair, and

$$p_j = \Pr(\text{fuel full} \mid \text{analyzer } j) = r_j\,G\,r_j^\dagger .$$

$G$ has four real coordinates, the four effects are independent, and the matrix that takes
$(G_{ee}, G_{e'e'}, \operatorname{Re} G_{ee'}, \operatorname{Im} G_{ee'})$ to $(p_1, \dots, p_4)$ has the
condition number 4.623. Four counts give $G$.

**Proposition 10 (the sectors from $G$).** In a state of independent orbitals, in the basis $(\phi_-, \phi_+)$,

$$p_2 = \det G,\qquad p_0 = \det(1 - G),\qquad p_X = G_{++} - \det G,\qquad p_{\rm ref} = G_{--} - \det G .$$

*Proof.* For such a state $\langle n_+ n_- \rangle = G_{++} G_{--} - |G_{+-}|^2 = \det G$. Then
$p_X = \langle n_+ \rangle - \langle n_+ n_- \rangle$, $p_{\rm ref} = \langle n_- \rangle - \langle n_+ n_- \rangle$ and
$p_0 = 1 - \langle n_+ \rangle - \langle n_- \rangle + \langle n_+ n_- \rangle$. ∎

The reader computes the four effects from its settings and from the coupling $\sigma$ that the dial gave. An
estimate is linear in the counts and may leave the range of a probability by its sampling error; it is reported
as it comes out.

### 6.4 The changed law

The changed law adds seam hopping round every port,

$$J = \sum_{b}\ \sum_{x, y \text{ consecutive round } b} \bigl(c_{bx}^\dagger c_{by} + c_{by}^\dagger c_{bx}\bigr),$$

with the coefficient one, in every step of a run: preparation, wait, analyzer and capture. It keeps the
alphabet, the routing, the number of objects and the symmetries of the cell, and the observer's reference keeps
its level under it. It does not keep the pair: the commutator of $J$ with the projector on the pair has the norm
$\sqrt2$. The start of a run under the changed law is the reference of the law.

## 7. The two operators on registers

The register form of the scalar channel ([PROOF.md](../../candidate_laws/scalar_channel/PROOF.md)) writes the number of a mode as a product
$B_i$ of $Z$ and the passage between two linked modes as a word $A_{ij}$, with
$c_i^\dagger c_j + c_j^\dagger c_i = \tfrac{i}{2} A_{ij} (B_i - B_j)$. For two modes that no link joins, $A_{ij}$ is
the product of the link words along a path with the phase $i^{\,(\text{links}) - 1}$, the value that a closed
loop of links with the sign $+1$ fixes. The four seams of the ring are joined through one tree: the two seams at
the port 0 through the fuel of that port, the two seams at the port 4 through the fuel of that port, and the two
fuels through the diagonal seam 0-4. Every link of the tree is a register of the cell. With the signs
$s = (+, -, +, -)$ of the profile,

$$Q_e = 4 - \sum_i B_i + \sum_{i<j} s_i s_j \bigl(i A_{ij} B_i - i A_{ij} B_j\bigr),\qquad n(e) = Q_e / 8 ,$$

$$K_e = \sum_i s_i \bigl(i A_{if} B_i - i A_{if} B_f\bigr),\qquad c(e)^\dagger f + f^\dagger c(e) = K_e / 4 ,$$

with $f$ the fuel of the port 0. The following are identities between polynomials of words with integer
coefficients on all 432 registers ([registers.py](registers.py)):

- $Q_e$ has 17 words of weight at most 15, and $Q_e^2 = 8\,Q_e$;
- $K_e$ has 8 words of weight at most 24, $K_e^2 = 2 Q_e + 8 R - 2 Q_e R$ with $R = 1 - B_f = 2 n_f$, and
  $[K_e, Q_e + 4R] = 0$; in the operators of the modes these read $G_e^2 = n(e) + n_f - 2\,n(e)\,n_f$ and
  $[G_e, n(e) + n_f] = 0$ for $G_e = K_e/4$;
- every word of $Q_e$ and of $K_e$ commutes with all 260 loops (4,420 and 2,080 pairs), $Q_e$ commutes with the
  number of objects and with every fuel and reference number, and $K_e$ with every fuel and reference number
  other than $f$;
- the words touch the registers of the ports 0, 2, 4 and 6, at most four of them in one word of $Q_e$ and three
  in one word of $K_e$.

One guard, the product of the code checks of those four ports, stands in front of each whole operator, the
identity term of $Q_e$ included. A guard on single words would not keep the identity $Q_e^2 = 8 Q_e$.

The signs of the profile are in the words. On seven modes joined as the tree (four seams, two fuels and the
diagonal seam), the words equal the operators of the modes under the map between registers and objects of the
register form, as dense matrices ([test_cage.py](test_cage.py),
`test_the_words_of_the_two_operators_against_dense_fermions`).

## 8. On a federation

Proposition 4 needs five seams at each port of the ring. On the federations of 20, 80 and 320 cells the rings
with that property are checked in integers, each with $Q e = 0$ and $C C^{\top} e = 10\,e$
([results.json](results.json), `construction.federations`). The ports with five seams and the seams between
them fall into one piece per cell, and the profiles with $Q e = 0$ on a piece number its seams minus the rank of
its incidence: 6 for a cell whose three ports with an external seam share no seam, 7 for a cell in which two of
them do. A ring through a port with an external seam has $C C^{\top} e \ne 10\,e$, and an object placed on its
profile leaves the pair. The instruments on registers exist for one cell; on a federation the statements are
those of the law of one object.

## 9. What is supplied, and where the statements end

Supplied: the law with its two coefficients $m$ and $c$; the prepared reference with its ramps; the two operators
on a profile with their guard; the pulses with their strengths and times; the fuel dial; the observer's
reference with the bounds on its step; the changed law as a control.

- The pulse supplies the phase and the work of the swing. The reference by itself does not swing, and a swing
  that starts with no supplied phase is work in progress.
- The line at $2E$ and the level difference from the dial are two routes to one number of the law: a record
  that mixes two levels swings at their difference. Their agreement tests the instruments and the description by
  two modes.
- The level of a cage is flat. A cage has no group speed, and what is written on it stays on its ring. The
  motion of a caged excitation is work in progress.
- The law has no term between two objects. The statements are statements of free evolution.
- The excitation leaves the number of objects as it is. It is a state of the singlet code of the register form,
  on which every one of the twelve currents acts as zero.
