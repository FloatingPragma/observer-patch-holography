# Positive transfer-energy windows from complete native records

The fair source fixes a positive transfer and a counted clock. These determine
finite positive spectral windows without fitting poles or selecting particle
names. Complete native record experiments reconstruct all their probabilities.
This does not implement that positive operator-valued measurement as one native
intervention on the working carrier.

## The source and its clock

Use the binary fair source with V ports, M seams, maximum degree Delta,
independent uniform seam choices, fresh fair coins and repair probability one.
Its fixed C/E mixture has population transfer

    P_q=I-G_q/(2M),       G_q=sum_e(I-S_e).

The [exchange construction](../native_exchange_clock/PROOF.md) proves positivity
and fixes counted time t=N/(2M). On primitive degrees at most r, it proves
0<=G<=2r Delta. Assume M>=2r Delta and set

    b=floor(M/(2r Delta)),     Q=P^b,     H=-2M log(P).

Bernoulli's inequality gives Q>=I/2 on these reducing blocks. Thus
H=-(2M/b)log(Q) is a unique nonnegative generator in this clock. The constant
transfer vector has energy zero. The local working density I/binomial(V,q)
is a different reference from this pure transfer ground.

For a degree-two contrast u=h/2 made from two disjoint native seams, the
[complete four-E identity](../native_correlation_sectors/PROOF.md) gives

    m_j=<u,Q^j u>=C(bj)/v_q,
    v_q=4 binomial(V-4,q-2)/binomial(V,q),     m_0=1.

It holds for 2<=q<=V-2. The outgoing contrasts are used in the first pair
and incoming contrasts in the second. All raw flags, zero outcomes and
measurement backaction are included. The stationary complete scalar score
has range [-1,1]. A generator or spectral calculation is not an expanded
execution of every native history.

## A complete positive spectral detector

For any integer n>=1 define n+1 effects on the transfer Hilbert space:

    B_k(Q)=binomial(n,k) Q^k(I-Q)^(n-k),    0<=k<=n.

Since 0<=Q<=I, every effect is positive and the binomial theorem gives
sum_k B_k=I. No eigenvalue, dark component or overflow bin is removed.
For a normalized transfer probe u, the exact probabilities are

    p_k=<u,B_k(Q)u>
       =sum_(j=k)^n (-1)^(j-k) binomial(n,k) binomial(n-k,j-k) m_j.   (1)

The complete native lag statistics therefore determine a positive normalized
histogram. The effect family is defined independently of the probe, while a
probe sees only its own spectral measure. A spectral component orthogonal to
that probe has zero weight; this does not eliminate it from the theory.

Equation (1) is a signed combination of different experiments. It supplies
neither an instrument realizing B_k on the working carrier nor its conditional
post-measurement state. Positive reconstructed probabilities alone do not prove
the operational premise of a quantum measurement theorem.

## Finite energy resolution

At a transfer eigenvalue x in [1/2,1], the detector output K has the binomial
law Bin(n,x). Assign k>0 the energy label

    E_k=-(2M/b)log(k/n),

and retain k=0 as an explicit overflow label. For 0<delta<1/2, Chebyshev's
inequality and x(1-x)<=1/4 give

    Pr(|K/n-x|>=delta)<=min(1,1/(4n delta^2)).                         (2)

On the complementary event k/n>=1/2-delta>0. The mean-value theorem gives

    |E_k-H(x)|<=(2M/b) delta/(1/2-delta)
              <=8r Delta delta/(1/2-delta).                          (3)

The same bounds hold after integrating any positive probe spectral measure.
The overflow probability is at most 2^-n. For a threshold a, the estimated
transfer cumulative law obeys

    mu([1/2,a-delta])-epsilon
      <=Pr(K/n<=a)<=mu([1/2,a+delta])+epsilon,
    epsilon=min(1,1/(4n delta^2)),                                    (4)

where intervals are intersected with [1/2,1]. This bounds spectral windows;
it does not isolate an atom without a separate spectral separation bound.

Energy labels have exact rational enclosures. For 0<y<=1,

    sum_(j=1)^s (1-y)^j/j <=-log(y)
      <=sum_(j=1)^s (1-y)^j/j +(1-y)^(s+1)/((s+1)y).                 (5)

The executable certificate uses these intervals, without rounded logarithms
as theorem inputs. On the same reducing degree blocks,

    ||H-G||<=r^2 Delta^2/(M-r Delta).

This adds to (3) when the label is compared with the local exchange energy.
It does not convert the exchange response into a Lorentzian dispersion law.

## Complete-vector errors and native preparation

The zeroth moment is fixed algebraically to one and needs no experiment.
For the coefficient table c_kj in (1), define A_k=sum_(j=1)^n |c_kj|.
The binomial identities give

    A=sum_k A_k=3^n-1.

If every normalized nonzero-lag score has error at most eta, the histogram
error is at most A eta/2 in half l1 distance. Approximate moments can give
negative bin estimates. The implementation retains them and this error bound;
it does not clip, normalize or label them exact probabilities. Their sum is
exactly one because sum_k c_kj is one for j=0 and zero otherwise.

For n lag settings and L independent repetitions per setting, each raw score
has variance at most one. Chebyshev's inequality and a union bound give half
l1 sampling error at most s, with failure probability at most p, whenever

    L>=n A^2/(4p v_q^2 s^2).                                         (6)

Fresh independent repetitions are a premise. Increasing n improves spectral
resolution while this sufficient statistical bound grows exponentially; it
is not a claim of practical spectral extraction. A separately proved stop
probability p_stop adds A p_stop/(2v_q) when stopped scores are set to zero.
Stops and all raw flags remain in the experiment, without conditioning on
success.

The finite fair twelve-port total-two source has v_q=2/33. The
[local preparation](../native_fair_local_preparation/PROOF.md) uses eleven
tree commands and b_p blocks of sixty E commands. With b_p=64 it has 3851
attempts and prepares I/66 with half trace error at most 2^-64 for arbitrary
input and untouched reference, after marginalizing its preparation records
for this statistic. On this global fixed sector the conditional reference
is I/66 tensor the untouched reference marginal. No preparation value is
fed back into the following score. The actual archive retains every value.
Thus the normalized score error is at most 2^-63/v_q.

For n=16, the source has M=30, Delta=5, r=2 and b=1. Each of the sixteen
setting programs includes the same preparation, four complete E probes and
its prescribed lag. The longest program has 3871 attempts and 3856 syntax
nodes, within the declared admission limits. Its lazy instrument is defined;
the expanded host's branch ceiling is retained. The exact ideal histogram
is evaluated from the native transfer. The certificate does not claim that
the displayed independent repetitions or all expanded histories were run.

## Capacity and physical identification

At fixed r and n, b grows proportionally to M, and the clock coefficient
2M/b stays at most 8r Delta. Consequently (2)-(3) have bounds independent
of capacity. The largest native lag bn grows with M; a fixed host attempt
ceiling cannot admit these experiments at all capacities. The executable
polynomial evaluator bounds n by 64. The theorem holds for every finite n,
without changing any native resource limit.

The twelve-port preparation proof does not prepare a whole growing canonical
source. The growing budgets retain their global-reference preparation
premise. The general all-sector source does not make the degree-two probe
nonzero at q=0,1,V-1,V. Those sectors remain in the source, outside this
particular normalized contrast domain.

Reading the simulator fixes P and the attempt count. Converting the counted
generator into joules requires one common physical time and Planck constant
identification. Under t_phys=a*t_count, E_phys=(hbar/a)H gives the same
dimensionless exponent. Knowledge of the source does not select a numerical
laboratory conversion. Dimensionless spectral ratios are invariant under
this common factor, but assigning those energies to rest masses requires
particle sectors and a common relativistic dispersion/intervention map.
The construction supplies no stress tensor, conserved physical energy,
interacting chiral continuum or measured particle mass.
