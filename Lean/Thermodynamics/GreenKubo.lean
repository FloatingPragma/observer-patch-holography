import FluctuationTheorems

namespace OPH.Thermodynamics

/-!
# Finite Green--Kubo transport from reversible repair dynamics

This module isolates a finite theorem package for equilibrium transport.
For a finite reversible stochastic kernel `K`, the positive generator is
`L = I - K`.  A Poisson solver `R` on centered currents obeys `L (R j) = j`.
The pairing `⟨j, R k⟩_pi` is then symmetric and its diagonal is nonnegative.
When `R` is linear and solves the Poisson equation on the whole centered
sector, every finite matrix of these coefficients has a nonnegative
quadratic form.  The positivity is derived from the exact Dirichlet-form
identity for `L`; it is not assumed for `R`.

At every finite cutoff, the same coefficient equals the discrete
Green--Kubo sum plus an exact propagated-resolvent remainder.  It reduces to
the finite sum only under literal finite-step extinction of that remainder.
Finite-step extinction is not used as a synonym for asymptotic decay or
infinite-time convergence.  No continuum limit, hydrodynamic scaling,
source attachment, physical clock, conductivity normalization, or measured
transport number is proved here.

The final section specializes the construction to the typed
conditional-resampling kernel.  A full-fibre heat bath is a projector.  A
current whose conditional weighted sum vanishes on every repaired-visible
fibre is killed after one step, so every positive-lag correlation vanishes.
This is a transport no-go for using one idempotent full-fibre resampling as a
long-memory dynamics; nontrivial transport requires a nonprojective local or
random-scan evolution, or another source-derived time generator.
-/

variable {Omega : Type*} [Fintype Omega] [DecidableEq Omega]
variable {B : Type*} [DecidableEq B]

/-! ## Centered equilibrium currents and finite dynamics -/

/-- Equilibrium mean with respect to a finite reference weight. -/
noncomputable def equilibriumMean (pi : Omega -> Real) (f : Omega -> Real) : Real :=
  Finset.sum Finset.univ fun x => pi x * f x

/-- A current with exactly zero equilibrium mean.  Using a subtype keeps the
Poisson and Green--Kubo interfaces restricted to the mean-zero sector where
the stationary constant mode has been removed. -/
def CenteredCurrent (pi : Omega -> Real) :=
  {f : Omega -> Real // equilibriumMean pi f = 0}

/-- Center an observable under a normalized reference weight. -/
noncomputable def centerCurrent (pi : Omega -> Real)
    (hpi : Finset.sum Finset.univ pi = 1) (f : Omega -> Real) :
    CenteredCurrent pi := by
  refine ⟨fun x => f x - equilibriumMean pi f, ?_⟩
  unfold equilibriumMean
  calc
    (Finset.sum Finset.univ fun x => pi x *
        (f x - Finset.sum Finset.univ fun y => pi y * f y))
        = (Finset.sum Finset.univ fun x => pi x * f x) -
          (Finset.sum Finset.univ fun x =>
            pi x * (Finset.sum Finset.univ fun y => pi y * f y)) := by
            rw [← Finset.sum_sub_distrib]
            apply Finset.sum_congr rfl
            intro x _
            ring
    _ = (Finset.sum Finset.univ fun x => pi x * f x) -
          (Finset.sum Finset.univ pi) *
            (Finset.sum Finset.univ fun y => pi y * f y) := by
            rw [Finset.sum_mul]
    _ = 0 := by rw [hpi]; ring

/-- Weighted real equilibrium pairing. -/
noncomputable def equilibriumPair (pi : Omega -> Real)
    (f g : Omega -> Real) : Real :=
  Finset.sum Finset.univ fun x => pi x * f x * g x

omit [DecidableEq Omega] in
theorem equilibriumPair_comm (pi : Omega -> Real) (f g : Omega -> Real) :
    equilibriumPair pi f g = equilibriumPair pi g f := by
  unfold equilibriumPair
  apply Finset.sum_congr rfl
  intro x _
  ring

omit [DecidableEq Omega] in
theorem equilibriumPair_add_right (pi : Omega -> Real)
    (f g h : Omega -> Real) :
    equilibriumPair pi f (fun x => g x + h x) =
      equilibriumPair pi f g + equilibriumPair pi f h := by
  unfold equilibriumPair
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro x _
  ring

omit [DecidableEq Omega] in
theorem equilibriumPair_sub_right (pi : Omega -> Real)
    (f g h : Omega -> Real) :
    equilibriumPair pi f (fun x => g x - h x) =
      equilibriumPair pi f g - equilibriumPair pi f h := by
  unfold equilibriumPair
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro x _
  ring

omit [DecidableEq Omega] in
theorem equilibriumPair_sub_left (pi : Omega -> Real)
    (f g h : Omega -> Real) :
    equilibriumPair pi (fun x => f x - g x) h =
      equilibriumPair pi f h - equilibriumPair pi g h := by
  unfold equilibriumPair
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro x _
  ring

omit [DecidableEq Omega] in
theorem equilibriumPair_smul_left (pi : Omega -> Real) (a : Real)
    (f g : Omega -> Real) :
    equilibriumPair pi (a • f) g = a * equilibriumPair pi f g := by
  unfold equilibriumPair
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro x _
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

omit [DecidableEq Omega] in
theorem equilibriumPair_smul_right (pi : Omega -> Real) (a : Real)
    (f g : Omega -> Real) :
    equilibriumPair pi f (a • g) = a * equilibriumPair pi f g := by
  unfold equilibriumPair
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro x _
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

omit [DecidableEq Omega] in
theorem equilibriumPair_sum_left {I : Type*} [Fintype I]
    (pi : Omega -> Real) (f : I -> Omega -> Real) (g : Omega -> Real) :
    equilibriumPair pi (Finset.sum Finset.univ f) g =
      Finset.sum Finset.univ fun i => equilibriumPair pi (f i) g := by
  unfold equilibriumPair
  simp only [Finset.sum_apply]
  calc
    (Finset.sum Finset.univ fun x =>
      pi x * (Finset.sum Finset.univ fun i => f i x) * g x) =
        Finset.sum Finset.univ (fun x =>
          Finset.sum Finset.univ fun i => pi x * f i x * g x) := by
            apply Finset.sum_congr rfl
            intro x _
            rw [Finset.mul_sum, Finset.sum_mul]
    _ = Finset.sum Finset.univ (fun i =>
        Finset.sum Finset.univ fun x => pi x * f i x * g x) :=
          Finset.sum_comm

omit [DecidableEq Omega] in
theorem equilibriumPair_sum_right {I : Type*} [Fintype I]
    (pi : Omega -> Real) (f : Omega -> Real) (g : I -> Omega -> Real) :
    equilibriumPair pi f (Finset.sum Finset.univ g) =
      Finset.sum Finset.univ fun i => equilibriumPair pi f (g i) := by
  unfold equilibriumPair
  simp only [Finset.sum_apply]
  calc
    (Finset.sum Finset.univ fun x =>
      pi x * f x * (Finset.sum Finset.univ fun i => g i x)) =
        Finset.sum Finset.univ (fun x =>
          Finset.sum Finset.univ fun i => pi x * f x * g i x) := by
            apply Finset.sum_congr rfl
            intro x _
            rw [Finset.mul_sum]
    _ = Finset.sum Finset.univ (fun i =>
        Finset.sum Finset.univ fun x => pi x * f x * g i x) :=
          Finset.sum_comm

/-- Backward action of a finite transition kernel on observables. -/
noncomputable def kernelAct (K : Omega -> Omega -> Real)
    (f : Omega -> Real) (x : Omega) : Real :=
  Finset.sum Finset.univ fun y => K x y * f y

/-- Positive discrete generator `L = I - K`. -/
noncomputable def dissipationOperator (K : Omega -> Omega -> Real)
    (f : Omega -> Real) (x : Omega) : Real :=
  f x - kernelAct K f x

/-- The `n`-step propagated observable. -/
noncomputable def kernelIterate (K : Omega -> Omega -> Real) :
    Nat -> (Omega -> Real) -> (Omega -> Real)
  | 0, f => f
  | n + 1, f => kernelAct K (kernelIterate K n f)

/-- Equilibrium lag-`n` correlation. -/
noncomputable def lagCorrelation (pi : Omega -> Real)
    (K : Omega -> Omega -> Real) (f g : Omega -> Real) (n : Nat) : Real :=
  equilibriumPair pi f (kernelIterate K n g)

/-- Finite geometric resolvent `I + K + ... + K^N`. -/
noncomputable def finiteResolvent (K : Omega -> Omega -> Real) :
    Nat -> (Omega -> Real) -> (Omega -> Real)
  | 0, f => f
  | n + 1, f => fun x =>
      finiteResolvent K n f x + kernelIterate K (n + 1) f x

/-- Discrete Green--Kubo correlation sum from lag zero through lag `N`. -/
noncomputable def integratedCorrelation (pi : Omega -> Real)
    (K : Omega -> Omega -> Real) (f g : Omega -> Real) : Nat -> Real
  | 0 => lagCorrelation pi K f g 0
  | n + 1 => integratedCorrelation pi K f g n +
      lagCorrelation pi K f g (n + 1)

omit [DecidableEq Omega] in
theorem equilibriumPair_kernelAct (pi : Omega -> Real)
    (K : Omega -> Omega -> Real) (f g : Omega -> Real) :
    equilibriumPair pi f (kernelAct K g) =
      Finset.sum Finset.univ fun x =>
        Finset.sum Finset.univ fun y => pi x * K x y * f x * g y := by
  unfold equilibriumPair kernelAct
  apply Finset.sum_congr rfl
  intro x _
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro y _
  ring

omit [DecidableEq Omega] in
/-- Detailed balance makes the transition operator self-adjoint in the
equilibrium pairing. -/
theorem kernelAct_pair_symm (pi : Omega -> Real)
    (K : Omega -> Omega -> Real)
    (hdb : forall x y, pi x * K x y = pi y * K y x)
    (f g : Omega -> Real) :
    equilibriumPair pi f (kernelAct K g) =
      equilibriumPair pi g (kernelAct K f) := by
  rw [equilibriumPair_kernelAct, equilibriumPair_kernelAct]
  exact correlation_symm pi K hdb f g

omit [DecidableEq Omega] in
/-- The positive generator `I-K` is self-adjoint under detailed balance. -/
theorem dissipation_selfAdjoint (pi : Omega -> Real)
    (K : Omega -> Omega -> Real)
    (hdb : forall x y, pi x * K x y = pi y * K y x)
    (f g : Omega -> Real) :
    equilibriumPair pi f (dissipationOperator K g) =
      equilibriumPair pi (dissipationOperator K f) g := by
  unfold dissipationOperator
  rw [equilibriumPair_sub_right, equilibriumPair_sub_left]
  rw [kernelAct_pair_symm pi K hdb f g]
  rw [equilibriumPair_comm pi g (kernelAct K f)]

/-! ## Exact Dirichlet identity and Green--Kubo reciprocity -/

/-- Dirichlet form of a finite transition kernel. -/
noncomputable def dirichletForm (pi : Omega -> Real)
    (K : Omega -> Omega -> Real) (f : Omega -> Real) : Real :=
  (1 / 2 : Real) * Finset.sum Finset.univ (fun x =>
    Finset.sum Finset.univ (fun y =>
      pi x * K x y * (f x - f y) ^ 2))

omit [DecidableEq Omega] in
/-- Exact reversible Markov identity
`<f,(I-K)f> = (1/2) sum pi_x K_xy (f_x-f_y)^2`. -/
theorem dissipation_eq_dirichlet (pi : Omega -> Real)
    (K : Omega -> Omega -> Real)
    (hK1 : forall x, Finset.sum Finset.univ (K x) = 1)
    (hdb : forall x y, pi x * K x y = pi y * K y x)
    (f : Omega -> Real) :
    equilibriumPair pi f (dissipationOperator K f) =
      dirichletForm pi K f := by
  let A : Real := equilibriumPair pi f f
  let C : Real := equilibriumPair pi f (kernelAct K f)
  let Sxx : Real := Finset.sum Finset.univ fun x =>
    Finset.sum Finset.univ fun y => pi x * K x y * (f x) ^ 2
  let Syy : Real := Finset.sum Finset.univ fun x =>
    Finset.sum Finset.univ fun y => pi x * K x y * (f y) ^ 2
  let Sxy : Real := Finset.sum Finset.univ fun x =>
    Finset.sum Finset.univ fun y => pi x * K x y * f x * f y
  have hxx : Sxx = A := by
    unfold Sxx A equilibriumPair
    apply Finset.sum_congr rfl
    intro x _
    calc
      (Finset.sum Finset.univ fun y => pi x * K x y * (f x) ^ 2) =
          (pi x * (f x) ^ 2) * Finset.sum Finset.univ (K x) := by
            rw [Finset.mul_sum]
            apply Finset.sum_congr rfl
            intro y _
            ring
      _ = pi x * f x * f x := by rw [hK1 x]; ring
  have hyy : Syy = A := by
    unfold Syy A equilibriumPair
    calc
      (Finset.sum Finset.univ fun x =>
          Finset.sum Finset.univ fun y => pi x * K x y * (f y) ^ 2) =
          Finset.sum Finset.univ (fun x =>
            Finset.sum Finset.univ fun y => pi y * K y x * (f y) ^ 2) := by
              apply Finset.sum_congr rfl
              intro x _
              apply Finset.sum_congr rfl
              intro y _
              rw [hdb x y]
      _ = Finset.sum Finset.univ (fun y =>
          Finset.sum Finset.univ fun x => pi y * K y x * (f y) ^ 2) :=
            Finset.sum_comm
      _ = Finset.sum Finset.univ (fun y => pi y * f y * f y) := by
            apply Finset.sum_congr rfl
            intro y _
            calc
              (Finset.sum Finset.univ fun x => pi y * K y x * (f y) ^ 2) =
                  (pi y * (f y) ^ 2) * Finset.sum Finset.univ (K y) := by
                    rw [Finset.mul_sum]
                    apply Finset.sum_congr rfl
                    intro x _
                    ring
              _ = pi y * f y * f y := by rw [hK1 y]; ring
  have hxy : Sxy = C := by
    unfold Sxy C
    exact (equilibriumPair_kernelAct pi K f f).symm
  have hsquare :
      (Finset.sum Finset.univ fun x =>
        Finset.sum Finset.univ fun y =>
          pi x * K x y * (f x - f y) ^ 2) =
        Sxx + Syy - 2 * Sxy := by
    unfold Sxx Syy Sxy
    calc
      (Finset.sum Finset.univ fun x =>
          Finset.sum Finset.univ fun y =>
            pi x * K x y * (f x - f y) ^ 2) =
          Finset.sum Finset.univ (fun x =>
            Finset.sum Finset.univ fun y =>
              (pi x * K x y * (f x) ^ 2 +
                pi x * K x y * (f y) ^ 2) -
                2 * (pi x * K x y * f x * f y)) := by
                  apply Finset.sum_congr rfl
                  intro x _
                  apply Finset.sum_congr rfl
                  intro y _
                  ring
      _ = (Finset.sum Finset.univ fun x =>
            Finset.sum Finset.univ fun y => pi x * K x y * (f x) ^ 2) +
          (Finset.sum Finset.univ fun x =>
            Finset.sum Finset.univ fun y => pi x * K x y * (f y) ^ 2) -
          2 * (Finset.sum Finset.univ fun x =>
            Finset.sum Finset.univ fun y => pi x * K x y * f x * f y) := by
              simp only [Finset.sum_sub_distrib, Finset.sum_add_distrib,
                Finset.mul_sum]
  unfold dissipationOperator
  rw [equilibriumPair_sub_right]
  change A - C = dirichletForm pi K f
  unfold dirichletForm
  rw [hsquare, hxx, hyy, hxy]
  ring

omit [DecidableEq Omega] in
/-- Positivity of the generator pairing follows from nonnegative reference
weights and transition probabilities. -/
theorem dissipation_nonneg (pi : Omega -> Real)
    (K : Omega -> Omega -> Real)
    (hpi0 : forall x, 0 <= pi x)
    (hK0 : forall x y, 0 <= K x y)
    (hK1 : forall x, Finset.sum Finset.univ (K x) = 1)
    (hdb : forall x y, pi x * K x y = pi y * K y x)
    (f : Omega -> Real) :
    0 <= equilibriumPair pi f (dissipationOperator K f) := by
  rw [dissipation_eq_dirichlet pi K hK1 hdb f]
  unfold dirichletForm
  apply mul_nonneg
  · norm_num
  · apply Finset.sum_nonneg
    intro x _
    apply Finset.sum_nonneg
    intro y _
    exact mul_nonneg (mul_nonneg (hpi0 x) (hK0 x y)) (sq_nonneg _)

/-- Green--Kubo pairing associated with a proposed Poisson solver on the
centered-current sector. -/
noncomputable def greenKuboPair (pi : Omega -> Real)
    (R : (Omega -> Real) -> (Omega -> Real))
    (j k : CenteredCurrent pi) : Real :=
  equilibriumPair pi j.1 (R k.1)

omit [DecidableEq Omega] in
/-- A centered current pairs to zero with every constant observable.  This
removes the additive-constant ambiguity of finite Poisson solutions. -/
theorem centeredCurrent_pair_const_zero (pi : Omega -> Real)
    (j : CenteredCurrent pi) (c : Real) :
    equilibriumPair pi j.1 (fun _ => c) = 0 := by
  unfold equilibriumPair equilibriumMean at *
  calc
    (Finset.sum Finset.univ fun x => pi x * j.1 x * c) =
        (Finset.sum Finset.univ fun x => pi x * j.1 x) * c := by
          rw [Finset.sum_mul]
    _ = 0 := by
      have hjmean : (Finset.sum Finset.univ fun x => pi x * j.1 x) = 0 := by
        simpa [equilibriumMean] using j.2
      rw [hjmean]
      simp

omit [DecidableEq Omega] in
/-- Green--Kubo coefficients are invariant under a current-dependent
constant shift of a Poisson solver. -/
theorem greenKuboPair_add_solver_constant (pi : Omega -> Real)
    (R : (Omega -> Real) -> (Omega -> Real))
    (c : (Omega -> Real) -> Real) (j k : CenteredCurrent pi) :
    greenKuboPair pi (fun f x => R f x + c f) j k =
      greenKuboPair pi R j k := by
  unfold greenKuboPair
  rw [equilibriumPair_add_right,
    centeredCurrent_pair_const_zero pi j (c k.1), add_zero]

omit [DecidableEq Omega] in
/-- Exact Onsager symmetry of the resolvent pairing.  The only inverse
premises are the two displayed Poisson equations for the currents used. -/
theorem greenKuboPair_symm_of_poisson (pi : Omega -> Real)
    (K : Omega -> Omega -> Real)
    (hdb : forall x y, pi x * K x y = pi y * K y x)
    (R : (Omega -> Real) -> (Omega -> Real))
    (j k : CenteredCurrent pi)
    (hj : dissipationOperator K (R j.1) = j.1)
    (hk : dissipationOperator K (R k.1) = k.1) :
    greenKuboPair pi R j k = greenKuboPair pi R k j := by
  unfold greenKuboPair
  calc
    equilibriumPair pi j.1 (R k.1) =
        equilibriumPair pi (dissipationOperator K (R j.1)) (R k.1) := by
          rw [hj]
    _ = equilibriumPair pi (R j.1) (dissipationOperator K (R k.1)) := by
          symm
          exact dissipation_selfAdjoint pi K hdb (R j.1) (R k.1)
    _ = equilibriumPair pi (R j.1) k.1 := by rw [hk]
    _ = equilibriumPair pi k.1 (R j.1) := equilibriumPair_comm pi _ _

omit [DecidableEq Omega] in
/-- Positive semidefiniteness of the Green--Kubo coefficient.  It is a
consequence of the reversible Markov Dirichlet form and the Poisson equation,
not an independent positivity assumption on the solver. -/
theorem greenKuboPair_nonneg_of_poisson (pi : Omega -> Real)
    (K : Omega -> Omega -> Real)
    (hpi0 : forall x, 0 <= pi x)
    (hK0 : forall x y, 0 <= K x y)
    (hK1 : forall x, Finset.sum Finset.univ (K x) = 1)
    (hdb : forall x y, pi x * K x y = pi y * K y x)
    (R : (Omega -> Real) -> (Omega -> Real))
    (j : CenteredCurrent pi)
    (hj : dissipationOperator K (R j.1) = j.1) :
    0 <= greenKuboPair pi R j j := by
  unfold greenKuboPair
  calc
    0 <= equilibriumPair pi (R j.1)
        (dissipationOperator K (R j.1)) :=
          dissipation_nonneg pi K hpi0 hK0 hK1 hdb (R j.1)
    _ = equilibriumPair pi (R j.1) j.1 := by rw [hj]
    _ = equilibriumPair pi j.1 (R j.1) := equilibriumPair_comm pi _ _

/-! ## Exact finite correlation representation -/

omit [DecidableEq Omega] in
theorem kernelAct_add (K : Omega -> Omega -> Real)
    (f g : Omega -> Real) :
    kernelAct K (fun x => f x + g x) =
      fun x => kernelAct K f x + kernelAct K g x := by
  funext x
  unfold kernelAct
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro y _
  ring

omit [DecidableEq Omega] in
theorem dissipation_add (K : Omega -> Omega -> Real)
    (f g : Omega -> Real) :
    dissipationOperator K (fun x => f x + g x) =
      fun x => dissipationOperator K f x + dissipationOperator K g x := by
  funext x
  unfold dissipationOperator
  rw [kernelAct_add]
  ring

omit [DecidableEq Omega] in
/-- Telescoping resolvent equation with an exact remainder:
`(I-K)(I+...+K^N)f = f-K^(N+1)f`. -/
theorem finiteResolvent_poisson (K : Omega -> Omega -> Real)
    (N : Nat) (f : Omega -> Real) :
    dissipationOperator K (finiteResolvent K N f) =
      fun x => f x - kernelIterate K (N + 1) f x := by
  induction N with
  | zero => rfl
  | succ N ih =>
      rw [finiteResolvent, dissipation_add, ih]
      funext x
      unfold dissipationOperator
      change (f x - kernelIterate K (N + 1) f x) +
          (kernelIterate K (N + 1) f x -
            kernelIterate K (N + 2) f x) =
        f x - kernelIterate K (N + 2) f x
      ring

omit [DecidableEq Omega] in
/-- The resolvent pairing is exactly the recursively displayed finite sum of
lagged equilibrium correlations. -/
theorem greenKuboPair_finiteResolvent_eq_integratedCorrelation
    (pi : Omega -> Real) (K : Omega -> Omega -> Real)
    (N : Nat) (j k : CenteredCurrent pi) :
    greenKuboPair pi (finiteResolvent K N) j k =
      integratedCorrelation pi K j.1 k.1 N := by
  unfold greenKuboPair
  induction N with
  | zero => rfl
  | succ N ih =>
      rw [finiteResolvent, integratedCorrelation,
        equilibriumPair_add_right, ih]
      rfl

omit [DecidableEq Omega] in
theorem finiteResolvent_solves_of_tail_zero
    (K : Omega -> Omega -> Real) (N : Nat) (f : Omega -> Real)
    (htail : kernelIterate K (N + 1) f = 0) :
    dissipationOperator K (finiteResolvent K N f) = f := by
  rw [finiteResolvent_poisson, htail]
  funext x
  simp

omit [DecidableEq Omega] in
/-- Exact symmetry of a finite time-integrated correlation when both
propagated currents vanish after the declared cutoff. -/
theorem integratedCorrelation_symm_of_tail_zero
    (pi : Omega -> Real) (K : Omega -> Omega -> Real)
    (hdb : forall x y, pi x * K x y = pi y * K y x)
    (N : Nat) (j k : CenteredCurrent pi)
    (hj : kernelIterate K (N + 1) j.1 = 0)
    (hk : kernelIterate K (N + 1) k.1 = 0) :
    integratedCorrelation pi K j.1 k.1 N =
      integratedCorrelation pi K k.1 j.1 N := by
  rw [← greenKuboPair_finiteResolvent_eq_integratedCorrelation,
    ← greenKuboPair_finiteResolvent_eq_integratedCorrelation]
  exact greenKuboPair_symm_of_poisson pi K hdb (finiteResolvent K N) j k
    (finiteResolvent_solves_of_tail_zero K N j.1 hj)
    (finiteResolvent_solves_of_tail_zero K N k.1 hk)

omit [DecidableEq Omega] in
/-- Exact positive semidefiniteness of a finite Green--Kubo sum under the
same explicit tail condition. -/
theorem integratedCorrelation_nonneg_of_tail_zero
    (pi : Omega -> Real) (K : Omega -> Omega -> Real)
    (hpi0 : forall x, 0 <= pi x)
    (hK0 : forall x y, 0 <= K x y)
    (hK1 : forall x, Finset.sum Finset.univ (K x) = 1)
    (hdb : forall x y, pi x * K x y = pi y * K y x)
    (N : Nat) (j : CenteredCurrent pi)
    (hj : kernelIterate K (N + 1) j.1 = 0) :
    0 <= integratedCorrelation pi K j.1 j.1 N := by
  rw [← greenKuboPair_finiteResolvent_eq_integratedCorrelation]
  exact greenKuboPair_nonneg_of_poisson pi K hpi0 hK0 hK1 hdb
    (finiteResolvent K N) j
    (finiteResolvent_solves_of_tail_zero K N j.1 hj)

/-! ## Conditional-resampling projector no-go -/

variable {pi : Omega -> Real} {b : Omega -> B}

omit [DecidableEq Omega] in
/-- Composition formula for the observable action of a finite kernel. -/
theorem kernelAct_comp (K : Omega -> Omega -> Real) (f : Omega -> Real)
    (x : Omega) :
    kernelAct K (kernelAct K f) x =
      Finset.sum Finset.univ fun z =>
        (Finset.sum Finset.univ fun y => K x y * K y z) * f z := by
  unfold kernelAct
  calc
    (Finset.sum Finset.univ fun y =>
        K x y * Finset.sum Finset.univ (fun z => K y z * f z)) =
      Finset.sum Finset.univ (fun y =>
        Finset.sum Finset.univ fun z => K x y * K y z * f z) := by
          apply Finset.sum_congr rfl
          intro y _
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro z _
          ring
    _ = Finset.sum Finset.univ (fun z =>
        Finset.sum Finset.univ fun y => K x y * K y z * f z) :=
          Finset.sum_comm
    _ = Finset.sum Finset.univ fun z =>
        (Finset.sum Finset.univ fun y => K x y * K y z) * f z := by
          apply Finset.sum_congr rfl
          intro z _
          rw [Finset.sum_mul]

/-- The finite conditional-resampling heat bath acts as a projector on
observables. -/
theorem heatBath_kernelAct_idempotent (hpi : forall x, 0 < pi x)
    (f : Omega -> Real) :
    kernelAct (heatBath pi b) (kernelAct (heatBath pi b) f) =
      kernelAct (heatBath pi b) f := by
  funext x
  rw [kernelAct_comp]
  unfold kernelAct
  apply Finset.sum_congr rfl
  intro z _
  rw [heatBath_idempotent hpi x z]

omit [DecidableEq Omega] in
/-- A fibre-weighted current with zero conditional sum is annihilated by
one full-fibre resampling step. -/
theorem heatBath_annihilates_of_fiber_sum_zero
    (j : Omega -> Real)
    (hfiber : forall x,
      Finset.sum (Finset.univ.filter fun y => b y = b x)
        (fun y => pi y * j y) = 0) :
    kernelAct (heatBath pi b) j = 0 := by
  funext x
  unfold kernelAct heatBath
  change (Finset.sum Finset.univ fun y =>
    (if b y = b x then pi y / fiberMass pi b x else 0) * j y) = 0
  simp only [ite_mul, zero_mul]
  rw [Finset.sum_ite, Finset.sum_const_zero, add_zero]
  calc
    (Finset.sum (Finset.univ.filter fun y => b y = b x)
        fun y => pi y / fiberMass pi b x * j y) =
      (Finset.sum (Finset.univ.filter fun y => b y = b x)
        fun y => pi y * j y) / fiberMass pi b x := by
          rw [Finset.sum_div]
          apply Finset.sum_congr rfl
          intro y _
          ring
    _ = 0 := by rw [hfiber x]; simp

/-- After its first step, an idempotent heat bath has no new time
dependence: every positive iterate equals the one-step conditional mean. -/
theorem heatBath_kernelIterate_succ (hpi : forall x, 0 < pi x)
    (f : Omega -> Real) (n : Nat) :
    kernelIterate (heatBath pi b) (n + 1) f =
      kernelAct (heatBath pi b) f := by
  induction n with
  | zero => rfl
  | succ n ih =>
      rw [kernelIterate, ih]
      exact heatBath_kernelAct_idempotent hpi f

/-- Every strictly positive-lag heat-bath correlation equals the one-step
correlation.  A projector has no decaying nonzero memory tail. -/
theorem heatBath_positive_lag_correlation_constant
    (hpi : forall x, 0 < pi x) (f g : Omega -> Real) (n : Nat) :
    lagCorrelation pi (heatBath pi b) f g (n + 1) =
      lagCorrelation pi (heatBath pi b) f g 1 := by
  unfold lagCorrelation
  rw [heatBath_kernelIterate_succ hpi g n]
  rfl

/-- Exact partial-sum formula for a projector dynamics.  Its one-sided
Green--Kubo sum is the equal-time correlation plus `N` copies of the
one-step correlation. -/
theorem heatBath_integratedCorrelation_formula
    (hpi : forall x, 0 < pi x) (f g : Omega -> Real) (N : Nat) :
    integratedCorrelation pi (heatBath pi b) f g N =
      equilibriumPair pi f g + (N : Real) *
        lagCorrelation pi (heatBath pi b) f g 1 := by
  induction N with
  | zero => simp [integratedCorrelation, lagCorrelation, kernelIterate]
  | succ N ih =>
      rw [integratedCorrelation, ih,
        heatBath_positive_lag_correlation_constant hpi f g N]
      push_cast
      ring

/-- If the one-step correlation is nonzero, consecutive projector
Green--Kubo partial sums never stabilize.  This is the complementary branch
of the full-fibre transport no-go. -/
theorem heatBath_integratedCorrelation_not_stable
    (hpi : forall x, 0 < pi x) (f g : Omega -> Real)
    (hstep : lagCorrelation pi (heatBath pi b) f g 1 ≠ 0) (N : Nat) :
    integratedCorrelation pi (heatBath pi b) f g (N + 1) ≠
      integratedCorrelation pi (heatBath pi b) f g N := by
  rw [integratedCorrelation,
    heatBath_positive_lag_correlation_constant hpi f g N]
  intro h
  apply hstep
  linarith

/-- Projector transport no-go: a fibre-centered current has zero
correlation at every strictly positive lag. -/
theorem heatBath_positive_lag_correlation_zero
    (hpi : forall x, 0 < pi x) (j : CenteredCurrent pi)
    (hfiber : forall x,
      Finset.sum (Finset.univ.filter fun y => b y = b x)
        (fun y => pi y * j.1 y) = 0)
    (n : Nat) :
    lagCorrelation pi (heatBath pi b) j.1 j.1 (n + 1) = 0 := by
  unfold lagCorrelation
  rw [heatBath_kernelIterate_succ hpi j.1 n,
    heatBath_annihilates_of_fiber_sum_zero j.1 hfiber]
  unfold equilibriumPair
  simp

/-- Therefore the finite Green--Kubo sum of an annihilated heat-bath current
contains only its equal-time variance, for every cutoff. -/
theorem heatBath_integratedCorrelation_eq_equalTime
    (hpi : forall x, 0 < pi x) (j : CenteredCurrent pi)
    (hfiber : forall x,
      Finset.sum (Finset.univ.filter fun y => b y = b x)
        (fun y => pi y * j.1 y) = 0)
    (N : Nat) :
    integratedCorrelation pi (heatBath pi b) j.1 j.1 N =
      equilibriumPair pi j.1 j.1 := by
  induction N with
  | zero => rfl
  | succ N ih =>
      rw [integratedCorrelation, ih,
        heatBath_positive_lag_correlation_zero hpi j hfiber N, add_zero]

/-- Nonvacuity control: the projector no-go removes memory, not equal-time
fluctuations.  A current of positive equilibrium variance retains a strictly
positive finite Green--Kubo coefficient at every cutoff. -/
theorem heatBath_integratedCorrelation_pos_of_variance_pos
    (hpi : forall x, 0 < pi x) (j : CenteredCurrent pi)
    (hfiber : forall x,
      Finset.sum (Finset.univ.filter fun y => b y = b x)
        (fun y => pi y * j.1 y) = 0)
    (hvar : 0 < equilibriumPair pi j.1 j.1) (N : Nat) :
    0 < integratedCorrelation pi (heatBath pi b) j.1 j.1 N := by
  rw [heatBath_integratedCorrelation_eq_equalTime hpi j hfiber N]
  exact hvar

/-! ### Explicit two-state nonvacuity witness -/

/-- Uniform two-state equilibrium weight. -/
noncomputable def binaryWeight (_ : Fin 2) : Real := 1 / 2

/-- One repaired-visible fibre containing both states. -/
def binaryVisible (_ : Fin 2) : Unit := ()

/-- Opposite signed current on the two states. -/
def binaryCurrentRaw (x : Fin 2) : Real := if x = 0 then 1 else -1

/-- The binary current is centered under the uniform equilibrium weight. -/
noncomputable def binaryCurrent : CenteredCurrent binaryWeight := by
  refine ⟨binaryCurrentRaw, ?_⟩
  unfold equilibriumMean
  rw [Fin.sum_univ_two]
  norm_num [binaryWeight, binaryCurrentRaw]

theorem binaryWeight_pos (x : Fin 2) : 0 < binaryWeight x := by
  norm_num [binaryWeight]

theorem binaryCurrent_fiber_sum_zero (x : Fin 2) :
    Finset.sum (Finset.univ.filter fun y => binaryVisible y = binaryVisible x)
      (fun y => binaryWeight y * binaryCurrent.1 y) = 0 := by
  simp only [binaryVisible, Finset.filter_true]
  rw [Fin.sum_univ_two]
  norm_num [binaryWeight, binaryCurrent, binaryCurrentRaw]

theorem binaryCurrent_variance :
    equilibriumPair binaryWeight binaryCurrent.1 binaryCurrent.1 = 1 := by
  unfold equilibriumPair
  rw [Fin.sum_univ_two]
  norm_num [binaryWeight, binaryCurrent, binaryCurrentRaw]

/-- Concrete witness: full-fibre resampling kills all positive-lag memory,
while its one-sided finite coefficient is exactly one at every cutoff. -/
theorem binary_heatBath_integratedCorrelation_eq_one (N : Nat) :
    integratedCorrelation binaryWeight
      (heatBath binaryWeight binaryVisible) binaryCurrent.1 binaryCurrent.1 N = 1 := by
  rw [heatBath_integratedCorrelation_eq_equalTime binaryWeight_pos binaryCurrent
    binaryCurrent_fiber_sum_zero N, binaryCurrent_variance]

/-! ## Negative control: nondissipating dynamics -/

/-- Identity transition kernel. -/
def identityKernel (x y : Omega) : Real := if y = x then 1 else 0

theorem identityKernel_act (f : Omega -> Real) :
    kernelAct (identityKernel : Omega -> Omega -> Real) f = f := by
  funext x
  simp [kernelAct, identityKernel]

/-- The identity dynamics has zero positive generator. -/
theorem identityKernel_dissipation (f : Omega -> Real) :
    dissipationOperator (identityKernel : Omega -> Omega -> Real) f = 0 := by
  funext x
  unfold dissipationOperator
  rw [identityKernel_act]
  simp

/-- Negative control: a nonzero centered current has no Poisson solution for
the nondissipating identity dynamics.  Thus the inverse/tail premise in the
Green--Kubo theorem cannot be erased. -/
theorem identityKernel_no_poisson_of_ne_zero (j : CenteredCurrent pi)
    (hj : j.1 ≠ 0) :
    ¬ (exists f : Omega -> Real,
      dissipationOperator (identityKernel : Omega -> Omega -> Real) f = j.1) := by
  rintro ⟨f, hf⟩
  rw [identityKernel_dissipation] at hf
  exact hj hf.symm

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.dissipation_eq_dirichlet
#print axioms OPH.Thermodynamics.greenKuboPair_symm_of_poisson
#print axioms OPH.Thermodynamics.greenKuboPair_nonneg_of_poisson
#print axioms OPH.Thermodynamics.integratedCorrelation_symm_of_tail_zero
#print axioms OPH.Thermodynamics.integratedCorrelation_nonneg_of_tail_zero
#print axioms OPH.Thermodynamics.heatBath_integratedCorrelation_formula
#print axioms OPH.Thermodynamics.heatBath_integratedCorrelation_not_stable
#print axioms OPH.Thermodynamics.heatBath_integratedCorrelation_eq_equalTime
#print axioms OPH.Thermodynamics.binary_heatBath_integratedCorrelation_eq_one
#print axioms OPH.Thermodynamics.identityKernel_no_poisson_of_ne_zero
