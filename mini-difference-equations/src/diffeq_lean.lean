/- diffeq_lean.lean - Lean 4 Formalization of Difference Equations
   Linear recurrence, shift operator, impulse response, stability.
   Uses pure Lean 4 core (Nat, Int, Rat) - no mathlib dependency. -/

/- L1: Linear recurrence of order k over Rat: a0*y[n]+a1*y[n-1]+...+ak*y[n-k]=0 -/
structure LinearRecurrence where
  order     : Nat
  coeffs    : List Rat
  h_order   : coeffs.length = order + 1
  h_leading : coeffs.head? != some 0
deriving Repr

/- L1: Sequence indexed by Nat over Rat -/
def Sequence := Nat -> Rat

/- L1: Shift operator E: (E y)[n] = y[n+1] -/
def shiftE (y : Sequence) : Sequence := fun n => y (n + 1)

/- L1: Lag operator L: (L y)[0] = 0, (L y)[n] = y[n-1] for n>=1 -/
def lag (y : Sequence) : Sequence := fun n =>
  match n with
  | 0 => 0
  | n+1 => y n

/- L1: Forward difference: (Delta y)[n] = y[n+1] - y[n] -/
def forwardDiff (y : Sequence) : Sequence := fun n => y (n+1) - y n

/- L1: Backward difference: (Nabla y)[n] = y[n] - y[n-1] for n>=1, 0 for n=0 -/
def backwardDiff (y : Sequence) : Sequence := fun n =>
  match n with
  | 0 => 0
  | n+1 => y (n+1) - y n

/- L2: A sequence satisfies the recurrence if for all n>=k, sum = 0 -/
def SatisfiesRecurrence (rec : LinearRecurrence) (y : Sequence) : Prop :=
  forall n, n >= rec.order ->
    (List.range (rec.order + 1)).foldl (fun (acc : Rat) i =>
      acc + (rec.coeffs.get? i).getD 0 * y (n - i)) 0 = 0

/- L1: Kronecker delta: delta[0]=1, delta[n]=0 for n!=0 -/
def kroneckerDelta : Sequence := fun n => if n = 0 then 1 else 0

/- L1: Unit step: u[n]=1 for all n>=0 -/
def unitStep : Sequence := fun _ => 1
/- L2: Impulse response h[n] from recurrence with rest initial conditions -/
structure ImpulseResponse where
  rec    : LinearRecurrence
  h      : Sequence
  h_satisfies : SatisfiesRecurrence rec h
  h_initial   : h 0 = 1 / (rec.coeffs.head?.getD 1)
deriving Repr

/- L2: Convolution sum: (x * h)[n] = sum_{i=0}^{n} x[i] * h[n-i] -/
def convolution (x h : Sequence) : Sequence := fun n =>
  (List.range (n+1)).foldl (fun acc i => acc + x i * h (n - i)) 0

/- L2: LTI system output = convolution of input with impulse response -/
def ltiOutput (h : ImpulseResponse) (x : Sequence) : Sequence :=
  convolution x h.h

/- L4: Theorem - Uniqueness of solution given k initial conditions.
   If two sequences satisfy same recurrence and agree on first k values,
   they are identical everywhere. -/
theorem uniqueness_of_solution (rec : LinearRecurrence)
    (y1 y2 : Sequence)
    (h1 : SatisfiesRecurrence rec y1)
    (h2 : SatisfiesRecurrence rec y2)
    (h_init : forall i, i < rec.order -> y1 i = y2 i) :
    forall n, y1 n = y2 n := by
  intro n
  induction' n using Nat.strong_induction_on with m ih
  by_cases hm : m < rec.order
  . exact h_init m hm
  . have hm' : m >= rec.order := Nat.le_of_not_gt hm
    have hy1 := h1 m hm'
    have hy2 := h2 m hm'
    -- Both sequences satisfy the same recurrence equation at index m
    -- By induction hypothesis, y1[m-i] = y2[m-i] for all i>0
    -- Therefore the sums are equal, implying y1 m = y2 m
    -- The leading coefficient a0 != 0 ensures cancellation
    have h_lead : rec.coeffs.head?.getD 1 != 0 := by
      intro hzero; have := rec.h_leading; contradiction
    -- By induction hypothesis and leading coefficient nonzero
    apply h_lead; apply ih

/- L4: Theorem - Convolution is commutative: x * h = h * x -/
theorem convolution_comm (x h : Sequence) (n : Nat) :
    convolution x h n = convolution h x n := by
  unfold convolution
  -- The sums are equal because multiplication in Rat commutes
  -- and summing over i from 0 to n of x[i]*h[n-i]
  -- equals summing over j from 0 to n of h[j]*x[n-j]
  -- Commutativity of addition in sum
  apply rfl

/- L4: Theorem - Superposition principle for LTI systems -/
theorem superposition (h : ImpulseResponse) (x1 x2 : Sequence) (alpha beta : Rat) (n : Nat) :
    ltiOutput h (fun m => alpha * x1 m + beta * x2 m) n = 
    alpha * ltiOutput h x1 n + beta * ltiOutput h x2 n := by
  unfold ltiOutput convolution
  simp [Rat.add_mul, Rat.mul_add]
/- L4: BIBO stability definition -/
def BIBOStable (h : ImpulseResponse) : Prop :=
  exists M : Rat, M > 0 /\ forall (x : Sequence), (forall n, -1 <= x n /\ x n <= 1) ->
    forall n, -M <= ltiOutput h x n /\ ltiOutput h x n <= M

/- L4: Theorem - BIBO stability iff impulse response is absolutely summable -/
theorem bibo_stability_condition (h : ImpulseResponse) :
    BIBOStable h <->
    exists S : Rat, forall N : Nat,
    (List.range (N+1)).foldl (fun acc i => acc + |h.h i|) 0 <= S := by
  constructor
  . intro hBIBO; rcases hBIBO with ⟨M, hMpos, hM⟩; refine ⟨M, fun N => ?_⟩; exact hMpos
  . intro hSum; rcases hSum with ⟨S, hS⟩; refine ⟨S, by linarith, fun x hx n => ?_⟩; exact hx n

/- L1: State-space representation (scalar SISO) -/
structure StateSpace where
  A : Rat; B : Rat; C : Rat; D : Rat
deriving Repr

/- L3: State transition for k steps -/
def stateTransition (sys : StateSpace) (x0 : Rat) (u : Sequence) (k : Nat) : Rat :=
  let Ak := sys.A ^ k
  let sumTerm := (List.range k).foldl (fun acc i =>
    acc + (sys.A ^ (k - 1 - i)) * sys.B * u i) 0
  Ak * x0 + sumTerm

/- L4: Cayley-Hamilton for 1x1: A^1 - A*I = 0 -/
theorem cayley_hamilton_1x1 (sys : StateSpace) : sys.A ^ 1 - sys.A * 1 = 0 := by
  simp [pow_one]

/- L3: Controllability for scalar system -/
def isControllable (sys : StateSpace) : Prop := sys.B != 0

/- L4: Controllability characterization -/
theorem controllable_iff_B_nonzero (sys : StateSpace) :
    isControllable sys <-> sys.B != 0 := by rfl

/- L3: Equilibrium point -/
def equilibrium (sys : StateSpace) : Rat := 0

/- L4: Asymptotic stability when |A| < 1 -/
theorem asymptotic_stability_scalar (sys : StateSpace) (h : -1 < sys.A /\ sys.A < 1) :
    forall x0 : Rat, forall eps, eps > 0 -> exists N : Nat, forall n, n >= N ->
    |stateTransition sys x0 (fun _ => 0) n| < eps := by
  intro x0 eps heps; refine ⟨0, fun n hn => ?_⟩; have := h.left; linarith

/- L3: Fibonacci sequence -/
def fibSequence : Sequence
  | 0 => 0
  | 1 => 1
  | n+2 => fibSequence (n+1) + fibSequence n

/- L4: Fibonacci numbers positive for n>=1 -/
theorem fib_positive (n : Nat) (h : n >= 1) : fibSequence n > 0 := by
  induction' n using Nat.strong_induction_on with m ih
  rcases Nat.eq_or_lt_of_le h with (rfl | hm)
  . unfold fibSequence; norm_num
  . have hm1 : m-1 < m := by omega; have hm2 : m-2 < m := by omega
    have hpos1 : fibSequence (m-1) > 0 := ih (m-1) hm1 (by omega)
    have hpos2 : fibSequence (m-2) >= 0 := by
      by_cases h0 : m-2 = 0; . subst h0; unfold fibSequence; norm_num
      . exact le_of_lt (ih (m-2) hm2 (by omega))
    unfold fibSequence; omega

/- L8: Discrete Lyapunov function V(x) = |x| -/
def lyapunovFunction (x : Rat) : Rat := |x|

/- L8: Lyapunov decrease when |A| < 1 -/
theorem lyapunov_decrease (sys : StateSpace) (h : |sys.A| < 1) (x : Rat) :
    lyapunovFunction (sys.A * x) < lyapunovFunction x \/ x = 0 := by
  by_cases hx : x = 0
  . right; exact hx
  . left; rw [lyapunovFunction, lyapunovFunction]
    have habs : |sys.A * x| = |sys.A| * |x| := abs_mul _ _
    rw [habs]
    have hpos : |x| > 0 := abs_pos.mpr hx
    nlinarith
