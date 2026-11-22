import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FunProp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# The elliptic endpoint identity for triangular spherical prisms

`R` is the observation radius, `s` is a radial endpoint, and `c` is a direction cosine.
The Legendre integrals use signed intervals from zero to the actual amplitude.
The theorems below concern exact real arithmetic.
-/

noncomputable section

open MeasureTheory

namespace TSP

def distance (R s c : ℝ) : ℝ :=
  Real.sqrt ((s - R) ^ 2 + 2 * R * s * (1 - c))

def radialPrimitive (R s c : ℝ) : ℝ :=
  let D := distance R s c
  D ^ 3 / 3 + R * c / 2 *
    ((s - R * c) * D + R ^ 2 * (1 - c ^ 2) * Real.log (s - R * c + D))

def numeratorPlus (R s h : ℝ) : ℝ :=
  (s - R) * ((s - R) ^ 2 * (2 * s + R) - 3 * R ^ 3 * (1 - h ^ 2)) / 12

def numeratorMinus (R s h : ℝ) : ℝ :=
  (s + R) * ((s + R) ^ 2 * (2 * s - R) + 3 * R ^ 3 * (1 - h ^ 2)) / 12

def directionCosine (h φ : ℝ) : ℝ := h * Real.cos φ

def ellipticRoot (m t : ℝ) : ℝ := Real.sqrt (1 - m * Real.sin t ^ 2)

def legendreF (m v : ℝ) : ℝ := ∫ t in (0 : ℝ)..v, (ellipticRoot m t)⁻¹

def legendreE (m v : ℝ) : ℝ := ∫ t in (0 : ℝ)..v, ellipticRoot m t

def legendrePi (ν m v : ℝ) : ℝ :=
  ∫ t in (0 : ℝ)..v, ((1 - ν * Real.sin t ^ 2) * ellipticRoot m t)⁻¹

def ellipticParameter (R s h : ℝ) : ℝ := -(4 * R * s * h) / distance R s h ^ 2

def characteristicPlus (h : ℝ) : ℝ := -(2 * h) / (1 - h)

def characteristicMinus (h : ℝ) : ℝ := 2 * h / (1 + h)

def ellipticEndpoint (R s h φ : ℝ) : ℝ :=
  let c := directionCosine h φ
  let d := distance R s h
  let m := ellipticParameter R s h
  R ^ 3 * h * Real.sin φ / 2 * Real.log (s - R * c + distance R s c) +
    R ^ 2 * d / 2 * legendreE m (φ / 2) -
    R ^ 2 * (3 * R ^ 2 + s ^ 2) / (6 * d) * legendreF m (φ / 2) +
    2 * numeratorPlus R s h / ((1 - h) * d) *
      legendrePi (characteristicPlus h) m (φ / 2) +
    2 * numeratorMinus R s h / ((1 + h) * d) *
      legendrePi (characteristicMinus h) m (φ / 2)

theorem distance_rad_sq {R s c : ℝ} (hR : 0 < R) (hs : 0 < s) (hc : c ≤ 1) :
    (distance R s c) ^ 2 = (s - R) ^ 2 + 2 * R * s * (1 - c) := by
  have hsub : 0 ≤ 1 - c := sub_nonneg.mpr hc
  have h : 0 ≤ (s - R) ^ 2 + 2 * R * s * (1 - c) := by positivity
  exact Real.sq_sqrt h

theorem distance_pos {R s c : ℝ} (hR : 0 < R) (hs : 0 < s) (hc : c < 1) :
    0 < distance R s c := by
  unfold distance
  apply Real.sqrt_pos.2
  have h : 0 < 2 * R * s * (1 - c) :=
    mul_pos (by positivity) (sub_pos.mpr hc)
  linarith [sq_nonneg (s - R)]

theorem distance_sq {R s c : ℝ} (hR : 0 < R) (hs : 0 < s) (hc : c ≤ 1) :
    (distance R s c) ^ 2 = R ^ 2 + s ^ 2 - 2 * R * s * c := by
  rw [distance_rad_sq hR hs hc]
  ring

theorem one_sub_sq_pos {c : ℝ} (hl : -1 < c) (hu : c < 1) : 0 < 1 - c ^ 2 := by
  have h : 0 < (1 - c) * (1 + c) := mul_pos (sub_pos.mpr hu) (by linarith)
  linarith only [h]

private theorem logArgument_pos {R s c : ℝ} (hR : 0 < R) (hs : 0 < s)
    (hl : -1 < c) (hu : c < 1) : 0 < s - R * c + distance R s c := by
  have hD := distance_pos hR hs hu
  have hsq := distance_sq hR hs hu.le
  have hprod : (s - R * c + distance R s c) * (distance R s c - (s - R * c)) =
      R ^ 2 * (1 - c ^ 2) := by linear_combination hsq
  by_contra! h
  have hdiff : 0 ≤ distance R s c - (s - R * c) := by linarith
  have hnonpos := mul_nonpos_of_nonpos_of_nonneg h hdiff
  rw [hprod] at hnonpos
  have hpos := mul_pos (pow_pos hR 2) (one_sub_sq_pos hl hu)
  linarith

theorem dir_bounds {h : ℝ} (hh : 0 < h) (h1 : h < 1) (φ : ℝ) :
    -1 < directionCosine h φ ∧ directionCosine h φ < 1 := by
  have hlo := mul_le_mul_of_nonneg_left (Real.neg_one_le_cos φ) hh.le
  have hhi := mul_le_mul_of_nonneg_left (Real.cos_le_one φ) hh.le
  unfold directionCosine
  constructor <;> linarith

private theorem ellipticRoot_pos {m : ℝ} (hm : m ≤ 0) (t : ℝ) :
    0 < ellipticRoot m t := by
  unfold ellipticRoot
  exact Real.sqrt_pos.2 (by linarith [mul_nonpos_of_nonpos_of_nonneg hm (sq_nonneg (Real.sin t))])

private theorem characteristicFactor_pos {ν : ℝ} (hν : ν < 1) (t : ℝ) :
    0 < 1 - ν * Real.sin t ^ 2 := by
  have hs : Real.sin t ^ 2 ≤ 1 := by linarith [Real.sin_sq_add_cos_sq t, sq_nonneg (Real.cos t)]
  by_cases hp : 0 ≤ ν
  · linarith [mul_le_mul_of_nonneg_left hs hp]
  · linarith [mul_nonpos_of_nonpos_of_nonneg (le_of_not_ge hp) (sq_nonneg (Real.sin t))]

private theorem hasDerivAt_intervalPrimitive {f : ℝ → ℝ} (hf : Continuous f) (v : ℝ) :
    HasDerivAt (fun u => ∫ t in (0 : ℝ)..u, f t) (f v) v :=
  intervalIntegral.integral_hasDerivAt_right (hf.intervalIntegrable 0 v)
    (hf.stronglyMeasurableAtFilter _ _) hf.continuousAt

private theorem continuous_ellipticRoot (m : ℝ) : Continuous (ellipticRoot m) := by
  unfold ellipticRoot
  fun_prop

private theorem hasDerivAt_legendreF {m : ℝ} (hm : m ≤ 0) (v : ℝ) :
    HasDerivAt (legendreF m) (ellipticRoot m v)⁻¹ v :=
  hasDerivAt_intervalPrimitive
    ((continuous_ellipticRoot m).inv₀ (fun t => (ellipticRoot_pos hm t).ne')) v

private theorem hasDerivAt_legendreE (m v : ℝ) :
    HasDerivAt (legendreE m) (ellipticRoot m v) v :=
  hasDerivAt_intervalPrimitive (continuous_ellipticRoot m) v

private theorem hasDerivAt_legendrePi {ν m : ℝ} (hν : ν < 1) (hm : m ≤ 0) (v : ℝ) :
    HasDerivAt (legendrePi ν m) ((1 - ν * Real.sin v ^ 2) * ellipticRoot m v)⁻¹ v := by
  have hfactor : Continuous (fun t => 1 - ν * Real.sin t ^ 2) := by fun_prop
  exact hasDerivAt_intervalPrimitive ((hfactor.mul (continuous_ellipticRoot m)).inv₀
    (fun t => (mul_pos (characteristicFactor_pos hν t) (ellipticRoot_pos hm t)).ne')) v

private theorem half_sin_squared (φ : ℝ) : 2 * Real.sin (φ / 2) ^ 2 = 1 - Real.cos φ := by
  have h : φ = 2 * (φ / 2) := by ring
  calc
    _ = 1 - (2 * Real.cos (φ / 2) ^ 2 - 1) := by
      linarith [Real.sin_sq_add_cos_sq (φ / 2)]
    _ = 1 - Real.cos φ := by rw [← Real.cos_two_mul, ← h]

private theorem distance_eq_scaled_root {R s h : ℝ}
    (hR : 0 < R) (hs : 0 < s) (h1 : h < 1) (φ : ℝ) :
    distance R s (directionCosine h φ) =
      distance R s h * ellipticRoot (ellipticParameter R s h) (φ / 2) := by
  have hd := distance_pos hR hs h1
  have hd_sq := distance_rad_sq hR hs h1.le
  have hsin := half_sin_squared φ
  have hrad : (s - R) ^ 2 + 2 * R * s * (1 - directionCosine h φ) =
      (distance R s h) ^ 2 * (1 - ellipticParameter R s h * Real.sin (φ / 2) ^ 2) := by
    unfold ellipticParameter directionCosine
    field_simp [hd.ne']
    linear_combination -hd_sq - 2 * R * s * h * hsin
  change Real.sqrt _ = distance R s h * Real.sqrt _
  rw [hrad, Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hd.le]

private theorem factorPlus {h : ℝ} (hh : h < 1) (φ : ℝ) :
    1 - characteristicPlus h * Real.sin (φ / 2) ^ 2 =
      (1 - directionCosine h φ) / (1 - h) := by
  unfold characteristicPlus directionCosine
  have hne : 1 - h ≠ 0 := sub_ne_zero.mpr hh.ne'
  field_simp [hne]
  linear_combination h * half_sin_squared φ

private theorem factorMinus {h : ℝ} (hh : 0 < h) (φ : ℝ) :
    1 - characteristicMinus h * Real.sin (φ / 2) ^ 2 =
      (1 + directionCosine h φ) / (1 + h) := by
  unfold characteristicMinus directionCosine
  have hne : 1 + h ≠ 0 := by positivity
  field_simp [hne]
  linear_combination -h * half_sin_squared φ

private theorem characteristicPlus_lt_one {h : ℝ} (hh : 0 < h) (h1 : h < 1) :
    characteristicPlus h < 1 := by
  have hden : 0 < 1 - h := sub_pos.mpr h1
  have hle : characteristicPlus h ≤ 0 := div_nonpos_of_nonpos_of_nonneg (by linarith) hden.le
  linarith

private theorem characteristicMinus_lt_one {h : ℝ} (hh : 0 < h) (h1 : h < 1) :
    characteristicMinus h < 1 := by
  unfold characteristicMinus
  apply (div_lt_iff₀ (by positivity : 0 < 1 + h)).2
  linarith

private theorem hasDerivAt_logTerm {R s h : ℝ}
    (hR : 0 < R) (hs : 0 < s) (hh : 0 < h) (h1 : h < 1) (φ : ℝ) :
    let c := directionCosine h φ
    let D := distance R s c
    HasDerivAt
      (fun t => R ^ 3 * h * Real.sin t / 2 *
        Real.log (s - R * directionCosine h t + distance R s (directionCosine h t)))
      (R ^ 3 * c / 2 * Real.log (s - R * c + D) +
        R ^ 4 * (h ^ 2 - c ^ 2) * (D + s) / (2 * D * (s - R * c + D))) φ := by
  dsimp only
  let c := directionCosine h φ
  let D := distance R s c
  have ⟨hcl, hcu⟩ := dir_bounds hh h1 φ
  have hD : 0 < D := distance_pos hR hs hcu
  have hlog := logArgument_pos hR hs hcl hcu
  have hcderiv : HasDerivAt (directionCosine h) (-(h * Real.sin φ)) φ := by
    simpa [directionCosine, mul_neg] using (Real.hasDerivAt_cos φ).const_mul h
  have hr := ((hcderiv.const_sub 1).const_mul (2 * R * s)).const_add ((s - R) ^ 2)
  have hrad : 0 < (s - R) ^ 2 + 2 * R * s * (1 - c) := by
    have h : 0 < 2 * R * s * (1 - c) :=
      mul_pos (by positivity) (sub_pos.mpr hcu)
    linarith [sq_nonneg (s - R)]
  have hdist : HasDerivAt (fun t => distance R s (directionCosine h t))
      (R * s * h * Real.sin φ / D) φ := by
    refine ((Real.hasDerivAt_sqrt hrad.ne').comp φ hr).congr_deriv ?_
    change 1 / (2 * D) * (2 * R * s * -( -(h * Real.sin φ))) = _
    field_simp [hD.ne']
    ring
  have ha := ((hcderiv.const_mul R).const_sub s).add hdist
  have hl := (Real.hasDerivAt_log hlog.ne').comp φ ha
  have hk := ((Real.hasDerivAt_sin φ).const_mul (R ^ 3 * h)).div_const 2
  refine (hk.mul hl).congr_deriv ?_
  have hsin : (h * Real.sin φ) ^ 2 = h ^ 2 - c ^ 2 := by
    dsimp [c, directionCosine]
    linear_combination h ^ 2 * Real.sin_sq_add_cos_sq φ
  change R ^ 3 * h * Real.cos φ / 2 * Real.log (s - R * c + D) +
      R ^ 3 * h * Real.sin φ / 2 *
        ((s - R * c + D)⁻¹ * (-(R * -(h * Real.sin φ)) + R * s * h * Real.sin φ / D)) = _
  have hcos : h * Real.cos φ = c := rfl
  rw [mul_assoc (R ^ 3) h (Real.cos φ), hcos, ← hsin]
  field_simp [hD.ne', hlog.ne']
  ring

private theorem endpoint_polynomial (R s c h : ℝ) :
    let q := R ^ 2 + s ^ 2 - 2 * R * s * c
    R ^ 3 * (h ^ 2 - c ^ 2) * (R - s * c) / 2 + R ^ 2 * q * (1 - c ^ 2) / 4 -
      R ^ 2 * (3 * R ^ 2 + s ^ 2) * (1 - c ^ 2) / 12 +
      numeratorPlus R s h * (1 + c) + numeratorMinus R s h * (1 - c) =
    q ^ 2 / 3 + R * c * (s - R * c) * q / 2 := by
  simp only [numeratorPlus, numeratorMinus]
  ring

private theorem rationalize_log_derivative {R s c : ℝ} (h : ℝ)
    (hR : 0 < R) (hs : 0 < s) (hl : -1 < c) (hu : c < 1) :
    let D := distance R s c
    R ^ 4 * (h ^ 2 - c ^ 2) * (D + s) / (2 * D * (s - R * c + D)) =
      (R ^ 3 * (h ^ 2 - c ^ 2) * (R - s * c) / 2) / (D * (1 - c ^ 2)) +
      R ^ 3 * c * (h ^ 2 - c ^ 2) / (2 * (1 - c ^ 2)) := by
  dsimp only
  let D := distance R s c
  let a := s - R * c
  have hD := distance_pos hR hs hu
  have hA := logArgument_pos hR hs hl hu
  have hC := one_sub_sq_pos hl hu
  have hsq : D ^ 2 = R ^ 2 + s ^ 2 - 2 * R * s * c := distance_sq hR hs hu.le
  have hprod : (a + D) * (D - a) = R ^ 2 * (1 - c ^ 2) := by
    dsimp [a]
    linear_combination hsq
  have hrel : (D + s) * (D - a) = R * (R - s * c) + R * c * D := by
    dsimp [a]
    linear_combination hsq
  change R ^ 4 * (h ^ 2 - c ^ 2) * (D + s) / (2 * D * (a + D)) = _
  calc
    _ = (R ^ 2 * (h ^ 2 - c ^ 2) * (D + s) * ((a + D) * (D - a))) /
          (2 * D * (a + D) * (1 - c ^ 2)) := by
          rw [hprod]
          field_simp [hD.ne', hA.ne', hC.ne']
          ring
    _ = R ^ 2 * (h ^ 2 - c ^ 2) * ((D + s) * (D - a)) /
          (2 * D * (1 - c ^ 2)) := by
          field_simp [hD.ne', hA.ne', hC.ne']
          ring
    _ = _ := by
          rw [hrel]
          field_simp [hD.ne', hC.ne']
          ring

private theorem endpoint_derivative_algebra {R s c : ℝ} (h : ℝ)
    (hR : 0 < R) (hs : 0 < s) (hl : -1 < c) (hu : c < 1) :
    let D := distance R s c
    (R ^ 3 * c / 2 * Real.log (s - R * c + D) +
      R ^ 4 * (h ^ 2 - c ^ 2) * (D + s) / (2 * D * (s - R * c + D))) +
      R ^ 2 * D / 4 - R ^ 2 * (3 * R ^ 2 + s ^ 2) / (12 * D) +
      numeratorPlus R s h / ((1 - c) * D) +
      numeratorMinus R s h / ((1 + c) * D) =
    radialPrimitive R s c / (1 - c ^ 2) + R ^ 3 * c * (h ^ 2 - c ^ 2) / (2 * (1 - c ^ 2)) := by
  dsimp only
  rw [rationalize_log_derivative h hR hs hl hu]
  let D := distance R s c
  have hD : 0 < D := distance_pos hR hs hu
  have hm : 0 < 1 - c := sub_pos.mpr hu
  have hp : 0 < 1 + c := by linarith
  have hC := one_sub_sq_pos hl hu
  have hq : D ^ 2 = R ^ 2 + s ^ 2 - 2 * R * s * c := distance_sq hR hs hu.le
  have hpoly := endpoint_polynomial R s c h
  dsimp only at hpoly
  rw [← hq] at hpoly
  have hrat :
      (R ^ 3 * (h ^ 2 - c ^ 2) * (R - s * c) / 2) / (D * (1 - c ^ 2)) +
        R ^ 2 * D / 4 - R ^ 2 * (3 * R ^ 2 + s ^ 2) / (12 * D) +
        numeratorPlus R s h / ((1 - c) * D) + numeratorMinus R s h / ((1 + c) * D) =
      (D ^ 3 / 3 + R * c / 2 * ((s - R * c) * D)) / (1 - c ^ 2) := by
    field_simp [hD.ne', hm.ne', hp.ne', hC.ne']
    linear_combination 576 * D ^ 3 * (1 - c ^ 2) ^ 2 * hpoly
  unfold radialPrimitive
  change _ = (D ^ 3 / 3 + R * c / 2 *
    ((s - R * c) * D + R ^ 2 * (1 - c ^ 2) * Real.log (s - R * c + D))) /
      (1 - c ^ 2) + _
  have hsplit :
      (D ^ 3 / 3 + R * c / 2 *
        ((s - R * c) * D + R ^ 2 * (1 - c ^ 2) * Real.log (s - R * c + D))) /
          (1 - c ^ 2) =
      (D ^ 3 / 3 + R * c / 2 * ((s - R * c) * D)) / (1 - c ^ 2) +
        R ^ 3 * c / 2 * Real.log (s - R * c + D) := by
    field_simp [hC.ne']
    ring
  rw [hsplit, ← hrat]
  ring

/-- The derivative of the single-radius primitive includes the radius-independent correction. -/
theorem hasDerivAt_ellipticEndpoint {R s h : ℝ}
    (hR : 0 < R) (hs : 0 < s) (hh : 0 < h) (h1 : h < 1) (φ : ℝ) :
    let c := directionCosine h φ
    HasDerivAt (ellipticEndpoint R s h)
      (radialPrimitive R s c / (1 - c ^ 2) +
        R ^ 3 * c * (h ^ 2 - c ^ 2) / (2 * (1 - c ^ 2))) φ := by
  let c := directionCosine h φ
  let D := distance R s c
  let d := distance R s h
  let m := ellipticParameter R s h
  let w := ellipticRoot m (φ / 2)
  have ⟨hcl, hcu⟩ := dir_bounds hh h1 φ
  have hd : 0 < d := distance_pos hR hs h1
  have hm : m ≤ 0 := div_nonpos_of_nonpos_of_nonneg
    (neg_nonpos.mpr (by positivity : 0 ≤ 4 * R * s * h)) (sq_nonneg d)
  have hw : 0 < w := ellipticRoot_pos hm _
  have hscale : D = d * w := distance_eq_scaled_root hR hs h1 φ
  have hminus : 0 < 1 - h := sub_pos.mpr h1
  have hplus : 0 < 1 + h := by positivity
  have hcminus : 0 < 1 - c := sub_pos.mpr hcu
  have hcplus : 0 < 1 + c := by linarith
  have hhalf := (hasDerivAt_id φ).div_const (2 : ℝ)
  have he := ((hasDerivAt_legendreE m (φ / 2)).comp φ hhalf).const_mul (R ^ 2 * d / 2)
  have hf := ((hasDerivAt_legendreF hm (φ / 2)).comp φ hhalf).const_mul
    (R ^ 2 * (3 * R ^ 2 + s ^ 2) / (6 * d))
  have ha := ((hasDerivAt_legendrePi (characteristicPlus_lt_one hh h1) hm (φ / 2)).comp φ
    hhalf).const_mul (2 * numeratorPlus R s h / ((1 - h) * d))
  have hb := ((hasDerivAt_legendrePi (characteristicMinus_lt_one hh h1) hm (φ / 2)).comp φ
    hhalf).const_mul (2 * numeratorMinus R s h / ((1 + h) * d))
  have he' := he.congr_deriv (show R ^ 2 * d / 2 * (w * (1 / 2)) = R ^ 2 * D / 4 by
    rw [hscale]
    ring)
  have hf' := hf.congr_deriv (show R ^ 2 * (3 * R ^ 2 + s ^ 2) / (6 * d) * (w⁻¹ * (1 / 2)) =
      R ^ 2 * (3 * R ^ 2 + s ^ 2) / (12 * D) by
    rw [hscale]
    field_simp [hd.ne', hw.ne']
    ring)
  have ha' := ha.congr_deriv (show 2 * numeratorPlus R s h / ((1 - h) * d) *
      (((1 - characteristicPlus h * Real.sin (φ / 2) ^ 2) * w)⁻¹ * (1 / 2)) =
      numeratorPlus R s h / ((1 - c) * D) by
    rw [factorPlus h1, hscale]
    change _ = numeratorPlus R s h / ((1 - directionCosine h φ) * (d * w))
    field_simp [hd.ne', hw.ne', hminus.ne', hcminus.ne']
    ring)
  have hb' := hb.congr_deriv (show 2 * numeratorMinus R s h / ((1 + h) * d) *
      (((1 - characteristicMinus h * Real.sin (φ / 2) ^ 2) * w)⁻¹ * (1 / 2)) =
      numeratorMinus R s h / ((1 + c) * D) by
    rw [factorMinus hh, hscale]
    change _ = numeratorMinus R s h / ((1 + directionCosine h φ) * (d * w))
    field_simp [hd.ne', hw.ne', hplus.ne', hcplus.ne']
    ring)
  exact (((((hasDerivAt_logTerm hR hs hh h1 φ).add he').sub hf').add ha').add hb').congr_deriv
    (endpoint_derivative_algebra h hR hs hcl hcu)


/-- The radial primitive integrates s times the source–observer distance. -/
theorem hasDerivAt_radialPrimitive {R s c : ℝ} (hR : 0 < R) (hs : 0 < s)
    (hl : -1 < c) (hu : c < 1) :
    HasDerivAt (fun t => radialPrimitive R t c) (s * distance R s c) s := by
  let D := distance R s c
  let a := s - R * c
  have hD : 0 < D := distance_pos hR hs hu
  have hA : 0 < a + D := logArgument_pos hR hs hl hu
  have hsq : D ^ 2 = a ^ 2 + R ^ 2 * (1 - c ^ 2) := by
    dsimp [a]
    linear_combination distance_sq hR hs hu.le
  have hid := hasDerivAt_id s
  have hr := ((hid.sub_const R).pow 2).add (hid.mul_const (2 * R * (1 - c)))
  have hrad : 0 < (s - R) ^ 2 + s * (2 * R * (1 - c)) := by
    have h : 0 < s * (2 * R * (1 - c)) :=
      mul_pos hs (mul_pos (by positivity) (sub_pos.mpr hu))
    linarith [sq_nonneg (s - R)]
  have hdist : HasDerivAt (fun t => distance R t c) (a / D) s := by
    have h := (Real.hasDerivAt_sqrt hrad.ne').comp s hr
    have hval : (s - R) ^ 2 + s * (2 * R * (1 - c)) =
        (s - R) ^ 2 + 2 * R * s * (1 - c) := by ring
    have heq : (fun t : ℝ => Real.sqrt ((t - R) ^ 2 + t * (2 * R * (1 - c)))) =
        (fun t => distance R t c) := by
      funext t
      unfold distance
      congr 1
      ring
    change HasDerivAt (fun t : ℝ => Real.sqrt ((t - R) ^ 2 + t * (2 * R * (1 - c)))) _ s at h
    rw [heq] at h
    refine h.congr_deriv ?_
    simp only [id_eq, hval]
    change 1 / (2 * D) * (2 * (s - R) ^ (2 - 1) * 1 + 1 * (2 * R * (1 - c))) = a / D
    dsimp [a]
    field_simp [hD.ne']
    ring
  have hlin := hid.sub_const (R * c)
  have hlog := (Real.hasDerivAt_log hA.ne').comp s (hlin.add hdist)
  have hlog' := hlog.congr_deriv (show (a + D)⁻¹ * (1 + a / D) = 1 / D by
    field_simp [hD.ne', hA.ne']
    ring)
  have hbracket :=
    ((hlin.mul hdist).add (hlog'.const_mul (R ^ 2 * (1 - c ^ 2)))).congr_deriv
      (show 1 * D + a * (a / D) + R ^ 2 * (1 - c ^ 2) * (1 / D) = 2 * D by
        field_simp [hD.ne']
        linear_combination -hsq)
  refine (((hdist.pow 3).div_const 3).add (hbracket.const_mul (R * c / 2))).congr_deriv ?_
  change (3 : ℝ) * D ^ (3 - 1) * (a / D) / 3 + R * c / 2 * (2 * D) = s * D
  dsimp [a]
  field_simp [hD.ne']
  ring

/-- Radial integration, at a direction strictly between the poles. -/
theorem radialIntegral_eq_sub {R r₁ r₂ c : ℝ} (hR : 0 < R) (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂)
    (hl : -1 < c) (hu : c < 1) :
    (∫ s in r₁..r₂, s * distance R s c) = radialPrimitive R r₂ c - radialPrimitive R r₁ c := by
  apply intervalIntegral.integral_eq_sub_of_hasDerivAt
  · intro s hs
    rw [Set.uIcc_of_le h₁₂] at hs
    exact hasDerivAt_radialPrimitive hR (lt_of_lt_of_le h₁ hs.1) hl hu
  · apply Continuous.intervalIntegrable
    unfold distance
    fun_prop

/-- Interpolation of the two polar values, with the north value at c = 1. -/
def polarInterpolation (Jplus Jminus c : ℝ) : ℝ :=
  (1 + c) / 2 * Jplus + (1 - c) / 2 * Jminus

def edgeKernel (R r₁ r₂ Jplus Jminus c : ℝ) : ℝ :=
  (radialPrimitive R r₂ c - radialPrimitive R r₁ c - polarInterpolation Jplus Jminus c) /
    (1 - c ^ 2)

/-- An exact special-function endpoint for the regularized edge integrand.
The two zero-parameter terms integrate the polar interpolation. -/
def regularizedEndpoint (R r₁ r₂ h Jplus Jminus φ : ℝ) : ℝ :=
  ellipticEndpoint R r₂ h φ - ellipticEndpoint R r₁ h φ -
    Jplus / (1 - h) * legendrePi (characteristicPlus h) 0 (φ / 2) -
    Jminus / (1 + h) * legendrePi (characteristicMinus h) 0 (φ / 2)

private theorem ellipticRoot_zero (t : ℝ) : ellipticRoot 0 t = 1 := by
  simp [ellipticRoot]

/-- The radius-independent correction cancels, and the pole subtractions have exact primitives. -/
theorem hasDerivAt_regularizedEndpoint {R r₁ r₂ h : ℝ}
    (hR : 0 < R) (h₁ : 0 < r₁) (h₂ : 0 < r₂) (hh : 0 < h) (h1 : h < 1)
    (Jplus Jminus φ : ℝ) :
    HasDerivAt (regularizedEndpoint R r₁ r₂ h Jplus Jminus)
      (edgeKernel R r₁ r₂ Jplus Jminus (directionCosine h φ)) φ := by
  let c := directionCosine h φ
  have ⟨hcl, hcu⟩ := dir_bounds hh h1 φ
  have hm : 0 < 1 - h := by linarith
  have hp : 0 < 1 + h := by positivity
  have hcminus : 0 < 1 - c := sub_pos.mpr hcu
  have hcplus : 0 < 1 + c := by linarith
  have hc := one_sub_sq_pos hcl hcu
  have hhalf := (hasDerivAt_id φ).div_const (2 : ℝ)
  have hz : (0 : ℝ) ≤ 0 := le_rfl
  have hn := ((hasDerivAt_legendrePi (characteristicPlus_lt_one hh h1) hz (φ / 2)).comp φ
    hhalf).const_mul (Jplus / (1 - h))
  have ht := ((hasDerivAt_legendrePi (characteristicMinus_lt_one hh h1) hz (φ / 2)).comp φ
    hhalf).const_mul (Jminus / (1 + h))
  have hn' := hn.congr_deriv (show Jplus / (1 - h) *
      (((1 - characteristicPlus h * Real.sin (φ / 2) ^ 2) * ellipticRoot 0 (φ / 2))⁻¹ *
        (1 / 2)) = Jplus / (2 * (1 - c)) by
    rw [factorPlus h1, ellipticRoot_zero, mul_one]
    change Jplus / (1 - h) * (((1 - c) / (1 - h))⁻¹ * (1 / 2)) = _
    rw [inv_div, ← mul_assoc, div_mul_div_cancel₀ hm.ne', ← mul_div_assoc, mul_one, div_div,
      mul_comm (1 - c)])
  have ht' := ht.congr_deriv (show Jminus / (1 + h) *
      (((1 - characteristicMinus h * Real.sin (φ / 2) ^ 2) * ellipticRoot 0 (φ / 2))⁻¹ *
        (1 / 2)) = Jminus / (2 * (1 + c)) by
    rw [factorMinus hh, ellipticRoot_zero, mul_one]
    change Jminus / (1 + h) * (((1 + c) / (1 + h))⁻¹ * (1 / 2)) = _
    rw [inv_div, ← mul_assoc, div_mul_div_cancel₀ hp.ne', ← mul_div_assoc, mul_one, div_div,
      mul_comm (1 + c)])
  have hQ₂ := hasDerivAt_ellipticEndpoint hR h₂ hh h1 φ
  have hQ₁ := hasDerivAt_ellipticEndpoint hR h₁ hh h1 φ
  refine (((hQ₂.sub hQ₁).sub hn').sub ht').congr_deriv ?_
  change (_ + _ - (_ + _) - Jplus / (2 * (1 - c)) - Jminus / (2 * (1 + c))) =
    edgeKernel R r₁ r₂ Jplus Jminus c
  unfold edgeKernel polarInterpolation
  field_simp [hc.ne', hcminus.ne', hcplus.ne']
  ring

/-- An oriented edge integral equals the difference of two explicit endpoint values. -/
theorem edgeIntegral_eq_endpoints {R r₁ r₂ h : ℝ}
    (hR : 0 < R) (h₁ : 0 < r₁) (h₂ : 0 < r₂) (hh : 0 < h) (h1 : h < 1)
    (Jplus Jminus α β : ℝ) :
    (∫ φ in α..β, edgeKernel R r₁ r₂ Jplus Jminus (directionCosine h φ)) =
      regularizedEndpoint R r₁ r₂ h Jplus Jminus β -
        regularizedEndpoint R r₁ r₂ h Jplus Jminus α := by
  apply intervalIntegral.integral_eq_sub_of_hasDerivAt
  · intro φ _
    exact hasDerivAt_regularizedEndpoint hR h₁ h₂ hh h1 Jplus Jminus φ
  · apply Continuous.intervalIntegrable
    have hc : Continuous (directionCosine h) := by
      unfold directionCosine
      fun_prop
    have hH (s : ℝ) (hs : 0 < s) :
        Continuous (fun φ => radialPrimitive R s (directionCosine h φ)) := by
      have hD : Continuous (fun φ => distance R s (directionCosine h φ)) := by
        unfold distance
        fun_prop
      have hlog : Continuous (fun φ =>
          Real.log (s - R * directionCosine h φ + distance R s (directionCosine h φ))) := by
        have harg : Continuous (fun φ =>
            s - R * directionCosine h φ + distance R s (directionCosine h φ)) := by fun_prop
        apply harg.log
        intro φ
        exact (logArgument_pos hR hs (dir_bounds hh h1 φ).1 (dir_bounds hh h1 φ).2).ne'
      unfold radialPrimitive
      dsimp only
      fun_prop
    unfold edgeKernel polarInterpolation
    refine Continuous.div ?_ (by fun_prop) ?_
    · have h₁' := hH r₁ h₁
      have h₂' := hH r₂ h₂
      fun_prop
    · intro φ
      exact (one_sub_sq_pos (dir_bounds hh h1 φ).1 (dir_bounds hh h1 φ).2).ne'

end TSP
