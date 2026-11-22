import TSPVolume

/-!
# The exterior shell kernel and the endpoint formula

The fixed-radius kernel uses a positive denominator at both poles when R > r₂.
The local curl calculation applies even when an observation pole lies inside the angular region.
-/

open MeasureTheory Set
open scoped BigOperators Topology

noncomputable section

namespace TSP

def fixedRadiusKernel (R s c : ℝ) : ℝ :=
  s ^ 3 / (distance R s c + (R - s * c))

def fixedRadiusKernelDerivative (R s c : ℝ) : ℝ :=
  s ^ 3 * (R * s / distance R s c + s) / (distance R s c + (R - s * c)) ^ 2

def shellKernel (R r₁ r₂ c : ℝ) : ℝ :=
  ∫ s in r₁..r₂, fixedRadiusKernel R s c

def northPolarMoment (R r₁ r₂ : ℝ) : ℝ :=
  R * (r₂ ^ 2 - r₁ ^ 2) / 2 - (r₂ ^ 3 - r₁ ^ 3) / 3

def southPolarMoment (R r₁ r₂ : ℝ) : ℝ :=
  R * (r₂ ^ 2 - r₁ ^ 2) / 2 + (r₂ ^ 3 - r₁ ^ 3) / 3

theorem exterior_radicand_pos {R s c : ℝ}
    (hR : 0 < R) (hs : 0 < s) (hsR : s < R) (hc : c ≤ 1) :
    0 < (s - R) ^ 2 + 2 * R * s * (1 - c) := by
  have hnon : 0 ≤ 2 * R * s * (1 - c) :=
    mul_nonneg (by positivity) (sub_nonneg.mpr hc)
  exact add_pos_of_pos_of_nonneg (sq_pos_of_ne_zero (sub_ne_zero.mpr hsR.ne)) hnon

theorem exterior_denominator_pos {R s c : ℝ} (hs : 0 < s) (hsR : s < R)
    (hc : c ≤ 1) : 0 < distance R s c + (R - s * c) := by
  have h : s * c ≤ s := by simpa using mul_le_mul_of_nonneg_left hc hs.le
  exact add_pos_of_nonneg_of_pos (Real.sqrt_nonneg _) (by linarith)

private theorem fixedRadiusKernel_mul {R s c : ℝ}
    (hq : 0 ≤ (s - R) ^ 2 + 2 * R * s * (1 - c))
    (hne : distance R s c + (R - s * c) ≠ 0) :
    (1 - c ^ 2) * fixedRadiusKernel R s c =
      s * distance R s c - s * (R - s * c) := by
  have hd : distance R s c ^ 2 = (s - R) ^ 2 + 2 * R * s * (1 - c) :=
    Real.sq_sqrt hq
  unfold fixedRadiusKernel
  field_simp [hne]
  linear_combination -s * hd

private theorem hasDerivAt_distance_direction {R s c : ℝ}
    (hq : 0 < (s - R) ^ 2 + 2 * R * s * (1 - c)) :
    HasDerivAt (distance R s) (-(R * s) / distance R s c) c := by
  have hD : 0 < distance R s c := Real.sqrt_pos.2 hq
  have hq' := (((hasDerivAt_id c).const_sub 1).const_mul (2 * R * s)).const_add ((s - R) ^ 2)
  refine ((Real.hasDerivAt_sqrt hq.ne').comp c hq').congr_deriv ?_
  change 1 / (2 * distance R s c) * (2 * R * s * (-1)) = _
  field_simp [hD.ne']
  ring

theorem hasDerivAt_fixedRadiusKernel {R s c : ℝ}
    (hq : 0 < (s - R) ^ 2 + 2 * R * s * (1 - c))
    (hden : 0 < distance R s c + (R - s * c)) :
    HasDerivAt (fixedRadiusKernel R s) (fixedRadiusKernelDerivative R s c) c := by
  have hd := hasDerivAt_distance_direction hq
  have ha := ((hasDerivAt_id c).const_mul s).const_sub R
  refine (((hd.add ha).inv hden.ne').const_mul (s ^ 3)).congr_deriv ?_
  simp only [fixedRadiusKernelDerivative, id_eq, neg_div, mul_one, neg_add_rev, neg_neg,
    mul_div_assoc, add_comm]

/-- Local scalar coefficient of the regularized exterior derivative on the unit sphere. -/
theorem fixedRadiusKernel_curl {R s c : ℝ} (hs : 0 < s) (hsR : s < R) (hc : c ≤ 1) :
    2 * c * fixedRadiusKernel R s c - (1 - c ^ 2) * fixedRadiusKernelDerivative R s c =
      R * (s ^ 2 / distance R s c) - s ^ 2 := by
  have hq := exterior_radicand_pos (hs.trans hsR) hs hsR hc
  have hden := exterior_denominator_pos hs hsR hc
  have hd := hasDerivAt_distance_direction hq
  have ha := ((hasDerivAt_id c).const_mul s).const_sub R
  have hk := hasDerivAt_fixedRadiusKernel hq hden
  have hproduct := (((hasDerivAt_id c).pow 2).const_sub 1).mul hk
  have hradContinuous : Continuous (fun z : ℝ => (s - R) ^ 2 + 2 * R * s * (1 - z)) := by
    fun_prop
  have hdenContinuous : Continuous (fun z : ℝ => distance R s z + (R - s * z)) := by
    unfold distance
    fun_prop
  have hradNear : ∀ᶠ z in 𝓝 c, 0 < (s - R) ^ 2 + 2 * R * s * (1 - z) :=
    (isOpen_lt continuous_const hradContinuous).mem_nhds hq
  have hdenNear : ∀ᶠ z in 𝓝 c, 0 < distance R s z + (R - s * z) :=
    (isOpen_lt continuous_const hdenContinuous).mem_nhds hden
  have hlocal := hproduct.congr_of_eventuallyEq (f₁ :=
      fun z => s * distance R s z - s * (R - s * z)) (by
    filter_upwards [hradNear, hdenNear] with z hz hzne
    exact (fixedRadiusKernel_mul hz.le hzne.ne').symm)
  have href := (hd.const_mul s).sub (ha.const_mul s)
  have hEq := hlocal.deriv.symm.trans href.deriv
  simp only [id_eq, Nat.cast_ofNat, show (2 - 1 : ℕ) = 1 by decide, pow_one, mul_one] at hEq
  simp only [div_eq_mul_inv] at hEq ⊢
  linear_combination -hEq

/-- Coordinate identity for the curl of K(c)(X dY - Y dX).
The two vectors are tangent to the unit sphere at (X,Y,c). -/
theorem sphericalCurl_coordinate_identity (X Y c Xp Yp cp Xq Yq cq k k' : ℝ)
    (hunit : X ^ 2 + Y ^ 2 + c ^ 2 = 1)
    (hp : X * Xp + Y * Yp + c * cp = 0)
    (hq : X * Xq + Y * Yq + c * cq = 0) :
    let det := X * (Yp * cq - cp * Yq) + Y * (cp * Xq - Xp * cq) +
      c * (Xp * Yq - Yp * Xq)
    k' * (cp * (X * Yq - Y * Xq) - cq * (X * Yp - Y * Xp)) +
      2 * k * (Xp * Yq - Yp * Xq) =
    (2 * c * k - (1 - c ^ 2) * k') * det := by
  dsimp only
  linear_combination
    (2 * k + c * k') * ((X * Yq - Y * Xq) * hp - (X * Yp - Y * Xp) * hq -
      (Xp * Yq - Yp * Xq) * hunit)

private theorem integral_linearized_radius (R r₁ r₂ c : ℝ) :
    (∫ s in r₁..r₂, s * (R - s * c)) =
      R * (r₂ ^ 2 - r₁ ^ 2) / 2 - c * (r₂ ^ 3 - r₁ ^ 3) / 3 := by
  have h :
      (∫ s in r₁..r₂, s * (R - s * c)) =
        (R * r₂ ^ 2 / 2 - c * r₂ ^ 3 / 3) - (R * r₁ ^ 2 / 2 - c * r₁ ^ 3 / 3) := by
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt
    · intro s _
      have h₂ := (((hasDerivAt_id s).pow 2).const_mul R).div_const 2
      have h₃ := (((hasDerivAt_id s).pow 3).const_mul c).div_const 3
      refine (h₂.sub h₃).congr_deriv ?_
      simp only [id_eq, Nat.cast_ofNat]
      ring
    · exact Continuous.intervalIntegrable (by fun_prop) _ _
  rw [h]
  ring

private theorem linearizedRadius_eq_interpolation (R r₁ r₂ c : ℝ) :
    R * (r₂ ^ 2 - r₁ ^ 2) / 2 - c * (r₂ ^ 3 - r₁ ^ 3) / 3 =
      polarInterpolation (northPolarMoment R r₁ r₂) (southPolarMoment R r₁ r₂) c := by
  unfold polarInterpolation northPolarMoment southPolarMoment
  ring

/-- The polar values of the radial moment, as finite polynomials in the exterior branch. -/
theorem polarMoments_eq_integrals {R r₁ r₂ : ℝ} (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂)
    (h₂R : r₂ < R) :
    (∫ s in r₁..r₂, s * distance R s 1) = northPolarMoment R r₁ r₂ ∧
    (∫ s in r₁..r₂, s * distance R s (-1)) = southPolarMoment R r₁ r₂ := by
  have hR : 0 < R := h₁.trans_le h₁₂ |>.trans h₂R
  constructor
  · calc
      _ = ∫ s in r₁..r₂, s * (R - s * 1) := by
        apply intervalIntegral.integral_congr
        intro s hs
        rw [uIcc_of_le h₁₂] at hs
        have hsq : (s - R) ^ 2 + 2 * R * s * (1 - 1) = (R - s) ^ 2 := by ring
        change s * distance R s 1 = s * (R - s * 1)
        rw [distance, hsq, Real.sqrt_sq (by linarith [hs.2]), mul_one]
      _ = _ := by simpa [northPolarMoment] using integral_linearized_radius R r₁ r₂ 1
  · calc
      _ = ∫ s in r₁..r₂, s * (R - s * (-1)) := by
        apply intervalIntegral.integral_congr
        intro s hs
        rw [uIcc_of_le h₁₂] at hs
        have hs0 : 0 < s := h₁.trans_le hs.1
        have hsq : (s - R) ^ 2 + 2 * R * s * (1 - (-1)) = (R + s) ^ 2 := by ring
        change s * distance R s (-1) = s * (R - s * (-1))
        rw [distance, hsq, Real.sqrt_sq (by positivity : 0 ≤ R + s)]
        ring
      _ = _ := by rw [integral_linearized_radius, southPolarMoment]; ring

/-- The nonsingular sphere kernel restricts to the regularized edge integrand. -/
theorem shellKernel_eq_edgeKernel {R r₁ r₂ c : ℝ}
    (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) (h₂R : r₂ < R) (hl : -1 < c) (hu : c < 1) :
    shellKernel R r₁ r₂ c =
      edgeKernel R r₁ r₂ (northPolarMoment R r₁ r₂) (southPolarMoment R r₁ r₂) c := by
  have hR : 0 < R := h₁.trans_le h₁₂ |>.trans h₂R
  have hC : 0 < 1 - c ^ 2 := one_sub_sq_pos hl hu
  have hcont₁ : Continuous (fun s => s * distance R s c) := by
    unfold distance
    fun_prop
  have hcont₂ : Continuous (fun s => s * (R - s * c)) := by fun_prop
  unfold shellKernel edgeKernel
  apply (eq_div_iff hC.ne').mpr
  rw [mul_comm, ← intervalIntegral.integral_const_mul,
    ← radialIntegral_eq_sub hR h₁ h₁₂ hl hu, ← linearizedRadius_eq_interpolation,
    ← integral_linearized_radius, ← intervalIntegral.integral_sub
      (hcont₁.intervalIntegrable _ _) (hcont₂.intervalIntegrable _ _)]
  apply intervalIntegral.integral_congr
  intro s hs
  rw [uIcc_of_le h₁₂] at hs
  have hs0 : 0 < s := h₁.trans_le hs.1
  have hsR : s < R := hs.2.trans_lt h₂R
  exact fixedRadiusKernel_mul (exterior_radicand_pos hR hs0 hsR hu.le).le
    (exterior_denominator_pos hs0 hsR hu.le).ne'


/-- A generic, oriented, great-circle edge in the observer's angular coordinates. -/
structure GenericEdge where
  h : ℝ
  start : ℝ
  stop : ℝ
  κ : ℝ
  h_pos : 0 < h
  h_lt_one : h < 1

def boundaryIntegral (R r₁ r₂ : ℝ) (e : GenericEdge) : ℝ :=
  e.κ * ∫ φ in e.start..e.stop, shellKernel R r₁ r₂ (directionCosine e.h φ)

def boundaryEndpoints (R r₁ r₂ : ℝ) (e : GenericEdge) : ℝ :=
  let P := regularizedEndpoint R r₁ r₂ e.h (northPolarMoment R r₁ r₂) (southPolarMoment R r₁ r₂)
  e.κ * (P e.stop - P e.start)

theorem boundaryIntegral_eq_endpoints {R r₁ r₂ : ℝ}
    (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) (h₂R : r₂ < R) (e : GenericEdge) :
    boundaryIntegral R r₁ r₂ e = boundaryEndpoints R r₁ r₂ e := by
  have h₂ : 0 < r₂ := h₁.trans_le h₁₂
  have hR : 0 < R := h₂.trans h₂R
  unfold boundaryIntegral boundaryEndpoints
  congr 1
  calc
    _ = ∫ φ in e.start..e.stop,
          edgeKernel R r₁ r₂ (northPolarMoment R r₁ r₂) (southPolarMoment R r₁ r₂)
            (directionCosine e.h φ) := by
      apply intervalIntegral.integral_congr
      intro φ _
      have ⟨hl, hu⟩ := dir_bounds e.h_pos e.h_lt_one φ
      exact shellKernel_eq_edgeKernel h₁ h₁₂ h₂R hl hu
    _ = _ := edgeIntegral_eq_endpoints hR h₁ h₂ e.h_pos e.h_lt_one _ _ _ _


/-- A Cartesian component of a unit-speed great-circle parametrization. -/
def greatCircleComponent (a b t : ℝ) : ℝ := a * Real.cos t + b * Real.sin t

private theorem deriv_greatCircleComponent (a b t : ℝ) :
    deriv (greatCircleComponent a b) t = -(a * Real.sin t) + b * Real.cos t := by
  simpa only [mul_neg] using
    (((Real.hasDerivAt_cos t).const_mul a).add ((Real.hasDerivAt_sin t).const_mul b)).deriv

theorem greatCircle_longitude_coefficient (ax ay bx bY t : ℝ) :
    greatCircleComponent ax bx t * deriv (greatCircleComponent ay bY) t -
      greatCircleComponent ay bY t * deriv (greatCircleComponent ax bx) t =
    ax * bY - ay * bx := by
  simp only [deriv_greatCircleComponent, greatCircleComponent]
  linear_combination (ax * bY - ay * bx) * Real.sin_sq_add_cos_sq t

/-- Reduction for one edge written in its own great-circle parameter.
The supplied coefficients specify its horizontal components and its phase. -/
theorem greatCircleIntegral_eq_endpoints {R r₁ r₂ h : ℝ}
    (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) (h₂R : r₂ < R) (hh : 0 < h) (h1 : h < 1)
    (ax ay bx bY δ α β : ℝ) :
    let X := greatCircleComponent ax bx
    let Y := greatCircleComponent ay bY
    let C := greatCircleComponent (h * Real.cos δ) (h * Real.sin δ)
    let P := regularizedEndpoint R r₁ r₂ h (northPolarMoment R r₁ r₂) (southPolarMoment R r₁ r₂)
    (∫ t in α..β, shellKernel R r₁ r₂ (C t) * (X t * deriv Y t - Y t * deriv X t)) =
      (ax * bY - ay * bx) * (P (β - δ) - P (α - δ)) := by
  let e : GenericEdge := ⟨h, α - δ, β - δ, ax * bY - ay * bx, hh, h1⟩
  have hPhase (t : ℝ) :
      greatCircleComponent (h * Real.cos δ) (h * Real.sin δ) t = directionCosine h (t - δ) := by
    unfold greatCircleComponent directionCosine
    rw [Real.cos_sub]
    ring
  dsimp only
  simp only [greatCircle_longitude_coefficient, hPhase]
  have h :
      (∫ t in α..β, shellKernel R r₁ r₂ (directionCosine h (t - δ)) * (ax * bY - ay * bx)) =
        boundaryIntegral R r₁ r₂ e := by
    rw [intervalIntegral.integral_mul_const]
    change (∫ t in α..β, (fun φ => shellKernel R r₁ r₂ (directionCosine h φ)) (t - δ)) *
      (ax * bY - ay * bx) = _
    rw [intervalIntegral.integral_comp_sub_right (a := α) (b := β)
      (fun φ => shellKernel R r₁ r₂ (directionCosine h φ)) δ]
    unfold boundaryIntegral e
    exact mul_comm _ _
  exact h.trans (boundaryIntegral_eq_endpoints h₁ h₁₂ h₂R e)

def areaCorrection (r₁ r₂ : ℝ) : ℝ := (r₂ ^ 3 - r₁ ^ 3) / 3

theorem areaCorrection_eq_polarDifference (R r₁ r₂ : ℝ) :
    areaCorrection r₁ r₂ = (southPolarMoment R r₁ r₂ - northPolarMoment R r₁ r₂) / 2 := by
  unfold areaCorrection northPolarMoment southPolarMoment
  ring

/-- The geometric Stokes identity for an angular region and three oriented edges.
Instantiating this premise for a chart from three vertices is separate from the shell and
special-function proof. -/
def SurfaceFluxIdentity (T : Set Direction) (r₁ r₂ : ℝ) (x : Space)
    (edges : Fin 3 → GenericEdge) : Prop :=
  (∫ u in T,
    (‖x‖ * radialPotential ‖x‖ r₁ r₂ (observationCosine x u) - areaCorrection r₁ r₂)
      ∂angularMeasure) =
    ∑ e : Fin 3, boundaryIntegral ‖x‖ r₁ r₂ (edges e)

/-- The Newton volume integral equals the polar area term plus six angular endpoint
evaluations, given the oriented surface-to-boundary identity. Multiplying by Gρ
gives the potential under the positive Newton-kernel convention of the manuscript. -/
theorem volumeIntegral_eq_endpoints_of_surfaceFlux {T : Set Direction} (hT : IsCompact T)
    {r₁ r₂ : ℝ} (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) {x : Space} (hx : r₂ < ‖x‖)
    {edges : Fin 3 → GenericEdge} (hStokes : SurfaceFluxIdentity T r₁ r₂ x edges) :
    (∫ y in prism T r₁ r₂, ‖x - y‖⁻¹) =
      (angularMeasure.real T * areaCorrection r₁ r₂ +
        ∑ e : Fin 3, boundaryEndpoints ‖x‖ r₁ r₂ (edges e)) / ‖x‖ := by
  have hR : 0 < ‖x‖ := (h₁.trans_le h₁₂).trans hx
  have hI := integrableOn_radialPotential hT h₁ h₁₂ hx
  have hsum :
      (∑ e : Fin 3, boundaryIntegral ‖x‖ r₁ r₂ (edges e)) =
        ∑ e : Fin 3, boundaryEndpoints ‖x‖ r₁ r₂ (edges e) := by
    apply Finset.sum_congr rfl
    intro e _
    exact boundaryIntegral_eq_endpoints h₁ h₁₂ hx (edges e)
  unfold SurfaceFluxIdentity at hStokes
  rw [integral_sub (hI.const_mul _) (integrable_const _), integral_const_mul, setIntegral_const,
    smul_eq_mul, hsum] at hStokes
  rw [volumeIntegral_eq_angularIntegral hT h₁ h₁₂ hx]
  apply (eq_div_iff hR.ne').mpr
  linarith

end TSP
