import TSPShellDiff
import Mathlib.Analysis.Calculus.FDeriv.Symmetric
import Mathlib.MeasureTheory.Integral.DivergenceTheorem

/-!
# Green's theorem for a sphere chart

A normalized triangle can be parametrized on the unit square with one collapsed side.
The geometric area and edge changes of variables are explicit premises of the final theorem.
-/

noncomputable section

open MeasureTheory Set
open scoped BigOperators Topology

namespace TSP

abbrev Plane := ℝ × ℝ

abbrev Square : Set Plane := Icc (0, 0) (1, 1)

def sphereAreaDensity (X Y C : Plane → ℝ) (z : Plane) : ℝ :=
  let u := (1, 0)
  let v := (0, 1)
  let dx := fderiv ℝ X z
  let dy := fderiv ℝ Y z
  let dc := fderiv ℝ C z
  X z * (dy u * dc v - dc u * dy v) +
    Y z * (dc u * dx v - dx u * dc v) +
    C z * (dx u * dy v - dy u * dx v)

def shellFormComponent (K : ℝ → ℝ) (X Y C : Plane → ℝ) (d : Plane) (z : Plane) : ℝ :=
  K (C z) * (X z * fderiv ℝ Y z d - Y z * fderiv ℝ X z d)

private theorem sphere_tangent {X Y C : Plane → ℝ} {z d : Plane}
    (hX : DifferentiableAt ℝ X z) (hY : DifferentiableAt ℝ Y z) (hC : DifferentiableAt ℝ C z)
    (hSphere : ∀ᶠ w in 𝓝 z, X w ^ 2 + Y w ^ 2 + C w ^ 2 = 1) :
    X z * fderiv ℝ X z d + Y z * fderiv ℝ Y z d + C z * fderiv ℝ C z d = 0 := by
  have hN := ((hX.hasFDerivAt.mul hX.hasFDerivAt).add
    (hY.hasFDerivAt.mul hY.hasFDerivAt)).add (hC.hasFDerivAt.mul hC.hasFDerivAt)
  have hN' := hN.congr_of_eventuallyEq (f₁ := fun _ => (1 : ℝ)) (by
    filter_upwards [hSphere] with w hw
    exact (by simpa only [pow_two] using hw.symm))
  have hh := congrArg (fun L : Plane →L[ℝ] ℝ => L d)
    (hN'.unique (hasFDerivAt_const (1 : ℝ) z))
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.smul_apply,
    ContinuousLinearMap.zero_apply, smul_eq_mul] at hh
  linarith only [hh]

private theorem shell_form_derivatives {R r₁ r₂ : ℝ} {X Y C : Plane → ℝ} {z : Plane}
    (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) (h₂R : r₂ < R)
    (hX : ContDiffAt ℝ 2 X z) (hY : ContDiffAt ℝ 2 Y z) (hC : ContDiffAt ℝ 2 C z)
    (hSphere : ∀ᶠ w in 𝓝 z, X w ^ 2 + Y w ^ 2 + C w ^ 2 = 1) :
    let K := shellKernel R r₁ r₂
    let P := shellFormComponent K X Y C (1, 0)
    let Q := shellFormComponent K X Y C (0, 1)
    DifferentiableAt ℝ P z ∧ DifferentiableAt ℝ Q z ∧
      fderiv ℝ Q z (1, 0) - fderiv ℝ P z (0, 1) =
        (R * radialPotential R r₁ r₂ (C z) - areaCorrection r₁ r₂) *
          sphereAreaDensity X Y C z := by
  let u : Plane := (1, 0)
  let v : Plane := (0, 1)
  have hSphere₀ := hSphere.self_of_nhds
  have hc : C z ≤ 1 := by nlinarith [sq_nonneg (X z), sq_nonneg (Y z)]
  have hX₁ := hX.differentiableAt (by norm_num)
  have hY₁ := hY.differentiableAt (by norm_num)
  have hC₁ := hC.differentiableAt (by norm_num)
  have hDX := (hX.fderiv_right (m := 1) (by norm_num)).differentiableAt le_rfl
  have hDY := (hY.fderiv_right (m := 1) (by norm_num)).differentiableAt le_rfl
  have hK := (hasDerivAt_shellKernel h₁ h₁₂ h₂R hc).hasFDerivAt.comp z hC₁.hasFDerivAt
  have hPart (d : Plane) :=
    (hK.mul ((hX₁.hasFDerivAt.mul (hDY.hasFDerivAt.clm_apply (hasFDerivAt_const d z))).sub
      (hY₁.hasFDerivAt.mul (hDX.hasFDerivAt.clm_apply (hasFDerivAt_const d z)))))
  have hP := hPart u
  have hQ := hPart v
  change HasFDerivAt (shellFormComponent (shellKernel R r₁ r₂) X Y C u) _ z at hP
  change HasFDerivAt (shellFormComponent (shellKernel R r₁ r₂) X Y C v) _ z at hQ
  refine ⟨hP.differentiableAt, hQ.differentiableAt, ?_⟩
  let dx := fderiv ℝ X z
  let dy := fderiv ℝ Y z
  let dc := fderiv ℝ C z
  let k := shellKernel R r₁ r₂ (C z)
  let k' := shellKernelDerivative R r₁ r₂ (C z)
  have hXmix := (hX.isSymmSndFDerivAt (by simp [minSmoothness])).eq u v
  have hYmix := (hY.isSymmSndFDerivAt (by simp [minSmoothness])).eq u v
  change fderiv ℝ (shellFormComponent (shellKernel R r₁ r₂) X Y C v) z u -
    fderiv ℝ (shellFormComponent (shellKernel R r₁ r₂) X Y C u) z v = _
  calc
    _ = k' * (dc u * (X z * dy v - Y z * dx v) - dc v * (X z * dy u - Y z * dx u)) +
        2 * k * (dx u * dy v - dy u * dx v) := by
      rw [hQ.fderiv, hP.fderiv]
      simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.sub_apply,
        ContinuousLinearMap.smul_apply, ContinuousLinearMap.comp_apply,
        ContinuousLinearMap.smulRight_apply, ContinuousLinearMap.one_apply,
        ContinuousLinearMap.zero_apply, ContinuousLinearMap.flip_apply, zero_add, map_zero,
        smul_eq_mul, Function.comp_apply]
      rw [hXmix, hYmix]
      dsimp only [dx, dy, dc, k, k']
      ring
    _ = (2 * C z * k - (1 - C z ^ 2) * k') * sphereAreaDensity X Y C z :=
      sphericalCurl_coordinate_identity _ _ _ _ _ _ _ _ _ _ _ hSphere₀
        (sphere_tangent hX₁ hY₁ hC₁ hSphere) (sphere_tangent hX₁ hY₁ hC₁ hSphere)
    _ = _ := by
      dsimp only [k, k']
      rw [shellKernel_curl h₁ h₁₂ h₂R hc]

/-- The positively oriented integral along the four sides of the unit square. -/
def squareBoundary (P Q : Plane → ℝ) : ℝ :=
  (∫ t in (0 : ℝ)..1, P (t, 0)) + (∫ t in (0 : ℝ)..1, Q (1, t)) -
    (∫ t in (0 : ℝ)..1, P (t, 1)) - (∫ t in (0 : ℝ)..1, Q (0, t))

/-- Green's theorem for the pulled-back shell form.
The coordinate chart takes values on the unit sphere near the square. -/
theorem square_shell_stokes {R r₁ r₂ : ℝ} {X Y C : Plane → ℝ}
    (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) (h₂R : r₂ < R)
    (hChart : ∀ z ∈ Square, ContDiffAt ℝ 2 X z ∧ ContDiffAt ℝ 2 Y z ∧
      ContDiffAt ℝ 2 C z ∧ (∀ᶠ w in 𝓝 z, X w ^ 2 + Y w ^ 2 + C w ^ 2 = 1))
    (hInt : IntegrableOn (fun z =>
      (R * radialPotential R r₁ r₂ (C z) - areaCorrection r₁ r₂) * sphereAreaDensity X Y C z)
      Square) :
    (∫ z in Square,
      (R * radialPotential R r₁ r₂ (C z) - areaCorrection r₁ r₂) * sphereAreaDensity X Y C z) =
      squareBoundary (shellFormComponent (shellKernel R r₁ r₂) X Y C (1, 0))
        (shellFormComponent (shellKernel R r₁ r₂) X Y C (0, 1)) := by
  let P := shellFormComponent (shellKernel R r₁ r₂) X Y C (1, 0)
  let Q := shellFormComponent (shellKernel R r₁ r₂) X Y C (0, 1)
  have hAt (z : Plane) (hz : z ∈ Square) :=
    shell_form_derivatives h₁ h₁₂ h₂R (hChart z hz).1 (hChart z hz).2.1
      (hChart z hz).2.2.1 (hChart z hz).2.2.2
  have hCurl : EqOn (fun z =>
      (R * radialPotential R r₁ r₂ (C z) - areaCorrection r₁ r₂) * sphereAreaDensity X Y C z)
      (fun z => fderiv ℝ Q z (1, 0) + (-fderiv ℝ P z) (0, 1)) Square := by
    intro z hz
    change _ = fderiv ℝ Q z (1, 0) + (-fderiv ℝ P z) (0, 1)
    rw [ContinuousLinearMap.neg_apply, ← sub_eq_add_neg]
    exact (hAt z hz).2.2.symm
  have hsubset : Ioo (0 : ℝ) 1 ×ˢ Ioo (0 : ℝ) 1 \ ∅ ⊆ Square := fun z hz =>
    ⟨⟨hz.1.1.1.le, hz.1.2.1.le⟩, ⟨hz.1.1.2.le, hz.1.2.2.le⟩⟩
  have hdiv := integral_divergence_prod_Icc_of_hasFDerivWithinAt_off_countable_of_le
      Q (-P) (fderiv ℝ Q) (fun z => -fderiv ℝ P z) (0, 0) (1, 1)
      (by exact ⟨zero_le_one, zero_le_one⟩) ∅ countable_empty
      (fun z hz => (hAt z hz).2.1.continuousAt.continuousWithinAt)
      (fun z hz => (hAt z hz).1.continuousAt.neg.continuousWithinAt)
      (fun z hz => (hAt z (hsubset hz)).2.1.hasFDerivAt)
      (fun z hz => (hAt z (hsubset hz)).1.hasFDerivAt.neg)
      (hInt.congr_fun hCurl measurableSet_Icc)
  rw [setIntegral_congr_fun measurableSet_Icc hCurl, hdiv]
  unfold squareBoundary
  simp only [Pi.neg_apply, intervalIntegral.integral_neg, Prod.fst, Prod.snd]
  ring

/-- A chart side with a differentiable phase has the same integral as the oriented edge.
The interval may be reversed, as on the top side of the positively oriented square. -/
theorem chartSide_eq_edgeIntegral {R r₁ r₂ : ℝ}
    (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) (h₂R : r₂ < R) (e : GenericEdge)
    {a b : ℝ} {trace phase speed : ℝ → ℝ}
    (hStart : phase a = e.start) (hStop : phase b = e.stop)
    (hDiff : ∀ t ∈ uIcc a b, HasDerivAt phase (speed t) t)
    (hTrace : ∀ t ∈ uIcc a b, trace t =
      e.κ * shellKernel R r₁ r₂ (directionCosine e.h (phase t)) * speed t)
    (hInt : IntervalIntegrable trace volume a b) :
    (∫ t in a..b, trace t) = boundaryIntegral R r₁ r₂ e := by
  let P := regularizedEndpoint R r₁ r₂ e.h (northPolarMoment R r₁ r₂) (southPolarMoment R r₁ r₂)
  have h₂ := h₁.trans_le h₁₂
  have hR := h₂.trans h₂R
  have hPrimitive :
      (∫ t in a..b, trace t) = e.κ * P (phase b) - e.κ * P (phase a) := by
    apply intervalIntegral.integral_eq_sub_of_hasDerivAt _ hInt
    intro t ht
    have hp := hasDerivAt_regularizedEndpoint hR h₁ h₂ e.h_pos e.h_lt_one
      (northPolarMoment R r₁ r₂) (southPolarMoment R r₁ r₂) (phase t)
    refine ((hp.comp t (hDiff t ht)).const_mul e.κ).congr_deriv ?_
    have ⟨hl, hu⟩ := dir_bounds e.h_pos e.h_lt_one (phase t)
    rw [hTrace t ht, shellKernel_eq_edgeKernel h₁ h₁₂ h₂R hl hu]
    ring
  rw [hPrimitive, hStart, hStop, boundaryIntegral_eq_endpoints h₁ h₁₂ h₂R]
  change _ = e.κ * (P e.stop - P e.start)
  ring

/-- End-to-end endpoint evaluation from a sphere chart with verified coordinate geometry.
The two change-of-variables premises can be supplied by a positively oriented triangle chart. -/
theorem volumeIntegral_eq_endpoints_of_chart {T : Set Direction} (hT : IsCompact T)
    {r₁ r₂ : ℝ} (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) {x : Space} (hx : r₂ < ‖x‖)
    {edges : Fin 3 → GenericEdge} {X Y C : Plane → ℝ} (u : Plane → Direction)
    (hChart : ∀ z ∈ Square, ContDiffAt ℝ 2 X z ∧ ContDiffAt ℝ 2 Y z ∧
      ContDiffAt ℝ 2 C z ∧ (∀ᶠ w in 𝓝 z, X w ^ 2 + Y w ^ 2 + C w ^ 2 = 1))
    (hCosine : ∀ z ∈ Square, observationCosine x (u z) = C z)
    (hArea : ∀ f : Direction → ℝ, IntegrableOn f T angularMeasure →
      IntegrableOn (fun z => f (u z) * sphereAreaDensity X Y C z) Square ∧
        (∫ y in T, f y ∂angularMeasure) = ∫ z in Square, f (u z) * sphereAreaDensity X Y C z)
    (hEdges :
      squareBoundary (shellFormComponent (shellKernel ‖x‖ r₁ r₂) X Y C (1, 0))
          (shellFormComponent (shellKernel ‖x‖ r₁ r₂) X Y C (0, 1)) =
        ∑ e : Fin 3, boundaryIntegral ‖x‖ r₁ r₂ (edges e)) :
    (∫ y in prism T r₁ r₂, ‖x - y‖⁻¹) =
      (angularMeasure.real T * areaCorrection r₁ r₂ +
        ∑ e : Fin 3, boundaryEndpoints ‖x‖ r₁ r₂ (edges e)) / ‖x‖ := by
  apply volumeIntegral_eq_endpoints_of_surfaceFlux hT h₁ h₁₂ hx
  have hI := integrableOn_radialPotential hT h₁ h₁₂ hx
  let f : Direction → ℝ := fun y =>
    ‖x‖ * radialPotential ‖x‖ r₁ r₂ (observationCosine x y) - areaCorrection r₁ r₂
  have hf : IntegrableOn f T angularMeasure := (hI.const_mul _).sub (integrable_const _)
  obtain ⟨hInt, hEq⟩ := hArea f hf
  have hcongr : EqOn (fun z => f (u z) * sphereAreaDensity X Y C z) (fun z =>
      (‖x‖ * radialPotential ‖x‖ r₁ r₂ (C z) - areaCorrection r₁ r₂) * sphereAreaDensity X Y C z)
      Square := by
    intro z hz
    simp only [f, hCosine z hz]
  unfold SurfaceFluxIdentity
  rw [hEq, setIntegral_congr_fun measurableSet_Icc hcongr, square_shell_stokes h₁ h₁₂ hx hChart
    (hInt.congr_fun hcongr measurableSet_Icc), hEdges]

end TSP
