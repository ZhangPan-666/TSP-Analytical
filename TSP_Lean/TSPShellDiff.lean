import TSPSpherical
import Mathlib.Analysis.Calculus.ParametricIntegral

/-! Differentiation of the regularized shell kernel, including the poles. -/

noncomputable section

open MeasureTheory Metric Set
open scoped Topology

namespace TSP

def shellKernelDerivative (R r₁ r₂ c : ℝ) : ℝ :=
  ∫ s in r₁..r₂, fixedRadiusKernelDerivative R s c

private theorem exterior_kernel_neighborhood {R r₁ r₂ c : ℝ}
    (h₁ : 0 < r₁) (h₂R : r₂ < R) (hc : c ≤ 1) :
    ∃ ε > 0, ∀ y ∈ closedBall c ε, ∀ s ∈ Icc r₁ r₂,
      0 < (s - R) ^ 2 + 2 * R * s * (1 - y) ∧ 0 < distance R s y + (R - s * y) := by
  let rad : ℝ × ℝ → ℝ := fun z => (z.2 - R) ^ 2 + 2 * R * z.2 * (1 - z.1)
  let den : ℝ × ℝ → ℝ := fun z => distance R z.2 z.1 + (R - z.2 * z.1)
  have hrad : Continuous rad := by unfold rad; fun_prop
  have hden : Continuous den := by unfold den distance; fun_prop
  let U := {z | 0 < rad z} ∩ {z | 0 < den z}
  have hU : IsOpen U := (isOpen_lt continuous_const hrad).inter (isOpen_lt continuous_const hden)
  have hsubset : {c} ×ˢ Icc r₁ r₂ ⊆ U := by
    rintro ⟨y, s⟩ ⟨rfl, hs⟩
    have hpos := h₁.trans_le hs.1
    have hlt := hs.2.trans_lt h₂R
    exact ⟨exterior_radicand_pos (hpos.trans hlt) hpos hlt hc, exterior_denominator_pos hpos hlt hc⟩
  obtain ⟨A, B, hA, _, hcA, hAB, hUAB⟩ :=
    generalized_tube_lemma isCompact_singleton isCompact_Icc hU hsubset
  obtain ⟨ε, hε, hεA⟩ := Metric.mem_nhds_iff.mp (hA.mem_nhds (hcA (mem_singleton c)))
  refine ⟨ε / 2, by positivity, fun y hy s hs => ?_⟩
  have hyA := hεA (closedBall_subset_ball (by linarith : ε / 2 < ε) hy)
  have hpair : (y, s) ∈ A ×ˢ B := ⟨hyA, hAB hs⟩
  exact hUAB hpair

set_option maxHeartbeats 400000 in
/-- Derivative under the radial integral, including c = ±1 in the exterior branch. -/
theorem hasDerivAt_shellKernel {R r₁ r₂ c : ℝ}
    (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) (h₂R : r₂ < R) (hc : c ≤ 1) :
    HasDerivAt (shellKernel R r₁ r₂) (shellKernelDerivative R r₁ r₂ c) c := by
  obtain ⟨ε, hε, hpos⟩ := exterior_kernel_neighborhood h₁ h₂R hc
  let A := Icc r₁ r₂
  let C := closedBall c ε ×ˢ A
  let μ : Measure ℝ := volume.restrict A
  have hA : IsCompact A := isCompact_Icc
  letI : IsFiniteMeasure μ := ⟨by simpa [μ] using hA.measure_lt_top (μ := volume)⟩
  have hC : IsCompact C := (isCompact_closedBall c ε).prod hA
  have hden : Continuous (fun z : ℝ × ℝ => distance R z.2 z.1 + (R - z.2 * z.1)) := by
    unfold distance
    fun_prop
  have hden_ne (z : ℝ × ℝ) (hz : z ∈ C) :
      distance R z.2 z.1 + (R - z.2 * z.1) ≠ 0 := (hpos z.1 hz.1 z.2 hz.2).2.ne'
  have hD_ne (z : ℝ × ℝ) (hz : z ∈ C) : distance R z.2 z.1 ≠ 0 :=
    (Real.sqrt_pos.2 (hpos z.1 hz.1 z.2 hz.2).1).ne'
  have hK : ContinuousOn (fun z : ℝ × ℝ => fixedRadiusKernel R z.2 z.1) C := by
    unfold fixedRadiusKernel
    exact ContinuousOn.div (by fun_prop) hden.continuousOn hden_ne
  have hKD : ContinuousOn (fun z : ℝ × ℝ => fixedRadiusKernelDerivative R z.2 z.1) C := by
    unfold fixedRadiusKernelDerivative
    have hD : ContinuousOn (fun z : ℝ × ℝ => distance R z.2 z.1) C := by
      unfold distance
      fun_prop
    have hfrac := ContinuousOn.div (f := fun z : ℝ × ℝ => R * z.2)
      (by fun_prop) hD hD_ne
    exact ContinuousOn.div
      ((continuous_snd.pow 3).continuousOn.mul (hfrac.add continuous_snd.continuousOn))
      (hden.continuousOn.pow 2) (fun z hz => pow_ne_zero _ (hden_ne z hz))
  have hslice (y : ℝ) : Continuous (fun s : ℝ => (y, s)) := by fun_prop
  have hmaps (y : ℝ) (hy : y ∈ closedBall c ε) : MapsTo (fun s : ℝ => (y, s)) A C :=
    fun _ hs => ⟨hy, hs⟩
  have hK_i (y : ℝ) (hy : y ∈ closedBall c ε) : Integrable (fun s => fixedRadiusKernel R s y) μ := by
    have hs := hK.comp (hslice y).continuousOn (hmaps y hy)
    exact ContinuousOn.integrableOn_compact (μ := volume) hA hs
  have hcC : c ∈ closedBall c ε := mem_closedBall_self hε.le
  have hK'_i : Integrable (fun s => fixedRadiusKernelDerivative R s c) μ := by
    have hs := hKD.comp (hslice c).continuousOn (hmaps c hcC)
    exact ContinuousOn.integrableOn_compact (μ := volume) hA hs
  obtain ⟨M, hM⟩ := (hC.image_of_continuousOn hKD).isBounded.exists_norm_le
  have hrad : ∀ᵐ s ∂μ, s ∈ A := ae_restrict_mem measurableSet_Icc
  have hmeas : ∀ᶠ y in 𝓝 c, AEStronglyMeasurable (fun s => fixedRadiusKernel R s y) μ := by
    filter_upwards [Metric.ball_mem_nhds c hε] with y hy
    exact (hK_i y (ball_subset_closedBall hy)).aestronglyMeasurable
  have hdombound : ∀ᵐ s ∂μ, ∀ y ∈ ball c ε, ‖fixedRadiusKernelDerivative R s y‖ ≤ M := by
    filter_upwards [hrad] with s hs
    intro y hy
    exact hM _ ⟨(y, s), ⟨ball_subset_closedBall hy, hs⟩, rfl⟩
  have hdiff : ∀ᵐ s ∂μ, ∀ y ∈ ball c ε,
      HasDerivAt (fixedRadiusKernel R s) (fixedRadiusKernelDerivative R s y) y := by
    filter_upwards [hrad] with s hs
    intro y hy
    have hp := hpos y (ball_subset_closedBall hy) s hs
    exact hasDerivAt_fixedRadiusKernel hp.1 hp.2
  have hd := (hasDerivAt_integral_of_dominated_loc_of_deriv_le (μ := μ)
    (F := fun y s => fixedRadiusKernel R s y) (F' := fun y s => fixedRadiusKernelDerivative R s y)
    (bound := fun _ => M) hε hmeas (hK_i c hcC) hK'_i.aestronglyMeasurable hdombound
    (integrable_const M) hdiff).2
  have hinterval (f : ℝ → ℝ) : (∫ s in r₁..r₂, f s) = ∫ s, f s ∂μ := by
    rw [intervalIntegral.integral_of_le h₁₂]
    exact setIntegral_congr_set Ioc_ae_eq_Icc
  have hfun : shellKernel R r₁ r₂ = fun y => ∫ s, fixedRadiusKernel R s y ∂μ := by
    funext y
    exact hinterval _
  rw [hfun, shellKernelDerivative, hinterval]
  exact hd

/-- The curl coefficient of the radial shell, with the exact constant area correction. -/
theorem shellKernel_curl {R r₁ r₂ c : ℝ}
    (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) (h₂R : r₂ < R) (hc : c ≤ 1) :
    2 * c * shellKernel R r₁ r₂ c - (1 - c ^ 2) * shellKernelDerivative R r₁ r₂ c =
      R * radialPotential R r₁ r₂ c - areaCorrection r₁ r₂ := by
  have hr (s : ℝ) (hs : s ∈ uIcc r₁ r₂) : 0 < s ∧ s < R := by
    rw [uIcc_of_le h₁₂] at hs
    exact ⟨h₁.trans_le hs.1, hs.2.trans_lt h₂R⟩
  have hD : Continuous (fun s => distance R s c) := by unfold distance; fun_prop
  have hD_pos (s : ℝ) (hs : s ∈ uIcc r₁ r₂) : 0 < distance R s c :=
    Real.sqrt_pos.2 (exterior_radicand_pos ((hr s hs).1.trans (hr s hs).2)
      (hr s hs).1 (hr s hs).2 hc)
  have hden : Continuous (fun s => distance R s c + (R - s * c)) := by fun_prop
  have hden_ne (s : ℝ) (hs : s ∈ uIcc r₁ r₂) : distance R s c + (R - s * c) ≠ 0 :=
    (exterior_denominator_pos (hr s hs).1 (hr s hs).2 hc).ne'
  have hK : IntervalIntegrable (fun s => fixedRadiusKernel R s c) volume r₁ r₂ := by
    unfold fixedRadiusKernel
    exact (ContinuousOn.div (by fun_prop) hden.continuousOn hden_ne).intervalIntegrable
  have hK' : IntervalIntegrable (fun s => fixedRadiusKernelDerivative R s c) volume r₁ r₂ := by
    unfold fixedRadiusKernelDerivative
    have hfrac := ContinuousOn.div (f := fun s : ℝ => R * s) (by fun_prop) hD.continuousOn
      (fun s hs => (hD_pos s hs).ne')
    exact (ContinuousOn.div
      ((continuous_id.pow 3).continuousOn.mul (hfrac.add continuous_id.continuousOn))
      (hden.continuousOn.pow 2) (fun s hs => pow_ne_zero _ (hden_ne s hs))).intervalIntegrable
  have hP : IntervalIntegrable (fun s => s ^ 2 / distance R s c) volume r₁ r₂ :=
    (ContinuousOn.div (by fun_prop) hD.continuousOn (fun s hs => (hD_pos s hs).ne')).intervalIntegrable
  have hsq : IntervalIntegrable (fun s : ℝ => s ^ 2) volume r₁ r₂ :=
    (continuous_id.pow 2).intervalIntegrable _ _
  have hi : (∫ s in r₁..r₂, s ^ 2) = areaCorrection r₁ r₂ := by
    have h : (∫ s in r₁..r₂, s ^ 2) = r₂ ^ 3 / 3 - r₁ ^ 3 / 3 := by
      apply intervalIntegral.integral_eq_sub_of_hasDerivAt _ hsq
      intro s _
      refine (((hasDerivAt_id s).pow 3).div_const 3).congr_deriv ?_
      simp only [id_eq, Nat.cast_ofNat]
      ring
    rw [h, areaCorrection]
    ring
  rw [shellKernel, shellKernelDerivative, radialPotential, ← hi, ← intervalIntegral.integral_const_mul,
    ← intervalIntegral.integral_const_mul, ← intervalIntegral.integral_const_mul,
    ← intervalIntegral.integral_sub (hK.const_mul _) (hK'.const_mul _),
    ← intervalIntegral.integral_sub (hP.const_mul _) hsq]
  apply intervalIntegral.integral_congr
  intro s hs
  exact fixedRadiusKernel_curl (hr s hs).1 (hr s hs).2 hc

end TSP
