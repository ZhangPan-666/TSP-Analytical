import TSPEndpoint
import Mathlib.MeasureTheory.Constructions.HaarToSphere
import Mathlib.MeasureTheory.Measure.Haar.InnerProductSpace

/-!
# Polar reduction of the three-dimensional Newton potential

The shell is cut out by an angular region on the unit sphere and a radial interval.
Mathlib's polar disintegration of Euclidean volume supplies the factor s².
-/

noncomputable section

open MeasureTheory Metric Set
open scoped RealInnerProductSpace

namespace TSP

abbrev Space := EuclideanSpace ℝ (Fin 3)

abbrev Direction := sphere (0 : Space) 1

abbrev Radius := Ioi (0 : ℝ)

def angularMeasure : Measure Direction := (volume : Measure Space).toSphere

instance : IsFiniteMeasure angularMeasure := by
  unfold angularMeasure
  infer_instance

def radialBand (r₁ r₂ : ℝ) : Set Radius := Subtype.val ⁻¹' Icc r₁ r₂

def prism (T : Set Direction) (r₁ r₂ : ℝ) : Set Space :=
  (fun p : Direction × Radius => (p.2 : ℝ) • (p.1 : Space)) '' (T ×ˢ radialBand r₁ r₂)

def observationCosine (x : Space) (u : Direction) : ℝ := ⟪x, (u : Space)⟫ / ‖x‖

def radialPotential (R r₁ r₂ c : ℝ) : ℝ :=
  ∫ s in r₁..r₂, s ^ 2 / distance R s c

private theorem direction_norm (u : Direction) : ‖(u : Space)‖ = 1 :=
  mem_sphere_zero_iff_norm.mp u.property

theorem observationCosine_bounds {x : Space} (hx : x ≠ 0) (u : Direction) :
    observationCosine x u ∈ Icc (-1) 1 := by
  have hR : 0 < ‖x‖ := norm_pos_iff.mpr hx
  have h := abs_real_inner_le_norm x (u : Space)
  rw [direction_norm, mul_one] at h
  apply abs_le.mp
  unfold observationCosine
  rw [abs_div, abs_of_pos hR, div_le_iff₀ hR, one_mul]
  exact h

private theorem norm_sub_smul_direction {x : Space} (hx : x ≠ 0)
    (u : Direction) {s : ℝ} (hs : 0 ≤ s) :
    ‖x - s • (u : Space)‖ = distance ‖x‖ s (observationCosine x u) := by
  have hR : 0 < ‖x‖ := norm_pos_iff.mpr hx
  have hc : ‖x‖ * observationCosine x u = ⟪x, (u : Space)⟫ := by
    unfold observationCosine
    field_simp [hR.ne']
  unfold distance
  rw [← Real.sqrt_sq (norm_nonneg (x - s • (u : Space)))]
  congr 1
  rw [norm_sub_sq_real, real_inner_smul_right, norm_smul, Real.norm_eq_abs, abs_of_nonneg hs,
    direction_norm, mul_one, ← hc]
  ring

private theorem radialBand_image {r₁ r₂ : ℝ} (h₁ : 0 < r₁) :
    (Subtype.val : Radius → ℝ) '' radialBand r₁ r₂ = Icc r₁ r₂ := by
  refine image_preimage_eq_of_subset ?_
  rw [Subtype.range_coe]
  intro s hs
  exact h₁.trans_le hs.1

private theorem setIntegral_positiveRadius {B : Set ℝ} (hB : B ⊆ Ioi (0 : ℝ))
    (f : ℝ → ℝ) :
    (∫ s in (Subtype.val : Radius → ℝ) ⁻¹' B, f s
      ∂(volume : Measure ℝ).comap Subtype.val) = ∫ s in B, f s := by
  have h := (MeasurableEmbedding.subtype_coe measurableSet_Ioi).setIntegral_map
    (μ := (volume : Measure ℝ).comap (Subtype.val : Radius → ℝ)) f B
  rw [map_comap_subtype_coe measurableSet_Ioi, Measure.restrict_restrict_of_subset hB] at h
  exact h.symm

private theorem setIntegral_radialBand_pow {r₁ r₂ : ℝ} (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂)
    (f : ℝ → ℝ) :
    (∫ s in radialBand r₁ r₂, f s ∂Measure.volumeIoiPow 2) =
      ∫ s in r₁..r₂, s ^ 2 * f s := by
  have hB : Icc r₁ r₂ ⊆ Ioi (0 : ℝ) := fun _ hs => h₁.trans_le hs.1
  have hmeas : MeasurableSet (radialBand r₁ r₂) := measurableSet_Icc.preimage measurable_subtype_coe
  have hdensity : Measurable (fun s : Radius => ENNReal.ofReal ((s : ℝ) ^ 2)) :=
    (measurable_subtype_coe.pow_const _).ennreal_ofReal
  rw [Measure.volumeIoiPow, setIntegral_withDensity_eq_setIntegral_toReal_smul hdensity
    (Filter.Eventually.of_forall (fun _ => ENNReal.ofReal_lt_top)) _ hmeas]
  simp only [ENNReal.toReal_ofReal (sq_nonneg _), smul_eq_mul]
  rw [radialBand, setIntegral_positiveRadius hB (fun s => s ^ 2 * f s),
    intervalIntegral.integral_of_le h₁₂]
  exact setIntegral_congr_set Ioc_ae_eq_Icc.symm

private theorem prism_preimage (T : Set Direction) (r₁ r₂ : ℝ) :
    (Subtype.val : ({(0 : Space)}ᶜ : Set Space) → Space) ⁻¹' prism T r₁ r₂ =
      homeomorphUnitSphereProd Space ⁻¹' (T ×ˢ radialBand r₁ r₂) := by
  let e := homeomorphUnitSphereProd Space
  ext y
  constructor
  · rintro ⟨p, hp, heq⟩
    have h : e.symm p = y := Subtype.ext heq
    change e y ∈ T ×ˢ radialBand r₁ r₂
    rw [← h, e.apply_symm_apply]
    exact hp
  · intro hy
    exact ⟨e y, hy, congrArg Subtype.val (e.symm_apply_apply y)⟩

private theorem exteriorPolarReduction {T : Set Direction} (hT : IsCompact T)
    {r₁ r₂ : ℝ} (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) {x : Space} (hx : r₂ < ‖x‖) :
    IntegrableOn (fun u => radialPotential ‖x‖ r₁ r₂ (observationCosine x u)) T angularMeasure ∧
    (∫ y in prism T r₁ r₂, ‖x - y‖⁻¹) =
      ∫ u in T, radialPotential ‖x‖ r₁ r₂ (observationCosine x u) ∂angularMeasure := by
  let e := homeomorphUnitSphereProd Space
  let B := T ×ˢ radialBand r₁ r₂
  let f : Space → ℝ := fun y => ‖x - y‖⁻¹
  let g : Direction × Radius → ℝ := fun p => f ((p.2 : ℝ) • (p.1 : Space))
  have hR : 0 < ‖x‖ := (h₁.trans_le h₁₂).trans hx
  have hx0 : x ≠ 0 := norm_pos_iff.mp hR
  have hband : IsCompact (radialBand r₁ r₂) := by
    rw [Topology.IsEmbedding.subtypeVal.isCompact_iff, radialBand_image h₁]
    exact isCompact_Icc
  have hB : IsCompact B := hT.prod hband
  have hp : MeasurePreserving e ((volume : Measure Space).comap Subtype.val)
      (angularMeasure.prod (Measure.volumeIoiPow 2)) := by
    simpa [e, angularMeasure, finrank_euclideanSpace_fin] using
      (volume : Measure Space).measurePreserving_homeomorphUnitSphereProd
  have hB0 : prism T r₁ r₂ ⊆ ({(0 : Space)}ᶜ : Set Space) := by
    rintro _ ⟨p, _, rfl⟩
    exact (e.symm p).property
  have hNonzero : MeasurableSet ({(0 : Space)}ᶜ : Set Space) :=
    (measurableSet_singleton _).compl
  have hset := (MeasurableEmbedding.subtype_coe hNonzero).setIntegral_map
    (μ := (volume : Measure Space).comap Subtype.val) f (prism T r₁ r₂)
  rw [map_comap_subtype_coe hNonzero, Measure.restrict_restrict_of_subset hB0, prism_preimage]
    at hset
  have hcomp : g ∘ e = (fun p : ({(0 : Space)}ᶜ : Set Space) => f (p : Space)) := by
    funext p
    exact congrArg f (congrArg Subtype.val (e.symm_apply_apply p))
  have hgi : IntegrableOn g B (angularMeasure.prod (Measure.volumeIoiPow 2)) := by
    apply (hp.integrableOn_comp_preimage e.measurableEmbedding).mp
    change IntegrableOn (g ∘ e) (homeomorphUnitSphereProd Space ⁻¹' (T ×ˢ radialBand r₁ r₂)) _
    rw [hcomp, ← prism_preimage]
    change IntegrableOn (f ∘ Subtype.val) (Subtype.val ⁻¹' prism T r₁ r₂)
      ((volume : Measure Space).comap Subtype.val)
    rw [← (MeasurableEmbedding.subtype_coe hNonzero).integrableOn_map_iff,
      map_comap_subtype_coe hNonzero, IntegrableOn, Measure.restrict_restrict_of_subset hB0]
    have hC : IsCompact (prism T r₁ r₂) := hB.image (by fun_prop)
    change IntegrableOn f (prism T r₁ r₂) volume
    apply ContinuousOn.integrableOn_compact hC
    have hc : Continuous (fun y : Space => ‖x - y‖) := by fun_prop
    refine hc.continuousOn.inv₀ ?_
    rintro _ ⟨p, ⟨_, hp₂⟩, rfl⟩
    apply norm_ne_zero_iff.mpr
    intro hz
    have heq : x = (p.2 : ℝ) • (p.1 : Space) := sub_eq_zero.mp hz
    have hn := congrArg norm heq
    rw [norm_smul, Real.norm_eq_abs, abs_of_pos p.2.property, direction_norm, mul_one] at hn
    exact (not_le_of_gt hx) (hn ▸ hp₂.2)
  have hinner (u : Direction) :
      (∫ s in radialBand r₁ r₂, g (u, s) ∂Measure.volumeIoiPow 2) =
        radialPotential ‖x‖ r₁ r₂ (observationCosine x u) := by
    change (∫ s in radialBand r₁ r₂, ‖x - (s : ℝ) • (u : Space)‖⁻¹
      ∂Measure.volumeIoiPow 2) = _
    rw [setIntegral_radialBand_pow h₁ h₁₂ (fun s => ‖x - s • (u : Space)‖⁻¹)]
    apply intervalIntegral.integral_congr
    intro s hs
    rw [uIcc_of_le h₁₂] at hs
    change s ^ 2 * ‖x - s • (u : Space)‖⁻¹ = s ^ 2 / distance ‖x‖ s (observationCosine x u)
    rw [norm_sub_smul_direction hx0 u (h₁.trans_le hs.1).le, div_eq_mul_inv]
  constructor
  · have hgp : Integrable g ((angularMeasure.restrict T).prod
        ((Measure.volumeIoiPow 2).restrict (radialBand r₁ r₂))) := by
      rw [Measure.prod_restrict]
      exact hgi
    exact hgp.integral_prod_left.congr (Filter.Eventually.of_forall hinner)
  calc
    (∫ y in prism T r₁ r₂, f y) =
        ∫ p in e ⁻¹' B, f (p : Space) ∂(volume : Measure Space).comap Subtype.val := hset
    _ = ∫ p in B, g p ∂angularMeasure.prod (Measure.volumeIoiPow 2) := by
      have hh := hp.setIntegral_preimage_emb e.measurableEmbedding g B
      change (∫ p in e ⁻¹' B, (g ∘ e) p ∂(volume : Measure Space).comap Subtype.val) = _ at hh
      rwa [hcomp] at hh
    _ = ∫ u in T, ∫ s in radialBand r₁ r₂, g (u, s) ∂Measure.volumeIoiPow 2 ∂angularMeasure :=
      setIntegral_prod g hgi
    _ = _ := by
      apply setIntegral_congr_fun hT.measurableSet
      exact fun u _ => hinner u

/-- Polar coordinates and Fubini for an exterior observation point.
The angular region may be any compact subset of the sphere. -/
theorem volumeIntegral_eq_angularIntegral {T : Set Direction} (hT : IsCompact T)
    {r₁ r₂ : ℝ} (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) {x : Space} (hx : r₂ < ‖x‖) :
    (∫ y in prism T r₁ r₂, ‖x - y‖⁻¹) =
      ∫ u in T, radialPotential ‖x‖ r₁ r₂ (observationCosine x u) ∂angularMeasure :=
  (exteriorPolarReduction hT h₁ h₁₂ hx).2

theorem integrableOn_radialPotential {T : Set Direction} (hT : IsCompact T)
    {r₁ r₂ : ℝ} (h₁ : 0 < r₁) (h₁₂ : r₁ ≤ r₂) {x : Space} (hx : r₂ < ‖x‖) :
    IntegrableOn (fun u => radialPotential ‖x‖ r₁ r₂ (observationCosine x u)) T angularMeasure :=
  (exteriorPolarReduction hT h₁ h₁₂ hx).1

end TSP
