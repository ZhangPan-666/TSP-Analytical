# Lean 4 proof of the TSP volume-to-endpoint reduction

This project formalizes the Newton potential of a homogeneous triangular spherical prism (TSP). It covers observations at radius $R=\lVert x\rVert>r_2\ge r_1>0$ and oriented edges with direction cosine $c=h_e\cos\phi$, where $0<h_e<1$. The main theorem gives the following identity in exact real arithmetic.

$$
U(x)=\int_\Omega \frac{\mathrm{d}y}{\lVert x-y\rVert}
=\frac{1}{R}\left[
\omega_{\mathcal T}\frac{r_2^3-r_1^3}{3}
+\sum_{e=1}^3\kappa_e\bigl(P_e(\phi_{e,1})-P_e(\phi_{e,0})\bigr)
\right].
$$

The gravitational potential is $V=G\rho U$. On each edge, $P_e$ takes the difference $\mathcal Q_{r_2}-\mathcal Q_{r_1}$ of the radial endpoint primitives and subtracts the antiderivatives of the two polar interpolation terms. The latter are expressed as third-kind Legendre integrals with zero elliptic parameter. Signed interval integrals define the Legendre functions at the full amplitudes, including amplitudes that cross periods.

Texts S1–S2 of the companion Supporting Information give the analytical derivation. The entry point for reading the Lean proof is [`volumeIntegral_eq_endpoints_of_chart`](TSPChart.lean#L183). Its declaration states the hypotheses and conclusion, followed by the proof after `by`.

## Proof structure

| Step | Lean declarations |
| --- | --- |
| Polar coordinates and Fubini's theorem give the radial factor $s^2$ in the three-dimensional volume integral. | [TSPVolume.lean](TSPVolume.lean) contains `volumeIntegral_eq_angularIntegral`. |
| Differentiate the radial primitive $H$ and the single-radius elliptic primitive $\mathcal Q_s$. | [TSPEndpoint.lean](TSPEndpoint.lean) contains `hasDerivAt_radialPrimitive` and `hasDerivAt_ellipticEndpoint`. |
| Regularize the polar terms, evaluate the polar moments, and derive the great-circle edge coefficient. | [TSPSpherical.lean](TSPSpherical.lean) contains `polarMoments_eq_integrals`, `shellKernel_eq_edgeKernel`, and `greatCircleIntegral_eq_endpoints`. |
| Differentiate the regularized radial integral, including at $c=\pm1$. | [TSPShellDiff.lean](TSPShellDiff.lean) contains `hasDerivAt_shellKernel` and `shellKernel_curl`. |
| Apply Green's theorem, reparametrize the sides by edge phase, and assemble the endpoint identity. | [TSPChart.lean](TSPChart.lean) contains `square_shell_stokes`, `chartSide_eq_edgeIntegral`, and `volumeIntegral_eq_endpoints_of_chart`. |

The regularized kernel is

$$
K(c)=\int_{r_1}^{r_2}\frac{s^3}{D_s(c)+R-sc}\,\mathrm{d}s,
\qquad D_s(c)=\sqrt{(s-R)^2+2Rs(1-c)}.
$$

Its denominator is positive at both poles in the exterior case. Taking the exterior derivative of the one-form $K(c)(X\,\mathrm{d}Y-Y\,\mathrm{d}X)$ gives the coefficient $R\int_{r_1}^{r_2}s^2/D_s(c)\,\mathrm{d}s-(r_2^3-r_1^3)/3$ multiplying the oriented area form. The angular derivative of $\mathcal Q_s$ contains an extra term independent of $s$. This term cancels in the difference between the two radial endpoints.

## Geometric hypotheses and scope

The main theorem takes an oriented sphere chart as input. The `hArea` hypothesis relates the spherical area measure to the Jacobian of the chart, and `hEdges` identifies the oriented boundary of the chart with the three great-circle edges. The theorem `square_shell_stokes` proves Green's theorem in the chart, and `chartSide_eq_edgeIntegral` reduces the edge correspondence to the phase derivative, its endpoint values, and the tangent coefficient.

For oriented unit vertices $a_1,a_2,a_3$ in a common open hemisphere, the geometric hypotheses can be addressed using the following parametrization of the unit square.

$$
p(\xi,\eta)=(1-\xi)a_1+\xi(1-\eta)a_2+\xi\eta a_3,
\qquad u(\xi,\eta)=p(\xi,\eta)/\lVert p(\xi,\eta)\rVert.
$$

The oriented area density is $\det(u,u_\xi,u_\eta)=\xi\det(a_1,a_2,a_3)/\lVert p\rVert^3$. It follows by comparing the Jacobian of the volume map $(s,\xi,\eta)\mapsto su(\xi,\eta)$ with the polar decomposition of volume. With positive orientation, the sides $\eta=0$, $\xi=1$, and the reversed side $\eta=1$ give the three edges in order. The side $\xi=0$ collapses to $a_1$. A theorem starting from the three vertices alone would need to establish these geometric hypotheses in Lean.

The formalized endpoint identity uses the Legendre representation. Further verification is required for its equivalence to the Carlson representation in the C++ implementation, floating-point evaluation, observation radii $R\le r_2$, and degenerate edge configurations.

## Build and check

Install [Elan](https://github.com/leanprover/elan) and run the following commands from this directory. The `lean-toolchain` file pins Lean 4.19.0, and `lakefile.toml` pins the mathlib revision.

```text
lake update
lake exe cache get
lake build
```

[TSPLean.lean](TSPLean.lean) runs `#print axioms` on seven principal theorems. All six source files were checked with the pinned toolchain. The seven theorems depend only on `propext`, `Classical.choice`, and `Quot.sound`.
