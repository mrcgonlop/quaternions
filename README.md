# Quaternions — Quaternionic Vacuum Electrodynamics Simulator

**A computational exploration of pre-Heaviside electrodynamics, vacuum polarizability, and longitudinal electromagnetic phenomena**

---

## Motivation

Modern electrodynamics, as taught and practiced, rests on the Heaviside-Gibbs reformulation of Maxwell's original theory. Between 1884 and 1887, Oliver Heaviside reduced Maxwell's 20 quaternion equations to the 4 vector equations we use today. This was presented as a simplification — a notational convenience that preserved all physical content while discarding the "unnecessary" complexity of Hamilton's quaternions.

This project begins from the premise that the simplification was not physically neutral. Specifically, three reductions occurred that may have eliminated real, observable phenomena from the theoretical framework:

1. **The quaternion scalar part was discarded.** Maxwell's quaternionic field potential Q = φ/c + Aₓi + Aᵧj + A_zk naturally couples the scalar potential φ and the vector potential **A** into a single algebraic object. When the quaternionic derivative acts on Q, it produces not just **E** and **B** but also a scalar field S = (1/c²)∂φ/∂t + ∇·**A**. In the Heaviside formulation, this scalar is set to zero by the Lorenz gauge condition. In the quaternion formulation, it is a dynamical degree of freedom — a propagating scalar longitudinal mode of the electromagnetic field.

2. **The potentials were demoted to mathematical fictions.** The Heaviside formulation treats **E** and **B** as the physically real quantities and (φ, **A**) as convenient but redundant intermediaries defined only up to a gauge transformation. The Aharonov-Bohm effect (1959, experimentally confirmed by Tonomura et al. in 1986) demonstrated unambiguously that the potentials have physical reality — a charged particle's quantum phase is affected by **A** even in regions where **E** = **B** = 0. If the potentials are primary, then gauge freedom is not mathematical redundancy but represents real physical degrees of freedom of the vacuum.

3. **Weber-type longitudinal interactions were rendered invisible.** Wilhelm Weber's electrodynamics (1846) predated Maxwell and described forces between charges that depend on their relative velocity and acceleration. These produce longitudinal forces along the line connecting charges — forces that the Lorentz force law (derived from the Heaviside formulation) cannot produce. André Assis showed that Weber's force law is recovered from the retarded potentials expanded to second order in v/c. The longitudinal forces are physically present in the full potential-based theory but vanish when gauge conditions are imposed and potentials are treated as secondary.

The convergence of these three lines of reasoning points toward a single unified framework: **quaternionic potential electrodynamics in a polarizable vacuum, without gauge fixing**. This project aims to simulate that framework, visualize its predictions, and identify experimental configurations where its predictions diverge from standard electrodynamics.

---

## Theoretical Foundation

### 1. Quaternionic Electrodynamics

#### 1.1 The Quaternion Potential

A Hamilton quaternion has the form:

```
q = w + xi + yj + zk
```

where i² = j² = k² = ijk = -1. The electromagnetic potential is naturally a quaternion:

```
Q = φ/c + Aₓi + Aᵧj + A_zk
```

The quaternionic nabla operator is:

```
∇_q = (1/c)∂/∂t + (∂/∂x)i + (∂/∂y)j + (∂/∂z)k
```

The quaternion product ∇_q Q (using Hamilton's non-commutative multiplication) produces:

```
∇_q Q = S + Eᵢi + Eⱼj + Eₖk + Bᵢi + Bⱼj + Bₖk
```

More precisely, decomposing the quaternion product:

- **Scalar part:** S = (1/c²)∂φ/∂t + ∇·A — this is the Lorenz gauge term, normally set to zero
- **Vector part (symmetric):** **E** = -∇φ - ∂**A**/∂t — the electric field
- **Vector part (antisymmetric):** **B** = ∇ × **A** — the magnetic field

The critical point: **S is a dynamical scalar field that propagates.** The wave equation for S, derived from the quaternionic Maxwell equations without gauge fixing, is:

```
□S = (1/c²)∂²S/∂t² - ∇²S = -ρ/ε₀ (source term)
```

This is a scalar longitudinal wave. It propagates at c in vacuum. It is driven by charge density fluctuations. And it is completely invisible in the standard Heaviside formulation because S is defined to be zero.

#### 1.2 Extended Maxwell Equations

When you refuse to impose gauge conditions, Maxwell's equations acquire additional terms. The extended system (following van Vlaenderen and Waser, 2001) is:

```
∇·E = ρ/ε₀ - ∂S/∂t          (modified Gauss's law)
∇·B = 0                        (no magnetic monopoles)
∇×E = -∂B/∂t                   (Faraday's law, unchanged)
∇×B = μ₀J + μ₀ε₀∂E/∂t + ∇S   (modified Ampère's law)
```

Plus the scalar field equation:

```
□S = -ρ/ε₀
```

The additional terms (-∂S/∂t in Gauss's law and ∇S in Ampère's law) couple the scalar field to the conventional vector fields. These terms are exactly what Heaviside eliminated. They represent longitudinal electromagnetic effects — compression and rarefaction of the potential field along the direction of propagation, rather than transverse oscillation perpendicular to it.

#### 1.3 Physical Interpretation

Transverse EM waves (light, radio) oscillate E and B perpendicular to the propagation direction. These are the only solutions when S = 0 (the gauge condition).

Longitudinal EM waves oscillate the potentials along the propagation direction. The scalar field S represents the "breathing mode" — a compression wave in the electromagnetic potential. Physically, this is a periodic fluctuation in the electrostatic potential that propagates through space, accompanied by a longitudinal oscillation of the vector potential.

This is directly analogous to the difference between transverse (shear) waves and longitudinal (compression) waves in an elastic solid. Standard EM theory only admits the transverse mode. The quaternionic theory admits both.

### 2. The Gauge Symmetry Debate

The physical status of the S field — and of the entire quaternionic extension — ultimately rests on a foundational question in physics that is less settled than textbooks suggest: **is gauge symmetry a redundancy in nature, or a redundancy in our description of nature?** This question runs through all of modern physics and is worth understanding in detail, because the answer determines whether this project is exploring a genuine physical degree of freedom or a mathematical artifact.

#### 2.1 What Gauge Symmetry Is

To make this concrete, consider an analogy. Imagine describing the topography of a landscape. You need to assign an altitude to every point. But altitude relative to what? You could measure from sea level, from the center of the Earth, from the top of a local hill — it doesn't matter, because the *differences* in altitude are what determine which way water flows. The choice of reference is the "gauge." Physics that depends only on altitude differences is "gauge invariant."

In electromagnetism, the potentials φ (voltage) and **A** (the vector potential) are like altitude. The electric and magnetic fields **E** and **B** are like the slopes — the differences. You can shift the potentials by adding any smooth function (a gauge transformation, φ → φ + ∂χ/∂t, **A** → **A** - ∇χ) without changing **E** and **B**, just as you can shift all altitudes by a constant without changing any slopes.

The standard interpretation says: since the physics only depends on the slopes (**E** and **B**), the absolute altitude (φ and **A**) is meaningless. The potentials are convenient mathematical fictions. Gauge freedom is redundancy in the description, not a feature of nature.

This seems perfectly reasonable. But there are cracks.

#### 2.2 The Aharonov-Bohm Crack

Consider a long solenoid — a coil of wire carrying current, creating a magnetic field **B** confined entirely inside the coil. Outside the coil, **B** = 0 everywhere. No magnetic field. But the vector potential **A** is *not* zero outside — it circulates around the solenoid like water flowing in a ring around a drain.

In 1959, Aharonov and Bohm predicted that an electron traveling outside the solenoid — in a region with zero field — would still be affected by the potential **A**. Its quantum phase would shift by an amount proportional to ∮**A**·d**l**, the integral of **A** around the solenoid. This was experimentally confirmed beyond any doubt by Tonomura in 1986 using electron holography.

The standard response is: "This doesn't break gauge invariance. The quantity ∮**A**·d**l** is gauge-invariant — it doesn't change when you perform a gauge transformation. So the electron is responding to a gauge-invariant quantity, not to the 'unphysical' potential itself."

This is mathematically correct. But the physics is uncomfortable. The electron is at point P, outside the solenoid. At point P, **E** = 0 and **B** = 0. There is no field at the electron's location. Yet the electron's behavior changes depending on what's happening inside the solenoid, which might be meters away. The information about the solenoid's state is carried to point P entirely by **A** — the supposedly unphysical potential.

This leaves three interpretive positions, each with real costs:

**Position 1 — Orthodox (gauge redundancy):** The potential **A** is not physical. The electron responds nonlocally to the flux ∫**B**·d**S** enclosed by its path, even though that flux is far away. Gauge invariance is exact, and the price is accepting a form of nonlocality in how electromagnetic effects propagate.

**Position 2 — Potential realism:** The potential **A** is the physically real thing. The electron responds locally to the **A** at its location. Gauge invariance is a mathematical property of the equations, not a statement about what's real. The price is that gauge freedom must be given physical meaning — if **A** is real, what does it mean that you can shift it by ∇χ?

**Position 3 — Fiber bundle (modern mathematical view):** Neither **A** nor **B** alone is the fundamental object. The real object is a geometric structure called a "connection on a fiber bundle" — a way of comparing quantum phases at different points in space. **A** is a local coordinate representation of this connection, and gauge transformations are changes of coordinates on the bundle. This is mathematically elegant and the basis for all modern gauge theories, but it sidesteps the physical question: what *is* this connection, physically, in the vacuum between the solenoid and the electron?

#### 2.3 A Deeper Analogy

To sharpen the stakes, think of a river flowing over terrain. The velocity field of the water (how fast it flows and in what direction at each point) is like **E** and **B** — it is what you directly measure and what exerts forces. The *depth* of the river at each point is like the potential — it is a scalar field that determines the velocity through its gradients.

The orthodox view says: only the velocity of the water (the current) is real. The depth is just a mathematical convenience for computing the velocity. If you uniformly raise the entire riverbed by one meter, nothing changes about the flow.

The potential-realist view says: the depth is physical. If you create a region where the riverbed is shaped so that the water surface is perfectly flat (no current, no flow — analogous to **E** = **B** = 0), but the depth varies in a topologically nontrivial way (imagine a hidden underwater mountain that the surface flow routes around), then a boat floating on the calm surface will still be affected by the underwater topography. It drifts, very slowly, in ways that the surface current (being zero) cannot explain. You have to look at the depth (the potential) to understand what is happening.

The Aharonov-Bohm effect is this situation. The electron is the boat. The zero-field region is the calm surface. The solenoid is the underwater mountain. And the electron drifts in ways that the surface fields cannot explain — but the potential can.

#### 2.4 Why This Debate Determines the S Field's Physical Status

The scalar field S = (1/c²)∂φ/∂t + ∇·**A** is analogous to asking whether the river can have a mode where the depth oscillates up and down without creating any surface current. A "breathing" mode of the potential that doesn't show up in **E** or **B**.

If gauge symmetry is exactly a redundancy (Position 1), then S = 0 is not a choice but a consequence. The depth cannot oscillate independently of the surface flow because the depth is not real. There are no scalar longitudinal modes. The Lorenz condition is automatic, not imposed.

If the potential is physical (Position 2), then S = 0 is an additional physical constraint — like a law that the river depth must stay constant. Removing that constraint means the river can support compression waves — the depth oscillating as a pulse propagates through — in addition to surface waves. That is the scalar longitudinal mode.

The honest answer is that this is not settled. QED, built on gauge invariance, is the most precisely tested theory in physics. But "it works as a calculational tool" and "it reflects the fundamental ontology of nature" are different claims, and the history of physics is full of cases where a framework that worked perfectly was later understood as an effective description of something deeper.

#### 2.5 Why the Debate Remains Open

The same foundational question arises in every gauge theory, not just electromagnetism. The Standard Model of particle physics is built on gauge theories (U(1) × SU(2) × SU(3)). The Yang-Mills problem — one of the seven Millennium Prize Problems with a million-dollar bounty — concerns the mathematical structure of non-abelian gauge theories. In gravity, the analogous question is whether spacetime coordinates are physically real or just labels (the substantivalism vs. relationalism debate that has run since Newton and Leibniz).

This project's approach — refuse to impose the gauge condition, let S be dynamical, simulate what happens, and look for observable differences — is a way to make progress on this foundational question through quantitative prediction rather than philosophical argument. If S produces measurable effects in simulated configurations, gauge symmetry is approximate at best. If S stays negligibly small in every realistic configuration even without being forced to zero, gauge symmetry is effectively exact and the question becomes philosophical rather than physical.

### 3. Identified Gaps in Standard Electrodynamics

Before proceeding with the extended framework, it is worth surveying where the standard Maxwell-Heaviside formulation has genuine unresolved issues or underexplored territory. These gaps motivate the exploration independently of any specific alternative theory.

#### 3.1 The Interpretation of Gauge Symmetry

As detailed in §2, the physical meaning of gauge invariance is an open question in the foundations of physics. The Aharonov-Bohm effect demonstrates that the potentials carry physical information that the fields do not encode locally. Whether this is best understood as nonlocality, potential realism, or fiber bundle geometry remains debated among working physicists and philosophers of physics. This is not a fringe concern — it touches the mathematical foundations of the entire Standard Model.

#### 3.2 The Near-Field Regime Is Underexplored Experimentally

Most precision tests of electrodynamics are in the radiation (far-field) regime: optics, radio, spectroscopy. The near-field regime — within a wavelength of sources, where retardation effects, longitudinal field components, and potential-dependent phenomena are strongest — is much less thoroughly tested. This matters because the differences between Weber-type electrodynamics and Lorentz electrodynamics are largest in the near field.

The standard response is "we have tested Maxwell's equations to extraordinary precision," which is true but almost entirely in the far-field/radiation regime. The near-field testing has gaps, particularly for transient (pulsed) phenomena where ∂/∂t terms are large. The simulation should specifically target near-field configurations where the quaternionic framework predicts the largest deviations.

#### 3.3 The Vacuum Has Structure That Classical EM Ignores

QED says the vacuum has structure: virtual pairs, zero-point fluctuations, the Euler-Heisenberg nonlinearity at high fields. The Casimir effect is real and measured. Vacuum birefringence (the vacuum acting as a polarization-dependent medium in strong magnetic fields) was measured by PVLAS and observed in neutron star spectra. These are small effects at laboratory scales, but they are real, and they mean the vacuum is not the inert stage that classical EM assumes.

The question of whether there is a useful classical effective description of these vacuum properties — which is essentially what the polarizable vacuum model attempts — is reasonable, even if the specific PV implementation has issues (see §5.3 for caveats).

#### 3.4 Topological Aspects of Classical EM Are Underexplored

The mathematical framework for topological solitons, knotted fields, and helicity conservation in EM is well-developed theoretically (Rañada, Irvine & Bouwmeester, and others), but the experimental exploration is in its infancy. The connection to plasma physics through spheromaks and Taylor relaxation is real and underappreciated outside the fusion community. The question "what happens to the field topology after plasma recombination?" appears to be genuinely underexplored — standard plasma physics experiments end measurement when the plasma is gone.

#### 3.5 Ball Lightning as an Empirical Anchor

Ball lightning is one of the few genuinely anomalous phenomena that mainstream physics acknowledges as real but unexplained. The observational database is large and includes credible witnesses (pilots, scientists, multiple independent observers of the same event). The key observational constraints are severe and collectively rule out all conventional models:

- Luminous spheres persisting for seconds to minutes (plasma at atmospheric pressure recombines in microseconds)
- Characteristic size 10-30 cm
- Occasionally passing through glass without damage (rules out particle-based models)
- Not obviously connected to any sustained energy source after formation
- Formed in association with lightning but not always at the strike point
- Mobile and coherent, maintaining spherical shape during motion

No conventional model satisfies all constraints simultaneously. Chemical models (burning silicon nanoparticles) can get seconds of luminosity but fail on glass penetration. Plasma models fail on lifetime. Microwave cavity models require an external source that is not observed. Ball lightning serves as a strong empirical anchor for this project: if the quaternionic framework can produce a model that satisfies all constraints simultaneously, that is significant regardless of what one thinks about the broader theoretical framework.

### 4. Weber Electrodynamics and Its Connection to Potentials

#### 4.1 Weber's Force Law

Weber's force between two charges q₁ and q₂ separated by distance r is:

```
F_W = (q₁q₂ / 4πε₀r²) · r̂ · [1 - (ṙ²/2c²) + (r·r̈/c²)]
```

where ṙ = dr/dt is the relative radial velocity and r̈ = d²r/dt² is the relative radial acceleration.

The first term (1) is Coulomb's law. The second term (-ṙ²/2c²) is a velocity-dependent correction. The third term (r·r̈/c²) is an acceleration-dependent correction. These corrections are always directed along the line connecting the charges — they are intrinsically longitudinal.

#### 4.2 Derivation from Retarded Potentials

The connection to the quaternionic potential framework was established by several authors (Assis 1994, Wesley 1990, and others). Start with the Liénard-Wiechert potentials for a moving charge:

```
φ(r,t) = q / 4πε₀ · 1 / (R - R·v/c)    evaluated at retarded time
A(r,t) = v/c² · φ(r,t)                    evaluated at retarded time
```

Expand these to second order in v/c (keeping terms up to v²/c² and a/c²), compute the electric field E = -∇φ - ∂A/∂t, and you recover Weber's force law.

This means Weber's "extra" forces are not ad hoc — they are the near-field, low-velocity content of the fully retarded potential theory. They were always in Maxwell's equations, hiding in the retardation of the potentials. The Heaviside simplification, by focusing on E and B and imposing gauge conditions, obscured these forces because they manifest as subtle correlations between the scalar and vector potentials that the gauge condition eliminates.

#### 4.3 Longitudinal Current Forces

The most experimentally accessible prediction unique to Weber's formulation is the existence of longitudinal forces in current-carrying conductors. When current flows through a wire, the drift velocity of electrons creates Weber-type forces along the wire axis. Standard Lorentz force only produces forces perpendicular to the current (the basis of motors and generators). Weber predicts additional forces parallel to the current.

Peter Graneau's experiments (1980s-2000s) demonstrated this: water jets carrying high current break into segments; mercury in U-tubes carrying current exhibits longitudinal forces; wires carrying high pulse current fragment in patterns inconsistent with Lorentz forces but consistent with Weber predictions.

The magnitude scales as I²/c², so these forces are small at ordinary currents but become significant at high pulse currents (kiloamps).

### 5. The Polarizable Vacuum

#### 5.1 Vacuum Permittivity as Dynamical Variable

The standard vacuum has constant permittivity ε₀ and permeability μ₀. The speed of light c = 1/√(ε₀μ₀) is therefore constant everywhere.

The Polarizable Vacuum (PV) model, developed by Harold Puthoff (1999, 2002) and building on earlier work by Dicke (1957) and Wilson (1921), treats these as the ground-state values of a dynamical medium. Near mass-energy concentrations, the vacuum polarizes:

```
ε(x) = ε₀ · K(x)
μ(x) = μ₀ · K(x)
c_local(x) = c₀ / K(x)
```

where K(x) ≥ 1 is the vacuum polarizability index. In the absence of mass-energy, K = 1 and we recover standard vacuum. Near a gravitating mass:

```
K(r) ≈ e^(2GM/rc²)
```

This reproduces the Schwarzschild metric predictions to post-Newtonian order: gravitational redshift (photon frequency shifts because c changes), light bending (refraction in a medium with variable refractive index n = K), and perihelion precession.

#### 5.2 Vacuum Polarization and Gauge Freedom

Here is the crucial synthesis: **if the vacuum is polarizable, then the gauge freedom of electrodynamics acquires physical meaning.**

A gauge transformation changes:

```
φ → φ + ∂χ/∂t
A → A - ∇χ
```

for arbitrary scalar function χ. In standard theory, this changes the potentials but not the fields — hence it is "unphysical." But if the potentials are primary (Aharonov-Bohm) and the vacuum is a dynamical medium, then a gauge transformation changes the state of the vacuum. Different gauges correspond to different vacuum polarization configurations.

The scalar field S = (1/c²)∂φ/∂t + ∇·**A** — which the Lorenz gauge sets to zero — can be reinterpreted as a measure of local vacuum polarization. Driving S away from zero means locally altering the vacuum state. This is not just a mathematical degree of freedom; it is the fundamental mechanism by which electromagnetic fields interact with the vacuum medium.

#### 5.3 Known Limitations and Caveats of the PV Model

**This section documents known issues with the PV framework that must be kept in mind during simulation design and interpretation of results.**

**The K field coupling is extraordinarily small at laboratory scales.** The Euler-Heisenberg effective Lagrangian, which describes QED vacuum nonlinearity, operates at field strengths approaching the Schwinger limit (~10¹⁸ V/m). The coupling constant η ≈ 2α/(45π) ≈ 10⁻⁴ multiplied by u_field/u_S means that for any laboratory-accessible field, the K deviation from 1 is astronomically small — of order 10⁻³⁰ or smaller for fields of a few Tesla. The simulation must track this quantitatively and report the actual predicted K deviation for every configuration studied. If the number is immeasurably small, that result is important and must be reported honestly.

**The plasma-K / vacuum-K analogy has limits.** The document's argument that "plasma effective-K" can substitute for vacuum K because the equations are "isomorphic" requires care. A plasma dielectric response and QED vacuum polarization are physically distinct phenomena that share mathematical structure in certain limits. The plasma response arises from real charged particles; the vacuum response arises from virtual pair polarization. They are governed by different physics at the microscopic level even if they produce similar macroscopic field equations. The simulation should clearly distinguish between vacuum-K effects (tiny, requiring Schwinger-scale fields) and plasma-K effects (large, achievable in laboratory plasmas) and not conflate the two.

**The Weber attraction threshold argument has a domain-of-validity issue.** The argument that in a polarized vacuum with K > 1, the Weber attraction threshold drops from v > c₀√2 to v > c₀√2/K appears to make Weber attraction achievable at sub-relativistic velocities. However, the Weber force law is itself a non-relativistic approximation — it is the v/c expansion of the retarded potential to second order. Pushing it into regimes where K-enhanced corrections become large means extrapolating a perturbative result outside its domain of validity. The self-consistent treatment would require returning to the full retarded potentials in a medium with variable K and re-deriving the force law, rather than simply substituting c → c/K into the approximate formula. The simulation should implement both the approximate (K-substituted Weber) and the full retarded potential calculation and compare them to quantify the error introduced by the approximation.

**The EVO/charge cluster self-binding argument has a bootstrapping problem.** The proposed mechanism — dense electrons → strong field → K >> 1 → enhanced Weber attraction → binding → denser electrons — is logically self-consistent as a description of a fixed point, but it requires overcoming Coulomb repulsion to reach the threshold density where the feedback loop engages. The document identifies this as a "bifurcation" with a threshold, which is correct, but the threshold density and the energy required to reach it from a dispersed initial state must be estimated quantitatively. Without that estimate, the self-sustaining pocket is a theoretical fixed point that may be dynamically inaccessible. The simulation should attempt to determine this threshold.

#### 5.4 Connection to Casimir Effect

The Casimir effect — the attractive force between closely spaced parallel conducting plates in vacuum — is conventionally explained as a consequence of vacuum zero-point energy: the plates restrict the modes of the vacuum electromagnetic field between them, creating a pressure imbalance.

In the PV model, the explanation is more direct: the plates constrain the vacuum polarization state. Between the plates, the vacuum is differently polarized than outside. The resulting gradient in ε(x) produces a force — essentially electrostatic attraction mediated by the vacuum as dielectric. The mathematical predictions converge, but the physical picture is different: the vacuum is not a sea of virtual photons but a polarizable medium responding to boundary conditions.

Schwinger, DeRaad, and Milton (1978) proved that the Casimir force is exactly equal to the sum of all pairwise retarded van der Waals interactions between atoms in the two plates — a derivation that does not invoke zero-point energy at all. In the virtual pair language: a pair forms coherently only if its coherence length `λ_pair ~ ℏ/(mₑc) = λ_C` fits in the available space. Between plates separated by `d < λ_C`, that pair mode is excluded. Outside the plates, no such constraint exists. The asymmetry — more coherent pair formation outside than inside, each pair attracted to its image charge in the nearby conductor — creates the net inward force. The mode-pressure picture and the van der Waals pair picture are mathematically identical; the latter makes the local physical mechanism explicit.

This distinction matters because it suggests engineering pathways: if you can create strong vacuum polarization gradients by other means (intense electromagnetic fields, specific geometrical configurations, topological field structures), you might produce Casimir-like forces without needing nanometer plate separations.

#### 5.5 The Quantum Vacuum as a Virtual Pair Plasma and K Field Dynamics

The classical K field has a concrete quantum mechanical interpretation. The QED vacuum contains virtual electron-positron pairs that spontaneously form, exist for time Δt ~ ℏ/(2mₑc²), and annihilate. During their brief existence, these pairs are displaced by an external field — exactly like a plasma's free electrons, but transient. The critical numbers:

```
Pair lifetime:          Δt ≈ 6.4 × 10⁻²² s
Maximum displacement:   Δx ~ c·Δt = λ_C ≈ 2.43 × 10⁻¹² m  (electron Compton wavelength)
Effective pair density: n_virt ~ λ_C⁻³ ≈ 1.7 × 10³⁷ m⁻³
Vacuum plasma frequency: ωₚ = mₑc²/ℏ ≈ 7.8 × 10²⁰ Hz
```

Below `ωₚ` (all laboratory and optical EM), the pairs polarize the vacuum — the Euler-Heisenberg effective Lagrangian is precisely this mean-field response. Above `ωₚ` (photon energies > 1.022 MeV), virtual pairs become real — the Schwinger pair production regime. The K field is the classical effective description of this polarization, valid at scales >> λ_C.

**Natural K field dynamics** follow from treating virtual pairs as driven oscillators restored by annihilation at rate ωₚ, driven by local field energy density u:

```
∂²K/∂t² = c²∇²K − ωₚ²(K − 1) + η · u_field / u_S
```

where η ≈ 2α/(45π) ≈ 10⁻⁴ is the Euler-Heisenberg coupling, and u_S = mₑ²c⁴/(ℏ²α) is the Schwinger energy density scale. The term ωₚ²(K − 1) is the pair annihilation — K restores to 1 at the Compton rate when the driving field is removed.

**Weber corrections enhanced by K:** Weber's force IS the retarded potential expanded to order (v/c)². In a region with K > 1, the local speed of light is c_local = c₀/K, so the retardation parameter becomes v/c_local = Kv/c₀. Weber corrections scale as K²(v/c₀)². The threshold for like-charge attraction in standard Weber requires v > c₀√2 (apparently superluminal). In a polarized vacuum this threshold becomes v > c₀√2/K — achievable at sub-relativistic velocities for K > √2. **Note the domain-of-validity caveat in §5.3: this K-substitution is an approximation whose accuracy must be validated against the full retarded potential calculation in the simulation.**

**Exotic Vacuum Objects as self-sustaining K pockets:** A dense electron cluster creates a strong self-field → local K >> 1 → K-enhanced Weber terms become attractive → cluster binds additional electrons → stronger self-field → higher K. This is a bifurcation: below a threshold electron density, thermal dispersion wins; above it, the self-consistent K pocket is self-sustaining. The S field provides the longitudinal coupling that transfers energy between the cluster's orbital dynamics and its self-field, completing the binding mechanism. **Note the bootstrapping caveat in §5.3: the threshold density for this bifurcation must be estimated quantitatively to determine whether it is dynamically accessible.**

**The unifying picture:** K > 1 in any confined geometry mediates effective attraction through three structurally identical mechanisms:

| Geometry | K source | Attraction mediated |
|----------|----------|---------------------|
| Parallel plates, gap ≪ λ_C | Restricted pair coherence length | Plates toward each other (Casimir) |
| Near rapidly moving charges | Retarded field energy density | Like-charge Weber attraction (K-enhanced) |
| Dense electron cluster | Cluster's own electrostatic self-field | Cluster self-binding (EVO/plasmoid) |

All three are the same phenomenon — vacuum polarization K > 1 in a constrained or self-generated region mediating an effective attractive force — in different spatial configurations.

### 6. Topological Solitons in QVED: Hopfions and Ball Lightning

#### 6.1 Q Field Topology

Unit quaternions form the group S³ — the 3-sphere. The third homotopy group π₃(S³) = ℤ means smooth maps from physical 3-space (compactified to S³) to the S³ value space are classified by an integer winding number. This is precisely the Skyrmion number — QVED has Skyrmionic solutions built into its mathematical structure without any additional assumptions. A configuration with winding number n ≠ 0 cannot be continuously deformed to n = 0 (the uniform vacuum): it is topologically protected against small perturbations.

#### 6.2 Hopfions

The Hopf fibration describes how the 3-sphere fibers over the 2-sphere as S¹ bundles, with a remarkable property: every pair of distinct fibres is linked exactly once. Applied to EM fields, if the field direction lives on S² (normalized **B** = **B**/|**B**|), then the preimage of any two distinct field orientations forms a pair of closed loops that link exactly once. This is the Hopfion configuration — a field whose field lines are all mutually Hopf-linked. Irvine & Bouwmeester (2008) derived exact Hopfion solutions of Maxwell's equations. In linear EM, these solutions are exact but not stable: they disperse at speed c. The Hopf topology is correct but there is no nonlinearity to sustain it.

#### 6.3 Ball Lightning as a Hopfion in (Q, K)

The (Q, K) system is nonlinear: K evolves under field pressure, and K > 1 feeds back into local wave speed. This nonlinearity provides the topological protection that free-space Hopfions lack:

1. **Topology**: E and B field lines form a Hopf fibration — every E field line links through every B field line exactly once. This requires simultaneous toroidal and poloidal B components, naturally created by a torus-knot (helical) discharge geometry.

2. **K boundary**: The self-consistent K > 1 region at the Hopfion core acts as a partial photon mirror (n = √K → Fresnel reflection at the K gradient boundary). This explains the observed glow of ball lightning: photons in the core decay slowly through the partially reflective K boundary, producing sustained luminosity without a persistent current source.

3. **K as dynamical Skyrme term**: Compression → increased local field energy density → K rises via η coupling → c_local = c₀/K decreases → effective pressure resists further compression. This is structurally analogous to the Skyrme stabilization term in nuclear physics — K provides a dynamically generated, field-energy-dependent stiffness. **Important caveat**: whether the quantitative feedback is strong enough to stabilize anything at accessible field strengths is exactly what the simulation needs to determine. This mechanism is an analogy supported by structural similarity, not a derivation from first principles. The simulation must answer the quantitative question.

4. **Glass penetration**: The Hopfion is a field configuration, not a charged particle plasma. A thin dielectric window is transparent to a slowly varying field topology — the field reorganizes around the obstacle. This is why ball lightning reportedly passes through glass without breaking it: there are no particles, only a topological field structure.

5. **Formation from lightning**: Multi-stroke lightning creates both azimuthal B (from the return stroke current) and B looping around the channel axis (from branching discharge geometry). This multi-stroke field topology is a natural seed for the Hopf-linked configuration.

**Topological charge density.** The integer topological charge of the Q field is computed from the normalized quaternion U = Q/|Q| (projection of Q onto S³):

```
ρ_topo = (1/24π²) ε^{ijk} Tr( ∂_i U U⁻¹ · ∂_j U U⁻¹ · ∂_k U U⁻¹ )
```

where ε^{ijk} is the Levi-Civita tensor and the trace is over quaternion components. Integrated over all space, this gives an integer n ∈ ℤ (the Hopf invariant / Baryon number). In simulation, this is monitored as a conserved diagnostic — an integer-valued quantity that cannot change smoothly, confirming topological stability.

#### 6.4 Derivable Geometry for Topological Discharge Synthesis

The geometry and parameters for creating a Hopf-linked EM configuration in the laboratory are derivable from the Hopfion structure and the (Q, K) field equations:

**The key geometric insight:** A conventional toroid (winding angle 0°) creates only toroidal B. A through-hole coil creates only poloidal B. A helix wound at **45°** on the torus surface creates equal toroidal and poloidal components from a single winding — no synchronization between two independent discharge banks required. This winding traces a (1,1) torus knot on the torus surface.

**Torus geometry:**
```
Major radius:  R = 5–10 cm  (sets ball lightning scale)
Minor radius:  r = 0.3–0.4 × R  (optimal r/R for Hopf fibration)
Winding angle: θ = 45°  (equal B_tor and B_pol → full Hopf symmetry)
Turns:         N = 20–40  (single helical winding, one bank)
```

**Discharge parameters (R = 8 cm, r = 3 cm, N = 30 turns at 45°):**
```
Toroid inductance:   L ~ μ₀ N² r² / (2R)  ≈  4 μH
Target B field:      B ~ 1–5 T             (energy density ~ 0.4–10 MJ/m³)
Peak current:        I ~ B × 2πR / (μ₀N)  ~ 10–50 kA
Bank energy:         E = ½ L I²            ~ 0.5–5 kJ
Capacitor bank:      C = 2E / V₀²; at V₀ = 10 kV: C ~ 10–100 μF
Rise time:           τ_rise ~ π√(LC)       ~ 0.3–1 μs
LC ring-down:        τ_LC = 2π√(LC)        ~ 1–5 μs  (sets topology dwell time)
```

**Cost estimate:** ~$500–2,000 in hardware (polycarbonate torus former, magnet wire, film/electrolytic capacitors rated >10 kV, spark gap or thyratron switch, oscilloscope with B-dot probes).

**Critical experiment:** Does the luminous toroidal plasma structure survive significantly longer than τ_LC after the discharge ends? The LC circuit stores energy only for τ_LC ~ 5 μs. If the EM structure persists for milliseconds, a self-sustaining K feedback loop is the only physical explanation within this framework.

**Success indicators:**
1. Luminous toroidal structure persisting > 100 × τ_LC after energy input ceases
2. B-dot probe decay slower than predicted LC ring-down (topological field energy storage)
3. Detachment of a mobile luminous structure from the discharge region

**The simulation's predictive role:** Phase 1.8 (K dynamics) can determine the minimum field energy density threshold for the K self-sustaining feedback loop to engage, and predict whether laboratory-accessible field strengths (~1–5 T) are sufficient to reach the critical K regime. This provides a theoretically grounded go/no-go prediction before hardware construction.

### 7. The Spheromak Connection: Fusion Research and Unexplored Territory

#### 7.1 Magnetic Helicity and Taylor Relaxation

The plasma physics community has been creating Hopf-topology field structures in the laboratory since the 1970s, under the name *spheromak*, without framing them topologically. The key quantity is **magnetic helicity**:

```
H_mag = ∫ A · B  dV
```

This integral measures the linking and knotting of magnetic field lines — it is the classical electromagnetic analogue of the topological charge from §6.1. In ideal magnetohydrodynamics (no resistivity), H_mag is exactly conserved. A plasma with nonzero magnetic helicity cannot evolve to a state with zero magnetic helicity; it is topologically constrained.

**J.B. Taylor (1974)** showed that when a turbulent, resistive plasma relaxes, it conserves global magnetic helicity even while dissipating energy. The minimum-energy state of a plasma with fixed magnetic helicity is the **spheromak** — a self-organized configuration with equal toroidal and poloidal magnetic fields, B_tor = B_pol, exactly the 45° field balance of the Hopfion. The plasma finds this configuration spontaneously through turbulent relaxation, without any engineered geometry. This is Taylor relaxation: the physical analogue of a dynamical system finding its attractor.

#### 7.2 The Spheromak as Hopfion-in-Plasma

A spheromak is a plasma-supported Hopfion. The plasma currents act as the source maintaining the Hopf-linked field configuration. Remove the plasma — either by cooling (recombination) or by the plasma escaping the confinement region — and the field configuration collapses in standard MHD. The QVED hypothesis introduces a new question that plasma physics has never asked: **if K feedback can sustain the field topology after the plasma recombines, does the spheromak transition to a field-only Hopfion?** The fusion community has always ended the experiment when the plasma is gone. QVED says: look at what happens in the microseconds after recombination.

#### 7.3 How Spheromaks Are Created Non-Destructively

A **coaxial plasma gun** (Marshall gun) injects magnetic helicity into a plasma without destroying the electrodes. Two coaxial cylindrical electrodes, a magnetic bias coil, and a gas puff valve. The gas ionizes, helicity is injected by J×B torques on the plasma, and the resulting spheromak detaches from the gun and travels into the chamber as a self-contained, self-organized field structure. This is repeatable thousands of times. Small-scale plasma guns have been built in university labs for $5,000–$20,000.

### 8. Non-Destructive Experimental Ladder

The following progression tests the Hopfion stabilization mechanism at increasing scales and field strengths, with each stage informing the next. The critical design principle is that each stage has a clear null hypothesis, a clear positive signal, and costs less than the next stage — so negative results save resources while positive results justify escalation.

**Stage 1 — Optical Hopfion in free space (topology creation, confirmed):**
Irvine & Bouwmeester (2008) demonstrated that an exact Hopfion solution of Maxwell's equations can be created in a laser beam using a spatial light modulator (SLM) with the appropriate holographic phase pattern. This is reproducible in any nonlinear optics laboratory. It confirms the topology can be created, measured, and verified to have the correct linking structure. Cost: access to an optics lab with a spatial light modulator (~$20k–$80k instrument, often already present in university labs). Result: topology confirmed, decay rate measured.

**Stage 2 — Optical Hopfion in a Kerr medium (K stabilization test):**
The Kerr nonlinearity of an optical medium is n = n₀ + n₂|E|², where the refractive index increases with field intensity. This is structurally identical to the K field: K = 1 + η·u/u_S, with higher field energy raising local n. Create the same optical Hopfion in a Kerr medium (CS₂ liquid, a photorefractive crystal, or a fiber with strong nonlinear coefficient). Measure whether the decay rate is slower than in air. If yes: the K-analogue mechanism stabilizes the Hopfion topology at optical scales. This is a direct, controlled, non-destructive test of the core stabilization physics — cleanly separated from all plasma and high-voltage complications. Cost: same optics lab plus a cuvette of CS₂ (~$50). This experiment could be published as a direct test of the QVED Hopfion model.

**Stage 3 — Superconducting helical torus (persistent knotted field):**
Wind the 45° helical torus coil in YBCO high-temperature superconductor tape (operates at 77 K, liquid nitrogen). Once the persistent current is established (by transformer coupling or direct charging), it circulates indefinitely without Ohmic loss. The knotted B field simply exists, stably. The coil can be quenched in a controlled way (warming it past T_c) to observe what happens to the field topology at the moment the persistent current ceases: does the topology survive the quench? This is the cleanest possible test of field-topology detachment from its current source. Cost: ~$10k–$30k (YBCO tape, small LN₂ dewar, cryostat). No high voltage, no ablation, fully repeatable.

**Stage 4 — Resonantly driven helical torus at low current (calibration):**
Wind the helical torus coil in heavy copper wire (4 AWG, safe up to ~1000 A for 10 ms pulses). Drive it at LC resonance, with energy recirculating rather than dumped. At 500–1000 A peak, B ~ 50–200 mT in the torus core — far below ablation threshold, far below K threshold. Purpose: calibrate the topology creation geometry, verify that B-dot probes correctly measure the expected toroidal/poloidal field ratio of 1:1 (confirming the helix creates the intended topology), and establish a measurement baseline for τ_LC ring-down. Cost: ~$500–$2,000. No plasma, no ablation, fully safe benchtop work.

**Stage 5 — Coaxial plasma gun / spheromak (helicity-injected topology, non-destructive):**
Build or acquire a small coaxial plasma gun. Establish spheromak formation (B-dot probe array, interferometry). Then focus measurement on the post-recombination phase: what does the field do in the 10–1000 μs after the plasma has cooled? Standard measurement stops before this window. A carefully arranged B-dot probe array with nanosecond-resolution ADC can capture any anomalous field persistence or detachment event. Cost: $5,000–$20,000. This is established plasma physics technique with a novel diagnostic focus.

**Stage 6 — Pulsed helical torus in dense gas (controlled ablation with topology seeding):**
The original proposal from §6.4, refined: rather than ablating the wire in air, fire the discharge in a controlled atmosphere — either a noble gas (argon, higher breakdown threshold, cleaner plasma chemistry) or in a partial vacuum. The ablating wire plasma is contained inside a quartz tube shaped as the torus. The knotted plasma can expand into the interior volume of the quartz torus, attempting to detach. The quartz tube serves as a controlled plasma vessel. B-dot probes are positioned outside the torus to detect any field that detaches outward. Cost: ~$2,000–$8,000 (existing capacitor bank from Stage 4 plus quartz torus tube plus argon supply). This is the highest-energy, highest-risk stage and should follow the calibration work of Stages 4 and 5.

### 9. The Fusion Energy Connection

#### 9.1 Fusion Confinement as a Topological Problem

The spheromak was invented not as a fundamental physics curiosity but as a potential path to fusion energy. The fundamental problem of fusion confinement — holding a 100-million-degree plasma long enough for net fusion power — is a topological problem: how do you create magnetic field structures that prevent plasma from leaking to the walls?

Tokamaks, stellarators, spheromaks, and field-reversed configurations (FRCs) are all variations on the same theme: create topologically nontrivial magnetic field configurations (nonzero helicity) that confine plasma by winding field lines around the plasma core. Decades of progress have improved confinement enormously. But the plasma still leaks, through turbulent cross-field transport that standard MHD theory cannot fully suppress.

#### 9.2 The Unexplored QVED Angle

In a fusion plasma at 10⁸ K with magnetic fields of 5–10 T, the direct QED vacuum K coupling (η ~ 10⁻⁴) is negligible — the Schwinger field E_S ~ 10¹⁸ V/m is unreachable. However, the **plasma itself** acts as a polarizable medium with its own effective dielectric response. For EM waves at frequencies below the plasma frequency ω_pe, the plasma has ε_r >> 1 — the plasma K is large and real. This plasma-K plays the same structural role as the vacuum K in the QVED equations: it modifies c_local, creates a refractive boundary, and alters the force balance. **See §5.3 for important caveats about the limits of this analogy.**

Nobody in fusion has modeled Taylor relaxation with a **variable effective K** — one that responds to local field energy density through the plasma dielectric. The QVED framework predicts that the minimum-energy state at fixed helicity changes when K is dynamical: the field equilibrium, the confinement geometry, and the transport properties are all modified. Specifically:

1. **K-gradient confinement**: A K > 1 region inside the plasma creates a refractive index boundary. Photons (radiation losses, a major energy drain in fusion) are partially reflected back into the plasma. This is a confinement mechanism with no analogue in standard MHD.

2. **K-modified Taylor state**: If K is allowed to vary, the minimum-energy state at fixed helicity is no longer the standard force-free spheromak (∇×B = λB). It is a modified state where the K gradient contributes to the force balance. This state might have better confinement properties than the standard Taylor state.

3. **S field energy channel**: The scalar longitudinal mode S, absent from standard MHD, provides an additional energy coupling pathway. In QVED extended mode, energy can flow through the S field on timescales and length scales different from transverse EM radiation. Whether this helps or hurts confinement depends on the geometry and requires simulation.

The simulation built in this project — once Phases 1.8 (dynamic K) and 2.1 (extended S) are complete — is directly applicable to modeling a small spheromak with QVED corrections. The grid resolution needed is modest (the spheromak physics is on scales >> λ_C). This is an untouched research direction that could yield novel confinement predictions testable in existing plasma physics facilities.

#### 9.3 Why the Fusion Community Has Not Explored This

The polarizable vacuum model sits outside mainstream theoretical physics — it is neither string theory nor lattice QCD nor standard model extensions, and it does not fit neatly into any funded research program. Plasma physicists do not read papers on quaternionic electrodynamics. The QED corrections to vacuum permittivity are universally dismissed as negligible in laboratory-accessible field regimes — correctly, for direct QED coupling. The key insight this framework offers — that the *structural mechanism* of K > 1 polarization and K-gradient confinement can be realized by plasma dielectric response at much lower field strengths than the Schwinger threshold — has not been articulated in the fusion literature.

This is a genuinely open theoretical territory. The QVED simulator, built to explore vacuum polarizability, turns out to be exactly the right tool for exploring whether dynamical-K confinement could contribute to solving the plasma confinement problem — because the field equations are isomorphic between vacuum K and plasma effective-K at the level of the continuum field theory.

### 10. Energy Considerations

This section addresses the question of energy generation and storage honestly, distinguishing between what is physically coherent, what is speculative, and what violates known physics.

#### 10.1 Topological Energy Storage

If a self-sustaining topological EM structure (a Hopfion stabilized by K feedback) can be created, it stores energy in the electromagnetic field configuration: E_stored = ∫ (ε₀E²/2 + B²/2μ₀) dV. For a 10 cm sphere with B ~ 1 T, this is roughly 1-2 kJ. This is intrinsically useful as an energy storage medium — like a battery made of pure field. The energy stored is the energy input during formation, minus formation losses. No thermodynamic laws are violated. The efficiency is determined by formation losses and extraction losses.

The energy budget for laboratory creation is accessible. Ball lightning observations suggest the phenomenon contains 10 J to 10 kJ. A capacitor bank of a few microfarads at 10-20 kV stores a few kilojoules — within range of benchtop equipment. The question is not whether enough energy can be delivered but whether it can be delivered in the right topology (see §6.4).

#### 10.2 Vacuum Energy Extraction — Honest Assessment

The proposition that energy could be extracted from the electromagnetic vacuum is a much harder claim. The vacuum does contain energy in the QED framework — zero-point fluctuations are real, and the Casimir effect demonstrates that boundary conditions can extract work from the vacuum. But the Casimir effect does not produce net free energy: work must be done to separate the plates again, and the round-trip is thermodynamically closed.

The PV framework suggests that driving K away from 1 and letting it relax back releases energy. But the energy to drive K away from 1 has to come from somewhere (the capacitor bank, the laser, etc.). If the vacuum is in its ground state (the lowest energy state), energy cannot be extracted from it by definition — there is no lower state to relax to. If the vacuum is in a metastable state (a local but not global energy minimum), then in principle a perturbation could trigger relaxation to the true ground state, releasing energy. Whether the QED vacuum is the true ground state or a metastable state is unknown and is connected to deep open questions about dark energy and the cosmological constant.

The honest framing: vacuum energy extraction is a possibility that the framework allows in principle but that requires specific conditions (metastable vacuum state, accessible transition path) that we do not currently know to exist. The simulation can model the energy balance of vacuum polarization and relaxation events. If the simulation shows energy being "created," the first conclusion should be a numerical error or an unphysical initial condition, not a discovery. Genuine energy extraction would require extraordinary evidence and extraordinary scrutiny.

#### 10.3 The Most Realistic Energy-Adjacent Outcome

The most exciting realistic outcome is not energy generation but rather the demonstration that electromagnetic fields in vacuum can form self-organized, topologically protected structures at macroscopic scales. That would be a fundamental discovery. Understanding the stabilization mechanism for such structures opens doors to engineered vacuum-field interactions whose applications cannot be fully predicted today.

Historical precedent: the discovery of superconductivity (1911) was a fundamental physics discovery about electron pairing and quantum coherence. It took decades to become technologically useful, but it eventually enabled MRI machines, particle accelerators, quantum computers, and other technologies that were unimaginable in 1911. Similarly, understanding self-sustaining topological EM structures — if they exist — would be a discovery whose applications unfold over time.

Storing and releasing conventionally sourced energy in topological field structures is a nearer-term and more defensible application than vacuum energy extraction. A compact, mobile, long-lived electromagnetic energy storage device is valuable even if the energy it stores comes from a conventional source.

### 11. On the Sociology of Unconventional EM Research

#### 11.1 Why Mainstream Academia Has Blind Spots Here

Academic physics has strong institutional incentives against exploring foundations of electrodynamics. Tenure committees do not reward it. Grant agencies do not fund it. The prevailing attitude is "Maxwell's equations work, QED is the most precisely tested theory in history, why would you question the foundations?" This is understandable but it does create genuine blind spots, particularly in the near-field, pulsed, and topological regimes where the standard framework has been least thoroughly tested.

The gap between "well-tested in the far-field regime" and "well-tested in all regimes" is real but sociologically invisible. Pointing this out is not anti-scientific — it is precisely the kind of observation that motivates new experiments.

#### 11.2 The Classification Question

It is historically documented that governments have funded research into unconventional electromagnetic phenomena. The US Department of Defense funded studies through the Advanced Aerospace Threat Identification Program (AATIP), which included contracted papers on topics like vacuum energy and polarizable vacuum models (Puthoff was involved in some of this work). These facts are documented through FOIA requests and public reporting.

However, "some of this research moved into classified programs" is a hypothesis that, by its nature, cannot be built upon. If the knowledge is locked away, it cannot be accessed, so it does not help this project. More importantly, the reasoning pattern of "the absence of mainstream research proves it was suppressed" is epistemically dangerous because it is unfalsifiable — it explains both the presence and absence of evidence equally well.

**The productive framing:** Mainstream academia has not explored these directions thoroughly, for a mix of sociological and theoretical reasons, and that creates an opportunity for anyone willing to do the work carefully and publish openly. This project operates entirely in the open, using published physics and publicly available information.

### 12. The Unified Framework

Bringing these threads together:

**Quaternionic Vacuum Electrodynamics** is an electromagnetic theory in which:

1. The fundamental quantity is the quaternionic potential Q = φ/c + **A**, representing both electromagnetic and vacuum states.
2. The vacuum is a polarizable medium with dynamical permittivity ε(x,t) and permeability μ(x,t).
3. No gauge condition is imposed. The scalar field S = (1/c²)∂φ/∂t + ∇·**A** is a physical degree of freedom representing the longitudinal/scalar mode of the EM field and the local vacuum polarization state.
4. The theory reduces to standard Maxwell electrodynamics when S = 0 and K = 1 (homogeneous, unpolarized vacuum with Lorenz gauge).
5. The theory reproduces Weber's longitudinal forces in the near-field, low-velocity limit for discrete charges.
6. The theory reproduces PV-model predictions for vacuum polarization in the presence of mass-energy.
7. Transverse EM waves (light) and longitudinal EM waves (scalar potential waves) are both solutions, with different coupling mechanisms and potentially different propagation characteristics in polarized vacuum regions.

The physical status of this framework depends on the resolution of the gauge symmetry question (§2). If gauge symmetry is exact, S = 0 is automatic and the framework reduces to standard EM with a variable-permittivity vacuum model. If gauge symmetry is approximate or emergent, S is a genuine physical field and the framework predicts new phenomena. The simulation is designed to determine, quantitatively, how large the predicted new effects are for experimentally accessible configurations — providing either measurable experimental targets or a demonstration that the effects are too small to detect with current technology.

---

## Implementation Status

A recurring question this README invites: "what does the simulator actually do today?" The catalog below documents the *physics targets*; this section documents how far down each target the implementation has reached. Status keys: ✅ implemented and tested · ⚠️ scenario implemented with deferred follow-on · 🚧 partial / pseudocode-only stub · ❌ not started.

| Scenario / Capability | Status | Notes |
|---|---|---|
| Standard FDTD wave propagation | ✅ | 32³ default grid, CFL-bounded leapfrog, CPML and Mur ABC available |
| Polarisable vacuum K-field | ✅ | Phase 1.8 K dynamics evolving via VacuumConfig |
| Topological charge diagnostic | ✅ | Skyrmion / Baryon-number invariant, 1.9 |
| Probe + FFT time-series | ✅ | Naive O(N²) DFT for FFT, ring-buffered probe history |
| Dipole radiation | ✅ | Validates standard EM in the K=1, S=0 limit |
| Bifilar coil + scalar wave | ✅ | Phase 2.3, demonstrates B-cancellation + S-mode excitation |
| Bifilar pair (Tx/Rx) | ✅ | Phase 2.4, scalar-wave propagation delay measurement |
| Vacuum K (polarisable demo) | ✅ | Phase 1.8 standalone scenario |
| Graneau wire (longitudinal Weber) | ✅ | Phase 3.3, end-peaked profile |
| Hairpin PCB trace | ✅ | Phase 7.3, antiparallel-strand modification verified |
| Aharonov-Bohm (solenoid) | ✅ | Phase 7.1, A ≠ 0 where B = 0; QVED ↔ standard EM equivalence verified |
| Toroidal AB (donut coil) | ✅ | Phase 7.2, linking flux ≈ NIr²/(2R); driven-AC mode deferred |
| Hopfion / ball lightning | ⚠️ | Phase 7.5, hedgehog Skyrmion ansatz; full Rañada-Irvine IC and η-sweep deferred |
| Brown / Biefeld asymmetric capacitor | 🚧 | Tier-3 stub with full pseudocode; needs conductor BC + Maxwell-stress integrator |
| Pulsed-circuit scalar excitation | 🚧 | Tier-3 stub; needs switched-current source variant |
| Charge cluster (EVO) stability | 🚧 | Tier-3 stub; needs seeded RNG + RMS-radius diagnostic |
| K-cycle resonator | ❌ | Phase 7.6, needs K-equation drive in `step_field_cpu` |
| Spheromak / Taylor with K | ❌ | Phase 7.7, needs `compute_magnetic_helicity` + Chandrasekhar-Kendall IC |
| Casimir | ❌ | Phase 4.2, needs conductor BC + Maxwell-stress integrator |
| Vacuum-modified wave propagation (lensing) | ❌ | Phase 4.3 |
| GPU migration | ❌ | Phase 5; CPU baseline performant enough at 32³ |
| Volume rendering (slabs) | ✅ | Phase 6.1 CPU prototype; GPU ray-march port deferred |
| RK4 streamlines (primary vector view) | ✅ | Phase 6.2 |
| Slice 2D inset | ✅ | Phase 6.0; 3D textured quad demoted, opt-in via `show_in_3d` |

**Recurring caveat.** Verifications throughout use **code units** (μ₀ = 1, ε₀ = 1). The qualitative physics — flux quantisation, end-peaked Graneau profile, S ≠ 0 in extended mode, etc. — is preserved; absolute force / flux magnitudes are unphysical without scaling factors. See per-scenario module headers for scaling notes.

**Verified equivalence claim (Phases 7.1, 7.2, 7.7).** For configurations with φ = 0 and ∇·A = 0 at t = 0 ("gauge-clean initial conditions"), QVED extended-mode evolution reduces to standard-mode evolution to floating-point precision. The QVED extension introduces no new physics for these configurations; the predictions diverge only when sources or initial conditions break gauge cleanliness.

---

## Phenomena Catalog

### Tier 1 — Experimentally Established, Theoretically Robust (Simulation Self-Validation)

These phenomena have solid experimental evidence and are explained by standard physics. Simulating them validates the simulation itself. The simulation must reproduce these results before any claims about extended framework predictions are credible.

#### 1.1 Standard EM Wave Propagation
- **What:** Transverse electromagnetic waves propagating at c, dipole radiation patterns, waveguide modes.
- **Status:** The bedrock of electrodynamics.
- **Simulation role:** If the simulation cannot reproduce these with high accuracy in the K = 1, S = 0 limit, it is broken. This is the first and most important validation target.

#### 1.2 Aharonov-Bohm Effect
- **What:** A charged particle's quantum phase is affected by the vector potential **A** even in field-free regions (**E** = **B** = 0).
- **Status:** Experimentally confirmed beyond doubt (Tonomura 1986, Webb 1985).
- **Quaternionic relevance:** Direct evidence that potentials are physically primary. The simulation should reproduce the phase accumulation ∮**A**·d**l** for specific geometries and show that the potentials carry information that the fields do not.
- **Simulation target:** Solenoid geometry with confined B field; demonstrate nonzero A outside solenoid; compute phase differences for different paths.

#### 1.3 Casimir Effect
- **What:** Attractive force between parallel conducting plates in vacuum, measured to high precision.
- **Status:** Experimentally confirmed; quantitative agreement with theory.
- **Quaternionic relevance:** In the PV framework, this is vacuum polarization between boundary conditions. The simulation should reproduce the force law F/A = -π²ℏc/240d⁴ from the ε(x) modification between plates.
- **Simulation target:** Two conducting boundary conditions in the vacuum grid; observe K(x) profile between plates; compute resulting force.

### Tier 2 — Experimentally Observed, Interpretation Contested

These phenomena have real experimental observations but their interpretation within either the standard or quaternionic framework is debated. Honest treatment requires acknowledging the contested nature of the evidence.

#### 2.1 Graneau's Exploding Wire/Water Arc Experiments
- **What:** High-current pulses through water jets and mercury pools produce longitudinal fragmentation inconsistent with Lorentz forces.
- **Status:** Experimentally reproduced by some groups; interpretation is genuinely contested. Mainstream explanation invokes magnetohydrodynamic instabilities rather than Weber forces. Independent replication has been limited. **This project treats the phenomena as real but the Weber interpretation as a hypothesis to be tested, not an established fact.**
- **Quaternionic relevance:** If the Weber interpretation is correct, this is a direct manifestation of longitudinal forces. The simulation should model high-current pulsed configurations and compare Weber and Lorentz force predictions, looking for configurations where the two frameworks make maximally different predictions.
- **Simulation target:** Model linear conductor with pulse current; compute both Weber and Lorentz force distributions; identify the geometry and current regime where the predictions diverge most — the "first divergence threshold."

#### 2.2 Scalar Longitudinal EM Waves
- **What:** Propagating compression waves in the electromagnetic potential, as distinct from transverse EM waves.
- **Status:** Claimed detection by several groups (Monstein & Wesley 2002; Meyl's experiments); not widely replicated or accepted.
- **Quaternionic relevance:** Core prediction of the non-gauge-fixed quaternionic theory. If S ≠ 0 propagating solutions exist in the simulation and correspond to physical degrees of freedom, this is the signature prediction of the framework.
- **Simulation target:** Drive oscillating charge configurations that preferentially excite the S field; observe propagation characteristics; identify antenna geometries that couple to S rather than to transverse modes.
- **Experimental target:** Bifilar coils (B fields cancel, A field survives) as transmitter/receiver pairs; look for coupling that cannot be explained by E, B radiation.

#### 2.3 Vacuum Polarization Gradients / Engineered Casimir Effects
- **What:** Using strong EM field configurations to locally modify vacuum permittivity and create forces from the resulting gradients.
- **Status:** Theoretical within PV model; no claimed experimental demonstration at accessible scales.
- **Quaternionic relevance:** If ε(x) is dynamical and couples to EM fields, then sufficiently intense field configurations should produce measurable permittivity changes. This is the "vacuum engineering" application.
- **Simulation target:** Model high-intensity EM field regions; compute induced K(x) profile; look for self-consistent solutions where the field modifies the vacuum which modifies the field propagation. **Critically: report the actual magnitude of predicted K deviations for laboratory-accessible field strengths.**

#### 2.4 Weber Force Anomalies in Circuit Design
- **What:** Circuit geometries where Weber longitudinal forces produce measurable effects — torques, stresses, or current distributions that differ from Lorentz predictions.
- **Status:** Predicted by Weber theory; partially supported by Graneau's work; not systematically explored in PCB-scale circuits.
- **Quaternionic relevance:** Most accessible experimental target. PCB traces carrying pulsed high current in specific geometries should experience forces that standard circuit theory doesn't predict.
- **Simulation target:** Design specific PCB trace geometries; model pulse current; compute Weber force map; identify geometries that maximize the anomalous longitudinal force.
- **Experimental target:** PCB with strain gauges or displacement sensors on traces designed to maximize Weber forces under pulsed current.

#### 2.5 Aharonov-Bohm Analogue Effects in Macroscopic Circuits
- **What:** Electromagnetic coupling between circuits mediated by the vector potential in regions where B = 0.
- **Status:** Quantum AB effect is confirmed; macroscopic classical analogue is predicted by the quaternionic framework but not systematically tested.
- **Quaternionic relevance:** If potentials are primary, then a toroidal solenoid (which confines B entirely inside the torus) should still produce measurable effects on external circuits through the external A field.
- **Simulation target:** Toroidal coil geometry; compute A field external to torus; model coupling to a secondary circuit; quantify predicted signal.
- **Experimental target:** Toroidal coil driven with AC; search for induced signal in an external loop that cannot be explained by stray B field leakage.

### Tier 3 — Speculative / Extraordinary Claims

These phenomena have limited or anecdotal evidence. They are simulation targets for after the framework is validated against Tier 1 and 2 phenomena. Results from Tier 3 simulations should be treated as predictions to be tested, not as evidence for the framework. **The simulation's primary value in Tier 3 is producing quantitative predictions — specific numbers that can be compared against specific measurements.**

#### 3.1 Ball Lightning as a Topological Hopfion
- **What:** Ball lightning as a self-sustaining topological soliton (Hopfion) in the (Q, K) field system.
- **Historical context:** Reliably reported by multiple independent witnesses including scientists. Duration up to minutes, glass penetration, spherical shape, spontaneous formation near lightning. No conventional plasma model produces stable, mobile, luminous spheres on these timescales.
- **Quaternionic connection:** See §6.3. The Q field lives on S³ with π₃(S³) = ℤ, giving Skyrmionic solutions. The Hopfion requires K nonlinearity for stability. K > 1 core creates partial photon reflection (luminosity), and K acts as stabilizing term (structural stability). Field soliton character (not particles) explains glass penetration.
- **Simulation target:** Hopfion initial condition in the (Q, K) grid; monitor topological charge integral over time; search for parameter regime where the Hopfion decays vs. self-sustains based on η and ωₚ; identify critical K threshold. Phase 1.9 (topological charge diagnostic) is a prerequisite.
- **Experimental target:** Helical (45°) torus coil discharge geometry; see §6.4 for full parameter derivation and §8 for the experimental ladder.

#### 3.2 T.T. Brown Electrogravitics
- **What:** Claimed thrust produced by high-voltage asymmetric capacitors, apparently in the direction of the positive plate, persisting (per Brown's claims) even in vacuum.
- **Historical context:** Brown filed patents in the 1920s-60s. Project Winterhaven (1952) proposed military applications. The mainstream explanation is ion wind in air; Brown claimed the effect persisted in vacuum.
- **Quaternionic connection:** In the PV model, a strong electric field could locally polarize the vacuum (modify K). An asymmetric capacitor would create an asymmetric K gradient, potentially producing a net force. The simulation can test whether realistic field strengths produce a meaningful K modification.
- **Simulation target:** Asymmetric capacitor geometry with realistic voltages (50-300kV); compute K(x) profile; determine if gradient produces a net force and its magnitude. **Expected outcome: the K modification at achievable voltages is likely immeasurably small (see §5.3). This is itself an important result.**

#### 3.3 "Cold Electricity" / Radiant Energy (Tesla, Bedini, Gray)
- **What:** Various claims of anomalous electrical phenomena in pulsed circuits: current without conventional electron flow, anomalous battery charging, COP > 1 energy systems.
- **Historical context:** Tesla described "radiant energy" from abruptly interrupted DC currents. Edwin Gray patented a pulsed motor. John Bedini built pulsed battery charging circuits.
- **Quaternionic connection:** If the scalar field S is a real propagating mode, then abruptly interrupted current (which produces large ∂φ/∂t and rapid changes in ∇·A) might preferentially excite S-mode propagation in the circuit. The "cold electricity" reports — current that doesn't heat wires, charges batteries in unusual ways — might be longitudinal potential waves rather than conventional transverse current flow.
- **Simulation target:** Model a pulsed circuit (sharp current interruption into an inductive load); analyze the S field excitation; determine whether longitudinal modes carry energy through the circuit differently than conventional current. **Note: every claimed COP > 1 energy system has either failed independent replication or been traced to measurement error. The simulation should focus on whether S-mode propagation exists and has distinctive signatures, not on whether it enables energy creation.**

#### 3.4 Exotic Vacuum Objects (EVOs) / Charge Clusters
- **What:** Claimed stable clusters of like-charged particles (electrons) that should repel but apparently form bound states.
- **Historical context:** Ken Shoulders' work on "charge clusters" (1980s-90s); Hal Puthoff's theoretical connection to vacuum engineering.
- **Quaternionic connection:** Weber acceleration-dependent term provides an attractive force; K-enhanced vacuum polarization could lower the Weber attraction threshold; topological stability (skyrmion-like configurations in the vacuum potential) could further stabilize such objects. **See §5.3 for the bootstrapping problem with this mechanism.**
- **Simulation target:** Model a dense cluster of like charges in the PV vacuum; include Weber acceleration terms; search for stable or metastable configurations; analyze the vacuum polarization profile around the cluster; **specifically estimate the threshold density required for the K feedback loop to engage and whether any known mechanism could produce it.**

#### 3.5 Vacuum Energy Extraction
- **What:** The proposition that the electromagnetic vacuum contains energy that could in principle be extracted.
- **Quaternionic connection:** See §10.2 for a detailed honest assessment. If the vacuum state can be driven away from equilibrium (K ≠ 1), the relaxation back to equilibrium would release energy — but only if the starting state is metastable, not if it is the true ground state.
- **Simulation target:** Model a strongly polarized vacuum region; allow it to relax; track energy flow through S field; determine whether energy is genuinely released or merely redistributed. **If the simulation shows net energy creation, treat this as a bug, not a discovery, until all numerical and initial-condition artifacts are exhaustively ruled out.**

---

## The First Divergence Threshold

The single most important output of the simulation is what we call the **first divergence threshold**: for each experimental configuration, the minimum field strength, current, geometry, or energy at which the quaternionic framework predicts a measurably different result from standard Maxwell electrodynamics.

For every simulation run, the output should include:

1. **The standard EM prediction** (S = 0, K = 1 limit of the same simulation)
2. **The quaternionic prediction** (full framework, no gauge fixing)
3. **The magnitude of the difference** in physically measurable quantities (force, voltage, field amplitude, frequency, phase)
4. **The experimental sensitivity required** to detect the difference

If the difference is 10⁻²⁰ of the standard prediction, that configuration is not experimentally accessible with current technology and should be deprioritized. If the difference is 10⁻³ or larger, it becomes a concrete experimental target. The simulation's job is to sort all proposed configurations by this ratio, from most accessible to least, creating a prioritized experimental program.

This approach protects against two failure modes: wasting resources on experiments that cannot possibly detect the predicted effect, and dismissing the framework without having tested it in the configurations where it makes its strongest predictions.

---

## Research Program

### Phase 1: Framework Validation
- Implement the quaternionic FDTD solver with PV vacuum model
- Reproduce standard EM results (transverse wave propagation, dipole radiation) as K=1, S=0 limit — **this must succeed before any extended framework simulation is trusted**
- Simulate Aharonov-Bohm geometry; verify potential structure
- Simulate Casimir geometry in PV model; verify force law
- Simulate Graneau-type configuration; compare Weber and Lorentz force predictions
- For every configuration: compute the first divergence threshold

### Phase 2: Novel Prediction Exploration
- Search for scalar wave propagation in specific antenna geometries
- Design and simulate bifilar coil transmitter/receiver for scalar mode coupling
- Model vacuum polarization feedback in high-field configurations
- Design PCB geometries optimized for Weber force detection
- Simulate toroidal AB-effect circuit coupling
- **Rank all configurations by first divergence threshold; publish the ranking**

### Phase 3: Experimental Design
- Select most promising simulation results (largest first divergence threshold) for experimental validation
- Design circuits and coil geometries using commercially available components
- Specify measurement requirements (sensitivity, bandwidth, shielding)
- Document experimental protocols for reproducibility
- Design experiments with clear null hypotheses

### Phase 4: Speculative Phenomena
- Model Hopfion stability in (Q, K) system; determine K threshold for self-sustenance
- Model asymmetric capacitor in PV vacuum (Brown electrogravitics)
- Simulate pulsed circuit scalar mode excitation
- Search for stable charge cluster configurations (EVOs)
- Explore vacuum energy dynamics
- **For all Tier 3 simulations: treat results as predictions, not confirmations. Report all numbers including disappointingly small ones.**

---

## Key References

### Primary Theoretical Sources
- Maxwell, J.C. (1865). "A Dynamical Theory of the Electromagnetic Field." — Original quaternion formulation
- Weber, W. (1846). "Elektrodynamische Maassbestimmungen" — Original Weber force law
- Assis, A.K.T. (1994). "Weber's Electrodynamics." Kluwer. — Modern treatment of Weber theory
- Puthoff, H.E. (2002). "Polarizable-Vacuum (PV) Approach to General Relativity." Found. Phys. 32, 927-943.
- van Vlaenderen, K.J. & Waser, A. (2001). "Generalization of Classical Electrodynamics to Admit a Scalar Field and Longitudinal Waves." Hadronic Journal 24, 609-628.

### Gauge Symmetry and Foundations
- Aharonov, Y. & Bohm, D. (1959). "Significance of Electromagnetic Potentials in the Quantum Theory." Phys. Rev. 115, 485-491. — Original AB effect prediction
- Healey, R. (2007). "Gauging What's Real." Oxford University Press. — Philosophical analysis of gauge symmetry and ontology
- Belot, G. (1998). "Understanding Electromagnetism." British Journal for the Philosophy of Science 49, 531-555. — Substantive treatment of the gauge interpretation debate

### Experimental References
- Tonomura, A. et al. (1986). "Evidence for Aharonov-Bohm Effect with Magnetic Field Completely Shielded from Electron Wave." Phys. Rev. Lett. 56, 792.
- Graneau, P. (1984). "Electromagnetic Jet-Propulsion in the Direction of Current Flow." Nature 295, 311.
- Graneau, P. & Graneau, N. (1996). "Newtonian Electrodynamics." World Scientific.
- Shoulders, K. (1991). "EV: A Tale of Discovery." — Charge cluster observations

### Plasma Physics and Fusion Confinement
- Taylor, J.B. (1974). "Relaxation of Toroidal Plasma and Generation of Reverse Magnetic Fields." Phys. Rev. Lett. 33, 1139. — Magnetic helicity conservation; derivation of the spheromak as minimum-energy state at fixed helicity
- Bellan, P.M. (2000). "Spheromaks: A Practical Application of Magnetohydrodynamic Dynamos and Plasma Self-Organization." Imperial College Press. — Comprehensive spheromak physics
- Marshall, J. (1960). "Performance of a Hydromagnetic Plasma Gun." Phys. Fluids 3, 134. — Original coaxial plasma gun design
- Woltjer, L. (1958). "A Theorem on Force-Free Magnetic Fields." Proc. Natl. Acad. Sci. 44, 489. — Conservation of magnetic helicity in ideal MHD

### Topological EM Theory
- Irvine, W.T.M. & Bouwmeester, D. (2008). "Linked and knotted beams of light." Nature Physics 4, 817-820. — Exact Hopfion solutions of Maxwell's equations
- Skyrme, T.H.R. (1962). "A Unified Field Theory of Mesons and Baryons." Nuclear Physics 31, 556-569. — Original Skyrmion/topological soliton theory
- Rañada, A.F. (1992). "Topological electromagnetism." J. Phys. A 25, 1621. — Topological charges in EM via Hopf fibration
- Bateman, H. (1915). "The Mathematical Analysis of Electrical and Optical Wave-Motion." Cambridge. — Foundation for Hopf-linked EM solutions

### Historical and Contextual
- Bearden, T.E. (1988). "Aids to the Establishment of a Far-Reaching Electromagnetic Technology." — Overview of unconventional EM phenomena (treat as historical guide, not rigorous source)
- LaViolette, P. (2008). "Secrets of Antigravity Propulsion." — T.T. Brown and electrogravitics history

---

## Epistemic Commitments

This project adopts the following epistemic standards:

1. **Quantitative predictions over qualitative arguments.** Every claimed effect must be accompanied by a number: how big, at what field strength, detectable with what sensitivity.
2. **Standard EM as the null hypothesis.** For every simulation, the standard Maxwell prediction is computed first. Deviations from standard EM are the signal; the signal-to-noise ratio determines experimental feasibility.
3. **Negative results are results.** If the simulation shows that K deviations are 10⁻³⁰ for all accessible field strengths, that is an important finding. It means the vacuum-K mechanism cannot explain the target phenomenon at those scales, and the search should move elsewhere.
4. **The tier system is binding.** Tier 3 phenomena are not explored until Tier 1 validation succeeds and Tier 2 predictions are quantified. Jumping to extraordinary claims without framework validation is not permitted.
5. **Reproducibility is mandatory.** All simulation code, parameters, and results are published openly. All experimental designs include sufficient detail for independent replication.
6. **Energy claims require extraordinary scrutiny.** Any simulation result suggesting net energy creation is treated as a numerical artifact until all other explanations are exhaustively eliminated.
7. **Honest reporting of all magnitudes.** Predicted effects that are too small to measure are reported as such, not hidden or rationalized away. The gap between theoretical possibility and experimental accessibility is always stated explicitly.

---

## License

This project is released under the MIT License. The theoretical framework synthesizes publicly available physics. No classified or restricted information is used or required.

---

## Contributing

This is an open research project. Contributions are welcome in:
- Theoretical refinement and mathematical rigor
- Simulation implementation and optimization
- Experimental design and execution
- Literature research and historical documentation
- Visualization and data analysis tools

See [ARCHITECTURE.md](ARCHITECTURE.md) for technical implementation details and [TODO.md](TODO.md) for the implementation task breakdown.