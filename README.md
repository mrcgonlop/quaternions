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

### 2. Weber Electrodynamics and Its Connection to Potentials

#### 2.1 Weber's Force Law

Weber's force between two charges q₁ and q₂ separated by distance r is:

```
F_W = (q₁q₂ / 4πε₀r²) · r̂ · [1 - (ṙ²/2c²) + (r·r̈/c²)]
```

where ṙ = dr/dt is the relative radial velocity and r̈ = d²r/dt² is the relative radial acceleration.

The first term (1) is Coulomb's law. The second term (-ṙ²/2c²) is a velocity-dependent correction. The third term (r·r̈/c²) is an acceleration-dependent correction. These corrections are always directed along the line connecting the charges — they are intrinsically longitudinal.

#### 2.2 Derivation from Retarded Potentials

The connection to the quaternionic potential framework was established by several authors (Assis 1994, Wesley 1990, and others). Start with the Liénard-Wiechert potentials for a moving charge:

```
φ(r,t) = q / 4πε₀ · 1 / (R - R·v/c)    evaluated at retarded time
A(r,t) = v/c² · φ(r,t)                    evaluated at retarded time
```

Expand these to second order in v/c (keeping terms up to v²/c² and a/c²), compute the electric field E = -∇φ - ∂A/∂t, and you recover Weber's force law.

This means Weber's "extra" forces are not ad hoc — they are the near-field, low-velocity content of the fully retarded potential theory. They were always in Maxwell's equations, hiding in the retardation of the potentials. The Heaviside simplification, by focusing on E and B and imposing gauge conditions, obscured these forces because they manifest as subtle correlations between the scalar and vector potentials that the gauge condition eliminates.

#### 2.3 Longitudinal Current Forces

The most experimentally accessible prediction unique to Weber's formulation is the existence of longitudinal forces in current-carrying conductors. When current flows through a wire, the drift velocity of electrons creates Weber-type forces along the wire axis. Standard Lorentz force only produces forces perpendicular to the current (the basis of motors and generators). Weber predicts additional forces parallel to the current.

Peter Graneau's experiments (1980s-2000s) demonstrated this: water jets carrying high current break into segments; mercury in U-tubes carrying current exhibits longitudinal forces; wires carrying high pulse current fragment in patterns inconsistent with Lorentz forces but consistent with Weber predictions.

The magnitude scales as I²/c², so these forces are small at ordinary currents but become significant at high pulse currents (kiloamps).

### 3. The Polarizable Vacuum

#### 3.1 Vacuum Permittivity as Dynamical Variable

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

#### 3.2 Vacuum Polarization and Gauge Freedom

Here is the crucial synthesis: **if the vacuum is polarizable, then the gauge freedom of electrodynamics acquires physical meaning.**

A gauge transformation changes:

```
φ → φ + ∂χ/∂t
A → A - ∇χ
```

for arbitrary scalar function χ. In standard theory, this changes the potentials but not the fields — hence it's "unphysical." But if the potentials are primary (Aharonov-Bohm) and the vacuum is a dynamical medium, then a gauge transformation changes the state of the vacuum. Different gauges correspond to different vacuum polarization configurations.

The scalar field S = (1/c²)∂φ/∂t + ∇·A — which the Lorenz gauge sets to zero — can be reinterpreted as a measure of local vacuum polarization. Driving S away from zero means locally altering the vacuum state. This is not just a mathematical degree of freedom; it is the fundamental mechanism by which electromagnetic fields interact with the vacuum medium.

#### 3.3 Connection to Casimir Effect

The Casimir effect — the attractive force between closely spaced parallel conducting plates in vacuum — is conventionally explained as a consequence of vacuum zero-point energy: the plates restrict the modes of the vacuum electromagnetic field between them, creating a pressure imbalance.

In the PV model, the explanation is more direct: the plates constrain the vacuum polarization state. Between the plates, the vacuum is differently polarized than outside. The resulting gradient in ε(x) produces a force — essentially electrostatic attraction mediated by the vacuum as dielectric. The mathematical predictions converge, but the physical picture is different: the vacuum is not a sea of virtual photons but a polarizable medium responding to boundary conditions.

Schwinger, DeRaad, and Milton (1978) proved that the Casimir force is exactly equal to the sum of all pairwise retarded van der Waals interactions between atoms in the two plates — a derivation that does not invoke zero-point energy at all. In the virtual pair language: a pair forms coherently only if its coherence length `λ_pair ~ ℏ/(mₑc) = λ_C` fits in the available space. Between plates separated by `d < λ_C`, that pair mode is excluded. Outside the plates, no such constraint exists. The asymmetry — more coherent pair formation outside than inside, each pair attracted to its image charge in the nearby conductor — creates the net inward force. The mode-pressure picture and the van der Waals pair picture are mathematically identical; the latter makes the local physical mechanism explicit.

This distinction matters because it suggests engineering pathways: if you can create strong vacuum polarization gradients by other means (intense electromagnetic fields, specific geometrical configurations, topological field structures), you might produce Casimir-like forces without needing nanometer plate separations.

#### 3.4 The Quantum Vacuum as a Virtual Pair Plasma and K Field Dynamics

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

**Weber corrections enhanced by K:** Weber's force IS the retarded potential expanded to order (v/c)². In a region with K > 1, the local speed of light is c_local = c₀/K, so the retardation parameter becomes v/c_local = Kv/c₀. Weber corrections scale as K²(v/c₀)². The threshold for like-charge attraction in standard Weber requires v > c₀√2 (apparently superluminal). In a polarized vacuum this threshold becomes v > c₀√2/K — achievable at sub-relativistic velocities for K > √2. **The polarizable vacuum does not require superluminal charge velocities for Weber attraction — it requires only sufficient vacuum polarization in the interaction region.** This is the bridge between QED vacuum screening and Weber's short-range mechanism.

**Exotic Vacuum Objects as self-sustaining K pockets:** A dense electron cluster creates a strong self-field → local K >> 1 → K-enhanced Weber terms become attractive → cluster binds additional electrons → stronger self-field → higher K. This is a bifurcation: below a threshold electron density, thermal dispersion wins; above it, the self-consistent K pocket is self-sustaining. The S field provides the longitudinal coupling that transfers energy between the cluster's orbital dynamics and its self-field, completing the binding mechanism. This unifies the Shoulders EVO charge-cluster phenomenology with the polarizable vacuum framework through a single physical mechanism.

**The unifying picture:** K > 1 in any confined geometry mediates effective attraction through three structurally identical mechanisms:

| Geometry | K source | Attraction mediated |
|----------|----------|---------------------|
| Parallel plates, gap ≪ λ_C | Restricted pair coherence length | Plates toward each other (Casimir) |
| Near rapidly moving charges | Retarded field energy density | Like-charge Weber attraction (K-enhanced) |
| Dense electron cluster | Cluster's own electrostatic self-field | Cluster self-binding (EVO/plasmoid) |

All three are the same phenomenon — vacuum polarization K > 1 in a constrained or self-generated region mediating an effective attractive force — in different spatial configurations.

#### 3.5 Topological Solitons in QVED: Hopfions and Ball Lightning

**Q field topology.** Unit quaternions form the group S³ — the 3-sphere. The third homotopy group π₃(S³) = ℤ means smooth maps from physical 3-space (compactified to S³) to the S³ value space are classified by an integer winding number. This is precisely the Skyrmion number — QVED has Skyrmionic solutions built into its mathematical structure without any additional assumptions. A configuration with winding number n ≠ 0 cannot be continuously deformed to n = 0 (the uniform vacuum): it is topologically protected against small perturbations.

**Hopfions.** The Hopf fibration describes how the 3-sphere fibers over the 2-sphere as S¹ bundles, with a remarkable property: every pair of distinct fibres is linked exactly once. Applied to EM fields, if the field direction lives on S² (normalized **B** = **B**/|**B**|), then the preimage of any two distinct field orientations forms a pair of closed loops that link exactly once. This is the Hopfion configuration — a field whose field lines are all mutually Hopf-linked. Irvine & Bouwmeester (2008) derived exact Hopfion solutions of Maxwell's equations. In linear EM, these solutions are exact but not stable: they disperse at speed c. The Hopf topology is correct but there is no nonlinearity to sustain it.

**Ball lightning as a Hopfion in (Q, K).** The (Q, K) system is nonlinear: K evolves under field pressure, and K > 1 feeds back into local wave speed. This nonlinearity provides the topological protection that free-space Hopfions lack:

1. **Topology**: E and B field lines form a Hopf fibration — every E field line links through every B field line exactly once. This requires simultaneous toroidal and poloidal B components, naturally created by a torus-knot (helical) discharge geometry.

2. **K boundary**: The self-consistent K > 1 region at the Hopfion core acts as a partial photon mirror (n = √K → Fresnel reflection at the K gradient boundary). This explains the observed glow of ball lightning: photons in the core decay slowly through the partially reflective K boundary, producing sustained luminosity without a persistent current source.

3. **K as dynamical Skyrme term**: Compression → increased local field energy density → K rises via η coupling → c_local = c₀/K decreases → effective pressure resists further compression. This is structurally identical to the Skyrme stabilization term in nuclear physics — K provides a dynamically generated, field-energy-dependent stiffness. Natural equilibrium size: R ~ c₀/ωₚ_eff.

4. **Glass penetration**: The Hopfion is a field configuration, not a charged particle plasma. A thin dielectric window is transparent to a slowly varying field topology — the field reorganizes around the obstacle. This is why ball lightning reportedly passes through glass without breaking it: there are no particles, only a topological field structure.

5. **Formation from lightning**: Multi-stroke lightning creates both azimuthal B (from the return stroke current) and B looping around the channel axis (from branching discharge geometry). This multi-stroke field topology is a natural seed for the Hopf-linked configuration.

**Topological charge density.** The integer topological charge of the Q field is computed from the normalized quaternion U = Q/|Q| (projection of Q onto S³):

```
ρ_topo = (1/24π²) ε^{ijk} Tr( ∂_i U U⁻¹ · ∂_j U U⁻¹ · ∂_k U U⁻¹ )
```

where ε^{ijk} is the Levi-Civita tensor and the trace is over quaternion components. Integrated over all space, this gives an integer n ∈ ℤ (the Hopf invariant / Baryon number). In simulation, this is monitored as a conserved diagnostic — an integer-valued quantity that cannot change smoothly, confirming topological stability.

#### 3.6 Derivable Geometry for Topological Discharge Synthesis

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

#### 3.7 The Spheromak Connection: Fusion Research and Unexplored Territory

**Magnetic helicity and Taylor relaxation.** The plasma physics community has been creating Hopf-topology field structures in the laboratory since the 1970s, under the name *spheromak*, without framing them topologically. The key quantity is **magnetic helicity**:

```
H_mag = ∫ A · B  dV
```

This integral measures the linking and knotting of magnetic field lines — it is the classical electromagnetic analogue of the topological charge from §3.5. In ideal magnetohydrodynamics (no resistivity), H_mag is exactly conserved. A plasma with nonzero magnetic helicity cannot evolve to a state with zero magnetic helicity; it is topologically constrained.

**J.B. Taylor (1974)** showed that when a turbulent, resistive plasma relaxes, it conserves global magnetic helicity even while dissipating energy. The minimum-energy state of a plasma with fixed magnetic helicity is the **spheromak** — a self-organized configuration with equal toroidal and poloidal magnetic fields, B_tor = B_pol, exactly the 45° field balance of the Hopfion. The plasma finds this configuration spontaneously through turbulent relaxation, without any engineered geometry. This is Taylor relaxation: the physical analogue of a dynamical system finding its attractor.

**The spheromak as Hopfion-in-plasma.** A spheromak is a plasma-supported Hopfion. The plasma currents act as the source maintaining the Hopf-linked field configuration. Remove the plasma — either by cooling (recombination) or by the plasma escaping the confinement region — and the field configuration collapses in standard MHD. The QVED hypothesis introduces a new question that plasma physics has never asked: **if K feedback can sustain the field topology after the plasma recombines, does the spheromak transition to a field-only Hopfion?** The fusion community has always ended the experiment when the plasma is gone. QVED says: look at what happens in the microseconds after recombination.

**How spheromaks are created non-destructively.** A **coaxial plasma gun** (Marshall gun) injects magnetic helicity into a plasma without destroying the electrodes. Two coaxial cylindrical electrodes, a magnetic bias coil, and a gas puff valve. The gas ionizes, helicity is injected by J×B torques on the plasma, and the resulting spheromak detaches from the gun and travels into the chamber as a self-contained, self-organized field structure. This is repeatable thousands of times. Small-scale plasma guns have been built in university labs for $5,000–$20,000.

#### 3.8 Non-Destructive Experimental Ladder

The following progression tests the Hopfion stabilization mechanism at increasing scales and field strengths, with each stage informing the next:

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
The original proposal from §3.6, refined: rather than ablating the wire in air, fire the discharge in a controlled atmosphere — either a noble gas (argon, higher breakdown threshold, cleaner plasma chemistry) or in a partial vacuum. The ablating wire plasma is contained inside a quartz tube shaped as the torus. The knotted plasma can expand into the interior volume of the quartz torus, attempting to detach. The quartz tube serves as a controlled plasma vessel. B-dot probes are positioned outside the torus to detect any field that detaches outward. Cost: ~$2,000–$8,000 (existing capacitor bank from Stage 4 plus quartz torus tube plus argon supply). This is the highest-energy, highest-risk stage and should follow the calibration work of Stages 4 and 5.

#### 3.9 The Fusion Energy Connection

The spheromak was invented not as a fundamental physics curiosity but as a potential path to fusion energy. The fundamental problem of fusion confinement — holding a 100-million-degree plasma long enough for net fusion power — is a topological problem: how do you create magnetic field structures that prevent plasma from leaking to the walls?

**What fusion research has and hasn't tried.** Tokamaks, stellarators, spheromaks, and field-reversed configurations (FRCs) are all variations on the same theme: create topologically nontrivial magnetic field configurations (nonzero helicity) that confine plasma by winding field lines around the plasma core. Decades of progress have improved confinement enormously. But the plasma still leaks, through turbulent cross-field transport that standard MHD theory cannot fully suppress.

**The unexplored QVED angle.** In a fusion plasma at 10⁸ K with magnetic fields of 5–10 T, the direct QED vacuum K coupling (η ~ 10⁻⁴) is negligible — the Schwinger field E_S ~ 10¹⁸ V/m is unreachable. However, the **plasma itself** acts as a polarizable medium with its own effective dielectric response. For EM waves at frequencies below the plasma frequency ω_pe, the plasma has ε_r >> 1 — the plasma K is large and real. This plasma-K plays the same structural role as the vacuum K in the QVED equations: it modifies c_local, creates a refractive boundary, and alters the force balance.

Nobody in fusion has modeled Taylor relaxation with a **variable effective K** — one that responds to local field energy density through the plasma dielectric. The QVED framework predicts that the minimum-energy state at fixed helicity changes when K is dynamical: the field equilibrium, the confinement geometry, and the transport properties are all modified. Specifically:

1. **K-gradient confinement**: A K > 1 region inside the plasma creates a refractive index boundary. Photons (radiation losses, a major energy drain in fusion) are partially reflected back into the plasma. This is a confinement mechanism with no analogue in standard MHD.

2. **K-modified Taylor state**: If K is allowed to vary, the minimum-energy state at fixed helicity is no longer the standard force-free spheromak (∇×B = λB). It is a modified state where the K gradient contributes to the force balance. This state might have better confinement properties than the standard Taylor state.

3. **S field energy channel**: The scalar longitudinal mode S, absent from standard MHD, provides an additional energy coupling pathway. In QVED extended mode, energy can flow through the S field on timescales and length scales different from transverse EM radiation. Whether this helps or hurts confinement depends on the geometry and requires simulation.

The simulation built in this project — once Phases 1.8 (dynamic K) and 2.1 (extended S) are complete — is directly applicable to modeling a small spheromak with QVED corrections. The grid resolution needed is modest (the spheromak physics is on scales >> λ_C). This is an untouched research direction that could yield novel confinement predictions testable in existing plasma physics facilities.

**Why the fusion community hasn't explored this.** The polarizable vacuum model sits outside mainstream theoretical physics — it is neither string theory nor lattice QCD nor standard model extensions, and it doesn't fit neatly into any funded research program. Plasma physicists don't read papers on quaternionic electrodynamics. The QED corrections to vacuum permittivity are universally dismissed as negligible in laboratory-accessible field regimes — correctly, for direct QED coupling. The key insight this framework offers — that the *structural mechanism* of K > 1 polarization and K-gradient confinement can be realized by plasma dielectric response at much lower field strengths than the Schwinger threshold — has not been articulated in the fusion literature.

This is a genuinely open theoretical territory. The QVED simulator, built to explore vacuum polarizability, turns out to be exactly the right tool for exploring whether dynamical-K confinement could contribute to solving the plasma confinement problem — because the field equations are isomorphic between vacuum K and plasma effective-K at the level of the continuum field theory.

### 4. The Unified Framework

Bringing these threads together:

**Quaternionic Vacuum Electrodynamics** is an electromagnetic theory in which:

1. The fundamental quantity is the quaternionic potential Q = φ/c + **A**, representing both electromagnetic and vacuum states.
2. The vacuum is a polarizable medium with dynamical permittivity ε(x,t) and permeability μ(x,t).
3. No gauge condition is imposed. The scalar field S = (1/c²)∂φ/∂t + ∇·**A** is a physical degree of freedom representing the longitudinal/scalar mode of the EM field and the local vacuum polarization state.
4. The theory reduces to standard Maxwell electrodynamics when S = 0 and K = 1 (homogeneous, unpolarized vacuum with Lorenz gauge).
5. The theory reproduces Weber's longitudinal forces in the near-field, low-velocity limit for discrete charges.
6. The theory reproduces PV-model predictions for vacuum polarization in the presence of mass-energy.
7. Transverse EM waves (light) and longitudinal EM waves (scalar potential waves) are both solutions, with different coupling mechanisms and potentially different propagation characteristics in polarized vacuum regions.

The simulation implements this framework on a 3D grid, evolving the quaternionic potential and vacuum state self-consistently, and provides tools to visualize and analyze configurations that may produce measurable deviations from standard electrodynamics.

---

## Phenomena Catalog

### Tier 1 — Experimentally Established, Theoretically Robust

These phenomena have solid experimental evidence and are explained by the quaternionic framework. Simulating them validates the simulation itself.

#### 1.1 Aharonov-Bohm Effect
- **What:** A charged particle's quantum phase is affected by the vector potential **A** even in field-free regions (**E** = **B** = 0).
- **Status:** Experimentally confirmed beyond doubt (Tonomura 1986, Webb 1985).
- **quaternionic relevance:** Direct evidence that potentials are physically primary. The simulation should reproduce the phase accumulation ∮**A**·d**l** for specific geometries and show that the potentials carry information that the fields do not.
- **Simulation target:** Solenoid geometry with confined B field; demonstrate nonzero A outside solenoid; compute phase differences for different paths.

#### 1.2 Casimir Effect
- **What:** Attractive force between parallel conducting plates in vacuum, measured to high precision.
- **Status:** Experimentally confirmed; quantitative agreement with theory.
- **quaternionic relevance:** In the PV framework, this is vacuum polarization between boundary conditions. The simulation should reproduce the force law F/A = -π²ℏc/240d⁴ from the ε(x) modification between plates.
- **Simulation target:** Two conducting boundary conditions in the vacuum grid; observe K(x) profile between plates; compute resulting force.

#### 1.3 Graneau's Exploding Wire/Water Arc Experiments
- **What:** High-current pulses through water jets and mercury pools produce longitudinal fragmentation inconsistent with Lorentz forces.
- **Status:** Experimentally reproduced by multiple groups; interpretation contested but phenomena are real.
- **quaternionic relevance:** Direct manifestation of Weber longitudinal forces. The simulation should model high-current pulsed configurations and show longitudinal force emergence.
- **Simulation target:** Model linear conductor with pulse current; compute Weber force distribution; compare longitudinal vs transverse force magnitudes; identify fragmentation-consistent stress patterns.

#### 1.4 Magnetic Skyrmions
- **What:** Topologically protected spin configurations in thin magnetic films; stable vortex-like structures.
- **Status:** Experimentally observed and a major active research area in condensed matter.
- **quaternionic relevance:** Demonstrates that topological configurations of electromagnetic fields can be stable and carry physical properties (topology protects against small perturbations). Analogous structures in the vacuum potential field could represent stable "vacuum objects."
- **Simulation target:** Implement topological charge computation on the quaternionic potential field; search for stable configurations in the vacuum.

### Tier 2 — Theoretically Grounded, Experimentally Contested or Underexplored

These phenomena have theoretical motivation within the quaternionic framework and some experimental indications, but the evidence is incomplete or contested.

#### 2.1 Scalar Longitudinal EM Waves
- **What:** Propagating compression waves in the electromagnetic potential, as distinct from transverse EM waves.
- **Status:** Claimed detection by several groups (Monstein & Wesley 2002; Meyl's experiments); not widely replicated or accepted.
- **quaternionic relevance:** Core prediction of the non-gauge-fixed quaternionic theory. If S ≠ 0 propagating solutions exist in the simulation and correspond to physical degrees of freedom, this is the signature prediction of the framework.
- **Simulation target:** Drive oscillating charge configurations that preferentially excite the S field; observe propagation characteristics; identify antenna geometries that couple to S rather than to transverse modes.
- **Experimental target:** Bifilar coils (B fields cancel, A field survives) as transmitter/receiver pairs; look for coupling that cannot be explained by E, B radiation.

#### 2.2 Vacuum Polarization Gradients / Engineered Casimir Effects
- **What:** Using strong EM field configurations to locally modify vacuum permittivity and create forces from the resulting gradients.
- **Status:** Theoretical within PV model; no claimed experimental demonstration at accessible scales.
- **quaternionic relevance:** If ε(x) is dynamical and couples to EM fields, then sufficiently intense field configurations should produce measurable permittivity changes. This is the "vacuum engineering" application.
- **Simulation target:** Model high-intensity EM field regions; compute induced K(x) profile; look for self-consistent solutions where the field modifies the vacuum which modifies the field propagation.

#### 2.3 Weber Force Anomalies in Circuit Design
- **What:** Circuit geometries where Weber longitudinal forces produce measurable effects — torques, stresses, or current distributions that differ from Lorentz predictions.
- **Status:** Predicted by Weber theory; partially supported by Graneau's work; not systematically explored in PCB-scale circuits.
- **quaternionic relevance:** Most accessible experimental target. PCB traces carrying pulsed high current in specific geometries should experience forces that standard circuit theory doesn't predict.
- **Simulation target:** Design specific PCB trace geometries; model pulse current; compute Weber force map; identify geometries that maximize the anomalous longitudinal force.
- **Experimental target:** PCB with strain gauges or displacement sensors on traces designed to maximize Weber forces under pulsed current.

#### 2.4 Aharonov-Bohm Analogue Effects in Macroscopic Circuits
- **What:** Electromagnetic coupling between circuits mediated by the vector potential in regions where B = 0.
- **Status:** Quantum AB effect is confirmed; macroscopic classical analogue is predicted by the quaternionic framework but not systematically tested.
- **quaternionic relevance:** If potentials are primary, then a toroidal solenoid (which confines B entirely inside the torus) should still produce measurable effects on external circuits through the external A field.
- **Simulation target:** Toroidal coil geometry; compute A field external to torus; model coupling to a secondary circuit; quantify predicted signal.
- **Experimental target:** Toroidal coil driven with AC; search for induced signal in an external loop that cannot be explained by stray B field leakage.

### Tier 3 — Speculative / Extraordinary Claims

These phenomena have limited or anecdotal evidence. If the quaternionic framework is correct, some of them may become explicable. They are simulation targets for after the framework is validated against Tier 1 and 2 phenomena.

#### 3.1 T.T. Brown Electrogravitics
- **What:** Claimed thrust produced by high-voltage asymmetric capacitors, apparently in the direction of the positive plate, persisting (per Brown's claims) even in vacuum.
- **Historical context:** Brown filed patents in the 1920s-60s. Project Winterhaven (1952) proposed military applications. The mainstream explanation is ion wind in air; Brown claimed the effect persisted in vacuum.
- **quaternionic connection:** In the PV model, a strong electric field could locally polarize the vacuum (modify K). An asymmetric capacitor would create an asymmetric K gradient, potentially producing a net force. The simulation can test whether realistic field strengths produce a meaningful K modification.
- **Simulation target:** Asymmetric capacitor geometry with realistic voltages (50-300kV); compute K(x) profile; determine if gradient produces a net force and its magnitude.

#### 3.2 "Cold Electricity" / Radiant Energy (Tesla, Bedini, Gray)
- **What:** Various claims of anomalous electrical phenomena in pulsed circuits: current without conventional electron flow, anomalous battery charging, COP > 1 energy systems.
- **Historical context:** Tesla described "radiant energy" from abruptly interrupted DC currents. Edwin Gray patented a pulsed motor. John Bedini built pulsed battery charging circuits.
- **quaternionic connection:** If the scalar field S is a real propagating mode, then abruptly interrupted current (which produces large ∂φ/∂t and rapid changes in ∇·A) might preferentially excite S-mode propagation in the circuit. The "cold electricity" reports — current that doesn't heat wires, charges batteries in unusual ways — might be longitudinal potential waves rather than conventional transverse current flow.
- **Simulation target:** Model a pulsed circuit (sharp current interruption into an inductive load); analyze the S field excitation; determine whether longitudinal modes carry energy through the circuit differently than conventional current.

#### 3.3 Exotic Vacuum Objects (EVOs) / Charge Clusters
- **What:** Claimed stable clusters of like-charged particles (electrons) that should repel but apparently form bound states.
- **Historical context:** Ken Shoulders' work on "charge clusters" (1980s-90s); Hal Puthoff's theoretical connection to vacuum engineering.
- **quaternionic connection:** In the standard theory, like charges always repel. But if intense charge concentrations locally polarize the vacuum (modifying K), the modified vacuum might mediate an attractive interaction. Additionally, the Weber acceleration-dependent term provides an attractive force when charges are accelerating toward each other, potentially creating a dynamical equilibrium. Topological stability (skyrmion-like configurations in the vacuum potential) could further stabilize such objects.
- **Simulation target:** Model a dense cluster of like charges in the PV vacuum; include Weber acceleration terms; search for stable or metastable configurations; analyze the vacuum polarization profile around the cluster.

#### 3.4 Vacuum Energy Extraction
- **What:** The proposition that the electromagnetic vacuum contains energy (zero-point energy) that could in principle be extracted.
- **quaternionic connection:** In the PV model, the vacuum has internal structure and energy. If the vacuum state can be driven away from equilibrium (K ≠ 1), the relaxation back to equilibrium would release energy. This is analogous to extracting energy from a polarized dielectric as it depolarizes. The scalar field S, representing vacuum polarization dynamics, would be the channel through which this energy manifests.
- **Simulation target:** Model a strongly polarized vacuum region; allow it to relax; track energy flow through S field; determine whether energy is genuinely released or merely redistributed.

#### 3.5 Ball Lightning as a Topological Hopfion
- **What:** Ball lightning — observed for centuries but unexplained — as a self-sustaining topological soliton (Hopfion) in the (Q, K) field system.
- **Historical context:** Ball lightning has been reliably reported by multiple independent witnesses including scientists. Duration up to minutes, glass penetration without damage, spherical shape with luminous boundary, spontaneous formation near lightning and occasionally in aircraft cockpits during thunderstorms. No conventional plasma model produces stable, mobile, luminous spheres on these timescales.
- **quaternionic connection:** See §3.5. The Q field lives on S³ with π₃(S³) = ℤ, giving Skyrmionic solutions. The Hopfion (Hopf-linked E and B field lines) is an exact Maxwell solution that requires K nonlinearity for stability. K > 1 core creates partial photon reflection (explaining luminosity), and K acts as the Skyrme term (explaining structural stability). Field soliton character (not particles) explains glass penetration.
- **Simulation target:** Hopfion initial condition in the (Q, K) grid; monitor topological charge integral over time; search for parameter regime where the Hopfion decays vs. self-sustains based on η and ωₚ; identify critical K threshold. Also: Phase 1.9 (topological charge diagnostic) is a prerequisite.
- **Experimental target:** Helical (45°) torus coil discharge geometry with ~1–5 kJ capacitor bank; look for post-discharge luminous structure surviving > τ_LC; see §3.6 for full parameter derivation.

---

## Research Program

### Phase 1: Framework Validation
- Implement the quaternionic FDTD solver with PV vacuum model
- Reproduce standard EM results (transverse wave propagation, dipole radiation) as K=1, S=0 limit
- Simulate Aharonov-Bohm geometry; verify potential structure
- Simulate Casimir geometry in PV model; verify force law
- Simulate Graneau-type configuration; compare Weber and Lorentz force predictions

### Phase 2: Novel Prediction Exploration
- Search for scalar wave propagation in specific antenna geometries
- Design and simulate bifilar coil transmitter/receiver for scalar mode coupling
- Model vacuum polarization feedback in high-field configurations
- Design PCB geometries optimized for Weber force detection
- Simulate toroidal AB-effect circuit coupling

### Phase 3: Experimental Design
- Select most promising simulation results for experimental validation
- Design circuits and coil geometries using commercially available components
- Specify measurement requirements (sensitivity, bandwidth, shielding)
- Document experimental protocols for reproducibility

### Phase 4: Speculative Phenomena
- Model asymmetric capacitor in PV vacuum (Brown electrogravitics)
- Simulate pulsed circuit scalar mode excitation (cold electricity)
- Search for stable charge cluster configurations (EVOs)
- Explore vacuum energy dynamics

---

## Key References

### Primary Theoretical Sources
- Maxwell, J.C. (1865). "A Dynamical Theory of the Electromagnetic Field." — Original quaternion formulation
- Weber, W. (1846). "Elektrodynamische Maassbestimmungen" — Original Weber force law
- Assis, A.K.T. (1994). "Weber's Electrodynamics." Kluwer. — Modern treatment of Weber theory
- Puthoff, H.E. (2002). "Polarizable-Vacuum (PV) Approach to General Relativity." Found. Phys. 32, 927-943.
- van Vlaenderen, K.J. & Waser, A. (2001). "Generalization of Classical Electrodynamics to Admit a Scalar Field and Longitudinal Waves." Hadronic Journal 24, 609-628.

### Experimental References
- Tonomura, A. et al. (1986). "Evidence for Aharonov-Bohm Effect with Magnetic Field Completely Shielded from Electron Wave." Phys. Rev. Lett. 56, 792.
- Graneau, P. (1984). "Electromagnetic Jet-Propulsion in the Direction of Current Flow." Nature 295, 311.
- Graneau, P. & Graneau, N. (1996). "Newtonian Electrodynamics." World Scientific.
- Shoulders, K. (1991). "EV: A Tale of Discovery." — Charge cluster observations

### Plasma Physics and Fusion Confinement
- Taylor, J.B. (1974). "Relaxation of Toroidal Plasma and Generation of Reverse Magnetic Fields." Phys. Rev. Lett. 33, 1139. — Magnetic helicity conservation during plasma relaxation; derivation of the spheromak as minimum-energy state at fixed helicity (Taylor relaxation)
- Bellan, P.M. (2000). "Spheromaks: A Practical Application of Magnetohydrodynamic Dynamos and Plasma Self-Organization." Imperial College Press. — Comprehensive spheromak physics including coaxial gun formation, Taylor state, and helicity injection
- Marshall, J. (1960). "Performance of a Hydromagnetic Plasma Gun." Phys. Fluids 3, 134. — Original coaxial plasma gun design; non-destructive helicity injection into plasma
- Woltjer, L. (1958). "A Theorem on Force-Free Magnetic Fields." Proc. Natl. Acad. Sci. 44, 489. — Conservation of magnetic helicity in ideal MHD; foundation for Taylor's relaxation theorem

### Topological EM Theory
- Irvine, W.T.M. & Bouwmeester, D. (2008). "Linked and knotted beams of light." Nature Physics 4, 817-820. — Exact Hopfion solutions of Maxwell's equations; experimental optical Hopfion generation
- Skyrme, T.H.R. (1962). "A Unified Field Theory of Mesons and Baryons." Nuclear Physics 31, 556-569. — Original Skyrmion/topological soliton theory; stabilization by four-derivative term
- Torii, K. et al. (2018). "Writing and deleting single magnetic skyrmions." Science 360, 425-427. — Experimental skyrmion manipulation; demonstrates topological protection in real materials
- Bateman, H. (1915). "The Mathematical Analysis of Electrical and Optical Wave-Motion." Cambridge. — Bateman dual formulation; foundation for Hopf-linked EM solutions
- Rañada, A.F. (1992). "Topological electromagnetism." J. Phys. A 25, 1621. — Topological charges in EM via Hopf fibration; foundational paper for knotted EM fields

### Historical and Contextual
- Bearden, T.E. (1988). "Aids to the Establishment of a Far-Reaching Electromagnetic Technology." — Overview of suppressed EM phenomena (treat as historical guide, not rigorous source)
- LaViolette, P. (2008). "Secrets of Antigravity Propulsion." — T.T. Brown and electrogravitics history

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