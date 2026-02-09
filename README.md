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

This distinction matters because it suggests engineering pathways: if you can create strong vacuum polarization gradients by other means (intense electromagnetic fields, specific geometrical configurations, topological field structures), you might produce Casimir-like forces without needing nanometer plate separations.

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