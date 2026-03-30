# Study Guide — Build Plan

## Format
- Single Markdown file with LaTeX math (`$...$` inline, `$$...$$` display)
- Renderable by: VS Code Markdown+Math, Obsidian, Typora, Pandoc → PDF
- Each chapter: historical context → key experiment → the equation derived from it → what it means physically
- No pre-interpreted "fields" or "particles" until the experiments that demanded those concepts

## Chapter Structure

### Part I: Static Electricity (600 BC – 1800)
- `[x]` **Ch 1 — The Amber Effect** — Thales of Miletus (~600 BC)
  - Amber rubbed with wool attracts light objects
  - The word "elektron" (ἤλεκτρον) = amber
  - No equations — pure phenomenology
  - *Key idea:* matter can acquire an invisible attractive property through friction

- `[x]` **Ch 2 — The Magnetic Earth** — William Gilbert (1600)
  - *De Magnete*: systematic distinction between electric vs magnetic attraction
  - The terrella model — Earth itself is a magnet
  - Invented the versorium (electroscope)
  - *Key idea:* magnetism and electricity are related but distinct phenomena

- `[x]` **Ch 3 — Generating and Storing Charge** — Otto von Guericke, Stephen Gray, du Fay, Leyden jar
  - Von Guericke's sulfur globe (1663) — first electrostatic generator
  - Gray (1729) — conductors vs insulators, charge flows through materials
  - Du Fay (1733) — two kinds of electricity: vitreous (+) and resinous (−)
  - Musschenbroek & von Kleist (1745) — the Leyden jar: first capacitor
  - *Key idea:* charge can be generated, conducted, classified, and stored

- `[x]` **Ch 4 — One Fluid, Lightning, and Grounding** — Benjamin Franklin (1752)
  - Single-fluid model: excess = positive, deficit = negative
  - Kite experiment: lightning is electrical discharge
  - Conservation of charge: charge is neither created nor destroyed, only transferred
  - The lightning rod — first applied electrostatics
  - *Key idea:* electricity is a single conserved substance; atmospheric lightning is the same stuff as Leyden jar sparks

- `[x]` **Ch 5 — The Inverse Square Law** — Coulomb (1785), Cavendish (1770s)
  - Cavendish: measured it first (unpublished), null experiment with concentric spheres
  - Coulomb's torsion balance: direct measurement of force vs distance and charge
  - $$F = \frac{1}{4\pi\varepsilon_0} \frac{q_1 q_2}{r^2} \hat{r}$$
  - Parallel to Newton's gravitational law — deliberate analogy
  - *Key idea:* the first quantitative law of electricity; force as a function of charge and distance

- `[x]` **Ch 6 — Mathematical Formalization of Potentials** — Poisson, Laplace, Green
  - Laplace's equation: $\nabla^2 \phi = 0$ (charge-free regions)
  - Poisson's equation: $\nabla^2 \phi = -\rho/\varepsilon_0$ (with charge)
  - Green's theorem: boundary value problems, the concept of a potential function
  - *Key idea:* the electric potential φ as a scalar field from which forces derive; the birth of field mathematics

### Part II: Current Electricity and Magnetism (1800 – 1831)
- `[x]` **Ch 7 — From Frogs to Batteries** — Galvani (1780s), Volta (1800)
  - Galvani's frog legs — "animal electricity"
  - Volta's rebuttal: the voltaic pile (first battery) — sustained current from dissimilar metals
  - The concept of electromotive force (EMF) and sustained current
  - *Key idea:* electricity can flow continuously, not just discharge; the birth of current electricity

- `[x]` **Ch 8 — The Magnetic Needle Deflects** — Ørsted (1820)
  - The lecture demonstration: current-carrying wire deflects a compass needle
  - Current produces magnetism — electricity and magnetism are linked
  - The effect is perpendicular to the wire (circular magnetic field)
  - *Key idea:* the first experimental proof that electricity causes magnetism; perpendicularity is fundamental

- `[x]` **Ch 9 — Quantifying the Magnetic Field of a Current** — Biot, Savart, Ampère (1820–1827)
  - Biot-Savart law: $$d\mathbf{B} = \frac{\mu_0}{4\pi} \frac{I \, d\boldsymbol{\ell} \times \hat{r}}{r^2}$$
  - Ampère's force law between current-carrying wires:
    $$\frac{F}{L} = \frac{\mu_0}{2\pi} \frac{I_1 I_2}{d}$$
  - Ampère's circuital law: $\oint \mathbf{B} \cdot d\boldsymbol{\ell} = \mu_0 I_{\text{enc}}$
  - Ampère's "swimming man" rule and right-hand rule
  - *Key idea:* magnetic effects of currents are fully quantifiable; forces between circuits obey precise mathematical laws

- `[x]` **Ch 10 — Resistance and Circuit Laws** — Ohm (1827), Kirchhoff (1845)
  - Ohm's law: $V = IR$, resistance as a property of materials
  - Kirchhoff's current law: $\sum I_{\text{in}} = \sum I_{\text{out}}$ (charge conservation at a node)
  - Kirchhoff's voltage law: $\sum V = 0$ around any closed loop (energy conservation)
  - Kirchhoff's discovery: electrical signals propagate at light speed in wires
  - *Key idea:* circuits as systems governed by conservation laws; the surprising speed of electrical signals

- `[x]` **Ch 11 — Thermoelectric Effects** — Seebeck (1821), Peltier (1834), Lord Kelvin
  - Seebeck effect: temperature difference → voltage in a bimetallic junction
  - Peltier effect: current → heating/cooling at a junction (reverse of Seebeck)
  - Thomson (Lord Kelvin) effect: current in a temperature gradient → heat absorption/emission
  - Thomson's thermodynamic analysis unifying all three
  - *Key idea:* electricity and heat are intimately coupled; energy conversion between thermal and electrical domains

- `[x]` **Ch 12 — Electromagnetic Induction** — Faraday (1831), Joseph Henry (1831), Lenz (1834)
  - Faraday's experiments: moving magnet → current in a coil
  - Joseph Henry: self-inductance, mutual inductance (discovered independently)
  - Faraday's law: $$\mathcal{E} = -\frac{d\Phi_B}{dt}$$
  - Lenz's law: induced current opposes the change that caused it
  - Faraday's lines of force — the birth of the field concept (no equations, pure geometric intuition)
  - *Key idea:* changing magnetism creates electricity (completing the circle); Faraday's field lines as the conceptual seed of field theory

### Part III: The Great Unifications (1840 – 1890)
- `[x]` **Ch 13 — Energy Conservation** — Joule (1840s), Lord Kelvin
  - Joule heating: $P = I^2 R$ — electrical energy → heat
  - Mechanical equivalent of heat
  - Lord Kelvin and the laws of thermodynamics
  - Energy as the unifying concept across electricity, magnetism, heat, and mechanics
  - *Key idea:* energy is conserved across all domains; this constrains what electromagnetic theories are possible

- `[x]` **Ch 14 — Weber's Electrodynamics** — Wilhelm Weber (1846)
  - Weber's force law: $$F_W = \frac{q_1 q_2}{4\pi\varepsilon_0 r^2} \hat{r} \left[1 - \frac{\dot{r}^2}{2c^2} + \frac{r\ddot{r}}{c^2}\right]$$
  - Velocity-dependent and acceleration-dependent corrections to Coulomb
  - Longitudinal forces along the line connecting charges
  - Weber-Kohlrausch experiment (1856): measured the ratio of electrostatic to electromagnetic units → obtained c (speed of light!) before Maxwell
  - *Key idea:* forces between charges depend on their relative motion; the speed of light appears naturally in electrodynamics; longitudinal forces exist in this framework

- `[x]` **Ch 15 — Maxwell's Quaternionic Equations** — James Clerk Maxwell (1865)
  - Hamilton's quaternions: $q = w + xi + yj + zk$, where $i^2 = j^2 = k^2 = ijk = -1$
  - Maxwell's original 20 equations in quaternion notation
  - The displacement current: $\mu_0 \varepsilon_0 \partial\mathbf{E}/\partial t$
  - The quaternion potential: $Q = \phi/c + A_x i + A_y j + A_z k$
  - The quaternionic nabla: $\nabla_q = \frac{1}{c}\frac{\partial}{\partial t} + \frac{\partial}{\partial x}i + \frac{\partial}{\partial y}j + \frac{\partial}{\partial z}k$
  - The product $\nabla_q Q$ produces: scalar part S, symmetric vector part E, antisymmetric vector part B
  - The scalar field: $S = \frac{1}{c^2}\frac{\partial\phi}{\partial t} + \nabla\cdot\mathbf{A}$
  - Prediction: electromagnetic waves propagate at $c = 1/\sqrt{\mu_0\varepsilon_0}$
  - *Key idea:* light is an electromagnetic wave; the quaternion formulation naturally contains a scalar field S that couples φ and A; this S is a dynamical degree of freedom

- `[x]` **Ch 16 — The Heaviside Reduction** — Oliver Heaviside (1884–1887)
  - Heaviside's 4 vector equations from Maxwell's 20:
    $$\nabla\cdot\mathbf{E} = \frac{\rho}{\varepsilon_0}$$
    $$\nabla\cdot\mathbf{B} = 0$$
    $$\nabla\times\mathbf{E} = -\frac{\partial\mathbf{B}}{\partial t}$$
    $$\nabla\times\mathbf{B} = \mu_0\mathbf{J} + \mu_0\varepsilon_0\frac{\partial\mathbf{E}}{\partial t}$$
  - What was gained: notational clarity, computational efficiency
  - What was lost: the scalar part S (set to zero by Lorenz gauge), quaternion structure, potential primacy
  - The Lorenz gauge condition: $\frac{1}{c^2}\frac{\partial\phi}{\partial t} + \nabla\cdot\mathbf{A} = 0$ (S ≡ 0 by definition)
  - Poynting vector (1884): $\mathbf{S}_P = \frac{1}{\mu_0}\mathbf{E}\times\mathbf{B}$ — energy flow in EM fields
  - *Key idea:* Heaviside's reformulation was presented as pure simplification but it eliminated a dynamical degree of freedom; the potentials were demoted from physical objects to mathematical conveniences

- `[x]` **Ch 17 — Experimental Confirmation of EM Waves** — Heinrich Hertz (1887)
  - Hertz's oscillator and resonator: generating and detecting radio waves
  - Measured wavelength, frequency, speed → confirmed $c = \lambda f$
  - Demonstrated reflection, refraction, polarization of EM waves
  - *Key idea:* Maxwell's prediction experimentally verified; but only transverse waves were tested — longitudinal modes (if they exist) would require different detection

- `[x]` **Ch 18 — Alternating Current and Wireless** — Tesla (1880s–1900s), Bell (1876)
  - Bell's telephone (1876): converting sound → electrical signal → sound (applied EM)
  - Tesla's rotating magnetic field (1882): the key insight behind AC motors
  - Polyphase AC system: generation, transmission, transformation
  - Tesla coil: resonant transformer, high-voltage high-frequency
  - Tesla's wireless transmission experiments — longitudinal vs transverse debate
  - *Key idea:* practical mastery of electromagnetic phenomena; Tesla's insistence on longitudinal wave effects connects to the quaternionic scalar mode

### Part IV: Electrons, Relativity, and the Lorentz Force (1890 – 1930)
- `[x]` **Ch 19 — The Electron** — J.J. Thomson (1897)
  - Cathode ray experiments: deflection by E and B fields
  - Measurement of charge-to-mass ratio: $e/m$
  - The electron as the first subatomic particle
  - *Key idea:* electric current is the flow of discrete charged particles; the microscopic basis of Ohm's law

- `[x]` **Ch 20 — The Ether Crisis** — Michelson & Morley (1887), FitzGerald, Lodge, Lorentz
  - The luminiferous ether: the supposed medium for EM waves
  - Michelson-Morley interferometer: null result — no ether wind detected
  - FitzGerald contraction: length contracts in direction of motion
  - Lodge's experiments: searching for ether drag
  - Lorentz transformations: $x' = \gamma(x - vt)$, $t' = \gamma(t - vx/c^2)$
  - The Lorentz force: $\mathbf{F} = q(\mathbf{E} + \mathbf{v}\times\mathbf{B})$
  - *Key idea:* no preferred frame; the speed of light is invariant; the Lorentz force replaces Weber's force in the standard framework — but note it has NO longitudinal component

- `[x]` **Ch 21 — Special Relativity and the Photon** — Einstein (1905)
  - Photoelectric effect: $E = hf - \phi$ — light as quanta
  - Special relativity: the speed of light is the same in all frames
  - $E = mc^2$: mass-energy equivalence
  - Electromagnetic field tensor: E and B unified into $F^{\mu\nu}$
  - Relativity as the successor to ether theory — but what about the vacuum structure?
  - *Key idea:* light is both wave and particle; spacetime replaces ether; but the vacuum may still have structure (Casimir, QED)

### Part V: Quantum Electrodynamics and the Vacuum (1920 – 1960)
- `[x]` **Ch 22 — Quantum Mechanics Meets EM** — Planck (1900), Dirac (1928)
  - Planck's blackbody radiation: $E = nhf$ — energy is quantized
  - Dirac equation: relativistic quantum mechanics, predicts antimatter
  - Magnetic monopoles: Dirac's argument from charge quantization
  - *Key idea:* the electromagnetic field is quantized; the vacuum has zero-point energy

- `[x]` **Ch 23 — QED: The Precise Theory** — Feynman, Schwinger, Tomonaga (1940s–1950s)
  - Path integral formulation (Feynman): sum over all possible photon exchanges
  - Renormalization: taming infinities in vacuum self-energy
  - Anomalous magnetic moment of the electron: theory matches experiment to 12 decimal places
  - Lamb shift: vacuum fluctuations shift atomic energy levels
  - *Key idea:* QED is the most precisely tested theory in physics; it is built on gauge invariance (the assumption that S = 0 is physical)

- `[x]` **Ch 24 — The Vacuum Has Structure** — Casimir (1948), Euler-Heisenberg
  - Casimir effect: $$F/A = -\frac{\pi^2 \hbar c}{240 d^4}$$ — attractive force between parallel plates from vacuum fluctuations
  - Euler-Heisenberg effective Lagrangian: vacuum becomes nonlinear in strong fields
  - Schwinger limit: $E_S = \frac{m_e^2 c^3}{e\hbar} \approx 1.3 \times 10^{18}$ V/m — pair production threshold
  - Vacuum birefringence: measured by PVLAS, observed in neutron star spectra
  - *Key idea:* the vacuum is not empty — it has measurable physical properties; this supports the idea of a polarizable vacuum

### Part V½: Nuclear Physics and Broken Symmetry
- `[x]` **Ch 25 — The Nuclear Force: How Binding Actually Works** — Rutherford (1911), Yukawa (1935)
  - Rutherford's nucleus: positive charge concentrated in 10⁻¹⁵ m
  - The binding problem: Coulomb repulsion between protons should blow nuclei apart
  - Yukawa's meson theory: massive mediator → short-range attraction
  - Yukawa potential: $V(r) = -g^2 e^{-r/r_0}/r$, range $r_0 = \hbar/(m_\pi c) \approx 1.4$ fm
  - Binding energy and mass deficit: $B = [\sum m_{\text{parts}} - m_{\text{nucleus}}]c^2$
  - Structural analogy to proposed EVO binding (Coulomb repulsion overcome by short-range attraction)
  - *Key idea:* like charges CAN be bound if a sufficiently strong short-range attraction exists; this is the pattern for K-enhanced Weber binding

- `[x]` **Ch 26 — Parity Violation and Broken Symmetry** — Lee & Yang (1956), Wu (1957), Prigogine
  - Parity (P): mirror symmetry of physical laws — assumed fundamental before 1957
  - The theta-tau puzzle: identical particles decaying with different parities
  - Lee & Yang: parity had never been TESTED in weak interactions — only assumed
  - Wu's experiment: cobalt-60 beta decay violates parity maximally
  - Every discrete symmetry tested has turned out broken (P, C, CP, T)
  - Prigogine: broken symmetry is the prerequisite for energy extraction from a reservoir
  - Connection to gauge symmetry: if "fundamental" symmetries like parity can break, gauge symmetry's exactness is an assumption, not a theorem
  - *Key idea:* symmetries are hypotheses to test, not axioms to assume; broken symmetry enables energy flow; the vacuum is an energy reservoir (Casimir); if gauge symmetry breaks (S ≠ 0), vacuum energy becomes accessible in principle

### Part VI: Potentials Are Physical (1959 – present)
- `[x]` **Ch 27 — The Aharonov-Bohm Effect** — Aharonov & Bohm (1959), Tonomura (1986)
  - Prediction: electron phase shifts from A even where B = 0
  - Tonomura's electron holography: definitive experimental confirmation
  - Three interpretive positions: gauge redundancy, potential realism, fiber bundles
  - *Key idea:* the vector potential A carries physical information that B does not encode locally; this undermines the claim that potentials are "just" mathematical

- `[x]` **Ch 28 — Topological Solitons** — Tony Skyrme (1961)
  - Skyrmions: topologically protected field configurations with integer winding number
  - $\pi_3(S^3) = \mathbb{Z}$: maps from space to the quaternion sphere are classified by integers
  - Topological charge cannot change continuously — stability without a confining potential
  - Connection to quaternions: unit quaternions form $S^3$, so quaternionic fields naturally support Skyrmions
  - *Key idea:* topology provides stability; field configurations with nonzero winding number are "knotted" and cannot unravel

### Part VIb: Plasma Physics and the Virtual Pair Picture (1920s – 1999)
- `[x]` **Ch 29 — Plasma Physics Fundamentals** — Langmuir (1920s), Debye, Tonks, Alfvén (1942)
  - Plasma = ionized gas: free electrons + ions responding to EM fields
  - Debye shielding: $\lambda_D = \sqrt{\varepsilon_0 k_B T/(n_e e^2)}$ — plasmas screen electric fields
  - Plasma frequency: $\omega_p = \sqrt{n_e e^2/(m_e \varepsilon_0)}$ — below this, waves are reflected
  - Plasma dielectric function: $\varepsilon(\omega) = \varepsilon_0(1 - \omega_p^2/\omega^2)$
  - Alfvén waves: transverse oscillations along magnetic field lines
  - Plasmas also support longitudinal modes (Langmuir waves, ion acoustic) — both transverse AND longitudinal, like QVED
  - Spheromaks and Taylor relaxation: plasma self-organizes into topologically linked field structures
  - *Key idea:* plasmas are the physical precedent for K field behavior; both transverse and longitudinal modes; self-organization into Hopf-linked topology is a real, observed plasma phenomenon

- `[x]` **Ch 30 — The Virtual Pair Plasma** — Schwinger (1951), Puthoff (1999)
  - QED vacuum as virtual $e^+e^-$ plasma: $n_{\text{virt}} \sim 10^{37}$ m⁻³, $\omega_{p,\text{vac}} \approx 7.8 \times 10^{20}$ Hz
  - Below $\omega_{p,\text{vac}}$ (all lab physics): virtual pairs polarize collectively → K field
  - Above $\omega_{p,\text{vac}}$: Schwinger pair production — vacuum "boils"
  - K field equation: Klein-Gordon with source term driven by field energy density
  - Plasma–vacuum isomorphism: mathematically exact, physically distinct (real vs virtual carriers)
  - Schwinger's source theory: vacuum as continuous polarizable medium (no virtual particles needed)
  - Vacuum energy density: $\sim 10^{23}$ J/m³ (Casimir proves it's real)
  - *Key idea:* K is the classical envelope of quantum vacuum polarization; the plasma analogy is exact in the simulator's regime; Schwinger provides rigorous QM foundation

### Part VII: The Unconventional Thread (1980s – present)
- `[x]` **Ch 31 — Longitudinal Forces Revisited** — Peter Graneau (1980s), André Assis (1994)
  - Graneau's experiments: exploding wires, mercury jets, water arc explosions
  - Forces along the current direction — inconsistent with Lorentz, consistent with Weber
  - Assis: Weber's force IS the retarded potential expanded to O(v²/c²)
  - The missing forces were always in Maxwell — hidden by gauge fixing
  - *Key idea:* experimental evidence for longitudinal EM forces that standard theory cannot produce; Weber and Maxwell agree when you keep the full potential structure

- `[x]` **Ch 32 — Electrogravitics and Anomalous Forces** — Thomas Townsend Brown (1920s–1960s)
  - Brown's experiments: asymmetric capacitors under high voltage show apparent thrust
  - The Biefeld-Brown effect
  - Controversy: ion wind vs genuine electrogravitic coupling
  - Connection to vacuum polarization: if K(x) varies, local c varies → effective refractive index gradient → radiation pressure asymmetry
  - *Key idea:* anomalous forces in strong-field EM remain unexplained; the vacuum polarization framework provides a candidate mechanism (but quantitative predictions are needed)

- `[x]` **Ch 33 — Exotic Vacuum Objects** — Ken Shoulders (1990s)
  - EVO (Exotic Vacuum Object) observations: coherent charge clusters (10⁶–10¹¹ electrons) that self-bind
  - Micron-scale, propagate through materials, leave characteristic tracks
  - Standard EM says they should explode from Coulomb repulsion
  - Proposed mechanism: K > 1 pocket → enhanced Weber attraction → self-binding
  - *Key idea:* experimental observations of charge clusters that defy standard EM; the QVED framework offers a binding mechanism through vacuum polarization feedback

- `[x]` **Ch 34 — Extended Maxwell Equations and the S Field** — van Vlaenderen & Waser (2001), Whittaker (1903)
  - Whittaker's scalar potential decomposition (1903): any EM field decomposable into scalar waves
  - Van Vlaenderen's extended Maxwell equations (the ones in your simulation):
    $$\nabla\cdot\mathbf{E} = \frac{\rho}{\varepsilon_0} - \frac{\partial S}{\partial t}$$
    $$\nabla\times\mathbf{B} = \mu_0\mathbf{J} + \mu_0\varepsilon_0\frac{\partial\mathbf{E}}{\partial t} + \nabla S$$
  - The scalar wave equation: $\Box S = -\rho/\varepsilon_0$
  - The polarizable vacuum: $K(x) \geq 1$, $c_{\text{local}} = c_0/K$
  - *Key idea:* removing the gauge constraint restores the scalar degree of freedom; this is what your simulation computes in extended_mode

- `[x]` **Ch 35 — Topological EM: Hopfions and Ball Lightning** — Rañada, Irvine & Bouwmeester (2008)
  - Hopf fibration: every pair of field lines links exactly once
  - Hopfion solutions to Maxwell's equations (exact but unstable in linear EM)
  - QVED nonlinearity (K field feedback) as stabilization mechanism
  - Ball lightning as a self-sustaining Hopfion in (Q, K)
  - *Key idea:* topology + nonlinearity = stable electromagnetic solitons; this is the ultimate target of your simulation

---

## Status: COMPLETE

All 35 chapters written. Document: `study_guide.md` (~27,000 words).

## Rendering to PDF

```bash
# With Pandoc + LaTeX installed:
pandoc study_guide.md -o study_guide.pdf --pdf-engine=xelatex -V geometry:margin=1in -V fontsize=11pt

# Or use Typora/Obsidian export
```
