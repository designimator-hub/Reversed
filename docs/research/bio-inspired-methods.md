# Methods survey: topology optimization, lattices, fractals, molecular and protein mechanics

This document collects the mathematics this software can build on. Each section covers the governing equations, what the method buys you in practice, and its limits. References are at the end and numbered [n].

Notation: ρ̄ = relative density (solid volume / total volume), E_s, ν_s = Young's modulus and Poisson's ratio of the solid, x_e = design density of element e.

---

## 1. Baseline: density-based topology optimization (the engine everything else plugs into)

### 1.1 SIMP compliance minimization [1][2][3]

$$
\min_{\mathbf{x}} \; c(\mathbf{x}) = \mathbf{U}^T\mathbf{K}(\mathbf{x})\mathbf{U}
\quad \text{s.t.} \quad \mathbf{K}\mathbf{U}=\mathbf{F},\;\; \frac{V(\mathbf{x})}{V_0}\le f,\;\; 0\le x_e\le 1
$$

Element stiffness interpolation, with penalty p ≈ 3:

$$
E_e(x_e) = E_{\min} + x_e^{\,p}\,(E_s - E_{\min})
$$

Sensitivity (the problem is self-adjoint, so no extra solve):

$$
\frac{\partial c}{\partial x_e} = -p\,x_e^{\,p-1}(E_s - E_{\min})\,\mathbf{u}_e^T\mathbf{k}_0\mathbf{u}_e
$$

### 1.2 Filtering and projection [4][5]

Density filter with radius r_min (sets minimum member size, removes checkerboards):

$$
\tilde{x}_e = \frac{\sum_i w_{ei}\,x_i}{\sum_i w_{ei}}, \qquad w_{ei}=\max(0,\,r_{\min}-\lVert \mathbf{c}_e-\mathbf{c}_i\rVert)
$$

Heaviside projection pushes the result toward crisp 0/1, with β increased in steps:

$$
\bar{x}_e = \frac{\tanh(\beta\eta)+\tanh(\beta(\tilde{x}_e-\eta))}{\tanh(\beta\eta)+\tanh(\beta(1-\eta))}
$$

Running eroded, intermediate and dilated projections (η = 0.6, 0.5, 0.4) together makes the design robust to manufacturing over- or under-size [5].

### 1.3 Optimizers

- Optimality criteria (OC): a fixed-point update, cheap, and good enough for one volume constraint.
- MMA, the method of moving asymptotes [6]: the standard for several constraints (stress, overhang, frequency).

### 1.4 Level-set alternative [7]

The boundary is the zero contour of φ, moved with normal velocity V_n taken from shape sensitivity:

$$
\frac{\partial \phi}{\partial t} + V_n\,\lvert\nabla\phi\rvert = 0
$$

Level sets give clean boundaries, and they map directly to an SDF kernel like PicoGK. Density methods converge more reliably. Recommendation: optimize with density, then convert to an SDF and finish with a few level-set steps for a clean surface.

### 1.5 Stress constraint (p-norm aggregation)

$$
\sigma_{PN} = \Big(\sum_e (\bar{x}_e^{\,q}\,\sigma^{vM}_e)^P\Big)^{1/P}\le \sigma_{\text{allow}}, \quad P\approx 8\text{–}12
$$

### 1.6 GPU scale

Wu, Dick and Westermann solved TO with GPU geometric multigrid at millions of voxels on one card [8]. Aage et al. reached over a billion voxels on a CPU cluster and produced bird-bone-like internal structure in an aircraft wing [9]. A desktop GPU with 16–24 GB handles roughly 10–50 million voxels with matrix-free storage. That is about 10–100x the resolution Fusion typically returns per outcome.

---

## 2. Bone: nature's own topology optimizer

### 2.1 Strain-energy remodeling [10][11]

Bone adds mass where the strain-energy density per unit mass exceeds a set point k, and removes it where it falls below:

$$
\frac{d\rho}{dt} = B\left(\frac{U}{\rho} - k\right), \qquad E = C\rho^{\gamma}
$$

Weinans et al. found this model converges to discrete, checkerboard-like patterns [11]. That is the same thing unfiltered SIMP does, and for the same mathematical reason: run with the OC update, the remodeling rule is nearly identical to SIMP. That is good news. The "bio-inspired" method and the engineering method agree, which is evidence the optimum is real.

Huiskes et al. added sensing cells (osteocytes) that signal over a distance, which produced trabecular architecture that aligns with load and re-aligns when the load changes [10]. In engineering terms that signal decay is a filter kernel.

### 2.2 Bone-like infill [12]

Wu, Aage, Westermann and Sigmund constrain the *local* volume fraction instead of the global one:

$$
\bar{\rho}_e^{\,\text{local}} = \frac{1}{N_e}\sum_{i\in\Omega_e}\bar{x}_i \;\le\; \alpha, \quad\text{aggregated by a p-norm}
$$

This forces porous, trabecular-like internal structure that is far more robust to unexpected loads and local damage than a solid-member Fusion outcome, at a small cost in stiffness for the design load case. It should be one switch in the UI: "Shell and trabecular core."

### 2.3 Mattheck's methods: CAO, SKO, tension triangles [13][14]

- CAO (computer-aided optimization): a fictitious thermal expansion in a surface layer, proportional to σ_vM − σ_ref, grows material where stress is high and shrinks it where stress is low. This is how trees grow at notches.
- SKO (soft kill option): set element modulus proportional to its stress and iterate. A crude topology optimizer.
- Tension triangles: a fillet you can draw by hand. Build a 45° isosceles triangle at the notch, take the midpoint of its hypotenuse, build a second triangle from there, repeat about three times, then round off. The resulting contour cuts notch stress about as well as a numerically optimized one.

Product use: an "organic fillet" tool that runs on any sharp internal corner in seconds, no FE needed. Many users would use this more than the full optimizer.

### 2.4 Principal stress lines and Michell structures [15][16]

Michell proved that minimum-weight trusses follow orthogonal families of principal stress lines [15]. Trabecular bone in the femoral head follows the same lines. Arora et al. built volumetric "Michell trusses" by tracing principal stress directions from one FE solve and placing beams along them [16].

Product use: a one-solve preview that traces stress lines from the initial solid within seconds, before full TO finishes. It also gives the lattice orientation field for Section 4.

---

## 3. Transport networks: slime mold, leaves, vessels

### 3.1 Physarum adaptive network [17]

Edges carry flux Q from a pressure difference, and conductivity D grows with use and decays without it:

$$
Q_{ij} = \frac{D_{ij}}{L_{ij}}(p_i-p_j), \qquad \sum_j Q_{ij} = s_i, \qquad \frac{dD_{ij}}{dt} = f(\lvert Q_{ij}\rvert) - \gamma D_{ij}
$$

With f(Q) = |Q|^μ/(1+|Q|^μ), Tero et al. reproduced the Tokyo rail network's cost, efficiency and fault tolerance [17]. μ > 1 gives tree-like networks (cheap). μ < 1 keeps loops (fault tolerant).

Product use: replace pressure with load and flux with force, and you get a fast load-path generator on a graph. It runs in milliseconds and makes a good interactive sketch layer for ribs and frames. Real stiffness still comes from the FE-based TO. Treat Physarum as a preview and as an initial guess for the optimizer.

### 3.2 Space colonization (leaf venation, tree branching) [18]

Scatter attraction points in the region to be served. Each vein node grows toward the mean direction of the points nearest it. A point is removed once a vein comes within a kill distance d_k.

$$
\mathbf{n}_v = \frac{\sum_{s\in S(v)} \frac{\mathbf{s}-\mathbf{v}}{\lVert\mathbf{s}-\mathbf{v}\rVert}}{\Big\lVert\sum_{s\in S(v)} \frac{\mathbf{s}-\mathbf{v}}{\lVert\mathbf{s}-\mathbf{v}\rVert}\Big\rVert}, \qquad \mathbf{v}' = \mathbf{v} + D\,\mathbf{n}_v
$$

Product use: tree supports for metal AM, branching cooling channels, cable routing.

### 3.3 Murray's law and constructal design [19][20]

For laminar flow at minimum pumping-plus-volume cost, parent and daughter channel radii satisfy:

$$
r_0^3 = \sum_i r_i^3
$$

Bejan's constructal law gives the same tree-like result for heat conduction paths [20].

Product use: size every branch of a cooling or fluid network automatically. This is a constraint you can apply to the Section 3.2 output at almost no cost.

---

## 4. Lattice mechanics: the rules that decide whether a lattice is good

### 4.1 Gibson–Ashby scaling [21]

$$
\frac{E^*}{E_s} = C_1\,\bar{\rho}^{\,n}, \qquad \frac{\sigma^*}{\sigma_{ys}} = C_2\,\bar{\rho}^{\,m}
$$

| Deformation mode | n (stiffness) | m (strength) | Example |
|---|---|---|---|
| Bending-dominated | 2 | 1.5 | open foam, BCC, Kelvin cell, diamond |
| Stretch-dominated | 1 | 1 | octet truss, plate lattices |

At ρ̄ = 0.1 a stretch-dominated lattice is roughly 10x stiffer than a bending-dominated one of the same weight. That single fact matters more than any choice of pretty unit cell.

### 4.2 Maxwell's rule and the Calladine extension [22][23][24]

For a pin-jointed frame with b struts and j joints:

$$
\text{2D: } M = b - 2j + 3, \qquad \text{3D: } M = b - 3j + 6, \qquad M = s - m
$$

where s is the number of states of self-stress and m the number of mechanisms [23]. Deshpande, Ashby and Fleck showed that a lattice is stretch-dominated when M ≥ 0 and every node has connectivity Z ≈ 12 in 3D (Z ≈ 6 in 2D) [22]. The octet truss (Z = 12) meets this, with E*/E_s ≈ ρ̄/9 [24].

Product use: compute M and Z for every lattice the software generates and show the user "stretch / bending" with the expected scaling exponent. No commercial tool shows this today.

### 4.3 The upper bound: Hashin–Shtrikman [25]

For an isotropic porous material the bulk modulus cannot exceed:

$$
K^* \le K_{HS} = \frac{4G_sK_s\,\bar{\rho}}{4G_s + 3K_s(1-\bar{\rho})}
$$

Truss lattices fall short of this bound. Berger, Wadley and McMeeking showed that closed-cell *plate* lattices (cubic plus octet plates) reach it [25].

Product use: plot every outcome's stiffness against the HS bound. It gives the user an absolute score, not only a ranking.

### 4.4 Triply periodic minimal surfaces (TPMS) [26][27]

Level-set approximations (cell size L, k = 2π/L, offset t controls density):

| Surface | Implicit function |
|---|---|
| Gyroid | sin kx cos ky + sin ky cos kz + sin kz cos kx = t |
| Schwarz P | cos kx + cos ky + cos kz = t |
| Schwarz D (diamond) | sin kx sin ky sin kz + sin kx cos ky cos kz + cos kx sin ky cos kz + cos kx cos ky sin kz = t |
| Neovius | 3(cos kx + cos ky + cos kz) + 4 cos kx cos ky cos kz = t |

Sheet TPMS (|f| ≤ t) are smooth, self-supporting for AM, free of stress concentrations at nodes, and drain powder well. Because they are implicit functions, they slot straight into an SDF kernel. Grading is just letting t and k vary in space: t(x) from the TO density field, and the local frame from principal stress directions (Section 2.4).

### 4.5 Homogenization and de-homogenization: the biggest lever [28][29]

Compute the effective elasticity of a unit cell Y once by solving the periodic cell problem:

$$
E^H_{ijkl} = \frac{1}{|Y|}\int_Y \big(\varepsilon^{0(ij)}_{pq} - \varepsilon^{(ij)}_{pq}(\chi)\big)\,E_{pqrs}\,\big(\varepsilon^{0(kl)}_{rs} - \varepsilon^{(kl)}_{rs}(\chi)\big)\,dY
$$

Tabulate E^H against each cell's parameters (density, aspect ratios). Then run TO on a coarse mesh where each element picks its own cell parameters and orientation. Finally, "de-homogenize": build a fine, continuous, oriented lattice that follows the optimized fields [29]. Groen and Sigmund got results within a few percent of full-resolution TO at a fraction of the cost [29].

Product use: this is how you beat Fusion at its own game. You get a macro shape plus a graded, stress-aligned lattice interior from one optimization, at a resolution no single-scale optimizer reaches.

---

## 5. Molecular mechanics applied to lattice generation

### 5.1 The force-field analogy

A classical molecular force field (AMBER/CHARMM form):

$$
E = \sum_{\text{bonds}} k_b(r-r_0)^2 + \sum_{\text{angles}} k_\theta(\theta-\theta_0)^2 + \sum_{\text{dihedrals}} \frac{V_n}{2}\big[1+\cos(n\varphi-\gamma)\big] + \sum_{i<j}\left[\frac{A_{ij}}{r_{ij}^{12}}-\frac{B_{ij}}{r_{ij}^{6}}\right] + \sum_{i<j}\frac{q_iq_j}{4\pi\varepsilon r_{ij}}
$$

| Molecular term | Lattice meaning |
|---|---|
| Bond stretch k_b | strut axial stiffness EA/L |
| Angle bend k_θ | joint rotational stiffness (pin vs rigid node), sets bending vs stretching behavior |
| Dihedral | torsion and out-of-plane stiffness of nodes |
| Lennard-Jones | node spacing and repulsion, used for *generating* the layout |
| Coulomb | long-range spacing control, used for density grading |

### 5.2 Crystal structures as a unit-cell library

Crystallography has already catalogued the stable, symmetric networks, and the Maxwell analysis carries over directly:

| Crystal | Node coordination Z | Lattice behavior |
|---|---|---|
| Diamond cubic | 4 | bending-dominated, very compliant, good for energy absorption |
| Simple cubic | 6 | stiff on axis, weak in shear (anisotropic) |
| BCC | 8 | bending-dominated |
| FCC (octet truss) | 12 | stretch-dominated, the benchmark |
| Kelvin (truncated octahedron) | 4 | bending-dominated, minimal surface area foam |
| Weaire–Phelan | 4 | lower surface area than Kelvin, near-isotropic foam |

Combining an FCC network with a simple cubic one (as in "octet + cube") gives a stiffness response that is close to isotropic.

### 5.3 Conformal lattice generation by particle relaxation [30][31]

The idea: treat lattice nodes as atoms in a box shaped like the part.

1. Seed N particles inside the design domain.
2. Give each particle a target spacing r_0(x) taken from the TO density or stress field (dense where stress is high).
3. Relax with a pair potential (Lennard-Jones or a soft "bubble" spring) plus a wall potential from the part's SDF that pins nodes onto the surface.
4. Minimize with FIRE [30], a molecular-dynamics minimizer with velocity mixing. When P = F·v > 0: v ← (1−α)v + α|v|F̂, and the time step grows. When P ≤ 0: v ← 0, the time step halves, and α resets.
5. Connect nodes by Delaunay tetrahedralization (gives Z ≈ 12–14, stretch-dominated) or by its dual Voronoi (Z = 4, bending-dominated, foam-like).

Shimada and Gossard used this same physics ("bubble mesh") to build high-quality meshes in 1995 [31], so the approach has a long, proven track record.

Why it is better than a periodic lattice clipped by the part boundary: no cut struts at the skin, nodes land on the surface, spacing follows stress, and you can turn the layout into a stretch- or bending-dominated lattice on purpose.

### 5.4 Form-finding with the force density method [32]

For a network with branch–node matrix C and force densities q (force/length) on the diagonal of Q, equilibrium is linear:

$$
\mathbf{D}=\mathbf{C}^T\mathbf{Q}\mathbf{C}, \qquad \mathbf{D}_{ff}\,\mathbf{x}_f = \mathbf{p}_f - \mathbf{D}_{fx}\,\mathbf{x}_x
$$

It takes one sparse linear solve per coordinate direction. Product use: tension-only and compression-only shells, cable nets, tensegrity, and funicular arches (pure compression, no bending, as in Gaudí's hanging chain models).

### 5.5 Auxetics [33]

Re-entrant and chiral cells give negative Poisson's ratio. Lakes made the first auxetic foam in 1987 [33]. Uses: indentation resistance, shells that curve both ways (synclastic), and energy absorption. In the force-field picture, auxetic behavior comes from angle terms that hinge open under tension.

---

## 6. Protein mechanics: what transfers and what does not

Proteins are chains of about 20 residue types that fold into one shape and then move in a few specific, low-energy ways. Four ideas transfer to engineering design. Direct "folding" of parts does not, except as an analogy for origami and kirigami.

### 6.1 Elastic network models: find a structure's soft modes [34][35]

The Anisotropic Network Model treats each residue as a node joined by identical springs to every neighbor within a cutoff r_c:

$$
V = \frac{\gamma}{2}\sum_{r^0_{ij}<r_c}\big(r_{ij}-r^0_{ij}\big)^2, \qquad \mathbf{H}_{ij} = -\frac{\gamma}{(r^0_{ij})^2}\,\mathbf{d}_{ij}\mathbf{d}_{ij}^T
$$

The lowest eigenvectors of the Hessian H match the functional motions measured in experiments, even though the model knows nothing of chemistry. The Gaussian Network Model gives fluctuations directly: ⟨ΔR_i·ΔR_j⟩ = (3k_BT/γ)[Γ⁻¹]_ij [35].

Product use:
- A "weak spot" view: run a fast modal analysis on any outcome or lattice and show its softest modes. You find buckling-prone and wobbly regions in seconds.
- An ENM is a coarse model of any lattice, so you can screen thousands of candidate topologies before running full FE.

### 6.2 Allostery-inspired mechanism design [36][37]

In allostery, a molecule binding at one site changes the shape of a distant site. Rocks et al. reproduced this in spring networks: removing a few chosen bonds from a disordered network makes a strain applied at a source produce a targeted strain at a distant target [36]. Goodrich, Liu and Nagel showed that because each bond contributes independently to bulk and shear moduli, pruning can set Poisson's ratio anywhere from about +0.5 to −1 [37].

Product use: compliant mechanisms and morphing parts, a market Fusion does not serve. The user marks an input region and an output region with the wanted motion, and the solver prunes a lattice to deliver it. Grippers, flexures, adaptive seating and shoe midsoles fit here.

### 6.3 Strain stiffening: the worm-like chain [38]

Collagen, skin and spider silk are soft at small strain and stiffen sharply as they stretch (a J-shaped curve). The worm-like chain model gives the force on a polymer of contour length L and persistence length P at extension x:

$$
\frac{FP}{k_BT} = \frac{1}{4(1-x/L)^2} - \frac14 + \frac{x}{L}
$$

Engineering version: wavy or helical struts that straighten under load. Small loads meet low stiffness (comfort, impact tolerance), and large loads meet high stiffness (protection). You choose the waviness amplitude and you place the transition strain. Product use: helmets, padding, orthotics, soft-robot skins.

### 6.4 Toughness: sacrificial bonds, hidden length, nanoconfinement [39][40][41]

- Bone: sacrificial bonds that break and re-form release "hidden length" and dissipate energy before the main structure fails [39].
- Spider silk: β-sheet crystals confined to a few nanometres are stiffer, stronger and tougher than larger ones [40].
- Size-dependent flaw tolerance: below a critical size, a structure becomes insensitive to flaws [41]:

$$
h^* \approx \alpha^2\,\frac{\gamma E}{\sigma_{th}^2}
$$

Engineering version: struts with built-in slack loops or frangible links that fail in sequence. The structure keeps load path redundancy while each link that fails absorbs energy. Product use: crash structures and one-time energy absorbers designed for a target force–displacement curve.

### 6.5 Energy landscapes and search: why folding is fast [42][43]

Proteins fold fast because evolution shaped their energy landscape into a funnel with little frustration [42]. Lattice and topology design problems are rugged, with many local optima. Replica exchange (parallel tempering), the standard sampler for rugged protein landscapes, transfers directly. Run R copies at temperatures T_1 < … < T_R and swap neighbors with probability:

$$
P_{\text{swap}} = \min\Big(1,\; \exp\big[(\beta_i-\beta_j)(E_i-E_j)\big]\Big)
$$

Product use: Fusion offers "multiple outcomes," but they come mostly from changing materials and manufacturing methods. Replica exchange over topology choices gives genuinely different, near-optimal designs, with a diversity measure the user can see.

### 6.6 Generative models: the AlphaFold and RFdiffusion lesson [44][45][46]

RFdiffusion generates new protein backbones by denoising, conditioned on a target function [45]. TopoDiff applied diffusion models to topology optimization and beat GANs on manufacturability and compliance [46].

Honest limits: these models need tens of thousands of solved examples. The Section 1 solver can produce that dataset itself overnight on one GPU. Use: instant first guesses (well under a second) that the physics solver then refines. The physics stays in charge, and the model only shortens the run.

---

## 7. Fractals and hierarchy

### 7.1 Generators

- L-systems [47]: string rewriting with a turtle interpreter. Example: F → F[+F]F[−F]F. Good for branching supports and trees.
- Iterated function systems [48]: a set of contractive affine maps w_i whose attractor is A = ∪ w_i(A). Good for self-similar lattices and textures.
- Phyllotaxis (Vogel) [49]: r_n = c√n, θ_n = n · 137.508°. Places points evenly on discs and surfaces, for perforations, pins and lattice nodes.
- Space-filling curves (Hilbert, Peano): continuous toolpaths and infill with local area control.

### 7.2 Measurement

Box-counting dimension, D = lim_{ε→0} log N(ε) / log(1/ε). Useful as a complexity metric to report, and to cap complexity for manufacturability.

### 7.3 What hierarchy does and does not buy [50][51][52]

- Stiffness: hierarchy does *not* improve linear stiffness per unit weight. Relative densities multiply across levels, so a bending-dominated hierarchy still scales as ρ̄² overall [50].
- Buckling-limited strength: hierarchy helps a great deal at low load per unit size. Rayneau-Kirkhope et al. showed hierarchical (fractal) frames from hollow tubes need much less material than conventional designs when buckling governs [51].
- Damage tolerance: Euplectella, the glass sponge, uses seven levels of hierarchy [52]. Its double diagonal bracing pattern has been shown to raise lattice buckling strength well above standard square lattices of the same mass [53].

Product use: apply fractal and hierarchical subdivision where analysis says *buckling* governs (long slender members, thin shells under compression), not as decoration everywhere.

---

## 8. Pattern formation: Turing reaction–diffusion [54][55]

Gray–Scott model:

$$
\frac{\partial u}{\partial t} = D_u\nabla^2u - uv^2 + F(1-u), \qquad \frac{\partial v}{\partial t} = D_v\nabla^2v + uv^2 - (F+k)v
$$

Changing (F, k) gives spots, stripes, labyrinths or holes [55]. Make the diffusion anisotropic, aligned with principal stress, and the stripes follow load paths. That gives fiber-like reinforcement patterns for shells, ribbing and skins. The cost is one GPU stencil per step, so it runs live.

---

## 9. Synthesis: where this beats Fusion generative design

| Capability | Fusion today | This project | Methods |
|---|---|---|---|
| Resolution | coarse, cloud-limited | 10–50M voxels on a desktop GPU | 1.6 |
| Interior | solid members | shell plus graded, stress-aligned lattice or trabecular core | 2.2, 4.4, 4.5 |
| Lattice quality | none in GD | Maxwell/Z check, HS-bound score, conformal nodes | 4.2, 4.3, 5.3 |
| Live feedback | wait for cloud | stress lines and Physarum preview in seconds | 2.4, 3.1 |
| Fillets | manual | tension-triangle tool | 2.3 |
| Diversity | material × process grid | replica exchange, diversity metric | 6.5 |
| Mechanisms | not offered | allostery-style pruning | 6.2 |
| Nonlinear targets | linear static | strain stiffening, energy absorption curves | 6.3, 6.4 |
| Fluid/thermal networks | limited | space colonization plus Murray sizing | 3.2, 3.3 |
| Robustness | single load case focus | eroded/dilated projection, local volume constraint | 1.2, 2.2 |

What makes it more intuitive:

1. The user paints intent instead of filling forms. Load, support, keep, avoid, "soft here," "stiff there," "this moves that."
2. Feedback within seconds (stress lines, Physarum paths), with full optimization filling in behind it.
3. Every outcome carries plain numbers: mass, stiffness, percent of theoretical bound, safety factor, stretch or bending, estimated print time.
4. Output is an SDF the user can keep editing with booleans, offsets and smoothing that never fail, in the PicoGK manner.

---

## 10. Build order (revised)

1. Solver core on the compute node: 3D voxel SIMP, GPU multigrid CG, filtering and projection, OC/MMA, job spec, web viewport. Validate against Liu and Tovar's top3d benchmarks [3].
2. SDF kernel and output: density to SDF, level-set cleanup, TPMS and Mattheck fillets, STL/3MF and OpenVDB export. Blender add-on and Fusion add-in as clients.
3. Lattice layer: unit-cell library with homogenized tables, Maxwell/Z/HS scoring, particle-relaxed conformal lattices, stress-aligned grading.
4. Multiscale: local volume constraint (bone-like infill), then homogenization-based TO with de-homogenization.
5. Previews and search: stress-line tracing, Physarum, space colonization with Murray sizing, replica-exchange outcome generation, ENM soft-mode view.
6. Nonlinear and mechanisms: allostery-style pruning, wavy-strut strain stiffening, energy-absorption targets (needs a nonlinear FE path).
7. Learned acceleration: train a diffusion model on the solver's own output for instant first guesses.

---

## References

1. Bendsøe, M.P., Kikuchi, N. (1988). Generating optimal topologies in structural design using a homogenization method. *Comput. Methods Appl. Mech. Eng.* 71, 197–224.
2. Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B.S., Sigmund, O. (2011). Efficient topology optimization in MATLAB using 88 lines of code. *Struct. Multidisc. Optim.* 43, 1–16.
3. Liu, K., Tovar, A. (2014). An efficient 3D topology optimization code written in Matlab. *Struct. Multidisc. Optim.* 50, 1175–1196.
4. Bourdin, B. (2001). Filters in topology optimization. *Int. J. Numer. Methods Eng.* 50, 2143–2158.
5. Wang, F., Lazarov, B.S., Sigmund, O. (2011). On projection methods, convergence and robust formulations in topology optimization. *Struct. Multidisc. Optim.* 43, 767–784.
6. Svanberg, K. (1987). The method of moving asymptotes: a new method for structural optimization. *Int. J. Numer. Methods Eng.* 24, 359–373.
7. Wang, M.Y., Wang, X., Guo, D. (2003). A level set method for structural topology optimization. *Comput. Methods Appl. Mech. Eng.* 192, 227–246.
8. Wu, J., Dick, C., Westermann, R. (2016). A system for high-resolution topology optimization. *IEEE Trans. Vis. Comput. Graph.* 22(3), 1195–1208.
9. Aage, N., Andreassen, E., Lazarov, B.S., Sigmund, O. (2017). Giga-voxel computational morphogenesis for structural design. *Nature* 550, 84–86.
10. Huiskes, R., Ruimerman, R., van Lenthe, G.H., Janssen, J.D. (2000). Effects of mechanical forces on maintenance and adaptation of form in trabecular bone. *Nature* 405, 704–706.
11. Weinans, H., Huiskes, R., Grootenboer, H.J. (1992). The behavior of adaptive bone-remodeling simulation models. *J. Biomech.* 25, 1425–1441.
12. Wu, J., Aage, N., Westermann, R., Sigmund, O. (2018). Infill optimization for additive manufacturing: approaching bone-like porous structures. *IEEE Trans. Vis. Comput. Graph.* 24(2), 1127–1140.
13. Mattheck, C., Burkhardt, S. (1990). A new method of structural shape optimization based on biological growth. *Int. J. Fatigue* 12, 185–190.
14. Mattheck, C. (1998). *Design in Nature: Learning from Trees.* Springer.
15. Michell, A.G.M. (1904). The limits of economy of material in frame-structures. *Phil. Mag.* 8, 589–597.
16. Arora, R., Jacobson, A., Langlois, T.R., Huang, Y., Mueller, C., Matusik, W., Shamir, A., Singh, K., Levin, D.I.W. (2019). Volumetric Michell trusses for parametric design & fabrication. *Proc. ACM Symp. Computational Fabrication (SCF '19).*
17. Tero, A., Takagi, S., Saigusa, T., Ito, K., Bebber, D.P., Fricker, M.D., Yumiki, K., Kobayashi, R., Nakagaki, T. (2010). Rules for biologically inspired adaptive network design. *Science* 327, 439–442.
18. Runions, A., Fuhrer, M., Lane, B., Federl, P., Rolland-Lagan, A.-G., Prusinkiewicz, P. (2005). Modeling and visualization of leaf venation patterns. *ACM Trans. Graph.* 24(3), 702–711.
19. Murray, C.D. (1926). The physiological principle of minimum work: I. The vascular system and the cost of blood volume. *PNAS* 12, 207–214.
20. Bejan, A. (1997). Constructal-theory network of conducting paths for cooling a heat generating volume. *Int. J. Heat Mass Transfer* 40, 799–816.
21. Gibson, L.J., Ashby, M.F. (1997). *Cellular Solids: Structure and Properties*, 2nd ed. Cambridge University Press.
22. Deshpande, V.S., Ashby, M.F., Fleck, N.A. (2001). Foam topology: bending versus stretching dominated architectures. *Acta Mater.* 49, 1035–1040.
23. Calladine, C.R. (1978). Buckminster Fuller's "tensegrity" structures and Clerk Maxwell's rules for the construction of stiff frames. *Int. J. Solids Struct.* 14, 161–172.
24. Deshpande, V.S., Fleck, N.A., Ashby, M.F. (2001). Effective properties of the octet-truss lattice material. *J. Mech. Phys. Solids* 49, 1747–1769.
25. Berger, J.B., Wadley, H.N.G., McMeeking, R.M. (2017). Mechanical metamaterials at the theoretical limit of isotropic elastic stiffness. *Nature* 543, 533–537.
26. Schoen, A.H. (1970). Infinite periodic minimal surfaces without self-intersections. NASA Technical Note D-5541.
27. Al-Ketan, O., Abu Al-Rub, R.K. (2019). Multifunctional mechanical metamaterials based on triply periodic minimal surface lattices. *Adv. Eng. Mater.* 21, 1900524.
28. Bendsøe, M.P., Sigmund, O. (2003). *Topology Optimization: Theory, Methods and Applications.* Springer.
29. Groen, J.P., Sigmund, O. (2018). Homogenization-based topology optimization for high-resolution manufacturable microstructures. *Int. J. Numer. Methods Eng.* 113, 1148–1163.
30. Bitzek, E., Koskinen, P., Gähler, F., Moseler, M., Gumbsch, P. (2006). Structural relaxation made simple. *Phys. Rev. Lett.* 97, 170201.
31. Shimada, K., Gossard, D.C. (1995). Bubble mesh: automated triangular meshing of non-manifold geometry by sphere packing. *Proc. 3rd ACM Symp. Solid Modeling and Applications*, 409–419.
32. Schek, H.-J. (1974). The force density method for form finding and computation of general networks. *Comput. Methods Appl. Mech. Eng.* 3, 115–134.
33. Lakes, R. (1987). Foam structures with a negative Poisson's ratio. *Science* 235, 1038–1040.
34. Atilgan, A.R., Durell, S.R., Jernigan, R.L., Demirel, M.C., Keskin, O., Bahar, I. (2001). Anisotropy of fluctuation dynamics of proteins with an elastic network model. *Biophys. J.* 80, 505–515.
35. Bahar, I., Atilgan, A.R., Erman, B. (1997). Direct evaluation of thermal fluctuations in proteins using a single-parameter harmonic potential. *Fold. Des.* 2, 173–181.
36. Rocks, J.W., Pashine, N., Bischofberger, I., Goodrich, C.P., Liu, A.J., Nagel, S.R. (2017). Designing allostery-inspired response in mechanical networks. *PNAS* 114, 2520–2525.
37. Goodrich, C.P., Liu, A.J., Nagel, S.R. (2015). The principle of independent bond-level response: tuning by pruning to exploit disorder for global behavior. *Phys. Rev. Lett.* 114, 225501.
38. Marko, J.F., Siggia, E.D. (1995). Stretching DNA. *Macromolecules* 28, 8759–8770.
39. Fantner, G.E., et al. (2005). Sacrificial bonds and hidden length dissipate energy as mineralized fibrils separate during bone fracture. *Nat. Mater.* 4, 612–616.
40. Keten, S., Xu, Z., Ihle, B., Buehler, M.J. (2010). Nanoconfinement controls stiffness, strength and mechanical toughness of β-sheet crystals in silk. *Nat. Mater.* 9, 359–367.
41. Gao, H., Ji, B., Jäger, I.L., Arzt, E., Fratzl, P. (2003). Materials become insensitive to flaws at nanoscale: lessons from nature. *PNAS* 100, 5597–5600.
42. Bryngelson, J.D., Wolynes, P.G. (1987). Spin glasses and the statistical mechanics of protein folding. *PNAS* 84, 7524–7528.
43. Hukushima, K., Nemoto, K. (1996). Exchange Monte Carlo method and application to spin glass simulations. *J. Phys. Soc. Jpn.* 65, 1604–1608.
44. Jumper, J., et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589.
45. Watson, J.L., et al. (2023). De novo design of protein structure and function with RFdiffusion. *Nature* 620, 1089–1100.
46. Mazé, F., Ahmed, F. (2023). Diffusion models beat GANs on topology optimization. *Proc. AAAI Conf. Artificial Intelligence* 37(8), 9108–9116.
47. Prusinkiewicz, P., Lindenmayer, A. (1990). *The Algorithmic Beauty of Plants.* Springer.
48. Barnsley, M.F. (1988). *Fractals Everywhere.* Academic Press.
49. Vogel, H. (1979). A better way to construct the sunflower head. *Math. Biosci.* 44, 179–189.
50. Lakes, R. (1993). Materials with structural hierarchy. *Nature* 361, 511–515.
51. Rayneau-Kirkhope, D., Mao, Y., Farr, R. (2012). Ultralight fractal structures from hollow tubes. *Phys. Rev. Lett.* 109, 204301.
52. Aizenberg, J., Weaver, J.C., Thanawala, M.S., Sundar, V.C., Morse, D.E., Fratzl, P. (2005). Skeleton of *Euplectella* sp.: structural hierarchy from the nanoscale to the macroscale. *Science* 309, 275–278.
53. Fernandes, M.C., Aizenberg, J., Weaver, J.C., Bertoldi, K. (2021). Mechanically robust lattices inspired by deep-sea glass sponges. *Nat. Mater.* 20, 237–241.
54. Turing, A.M. (1952). The chemical basis of morphogenesis. *Phil. Trans. R. Soc. Lond. B* 237, 37–72.
55. Pearson, J.E. (1993). Complex patterns in a simple system. *Science* 261, 189–192.
56. Fratzl, P., Weinkamer, R. (2007). Nature's hierarchical materials. *Prog. Mater. Sci.* 52, 1263–1334.
