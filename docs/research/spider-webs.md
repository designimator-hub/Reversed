# Spider webs as a source of structural methods

This document surveys what spider webs and spider silk can teach a generative design engine, with particular attention to funnel webs. Each section gives the governing equations where they exist, what is actually known (and from which paper), and where the evidence runs out. References are at the end and numbered [n]. Companion documents: [bio-inspired-methods.md](bio-inspired-methods.md) (the methods survey this builds on) and [patent-landscape.md](patent-landscape.md).

How to read the citations. Where volume and pages are given, they were confirmed by a search result during this research. Where a reference says "details unverified", the paper was found but its bibliographic details or its specific numbers could not be confirmed. Statements marked **speculative** are this document's own hypotheses, not published findings.

Notation: T = thread tension, μ = mass per unit length, L = thread length, q = T/L = force density, σ, ε = engineering stress and strain, E = Young's modulus, ρ = density.

---

## 1. Terminology: two very different "funnel webs"

"Funnel web" names two unrelated groups of spiders that build structurally different things. The software lessons differ, so keep them apart.

| | (a) Funnel-web spiders | (b) Funnel weavers |
|---|---|---|
| Family | Atracidae (mygalomorphs), e.g. *Atrax robustus*, *Hadronyche* spp. | Agelenidae (araneomorphs), e.g. *Agelena*, *Agelenopsis*, *Eratigena*, *Tegenaria* |
| Where | Eastern and southern Australia | Worldwide; grass spiders, house spiders |
| Structure | Silk-lined burrow in soil, with silk trip lines radiating over the ground from the entrance | Horizontal, non-sticky silk sheet joined to a tubular retreat, usually with a scaffold of threads above |
| Prey capture | Sensing only: trip lines report prey, the spider runs out and bites | Sheet slows and entangles prey; vibration tells the spider where it is |
| Structural interest | Low. A sensor array, not a load-bearing structure | Moderate. A prestressed membrane-like net plus a tube |

### 1.1 Atracidae: burrow and trip lines [1][2][3]

What the sources describe (Australian Museum and related pages; no primary mechanics paper was found):

- The burrow is a silk-lined tube in a cool, humid, sheltered place, typically under logs and rocks. The Australian Museum describes the entrance of an *Atrax robustus* burrow as a funnel-shaped silk "vestibule" that leads into a collapsed, tunnel-like section with one or two slit-like openings, then a short chamber leading down into a burrow that is often only weakly silk-lined and rarely more than about 30 cm deep [1].
- Irregular silk trip lines radiate from the entrance over the ground. The spider sits just inside with its front legs on the lines, and rushes out when prey (beetles, small skinks) or a mate crosses them [1][2].
- Not all burrows have trip lines. Some sources report burrows with a Y- or T-shaped silk entrance instead [2]. For *Hadronyche flindersi*, Gray (1984) describes a side chamber and an entrance with a soil door in a collapsible silk collar [3] (cited second-hand via Wikipedia; details unverified).

Construction sequence: no published step-by-step description was found. The general picture (dig, line with silk, extend a silk collar at the entrance, lay trip lines outward) is inferred from the finished structure and is **unverified**.

Mechanics: no quantitative study of atracid trip-line mechanics was found. The physics is simple and well understood for any taut thread (Section 4.1): a trip line is a one-dimensional waveguide, the spider's legs are the receivers, and an array of radiating lines gives a crude direction estimate because the line that vibrates identifies the bearing. That last point is a reasonable inference, **not** a published measurement on Atracidae.

### 1.2 Agelenidae: sheet and tubular retreat [4][5][6][7][8]

- The capture structure is a flat, slightly concave, non-adhesive silk sheet joined to a funnel-shaped retreat. For *Agelena* the sheet is reported at about 40–60 cm wide, usually 0–1.5 m above ground and most often around 60 cm [4][5]. The web stays in one place and the spider repairs and enlarges it as it grows [4].
- *Agelena labyrinthica* takes its name from the labyrinth of tunnels beyond the funnel; the outer wall of the retreat may be camouflaged with grass and leaves [5]. Sources disagree on whether the retreat is central or at one edge; it varies between genera.
- The sheet is not sticky. It is full of entangling filaments that the spider keeps adding as it walks across it [6]. Prey is detected by vibration; the spider runs out from the retreat, bites, and drags the prey back [6].
- Construction sequence: a detailed, primary description of agelenid sheet construction was not found in this session. Eberhard's monograph (2020) covers building behaviour of non-orb weavers (Ch. 5) and defines sheet webs (Ch. 9) [9]; it is the place to look. In outline (from general descriptions, **unverified** as a sequence): anchor and scaffold lines are laid first, the tube is spun at the retreat site, and the sheet is thickened incrementally over many nights by the spider's repeated walking and spinning.

Mechanics, what is measured:

- Vibration transmission through sheet webs of the hobo spider *Eratigena agrestis* (Agelenidae) was measured by laser Doppler vibrometry and compared with tangle webs of the black widow *Latrodectus hesperus*. Transmission efficiency varied more *within* a web than *between* webs, and longitudinal and transverse vibrations were transmitted about equally well [7] (from the abstract; full text not read).
- Pessman and Hebets (2025) had *Agelenopsis pennsylvanica* build webs under quiet or loud vibratory noise. Web transmission properties depended on both the spider's prior (urban vs rural origin) and current noise exposure; rural webs kept more energy in longer-range vibrations under noise, while urban webs damped short-range vibrations. Whether this comes from structure, tension or anchoring is not yet known [8].

What is **not** established:

- Sheet tension values for agelenid webs: none found.
- That the funnel *concentrates* vibrations: no published measurement found. A plausible mechanism (**speculative**) is horn-like focusing: energy flux carried by a wave in a converging membrane is conserved along the convergence, so by the WKB argument amplitude grows as the width shrinks, roughly

$$
A(s) \propto \frac{1}{\sqrt{w(s)\,Z(s)}}, \qquad Z = \sqrt{N\,\rho_s h}
$$

where w is the local width of the sheet section feeding the tube, N the membrane tension per unit width and ρ_s h the areal mass. This would make the spider at the tube mouth an amplified receiver. It is a testable hypothesis, not a finding. Damping and the irregular sheet structure could easily dominate.

### 1.3 Lesson for structural design

Neither funnel web is a high-performance structure in the way an orb web is. The transferable ideas are (i) a sensing topology (radial trip lines, a sheet that funnels to a receiver) and (ii) incremental, maintenance-based construction (agelenids keep adding silk where they walk, so the sheet densifies along used paths). Idea (ii) is the same "reinforce where used" logic as bone remodeling and Physarum in bio-inspired-methods.md Sections 2.1 and 3.1. For load-bearing design, the orb web and the 3D webs in Sections 2 and 5 carry more weight.

---

## 2. Orb webs: the best-studied baseline

### 2.1 Thread roles and prestress [10][11][12][13]

An orb web has frame threads around the edge, radial threads running from the hub to the frame, and a sticky capture spiral laid last. Radials and frame are major ampullate (dragline) silk: stiff and strong. The capture spiral is flagelliform silk coated with aqueous glue: compliant and very extensible [10][11].

The web is a prestressed tension network. Wirth and Barth (1992) measured forces in *Araneus* orb webs [12]; later work cites the finding that forces decrease from mooring threads to frame threads to radii [13], and that abrupt changes in radius tension are rare in finished webs [12] (second-hand; details unverified). Spiders also adjust tension actively: *Octonoba sybotides* appears to raise radial tension when hungry, so that it can sense smaller prey, possibly by changing the distortion of the radials at the hub [14].

### 2.2 The web as a tension network: force-density form finding [15][16][17]

Every thread in a finished web is in tension (slack threads carry no load). Equilibrium of a pin-jointed tension network is the force density method (FDM) of bio-inspired-methods.md Section 5.4. With branch–node matrix C and force densities q_i = T_i / L_i:

$$
\mathbf{D} = \mathbf{C}^T \mathrm{diag}(\mathbf{q})\, \mathbf{C}, \qquad \mathbf{D}_{ff}\,\mathbf{x}_f = \mathbf{p}_f - \mathbf{D}_{fx}\,\mathbf{x}_x
$$

one sparse solve per coordinate. Masmeijer et al. used exactly this (via the COMPAS library) to form-find spider-web-inspired cable nets, assigning force density 1 at the hub rising linearly to 4 at the boundary to mimic the tension gradient in natural webs [16][17]. Their later work 3D-prints networks with programmed tension [17].

A tension network is not rigid in the Maxwell sense (bio-inspired-methods.md Section 4.2): a planar orb web has many mechanisms out of plane. It is stabilized by prestress. The tangent stiffness of a thread of tension T, axial stiffness EA and unit direction e is

$$
\mathbf{k}_t = \frac{EA}{L}\,\mathbf{e}\mathbf{e}^T + \frac{T}{L}\,(\mathbf{I} - \mathbf{e}\mathbf{e}^T)
$$

The second term, the geometric stiffness, is what gives the web out-of-plane stiffness. This is the same mechanism as a tensegrity or a cable net, and it is why spiders can tune web response with tension alone (Section 4).

### 2.3 Nonlinear silk localizes damage [18][19][20][21]

Cranford, Tarakanova, Pugno and Buehler (2012) combined web loading experiments with simulations using a silk law derived from atomistic models of dragline silk, in a web of radial (dragline) and spiral threads [18]. Their central result: the silk's nonlinear stress–strain response (stiff at first, then a softening "yield" plateau, then strong stiffening before failure), together with the web's discrete geometry, makes the web robust. Under a local load the loaded thread yields and stretches, so it fails locally without shedding much load to its neighbours, and the web as a whole remains usable. Control simulations with linear-elastic and with elastic–perfectly-plastic threads were less damage tolerant [18].

A generic four-regime law with that shape, suitable for a strut model (the functional form is written here for engineering use; the parameters in [18] come from atomistic fits and should be taken from the paper, not from this document):

$$
\sigma(\varepsilon) =
\begin{cases}
E_1\,\varepsilon, & \varepsilon \le \varepsilon_y \quad \text{(stiff, glassy)}\\
\sigma_y + E_2\,(\varepsilon-\varepsilon_y), & \varepsilon_y < \varepsilon \le \varepsilon_s \quad \text{(softened plateau, } E_2 \ll E_1)\\
\sigma_s\,\exp\!\big[\alpha(\varepsilon-\varepsilon_s)\big], & \varepsilon_s < \varepsilon \le \varepsilon_f \quad \text{(stiffening)}\\
0, & \varepsilon > \varepsilon_f \quad \text{(failure)}
\end{cases}
$$

Why it localizes damage, in one line: in the plateau the tangent stiffness drops, so the overloaded member stops attracting load (load redistribution in a network is proportional to relative stiffness), while the stiffening branch stops it from simply running away.

Earlier and related work:

- Alam and Jenkins (2005) were the first to apply damage-tolerance analysis to orb webs [19] (Int. J. Damage Mech. 14(4), 365–384; details as cited in a reference list).
- Aoyanagi and Okumura (2010) gave a simple analytical model of an orb web, with a formal solution when no threads are broken [20]. A result cited from this line of work: when the radials are sufficiently strong relative to the spirals, the web is free of stress concentrations even with a few spiral threads broken [20] (second-hand).
- Meyer, Pugno and Cranford (2014) showed that a *compliant* anchor line greatly raises the effective adhesion strength of a silk-on-silk joint compared with a rigid substrate, through a nonlinear change of the peeling angle [21]. Connection design matters, not only thread design.

### 2.4 Energy absorption: radials, spiral and air [22][23][24][25]

When an insect of mass m hits at speed v, the web must absorb

$$
\tfrac12 m v^2 = \sum_i V_i \int_0^{\varepsilon_i^{\max}} \sigma_i\, d\varepsilon \;+\; W_{\text{aero}} \;+\; W_{\text{visco}}
$$

with the elastic energy partly returned (rebound) unless hysteresis dissipates it. Which term dominates is debated:

- Sensenig, Lorentz, Kelly and Blackledge (2012), using high-speed video of prey impacts on 45 webs plus silk testing, found that radial threads dominated energy absorption in all webs, potentially accounting for roughly 100% of the work of stopping prey in larger webs. They judged aerodynamic dissipation small [22].
- Lin, Edmonds and Vollrath (1995) found in a full-web DYNA3D model, and then in experiments, that aerodynamic damping matters for prey capture [23].
- Zaera, Soler and Teus (2014) revisited the drag-per-unit-length assumption with standard fluid mechanics, bringing in a dependence on thread diameter. Their simulations favour air drag as a meaningful factor that limits web damage, and they argue that it favours smaller, less dense webs under wind and larger, denser webs under prey impact [24].

Drag per unit length on a thread of diameter d moving at speed v normal to its axis:

$$
f_d = \tfrac12 \rho_{\text{air}}\, C_D(Re)\, d\, v^2, \qquad Re = \frac{\rho_{\text{air}}\, v\, d}{\eta_{\text{air}}}
$$

Silk threads are micrometres thick, so Re is well below 1 and C_D is large and strongly Re-dependent (the low-Re regime is where [24] corrects earlier work). For engineered structures at millimetre-to-metre scale, aerodynamic damping is negligible except for very light cable nets and membranes; do not carry this mechanism over by default.

The capture spiral's job is less to absorb energy than to retain prey: it is extensible enough not to break and sticky enough to hold. Sensenig, Agnarsson and Blackledge (2010) found across many orb weavers that silk quality and web architecture co-evolve: larger spiders make better silk and webs with greater stopping potential, and, after correcting for size, species with higher-quality silk use less of it [25]. That is a material-versus-geometry trade-off, which is exactly what a topology optimizer exploits.

### 2.5 Material distribution depends on load type [26]

Qin, Compton, Lewis and Buehler (2015) 3D-printed elastomeric web mimics and modelled them [26]. Findings directly usable as a design rule:

- The loading type dictates the optimal material distribution. A *homogeneous* thread diameter is better for localized (point) loads; *stronger radials with weaker spirals* is better for distributed loads.
- There is an asymptotic prey (load-patch) size beyond which web strength saturates.
- Failure mode and optimal strength depend on the radial-to-spiral diameter ratio [26].

---

## 3. Silk as a material

### 3.1 Dragline versus capture silk: numbers [10][11][27][28][29]

| Property | Major ampullate (dragline) | Flagelliform (capture spiral, *A. diadematus*) | Kevlar 49 (for scale) |
|---|---|---|---|
| Initial modulus | ~10 GPa (unverified) | ~0.003 GPa | ~130 GPa |
| Tensile strength | ~1.1 GPa (Gosline table); 1.1–1.8 GPa across reviews | ~0.5 GPa | ~3.6 GPa |
| Strain at break | ~0.20–0.27 | ~2.7 (270%) | ~0.03 |
| Toughness (work to fracture per volume) | ~160 MJ/m³ (Gosline table); 100–400 MJ/m³ across reviews | ~150 MJ/m³ | ~50 MJ/m³ |

Sources and caveats. The flagelliform row and the Kevlar row appear in a reproduced property table attributed to Gosline et al. (1999) [10] as surfaced in a search result [28]; the dragline ranges come from a review table [28][29]. The exact Gosline values for dragline (10 GPa, 1.1 GPa, 0.27, 160 MJ/m³) are the commonly quoted ones but were **not** re-verified against the original table in this session. Silk properties vary widely with species, spinning speed, humidity and test method, so treat all of these as order-of-magnitude design inputs.

Two points matter more than the exact numbers:

1. Dragline and capture silk are about equally *tough*, but get there differently: dragline by strength, capture silk by extensibility. Their stiffnesses differ by more than three orders of magnitude [11] (Denny 1976). Denny argued that viscid silk acts as the shock absorber and frame silk as the structural element [11].
2. Within one spider, the five fibrous silks are tuned for different jobs. In *Argiope argentata*, the glue-coated capture spiral is about an order of magnitude more extensible than the other silks, major ampullate is about 50% stronger but less extensible, and aciniform (wrapping) silk is the toughest, nearly twice the others [27].

Extremes: dragline from Darwin's bark spider, *Caerostris darwini*, which spans rivers, averaged about 350 MJ/m³ with some samples reaching 520 MJ/m³, more than twice any previously described silk [30].

### 3.2 Hierarchical structure: why silk is tough [31][32]

Dragline silk is a semicrystalline polymer: β-sheet nanocrystals (a few nanometres) cross-link a matrix of less ordered chains. Keten, Xu, Ihle and Buehler (2010) showed by simulation that β-sheet crystals confined to a few nanometres are stiffer, stronger and tougher than larger ones, because small crystals let the weak hydrogen bonds act cooperatively and fail by stick–slip rather than by brittle cracking [31]. Nova, Keten, Pugno, Redaelli and Buehler (2010) carried this to the fibril level: strength and toughness depend on nanocrystal size, with small crystals giving high performance [32].

The design reading, already in bio-inspired-methods.md Section 6.4: many small, weak, redundant links that fail in sequence beat a few strong ones. At the structural scale this is "many thin struts with a softening law" (Section 2.3).

### 3.3 Supercontraction [33][34]

Unrestrained dragline silk shrinks to roughly half its length when wetted or exposed to very high humidity; restrained, it develops substantial stress instead [33]. Argiope silk supercontracts more than Nephila silk, and major ampullate more than minor ampullate [33]. Mechanistically it is attributed to water plasticizing the amorphous matrix into a rubbery, entropy-elastic state [33]. Mortimer et al. (2016) suggest supercontraction may have evolved partly as a way to control web tension, and hence vibration transmission [34]. Work (1977) is often cited as the classic early study; it was not found in this session (details unverified).

Engineering analogue: humidity- or temperature-activated prestress. In a printed tension network, a shape-memory or hygroscopic material that contracts on activation could set or reset prestress after assembly (**speculative** for our software; it is a fabrication concern, not a design-solver feature).

### 3.4 Glue droplets as distributed dampers and stored length [35][36][37]

The capture spiral is a core filament coated with aqueous glue that breaks up into regularly spaced droplets (a Plateau–Rayleigh instability of the liquid coat) [35]. Vollrath and Edmonds (1989) proposed that the droplets act as microscopic *windlasses*: when the thread is slackened, the core filament coils up inside the droplets, keeping the thread taut [35][36] (original paper not found; described second-hand).

Elettro, Neukirch, Vollrath and Antkowiak (2016) confirmed the mechanism: the capture thread behaves as a solid in tension and as a liquid in compression, because the core filament buckles and spools inside the droplets along a force plateau set by the competition of elasticity and capillarity [37]. They built a synthetic version from oil droplets on a polymer filament [37].

An order-of-magnitude statement of the plateau (a scaling argument written for this document, not an equation quoted from [37]): spooling starts when the capillary pull on the filament exceeds what is needed to bend it inside the drop, and in the plateau the thread tension is held near

$$
T_{\text{plateau}} \sim 2\pi a\,\gamma
$$

for filament radius a and glue surface tension γ, independent of how much length is stored. Mechanically this is a tension-only element with a reservoir of slack and a near-constant small tension: it never goes slack and never buckles.

---

## 4. Vibration sensing: the web as an antenna

### 4.1 Wave mechanics of a thread and a sheet

A taut thread carries two wave types that matter to spiders:

$$
c_T = \sqrt{\frac{T}{\mu}} = \sqrt{\frac{\sigma}{\rho}} \quad \text{(transverse)}, \qquad c_L = \sqrt{\frac{E}{\rho}} \quad \text{(longitudinal)}
$$

Transverse speed is set by prestress, longitudinal speed by stiffness. Since σ ≪ E for a working thread, longitudinal waves are much faster. For a membrane-like sheet with tension N per unit width and areal mass m_s, transverse waves obey m_s ü = N∇²u, with speed √(N/m_s). Attenuation along a thread can be summarized as A(x) = A_0 e^{−αx}, with α depending on silk viscoelasticity, air damping and junction losses.

### 4.2 What Mortimer, Vollrath and colleagues showed [34][38][39][40][41]

- Mortimer, Soler, Siviour, Zaera and Vollrath (2016), with a finite-element model of *Araneus diadematus* webs, showed that both transverse and longitudinal wave amplitudes can be adjusted through web tension and dragline stiffness, both of which the spider controls [34].
- A review-style paper (Mortimer et al. 2019, *Integr. Comp. Biol.* 59, 1636–1645) summarizes energy loss, distortion and filtering of vibrations in orb webs, and notes that spider leg vibration sensors (lyriform slit organs) are least sensitive below about 30 Hz and most sensitive near 1 kHz [38].
- Remote versus direct monitoring: *Zygiella x-notata* monitors its web from a retreat via a signal thread. The thread adds no biologically relevant delay but attenuates the signal [39]. This is the closest orb-web analogue of a funnel weaver's retreat.
- Locating a source: FE models show that comparing longitudinal and transverse amplitudes across the eight legs can give the direction and distance of a vibration source. A signal thread hides that location information unless the input angle is known [40].
- In black widow (theridiid) gumfoot webs, high-speed video of prey capture shows dominant frequencies of 40–60 Hz [41] (conference abstract; details unverified).

### 4.3 Engineering translations of web sensing [16][42][43]

- Masmeijer, Méndez Echenagucia, Slavič, Loendersloot and Habtour study spider-web-inspired cable nets as structural-health-monitoring sensors, including the effect of hub eccentricity on sensing [16].
- Kang et al. (2014, *Nature* 516, 222–226) built an ultrasensitive crack-based strain and vibration sensor inspired by the spider's slit organ (not the web): nanoscale crack junctions in a platinum film on a polymer substrate that disconnect and reconnect under strain [42].
- Many flexible "spider-web" strain sensors exist, for example a laser-induced graphene sensor with a web layout reporting a gauge factor of 36.8 [43]. These borrow the web's appearance more than its mechanics.

For our software the sensing story is secondary, but one idea transfers: a structure's dynamic response (which modes are excited, how much signal reaches a point) can be shaped with prestress and stiffness distribution, which are design variables we already have (Section 8, P5).

---

## 5. Sheet webs and 3D webs: space frames and tension networks

### 5.1 Sheet webs (Linyphiidae) [44][45]

Linyphiid sheet webs are near-horizontal sheets, often with a tangle of "knockdown" threads above that intercepts flying insects and drops them onto the sheet, where the spider waits underneath [44]. Hormiga and Eberhard's survey of about 113 species found that nearly all build a more or less horizontal, continuous sheet with an open space just below it, and that the function of many architectural features remains obscure [44]. Benjamin and colleagues report that these sheets use no aqueous capture glue; aggregate-gland secretions instead cement the layers and threads together [44] (second-hand; details unverified). Structurally this is a two-layer system: a compliant, sacrificial interception layer over a load-carrying, sensing layer.

### 5.2 Cobwebs and gumfoot lines (Theridiidae) [41][46][47]

Benjamin and Zschokke (2003) describe construction in several theridiids and classify their webs into four types: gumfooted webs with a central retreat (*Achaearanea*-type) or a peripheral retreat (*Latrodectus*-type), and webs without gumfoot lines that are either viscid (*Theridion*-type) or sheet-like (*Coleosoma*-type) [46]. They suggest the cost of viscid silk drove web reduction in the family [46].

The gumfoot line is a preloaded, spring-loaded trap. A vertical thread is attached to the ground under tension by a weakly adhering attachment disc with glue near the bottom. When a walking insect touches the glue, the weak disc releases and the stored tension hoists the insect off the ground [41][47]. Cobweb spiders make two kinds of attachment disc: weak ones for gumfoot lines (designed to release) and strong ones for scaffolding (the thread breaks before the disc lets go) [47]. In engineering terms: a deliberately weak link that converts stored strain energy into motion, next to deliberately strong links in the same network. Kinematic measurements suggest the gumfoot threads operate past yield, about 5% shorter after the trap is sprung [41] (conference abstract; details unverified).

### 5.3 3D web architecture: scans and mechanics [48][49]

Su, Qin, Saraceno, Krell, Muehlethaler, Bisshop and Buehler (2018) built a sheet-laser scanning platform that images a web slice by slice and reconstructs a complete 3D web of the tent-web spider *Cyrtophora citricola* at micron resolution [48]. Su and Buehler (2020) then simulated that architecture with coarse-grained bead-spring models under stretching and projectile impact, and concluded that the interplay of silk nonlinearity and the redundancy of the complex 3D structure is crucial for robustness and resilience [49].

The 3D web is an irregular, highly redundant tension network with no compression members. In bio-inspired-methods.md terms it is closer to a prestressed cable net or a tensegrity without struts (supports and the spider's anchor points take the compression) than to a space frame. The lesson for lattices: redundancy plus a softening member law gives graceful failure even in an irregular network.

### 5.4 Web construction as an algorithm [50][51][52]

- Krink and Vollrath (1997) modelled orb-web construction by an "ideal" *Araneus diadematus* with a rule-based simulation and optimized its rules with a genetic algorithm [50] (*J. Theor. Biol.* 185, 321–331, as listed on a user-uploaded copy; details unverified).
- Krink and Vollrath (1998) used a virtual spider robot to study emergent properties of web-building behaviour, including how leg loss changes capture-spiral geometry [51] (*Proc. R. Soc. B* 265, 2051–2055; title as found, exact wording unverified).
- Vollrath and Krink (2020) review this work, combining the rule-controlled virtual spider ("Theseus", building a 2D network of rigid segments) with an evolutionary algorithm using cost–benefit evaluation, as inspiration for soft robotics [52].

These are behavioural models: local rules (follow a thread, measure the gap to the previous spiral turn, attach) produce a global structure. They are elegant but they do not optimize mechanics directly. For a design tool, the same structures can be reached faster by FDM plus parameter optimization (Section 8).

---

## 6. Engineering translations already published

| Work | What it does | Relevance |
|---|---|---|
| Qin et al. 2015 [26] | 3D-printed elastomeric orb webs; load type sets optimal radial/spiral diameter ratio | Direct design rule for polar lattices |
| Miniaci, Krushynska, Movchan, Bosia, Pugno 2016 [53] | *Nephila*-web-based lattice metamaterial with sub-wavelength band gaps for elastic waves; parametric optimization | Vibration isolation lattices |
| Krushynska, Bosia, Miniaci, Pugno 2017 [53] | Spider-web-structured labyrinthine acoustic metamaterial for low-frequency sound (*New J. Phys.*; details unverified) | Acoustic panels |
| Masmeijer et al. 2023, 2025 [16][17] | FDM form-found web-inspired cable nets for sensing; 3D-printed networks with programmed tension | Closest to our FDM module |
| "Energy absorber inspired by spider webs" 2025 [54] | 3D orb-web lattice cell, size-optimized with beam elements for energy absorption, printed and tested; reported gains of 108.9%, 90.4% and 27.5% for three variants | Crash/impact lattices (authors not captured; details unverified) |
| "Dynamic behaviour of spiderweb-inspired lattice structures" 2026 [55] | Bézier radial/spiral lattice, ABAQUS/Explicit impact simulations | Impact; simulation-only claims |
| "Bionic spiderweb lattice metamaterials" 2025 [56] | Energy absorption plus vibration isolation reported from about 25 Hz to 2000 Hz | Vibration isolation (details unverified) |
| "Investigation on the energy-absorbing properties of bionic spider web structure" 2023 [57] | Radial line length and spiral cross-section affect energy absorption more than radial diameter, spiral count and spacing | Parameter sensitivity |
| Elleithy and Acar 2026 [58] | "Stochastic spider webs": bio-inspired metamaterial design with structural optimization, validated by compression and three-point bending; web architecture favoured in bending | Closest to TO; found once, not confirmed on a second search (details unverified) |
| Elettro et al. 2016 [37] | Synthetic "liquid wire" from oil droplets on a filament | Stored-length tension elements |
| Kang et al. 2014 [42] | Slit-organ-inspired crack sensor | Sensors, not structures |

Gap: no paper was found that runs density-based topology optimization (SIMP) with an explicit spider-web principle, such as a softening member law for damage tolerance or a prestressed tension-network ground structure. The bio-inspired stiffener-growth method of Ding and Yamazaki (2004) draws on tree and leaf branching, not webs [59]. Fail-safe topology optimization (optimizing against the removal of local patches of material) exists independently of spiders, e.g. Jansen, Lombaert, Schevenels and Sigmund (2014), "Topology optimization of fail-safe structures using a simplified local damage model" (details unverified), and is the natural engineering home for the Cranford result.

---

## 7. Patent check (via WebSearch only; everything here is unverified)

Searches were limited to the WebSearch tool; patent databases could not be opened directly. Nothing below has been checked for claims, status or expiry.

| Document | What the search snippet says | Status known? |
|---|---|---|
| CN101950316A (Beihang University; inventors listed as 袁松梅 and 吕茂斌) | "Design method of bionic spider web composite material structure": ring-shaped and arc-shaped composite web-like structures built from core rings, with FE stress analysis used to place radial ("longitude") and circumferential ("latitude") elements. Filed 2010-08-19, published 2011-01-19 [60] | No |
| CN115561484B | MEMS bionic acceleration sensor based on a spider-web structure with cantilever beams [61] | Listed as granted (B), unverified |
| Unnumbered CN patent | Bionic cobweb topology for laying out fibre-grating sensors on a tunnel boring machine cutterhead [61] | No |
| Unnumbered CN patent | MXene/fibroin "nano-cobweb" sensors in clothing [61] | No |
| US patents | No granted US patent explicitly claiming a spider-web-inspired lattice, metamaterial or structural design method was found in two searches [62] | n/a |

Assessment (unverified): the concern for our software is low. Radial-plus-circumferential ("web") lattices are widely published in the academic literature listed in Section 6, which is prior art against broad claims. The one document that touches a *design method* is CN101950316A, which is a Chinese application (published 2011) for web-like composite layouts placed by FE stress analysis; its scope, grant status and term are unknown, and it would only matter for products sold in China. A proper search on Google Patents, Espacenet and CNIPA for "spider web" / "orb web" / 蜘蛛网 / 仿生 combined with lattice, cable net and topology optimization should be done before shipping any web-named feature.

---

## 8. Proposals for the software

Context: today the repository has the voxel SIMP solver (src/reversed_core: simp.py, fe.py, filters.py, multigrid.py, optimizers.py). The lattice layer, force-density form finding and nonlinear materials are planned (bio-inspired-methods.md Section 10, steps 3–6) but not yet built. Effort ratings below assume those modules exist when the proposal is built, and say so when they do not.

Value and effort are rated Low / Medium / High.

### P1. Damage-tolerant design via a softening member law (Cranford principle) — top proposal

What: give the user a "damage tolerant" option. In the lattice layer, struts use the four-regime law of Section 2.3 (stiff, soft plateau, stiffening). In the voxel solver, the equivalent is fail-safe TO: optimize the worst case over local damage patches.

Math. For lattices: nonlinear static or explicit analysis with the piecewise σ(ε) of Section 2.3, objective = energy absorbed to first loss of function, or residual stiffness after removing the most-loaded member. For SIMP: the fail-safe formulation

$$
\min_{\mathbf{x}} \; \max_{k=1..K} \; c_k(\mathbf{x}) = \mathbf{U}_k^T \mathbf{K}_k(\mathbf{x}) \mathbf{U}_k, \qquad E_e^{(k)} = \begin{cases} \varepsilon E_e(x_e) & e \in \Omega_k \\ E_e(x_e) & \text{otherwise} \end{cases}
$$

aggregated with a p-norm or KS function over K damage patches Ω_k, with the usual volume constraint and adjoint sensitivities per patch.

Fits: SIMP solver (fail-safe variant, available now) and nonlinear materials plus lattice layer (softening struts, later).

Value: High. Damage tolerance is a real customer need (aerospace fail-safe rules, drones, protective gear), Fusion does not offer it, and the evidence is solid (Cranford 2012 [18]; Su and Buehler 2020 [49]; fail-safe TO literature).
Effort: Medium for fail-safe SIMP (K extra solves per iteration; patches can be sampled stochastically to keep cost down). High for the full nonlinear lattice version, which needs the nonlinear FE path.

### P2. Prestressed tension networks with web-like force-density gradients — top proposal

What: a cable-net / tension-network designer. The user draws a frame (or picks boundary nodes on a part) and a hub or load region; the tool generates radial-plus-spiral (or arbitrary) topology, form-finds it by FDM with a force-density field, and reports thread tensions, out-of-plane stiffness and natural frequencies. Use cases: cable nets, tensile façades, trampoline-like membranes, sensor nets, tensioned lattices inside printed parts, drone arresting nets.

Math. FDM (Section 2.2) with q as design variables, for example q_i = q_hub + (q_rim − q_hub)·r_i/R as in [16]. Then optimize q (and member areas A_i) for an objective such as maximum out-of-plane stiffness under a point load at a given prestress budget, using the tangent stiffness of Section 2.2:

$$
\mathbf{K}_t(\mathbf{q}, \mathbf{A}) = \sum_i \left[\frac{E A_i}{L_i}\,\mathbf{e}_i\mathbf{e}_i^T + q_i\,(\mathbf{I}-\mathbf{e}_i\mathbf{e}_i^T)\right]
$$

with constraint T_i = q_i L_i ≤ σ_allow A_i and q_i ≥ q_min > 0 (all threads stay taut). Sensitivities of K_t with respect to q are explicit. Topology can be pruned in the same loop (force-density topology optimization exists in the cable literature [15]).

Fits: force-density form finding (planned, bio-inspired-methods.md Section 5.4); output can be meshed into the voxel/SDF pipeline as struts.

Value: Medium–High. It opens a market segment (tensile and cable structures) that the voxel solver cannot serve, and the math is linear and fast. Evidence for the web-specific gradients is moderate (measured tension decrease toward the hub [12][13]; Masmeijer et al. [16][17]); the FDM itself is mature.
Effort: Low–Medium once the FDM module exists (one sparse solve per form, small gradient-based loop).

### P3. Load-adaptive polar ("web") lattice infill — top proposal

What: a lattice family for plates, discs, hubs and pads with a radial-plus-circumferential topology, centred on a load point, whose radial/circumferential stiffness ratio is set from the load type, following Qin et al. [26]: uniform member sizes for concentrated loads, stronger radials for distributed loads.

Math. Node layout in polar coordinates r_j = R (j/J)^β (β grades spacing), θ_k = 2πk/K, struts radial and circumferential (optionally a logarithmic spiral r = a e^{bθ} instead of rings). Areas A_rad, A_circ chosen by a two-parameter optimization or by the rule above; check Maxwell's count (bio-inspired-methods.md Section 4.2): a pure radial-ring net is a quad mesh, bending-dominated in plane (each interior node has Z = 4), so add a diagonal or spiral family where in-plane stiffness matters. Homogenized properties are anisotropic and position-dependent; within the de-homogenization step (bio-inspired-methods.md Section 4.5) the orientation field can simply be the polar frame around the load.

Fits: lattice layer (as a unit-cell family and a conformal generator).

Value: Medium. It is a recognizable, explainable pattern that users will like and that has a measured design rule behind it, but for pure stiffness per weight it will not beat a stress-aligned octet or plate lattice. Its real advantage is graceful failure under point impact, which ties to P1.
Effort: Low (geometry generator plus two-parameter sizing).

### P4. Stored-length tension elements and sacrificial weak links (windlass and gumfoot ideas)

What: two special elements for the lattice/cable layer. (i) A "windlass" strut: tension-only with a slack reservoir and a constant low tension, from capture-thread droplets [37]. (ii) A designed weak link that releases at a set force (gumfoot attachment disc [47]) or a frangible link that fails in sequence (bio-inspired-methods.md Section 6.4).

Math. (i) T(δ) = T_plateau for δ < δ_stored, then the normal elastic law; never in compression. (ii) A cohesive-zone or force-limited link, F ≤ F_release, after which the member is removed (or its stiffness set to ε) and the network re-equilibrated.

Fits: nonlinear materials, lattice layer.

Value: Low–Medium. Useful for deployable structures, energy absorbers with a target force–displacement curve, and fuses. **Speculative** as a mainstream feature; the windlass in particular depends on a fabrication route (liquid-in-drop or a printed coiled spring) that we do not control.
Effort: Medium (element definitions are easy; the event-driven nonlinear solve is the cost).

### P5. Vibration-tuned structures: tension and stiffness as dynamic design variables

What: let the user shape dynamic response: raise or lower a target natural frequency, maximize transmission from a sensor site to a receiver (sensing), or minimize it (isolation), using prestress and member stiffness as the variables, as spiders do [34][38].

Math. Eigenproblem with prestress: (K_t(q, A) − ω² M) φ = 0, sensitivities dω²/dq_i = φ^T (∂K_t/∂q_i) φ / (φ^T M φ); or harmonic response u(ω) = (K_t − ω²M + iωC)^{-1} f with an objective on |u_receiver|. In SIMP, the corresponding objective is standard eigenfrequency or dynamic-compliance TO.

Fits: SIMP solver (eigenfrequency TO, independent of spiders) and the FDM/tension-network module of P2.

Value: Medium for the eigenfrequency part (which is a standard TO feature anyway); Low for sensing-transmission objectives, which serve a niche. The web literature gives the motivation, not new math.
Effort: Medium.

### P6. Rule-based construction as a generator (Krink–Vollrath style)

What: generate web-like cable nets inside arbitrary frames with local construction rules (lay radials to the frame, then a spiral with spacing tied to local radial spacing), then hand off to P2 for form finding.

Value: Low. **Speculative** as a product feature: the rules reproduce spider webs, not optimal structures, and FDM plus parameter optimization reaches better designs directly. It may be worth a weekend as a "sketch" generator with good visual appeal.
Effort: Low–Medium.

### P7. Funnel concentrator for sensing

What: shape a membrane or net so that it funnels vibration energy to one sensor (the agelenid hypothesis of Section 1.2).

Value: Low and **speculative**: there is no published measurement that the agelenid funnel concentrates vibration. Revisit only if a customer asks for passive vibration-collecting structures.
Effort: Medium (needs dynamic membrane analysis).

---

## 9. Synthesis

| Spider principle | Evidence strength | Engineering form | Our module | Proposal | Value / effort |
|---|---|---|---|---|---|
| Nonlinear silk (stiff, yield, stiffen) localizes damage | Strong: experiments plus simulation [18][49] | Softening member law; fail-safe TO | SIMP (fail-safe), nonlinear lattice | P1 | High / Medium |
| Prestressed tension network, tension gradient hub to frame | Moderate: measurements [12][13], engineering demos [16][17] | FDM form finding with optimized force densities | Force-density module | P2 | Medium–High / Low–Medium |
| Load type sets radial/spiral material ratio | Moderate: printed and simulated webs [26] | Load-adaptive polar lattice | Lattice layer | P3 | Medium / Low |
| Radials absorb most impact energy | Moderate, debated vs air drag [22][23][24] | Size the load-path members, not the infill, for impact | Lattice, nonlinear | Inputs to P1, P3 | — |
| Silk quality and geometry co-evolve | Moderate, comparative [25] | Material–geometry trade-off is what TO already does | SIMP | Confirms approach | — |
| Glue droplet windlass: solid in tension, liquid in compression | Strong for the mechanism [37] | Stored-length tension element | Nonlinear materials | P4 | Low–Medium / Medium |
| Gumfoot: weak release link next to strong links | Moderate [46][47] | Force-limited fuse links | Nonlinear materials | P4 | Low–Medium / Medium |
| Web tension and stiffness tune vibration | Strong for orb webs [34][38][39][40] | Dynamic TO with prestress variables | SIMP, FDM | P5 | Medium (eigen) / Medium |
| 3D redundancy plus nonlinearity gives resilience | Moderate, simulation of scanned web [48][49] | Redundant irregular tension networks | FDM, lattice | Part of P1/P2 | — |
| Funnel-web trip lines as a radial sensor array | Descriptive only [1][2] | Radial sensor layout | — | none | Low |
| Funnel concentrates vibration | None found | Horn-like membrane | — | P7 | Low, speculative |
| Construction rules as algorithm | Behavioural models [50][51][52] | Rule-based generator | FDM front end | P6 | Low, speculative |
| β-sheet nanoconfinement | Strong at molecular scale [31][32] | Many small redundant links (already in bio-inspired-methods.md 6.4) | Lattice | — | — |

Recommended order: P1 (fail-safe SIMP first, because the solver exists today), then P3 once the lattice layer exists, then P2 with the FDM module. P4–P7 wait for customer pull.

---

## References

1. Australian Museum. Funnel-web spiders (group) and Sydney funnel-web spider pages. https://australian.museum/learn/animals/spiders/funnel-web-spiders-group/ ; https://australian.museum/sydney-funnel-web-spider
2. Wikipedia. Sydney funnel-web spider; *Atrax*. https://en.wikipedia.org/wiki/Sydney_funnel-web_spider ; https://en.wikipedia.org/wiki/Atrax
3. Gray, M.R. (1984). Notes on the *Atrax adelaidensis* species group, including burrowing behaviour. *Records of the South Australian Museum* 18, 441–452 (as cited on Wikipedia; details unverified). https://en.wikipedia.org/wiki/Hadronyche_flindersi
4. Wikipedia. *Agelena*; Agelenidae. https://en.wikipedia.org/wiki/Agelena ; https://en.wikipedia.org/wiki/Agelenidae
5. Wikipedia. *Agelena labyrinthica*. https://en.wikipedia.org/wiki/Agelena_labyrinthica
6. Arizona-Sonora Desert Museum. Funnel-web spider (grass spider) natural history; Encyclopedia of Life, funnel weavers. https://desertmuseum.org/books/nhsd_funnelweb_spider.php ; https://eol.org/pages/8826/articles
7. (2016). Vibration transmission through sheet webs of hobo spiders (*Eratigena agrestis*) and tangle webs of western black widow spiders (*Latrodectus hesperus*). *J. Comp. Physiol. A* (authors not confirmed, possibly the Vibert/Gries group; details unverified). https://link.springer.com/article/10.1007/s00359-016-1113-0
8. Pessman, B.J., Hebets, E.A. (2025). Web transmission properties vary with a spider's past and current noise exposure. *Current Biology* 35(8), 1706–1715. https://www.sciencedirect.com/science/article/abs/pii/S0960982225002052 ; https://news.unl.edu/article/novel-spider-research-featured-in-new-york-times
9. Eberhard, W.G. (2020). *Spider Webs: Behavior, Function, and Evolution.* University of Chicago Press. https://www.press.uchicago.edu/ucp/books/book/isbn/9780226534602.html
10. Gosline, J.M., Guerette, P.A., Ortlepp, C.S., Savage, K.N. (1999). The mechanical design of spider silks: from fibroin sequence to mechanical function. *J. Exp. Biol.* 202, 3295–3303. https://citedrive.com/en/discovery/the-mechanical-design-of-spider-silks-from-fibroin-sequence-to-mechanical-function
11. Denny, M. (1976). The physical properties of spider's silk and their role in the design of orb-webs. *J. Exp. Biol.* 65, 483–506. https://cob.silverchair.com/jeb/article/65/2/483/22346/THE-PHYSICAL-PROPERTIES-OF-SPIDER-S-SILK-AND-THEIR
12. Wirth, E., Barth, F.G. (1992). Forces in the spider orb web. *J. Comp. Physiol. A* 171(3), 359–371 (as cited; details unverified). https://www.ibb.uni-stuttgart.de/dokumente/abschlussarbeiten/poster/2021/2021-06_Habrich_Poster.pdf ; https://www.americanarachnology.org/journal-joa/joa-all-volumes/detail/article/download/arac-53-01-62.pdf/?no_cache=1
13. Statement on tension decreasing from mooring threads to frame and radii in *Araneus diadematus*, as cited in [17]. https://arxiv.org/pdf/2509.05855
14. Watanabe, T. (2000). Web tuning of an orb-web spider, *Octonoba sybotides*, regulates prey-catching behaviour. *Proc. R. Soc. Lond. B* 267, 565–569. https://pmc.ncbi.nlm.nih.gov/articles/PMC1690565 ; https://www.sciencenews.org/article/hungry-spiders-tune-web-jiggliness
15. Topology optimization of force densities for form finding of cable structures (2026 preprint; authors not captured, details unverified). https://arxiv.org/html/2608.24177v1
16. Masmeijer, T.C.P., Méndez Echenagucia, T., Slavič, J., Loendersloot, R., Habtour, E. (2023). Effect of eccentricity on sensing in spider web inspired cable nets. *J. Phys.: Conf. Ser.* 2647(19), 192012. https://research.utwente.nl/en/publications/effect-of-eccentricity-on-sensing-in-spider-web-inspired-cable-ne/ ; https://www.researchgate.net/publication/374004391_Effect_of_eccentricity_on_sensing_in_spider_web_inspired_cable_nets
17. Masmeijer, T., Swain, C., Hill, J., Habtour, E. (2025). Programming tension in 3D printed networks inspired by spiderwebs. *Materials & Design*, doi 10.1016/j.matdes.2025.115281 (volume/pages unverified). https://arxiv.org/pdf/2509.05855 ; https://digital.lib.washington.edu/researchworks/items/b7f7ab5e-8a12-48b3-858c-d70fb58eb49f
18. Cranford, S.W., Tarakanova, A., Pugno, N.M., Buehler, M.J. (2012). Nonlinear material behaviour of spider silk yields robust webs. *Nature* 482, 72–76. https://ideas.repec.org/a/nat/nature/v482y2012i7383d10.1038_nature10739.html ; https://pugno.dicam.unitn.it/NP_PDF/175-NATURE12-spidersilk-web.pdf
19. Alam, M.S., Jenkins, C.H. (2005). Damage tolerance in naturally compliant structures. *Int. J. Damage Mech.* 14(4), 365–384 (as cited in a reference list). https://bhxb.buaa.edu.cn/bhzk/en/article/id/13776
20. Aoyanagi, Y., Okumura, K. (2010). Simple model for the mechanics of spider webs. *Phys. Rev. Lett.* 104, 038102. https://link.aps.org/doi/10.1103/PhysRevLett.104.038102
21. Meyer, A., Pugno, N.M., Cranford, S.W. (2014). Compliant threads maximize spider silk connection strength and toughness. *J. R. Soc. Interface* (volume/pages unverified). https://pmc.ncbi.nlm.nih.gov/articles/PMC4233707
22. Sensenig, A.T., Lorentz, K.A., Kelly, S.P., Blackledge, T.A. (2012). Spider orb webs rely on radial threads to absorb prey kinetic energy. *J. R. Soc. Interface* 9(73), 1880–1891. https://pmc.ncbi.nlm.nih.gov/articles/PMC3385755
23. Lin, L.H., Edmonds, D.T., Vollrath, F. (1995). Structural engineering of an orb-spider's web. *Nature* 373, 146–148. https://www.nature.com/articles/373146a0
24. Zaera, R., Soler, A., Teus, J. (2014). Uncovering changes in spider orb-web topology owing to aerodynamic effects. *J. R. Soc. Interface* 11(98), 20140484. https://e-archivo.uc3m.es/entities/publication/2225e060-0fa9-4c17-b7aa-fb26deb37ab0/full
25. Sensenig, A., Agnarsson, I., Blackledge, T.A. (2010). Behavioural and biomaterial coevolution in spider orb webs. *J. Evol. Biol.* 23(9), 1839–1856. https://datadryad.org/dataset/doi:10.5061/dryad.1827
26. Qin, Z., Compton, B.G., Lewis, J.A., Buehler, M.J. (2015). Structural optimization of 3D-printed synthetic spider webs for high strength. *Nat. Commun.* 6, 7038. https://pmc.ncbi.nlm.nih.gov/articles/PMC4479035 ; https://news.mit.edu/2015/3d-printed-webs-0515
27. Blackledge, T.A., Hayashi, C.Y. (2006). Silken toolkits: biomechanics of silk fibers spun by the orb web spider *Argiope argentata* (Fabricius 1775). *J. Exp. Biol.* 209(13), 2452–2461. https://ideaexchange.uakron.edu/biology_ideas/156
28. Silk property tables reproduced in reviews (dragline 1.1–1.8 GPa, 100–400 MJ/m³; flagelliform *A. diadematus* row; Kevlar 49 row). https://pmc.ncbi.nlm.nih.gov/articles/PMC7466511 ; https://www.americanarachnology.org/journal-joa/joa-all-articles/article/download/arac-033-02-0629.pdf/?no_cache=1
29. Wikipedia. Spider silk. https://en.wikipedia.org/wiki/Spider_silk
30. Agnarsson, I., Kuntner, M., Blackledge, T.A. (2010). Bioprospecting finds the toughest biological material: extraordinary silk from a giant riverine orb spider. *PLoS ONE* 5(9), e11234. https://www.plosone.org/article/info:doi/10.1371/journal.pone.0011234
31. Keten, S., Xu, Z., Ihle, B., Buehler, M.J. (2010). Nanoconfinement controls stiffness, strength and mechanical toughness of β-sheet crystals in silk. *Nat. Mater.* 9, 359–367. https://web.mit.edu/mbuehler/www/papers/Nature_Materials_2010.pdf
32. Nova, A., Keten, S., Pugno, N.M., Redaelli, A., Buehler, M.J. (2010). Molecular and nanostructural mechanisms of deformation, strength and toughness of spider silk fibrils. *Nano Lett.* 10, 2626–2634. https://pugno.dicam.unitn.it/NP_PDF/146-NL10-spidersilk.pdf
33. Supercontraction sources surfaced in search (various authors; Pugno group, Blackledge group, PRL 2022 model). https://iris.unitn.it/handle/11572/304122 ; https://works.bepress.com/todd_blackledge/2/ ; https://link.aps.org/accepted/10.1103/PhysRevLett.128.098101
34. Mortimer, B., Soler, A., Siviour, C.R., Zaera, R., Vollrath, F. (2016). Tuning the instrument: sonic properties in the spider's web. *J. R. Soc. Interface* 13(122), 20160341. https://pmc.ncbi.nlm.nih.gov/articles/PMC5046944
35. Elettro, H., et al. (2015 preprint). In-drop capillary spooling of spider capture thread inspires highly extensible fibres (describes Vollrath and Edmonds 1989 windlass hypothesis). https://ar5iv.arxiv.org/html/1501.00962
36. Vollrath, F., Edmonds, D.T. (1989). Modulation of the mechanical properties of spider silk by coating with water. *Nature* (details unverified; described second-hand in [35]).
37. Elettro, H., Neukirch, S., Vollrath, F., Antkowiak, A. (2016). In-drop capillary spooling of spider capture thread inspires hybrid fibers with mixed solid–liquid mechanical properties. *PNAS* 113(22), 6143–6147. https://pmc.ncbi.nlm.nih.gov/articles/PMC4896710
38. Mortimer, B., et al. (2019). A spider's vibration landscape: adaptations to promote vibrational information transfer in orb webs. *Integr. Comp. Biol.* 59, 1636–1645 (author list unverified). https://ora.ox.ac.uk/objects/uuid:6fa55b00-3f1a-4889-a231-45b1c45003d4
39. Mortimer, B., et al. (2018). Remote monitoring of vibrational information in spider webs. *Sci. Nat.* (details unverified). https://e-archivo.uc3m.es/entities/publication/2c98d0f4-59d9-4b2c-a003-b5000e0a6caa/full
40. Mortimer, B., et al. (2019). Decoding the locational information in the orb web vibrations of *Araneus diadematus* and *Zygiella x-notata* (details unverified). https://research-information.bris.ac.uk/en/publications/4bb31410-901f-405d-8555-93882edbcc4e
41. Kinematics of prey capture by gumfoot threads of the black widow (*Latrodectus hesperus*) (SICB abstract; authors and details unverified); The efficiency of black widow gumfoot lines (student abstract). https://sicb.org/?p=33925 ; https://scholarlycommons.pacific.edu/purcc/2002/events/29
42. Kang, D., et al. (2014). Ultrasensitive mechanical crack-based sensor inspired by the spider sensory system. *Nature* 516, 222–226. https://ideas.repec.org/a/nat/nature/v516y2014i7530d10.1038_nature14002.html
43. Wearable LIG flexible stress sensor based on spider web bionic structure (2023). *Coatings* 13(1), 155 (authors not captured). https://www.mdpi.com/2079-6412/13/1/155
44. Hormiga, G., Eberhard, W.G. (2023). Survey of linyphiid web architecture. *Bull. Mus. Comp. Zool.* 163 (title and details unverified); Encyclopedia of Life, Linyphiidae. https://bioone.org/journals/bulletin-of-the-museum-of-comparative-zoology/current ; https://www.eol.org/pages/8781/articles
45. Benjamin, S.P., et al. (2002). On linyphiid sheet webs lacking aqueous glue. *Acta Zoologica* (title and details unverified). https://edoc.unibas.ch/12897
46. Benjamin, S.P., Zschokke, S. (2003). Webs of theridiid spiders: construction, structure and evolution. *Biol. J. Linn. Soc.* 78, 293–305. https://edoc.unibas.ch/12740
47. Sahni, V., Harris, J., Blackledge, T.A., Dhinojwala, A. (2012). Cobweb-weaving spiders produce different attachment discs for locomotion and prey capture. *Nat. Commun.* (author list and details unverified). https://ideaexchange.uakron.edu/biology_ideas/184
48. Su, I., Qin, Z., Saraceno, T., Krell, A., Muehlethaler, R., Bisshop, A., Buehler, M.J. (2018). Imaging and analysis of a three-dimensional spider web architecture. *J. R. Soc. Interface* 15(146), 20180193. https://research-repository.griffith.edu.au/items/34c5bd03-fe54-4b0f-88d1-d4d12f44f5e3 ; https://cee.mit.edu/professor-markus-buehler-and-his-laboratory-for-atomistic-and-molecular-mechanics-published-research-that-shows-the-3d-architecture-of-spider-webs
49. Su, I., Buehler, M.J. (2020). Mesomechanics of a three-dimensional spider web. *J. Mech. Phys. Solids* 144, 104096. https://lammps.org/abstracts/abstract.20401.html
50. Krink, T., Vollrath, F. (1997). Analysing spider web-building behaviour with rule-based simulations and genetic algorithms. *J. Theor. Biol.* 185, 321–331 (details unverified). https://www.scribd.com/document/259242728/Spider-web-algorithm ; https://minds.wisconsin.edu/handle/1793/77342
51. Krink, T., Vollrath, F. (1998). Emergent properties in the behaviour of a virtual spider robot. *Proc. R. Soc. Lond. B* 265, 2051–2055 (title wording unverified). https://www.neuroscience.ox.ac.uk/publications/1311360
52. Vollrath, F., Krink, T. (2020). Spider webs inspiring soft robotics. *J. R. Soc. Interface* 17(172), 20200569. https://pmc.ncbi.nlm.nih.gov/articles/PMC7729045
53. Miniaci, M., Krushynska, A., Movchan, A.B., Bosia, F., Pugno, N.M. (2016). Spider web-inspired acoustic metamaterials. *Appl. Phys. Lett.* 109(7), 071905; Krushynska, A.O., Bosia, F., Miniaci, M., Pugno, N.M. (2017). Spider web-structured labyrinthine acoustic metamaterials for low-frequency sound control. *New J. Phys.* (details unverified). https://pugno.dicam.unitn.it/NP_PDF/299-APL16-spiderwebmetamaterials.pdf ; https://arxiv.org/pdf/1607.01407 ; https://doaj.org/article/6dd02d3e307443089bcc255832dd29b4
54. Energy absorber inspired by spider webs (2025). *Composite Structures* (authors and details unverified). https://www.sciencedirect.com/science/article/abs/pii/S0263822325003253
55. Dynamic behaviour of spiderweb-inspired lattice structures for enhanced crashworthiness and energy absorption under impact loading (2026) (authors and details unverified). https://www.sciencedirect.com/science/article/pii/S2590123026006432
56. Bionic spiderweb lattice metamaterials for energy absorption and vibration isolation (2025). *Int. J. Mech. Sci.* (authors and details unverified). https://www.sciencedirect.com/science/article/abs/pii/S0020740325004722
57. Investigation on the energy-absorbing properties of bionic spider web structure (2023). *Biomimetics* (authors not captured). https://pmc.ncbi.nlm.nih.gov/articles/PMC10669611/
58. Elleithy, ?, Acar, ? (2026). Stochastic spider webs: bio-inspired metamaterials design and structural optimization. *Smart Mater. Struct.*, doi 10.1088/1361-665X/aeb83d (surfaced once in search; not confirmed; details unverified). https://iopscience.iop.org/journal/0964-1726
59. Ding, X., Yamazaki, K. (2004). Adaptive growth technique of stiffener layout pattern for plate and shell structures (as described in later reviews; details unverified). https://www.frontiersin.org/articles/10.3389/fmats.2023.1277421/full
60. CN101950316A. Design method of bionic spider web composite material structure (Google Patents listing; unverified). https://patents.google.com/patent/CN101950316A/en
61. PatSnap Eureka, "spider web" patent topic listing (CN115561484B and others; unverified). https://eureka.patsnap.com/topic-patents-spider-web
62. Searches for US spider-web-inspired lattice/metamaterial/sensor patents returned no matching grants (patents surfaced, e.g. US8579073B2, are unrelated to webs). https://patents.google.com/patent/US8579073B2/en
