# Patent landscape (preliminary)

Date of search: 2026-10-28.

## Scope and limits of this search

- Method: web search only. This session's network policy blocked direct access to Google Patents, USPTO full text, Justia, OSTI and Espacenet, so claim text was **not** read for any patent below. Titles, numbers and summaries come from search-engine excerpts and should be confirmed at the source.
- Not legal advice. Freedom to operate needs a patent attorney reading the claims, legal status and expiry of each document in each country where the software will be sold or used.
- To deepen this search, allow these hosts in the environment's network settings: `patents.google.com`, `ppubs.uspto.gov`, `image-ppubs.uspto.gov`, `patents.justia.com`, `worldwide.espacenet.com`, `osti.gov`.

## Summary

| Area (methods doc section) | Patent activity found | Risk to this project |
|---|---|---|
| Molecular-dynamics / force-balanced graded lattices (5.3) | US 11738515 and continuation US 12233603 (Oak Ridge National Laboratory inventors) | **High: read claims before building 5.3** |
| Level-set TO with machining constraints (1.4, constraints) | Autodesk US 11734471, US 12210804; publication US 2018/0349531 | Medium if we add CNC constraints in the same way |
| Generative design with lattice infill, hollow topology | Autodesk US 10635088 and a family sharing the same specification (US 11321508, 11675333, 11663379, 12147740, 12223238) | Medium: claim scope unread |
| Stress-graded gyroid / TPMS | US 2020/0171753 A1 (reported abandoned); CN 111291512 A (thermal grading, China) | Low to medium |
| Lattice TO preserving structure, implicit field | WO 2015/106021 A1 (assignee not confirmed) | Unknown |
| Microstructure-based TO for AM | US 2020/0180228 (application; status unknown) | Unknown |
| Tree supports | US 9688024 "Adaptive supports for 3D printing"; applications US 2018/0065311, US 2018/0162058 | Medium for support generation (3.2) |
| Bone-derived porous metal structures | US 9993341 (implant-specific) | Low for general engineering use |
| Protective lattices | US 11864610 (multilayer lattice protective equipment), US 12303766 B2 (AM helmet energy attenuators) | Medium for helmet/padding products, low for the software |
| Allostery-inspired pruning (6.2) | None found | Low, published academic prior art since 2015–2017 |
| Physarum networks (3.1) | None found | Low |
| Diffusion-model initial guess for TO (6.6) | None found granted; extensive academic prior art from 2022 | Low to medium (fast-moving area) |
| Elastic network soft-mode analysis of lattices (6.1) | None found | Low |
| Replica exchange for design diversity (6.5) | None found | Low |
| Mattheck CAO/SKO, tension triangles (2.3) | Not searched in depth; methods published 1990–1998 | Low (any patents from that era have expired) |
| Homogenization / de-homogenization (4.5) | No patent found in search; dense academic literature | Low to medium |

## Details

### 1. Molecular-dynamics graded lattices: closest overlap

- US 11738515, "Systems and methods for generating a molecular dynamic graded lattice structure and their application to additive manufacturing." Inventors reported as Seokpum Kim, Ahmed Arabi Hassen, Lonnie J. Love, Vlastimil Kunc (Oak Ridge National Laboratory). Listed by the DOE patent database (OSTI 2222257).
- US 12233603 B2, "Slicer systems for generating a molecular dynamic graded lattice structure ..." (continuation).
- Reported method: node positions are force-balanced using a field (stress or temperature) to set node spacing. High-intensity regions get smaller cells. Triangulation and its dual graph produce a non-uniform polygonal (typically honeycomb) structure, then AM instructions.
- Overlap with our methods doc section 5.3 (particle relaxation with spacing from the stress field, then Delaunay or Voronoi) is direct.
- Options:
  1. Read the claims. Claims may be limited, for example to 2D honeycomb, to slicer output, or to a specific force law.
  2. Ask ORNL about licensing. DOE labs license routinely.
  3. Rely on prior art for anything outside the claims. Shimada and Gossard's bubble mesh (1995) and centroidal Voronoi tessellation (Lloyd 1982; Du, Faber and Gunzburger 1999) predate the filing.
  4. Design around it, for example by generating lattices by de-homogenization (4.5) instead of particle relaxation.
- Recommendation: do not build 5.3 until the claims are read.

### 2. Autodesk generative design patents

- US 11734471 and US 12210804: level-set TO that moves the boundary into a feasible region for given tool bits, tool heads and machining directions. Matches Morris, Butscher and Iorio (2019), "A subtractive manufacturing constraint for level set topology optimization" (Autodesk Research).
- US 10635088 "Hollow topology generation with lattices," plus a family sharing one specification covering damage prevention over load cycles, singularity and disconnection prevention, build-material strength models and controlled convergence.
- Implication: when we add milling constraints, use a formulation from the open literature that differs from the claimed one. Examples: Langelaar's AM filter (2016), or Gersborg and Andreasen's casting filter (2011). Have the claims reviewed first.

### 3. Graded TPMS

- US 2020/0171753 A1 "Additively-manufactured gradient gyroid lattice structures": filed 2018, reported abandoned. Abandoned applications become prior art that anyone can use.
- CN 111291512 A: variable-density gyroid for thermal insulation (China only).
- Academic prior art for stress-driven graded gyroids: Li et al., Computer-Aided Design (2018).

### 4. Tree supports

- US 9688024 "Adaptive supports for 3D printing" covers tree-like supports with branches ending at contact points.
- Autodesk Research, "Branching support structures for 3D printing" (Vanek, Galicia and Benes 2014 is a separate academic precursor). Schmidt and Umetani, "Branching support structures for 3D printing," SIGGRAPH Talks 2014.
- Implication: before shipping a support generator (3.2), read the claims of US 9688024.

### 5. Areas with no patents found

No patents turned up for allostery-style pruning, Physarum-based structural design, elastic-network soft-mode analysis of lattices, replica exchange for design diversity, or diffusion-model warm starts for TO. Each area has dated academic prior art, which makes later patents in those areas harder to obtain and easier to challenge. This is where the project has the most room, and it lines up with the features that most set it apart from Fusion (methods doc sections 6.1, 6.2, 6.5).

## Sources

Search results consulted (accessed 2026-10-28):

- US 11738515 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11738515
- US 12233603 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12233603
- OSTI DOE patent record 2222257: https://osti.gov/doepatents/biblio/2222257
- US 11734471 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11734471
- US 2018/0349531 (FreePatentsOnline): https://www.freepatentsonline.com/y2018/0349531.html
- Morris, Butscher, Iorio (2019), Autodesk Research PDF: https://damassets.autodesk.net/content/dam/autodesk/research/publications-assets/pdf/a-subtractive-manufacturing-constraint.pdf
- US 10635088 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10635088
- US 11321508 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11321508
- US 11675333 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11675333
- US 11663379 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11663379
- US 12147740 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12147740
- US 12223238 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12223238
- US 2020/0171753 A1 (Google Patents): https://patents.google.com/patent/US20200171753A1/en
- CN 111291512 A (Google Patents): https://patents.google.com/patent/CN111291512A/en
- WO 2015/106021 A1 (Google Patents): https://patents.google.com/patent/WO2015106021A1/en
- US 2020/0180228 (Justia): https://patents.justia.com/patent/20200180228
- US 9688024 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/9688024
- US 9993341 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/9993341
- US 11864610 (USPTO PDF): https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11864610
- US 12303766 B2 (USPTO Official Gazette): https://patentsgazette.uspto.gov/week20/OG/html/1534-3/US12303766-20250520.html
- Rocks et al. (2017) PNAS, allostery-inspired networks: https://www.pnas.org/doi/10.1073/pnas.1612139114
- Stress-field driven conformal lattice by circle packing (Heliyon 2023): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10031376/
- Mazé and Ahmed, Diffusion models beat GANs on topology optimization: https://arxiv.org/pdf/2208.09591
- LEAP 71 PicoGK repository (no patents referenced): https://github.com/leap71/PicoGK
