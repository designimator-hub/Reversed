# Architecture: hosted interface, local GPU

## Goal

The interface and accounts run on a server. The heavy compute runs on a GPU the user owns, or on a rented GPU when they have none. You build one solver, and it runs in both places unchanged.

## Components

```
 ┌────────────────────┐      ┌────────────────────┐      ┌────────────────────┐
 │  Front ends        │      │  Hosted server     │      │  Compute node      │
 │  - Web app         │◄────►│  - accounts        │◄────►│  (user's GPU or    │
 │  - Blender add-on  │ WSS  │  - job queue       │ WSS  │   rented GPU)      │
 │  - Fusion add-in   │      │  - relay           │ out- │  - FE solver       │
 └────────────────────┘      │  - result storage  │ bound│  - optimizer       │
          ▲                  └────────────────────┘      │  - SDF/voxel kernel│
          │                                              │  - meshing         │
          └──────── direct localhost (optional) ────────►└────────────────────┘
```

### Compute node ("solver daemon")

- Language: Python on top of NVIDIA Warp (Apache 2.0). Warp gives you CUDA kernels written in Python, sparse NanoVDB volumes, and a finite-element module (`warp.fem`). Hot loops that outgrow Warp move to C++/CUDA later.
- FE solver: matrix-free, on a voxel grid, with conjugate gradient and a geometric multigrid preconditioner. Wu, Dick and Westermann (2016) showed this solving millions of elements per iteration on a single GPU.
- Geometry: sparse signed-distance and density fields held as NanoVDB/OpenVDB grids. This is the same representation PicoGK uses, which keeps the Leap 71 libraries in reach.
- Fallback: a CPU path, so the node still runs on a machine without an NVIDIA GPU, only slower. An AMD/Apple path through WebGPU-native (wgpu) or Vulkan compute is a later option.

### Connection model

- The compute node opens an outbound WebSocket to the hosted server. This avoids port forwarding, firewall rules and the browser restrictions on public sites calling `localhost` (Chrome's Private Network Access rules).
- When the web app and the node sit on the same machine, a direct localhost connection can skip the relay for live previews.
- Each job is a JSON spec (domain, loads, fixtures, materials, constraints, objective) plus binary field data. Because the spec is the whole input, any job can be re-run, compared or handed to another node.

### Front ends

- Web app: Three.js viewport, a setup panel (loads, supports, keep/avoid regions), live density preview streamed from the node, and an outcome gallery.
- Blender add-on: Blender already ships OpenVDB, volume objects, and mesh-to-volume and volume-to-mesh nodes, so results drop in as native volumes or meshes. Licensing: any add-on that imports `bpy` must be GPL-compatible. Keep the add-on as a thin GPL client that talks to the solver over a socket, and keep the solver as a separate process under its own license.
- Fusion add-in: Fusion's Python/C++ API lets you pull geometry and setup out of a Fusion design and push the result back as a mesh. You keep Fusion for precise B-rep work and use this tool for generation. That plays to what each one does well.

## Why not Blender alone

Blender is a strong mesh and volume tool, but it lacks a units-exact B-rep kernel, constraint sketches and engineering drawings. As a front end it's a good choice; as the core it would limit the product.
