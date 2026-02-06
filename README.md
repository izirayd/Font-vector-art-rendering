# Font & Vector Art Rendering with Mesh Shaders

A DirectX 12 implementation of resolution-independent font and vector art rendering using **mesh shaders** and **per-primitive attributes**, based on the AMD GPUOpen article:

> [Font- and vector-art rendering with mesh shaders](https://gpuopen.com/learn/mesh_shaders/mesh_shaders-font_and_vector_art_rendering_with_mesh_shaders/) — Quirin Meyer, Bastian Kuth, Max Oberberger (AMD, 2024)

The core idea originates from the classic **Loop & Blinn (SIGGRAPH 2005)** technique for resolution-independent curve rendering, extended here with modern mesh shader capabilities available on AMD RDNA 2+ and NVIDIA Turing+ GPUs.

## Overview

This project was implemented with the assistance of an **LLM (Large Language Model)** following the concepts described in the AMD GPUOpen article. The LLM translated the theoretical approach from the blog post into a working Direct3D 12 application with:

- **TTF font parsing** via `stb_truetype`
- **Constrained triangulation** via `mapbox/earcut`
- **Mesh shader pipeline** with per-primitive attributes (HLSL SM 6.5)
- **Single dispatch per string** — each thread group renders one glyph

### How It Works

1. **Glyph Processing (CPU):** TrueType font glyphs are decomposed into contours. Each contour consists of linear segments and quadratic Bezier curves. A constrained Delaunay-style triangulation (earcut) produces solid triangles. Quadratic curve triangles are classified as **convex** or **concave** based on control point orientation.

2. **GPU Buffers:** All glyphlet data (vertices, indices, per-primitive attributes, glyphlet metadata) is packed into structured buffers and uploaded to the GPU once.

3. **Mesh Shader:** A single `DispatchMesh(numCharacters, 1, 1)` call renders the entire string. Each thread group reads its character from a text buffer, fetches the corresponding glyphlet, and emits vertices with transformed positions and per-primitive triangle type attributes.

4. **Pixel Shader:** Uses `SV_Barycentrics` to evaluate the canonical quadratic Bezier function `u² - v`. Fragments are discarded based on the per-primitive attribute:
   - **SOLID** — always filled
   - **CONVEX** — discard where `u² - v > 0`
   - **CONCAVE** — discard where `u² - v < 0`

This produces infinitely smooth glyph silhouettes at any zoom level without tessellation.

## Project Structure

```
Font-vector-art-rendering/
├── CMakeLists.txt              # Build configuration
├── build.bat                   # One-click build script (MSVC + Ninja)
├── shaders/
│   └── font_shader.hlsl        # Mesh shader + pixel shader (SM 6.5)
├── src/
│   ├── main.cpp                # Win32 window, D3D12 init, render loop
│   ├── font_processor.h        # Font processing interface
│   └── font_processor.cpp      # TTF parsing, triangulation, GPU data packing
├── third_party/
│   ├── stb_truetype.h          # TTF font loader (single-header)
│   └── mapbox/earcut.hpp       # Polygon triangulation (single-header)
└── test_font.cpp               # Console test for glyph processing
```

## Requirements

- **OS:** Windows 10/11
- **GPU:** Mesh shader support required — AMD RDNA 2 (RX 6000+) or NVIDIA Turing (RTX 20xx+)
- **Compiler:** MSVC 2019+ (C++17)
- **CMake:** 3.20+
- **Windows SDK:** 10.0.22000.0+ (includes DXC shader compiler)

## Building

### Option 1: Using build.bat

```bat
build.bat
```

This sets up MSVC environment, runs CMake with Ninja, and compiles everything including shader compilation via DXC.

### Option 2: Manual CMake build

```bat
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022" -A x64
cmake --build . --config Release
```

The build system automatically compiles `shaders/font_shader.hlsl` into header files (`font_ms.h`, `font_ps.h`) using the DXC compiler:
- Mesh shader: `-T ms_6_5 -E MSMain`
- Pixel shader: `-T ps_6_5 -E PSMain`

### Running

```bat
build\Release\MeshShaderFont.exe
```

or in Debug:

```bat
build\Debug\MeshShaderFont.exe
```

The application opens a 1280x720 window and renders the text `"QWERTY18976BP"` using Arial. Press **Escape** to exit.

## Key Technical Details

| Feature | Implementation |
|---|---|
| Font format | TrueType (.ttf) via `stb_truetype` |
| Curve types | Linear segments + quadratic Bezier |
| Triangulation | `mapbox::earcut` with constrained edges |
| Shader model | SM 6.5 (mesh + pixel shader) |
| Per-primitive attributes | Triangle type (SOLID / CONVEX / CONCAVE) via `BLENDINDICES0` |
| Barycentric evaluation | `SV_Barycentrics` maps to canonical Bezier `u²-v` |
| Rendering | Single `DispatchMesh()` call per string |
| Projection | Orthographic, auto-fitted to text bounding box |

## References

- [Font- and vector-art rendering with mesh shaders — AMD GPUOpen](https://gpuopen.com/learn/mesh_shaders/mesh_shaders-font_and_vector_art_rendering_with_mesh_shaders/)
- Loop, C. and Blinn, J. *Resolution Independent Curve Rendering using Programmable Graphics Hardware.* SIGGRAPH 2005.
- [stb_truetype](https://github.com/nothings/stb) — Sean Barrett
- [earcut.hpp](https://github.com/mapbox/earcut.hpp) — Mapbox

## License

This project is provided as-is for educational and research purposes.
