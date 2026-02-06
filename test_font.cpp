#include <cstdio>
#include "font_processor.h"

int main() {
    FontProcessor fp;
    if (!fp.loadFont("C:\\Windows\\Fonts\\arial.ttf", 100.0f)) {
        printf("Failed to load font\n");
        return 1;
    }
    printf("Font loaded successfully.\n");

    const char* text = "QWERTY18976BP";
    for (int i = 0; text[i]; ++i) {
        Glyphlet gl = fp.processGlyph(text[i]);
        printf("  '%c': %zu vertices, %zu triangles (solid=%d, convex=%d, concave=%d), advance=%.2f\n",
               text[i],
               gl.vertices.size(), gl.triangles.size(),
               (int)std::count(gl.primAttrs.begin(), gl.primAttrs.end(), PRIM_SOLID),
               (int)std::count(gl.primAttrs.begin(), gl.primAttrs.end(), PRIM_CONVEX),
               (int)std::count(gl.primAttrs.begin(), gl.primAttrs.end(), PRIM_CONCAVE),
               gl.advanceWidth);
    }

    FontGPUData gpuData = fp.buildGPUData(text);
    printf("\nGPU Data:\n");
    printf("  Total vertices: %zu\n", gpuData.allVertices.size());
    printf("  Total triangles: %zu\n", gpuData.allTriangles.size());
    printf("  Total prim attrs: %zu\n", gpuData.allPrimAttrs.size());

    auto layout = fp.layoutString(text, 0.0f, 0.0f);
    printf("\nLayout:\n");
    for (auto& ci : layout) {
        printf("  '%c' at (%.2f, %.2f)\n", (char)ci.character, ci.posX, ci.posY);
    }
    return 0;
}
