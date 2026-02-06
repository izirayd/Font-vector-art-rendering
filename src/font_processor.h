#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <string>
#include <filesystem>
#include <unordered_map>

// Triangle type constants matching HLSL
static constexpr uint32_t PRIM_SOLID   = 0;
static constexpr uint32_t PRIM_CONVEX  = 1;
static constexpr uint32_t PRIM_CONCAVE = 2;

struct Glyphlet {
    std::vector<std::array<float, 2>>      vertices;   // 2D positions
    std::vector<std::array<uint32_t, 3>>   triangles;  // index triples into vertices
    std::vector<uint32_t>                  primAttrs;   // SOLID/CONVEX/CONCAVE per triangle
    float advanceWidth = 0.0f;                          // horizontal advance in normalized units
    float lsb          = 0.0f;                          // left side bearing
};

// GPU-side struct (must match HLSL)
struct GlyphletInfo {
    uint32_t vertexBaseIndex;
    uint32_t triangleBaseIndex;
    uint32_t vertexCount;
    uint32_t primitiveCount;
};

// GPU-side struct (must match HLSL)
struct CharacterRenderInfo {
    float    posX;
    float    posY;
    uint32_t character;
    uint32_t padding;
};

// Packed font data ready for GPU upload
struct FontGPUData {
    std::vector<std::array<float, 2>>      allVertices;
    std::vector<std::array<uint32_t, 3>>   allTriangles;
    std::vector<uint32_t>                  allPrimAttrs;
    std::vector<GlyphletInfo>              glyphletInfos;  // indexed by ASCII code (256 entries)
};

class FontProcessor {
public:
    FontProcessor() = default;
    ~FontProcessor();

    FontProcessor(const FontProcessor&) = delete;
    FontProcessor& operator=(const FontProcessor&) = delete;
    FontProcessor(FontProcessor&& other) noexcept;
    FontProcessor& operator=(FontProcessor&& other) noexcept;

    // Load a font file (TTF/OTF/TTC). Returns false on failure.
    bool loadFont(const std::filesystem::path& fontPath, float pixelHeight);

    // Process a specific character and return its glyphlet.
    Glyphlet processGlyph(int codepoint);

    // Build packed GPU data for a set of characters.
    FontGPUData buildGPUData(const std::string& characters);

    // Build CharacterRenderInfo array to lay out a string.
    std::vector<CharacterRenderInfo> layoutString(const std::string& text, float startX, float startY);

private:
    std::vector<uint8_t>                         m_fontFileData;
    void*                                        m_fontInfo = nullptr;  // stbtt_fontinfo*
    float                                        m_scale    = 1.0f;
    std::unordered_map<int, Glyphlet>            m_cache;
};
