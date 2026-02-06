#define STB_TRUETYPE_IMPLEMENTATION
#include "stb_truetype.h"
#include "font_processor.h"

#include <mapbox/earcut.hpp>

#include <algorithm>
#include <cmath>
#include <cassert>
#include <fstream>
#include <numeric>
#include <set>
#include <utility>

// ---------------------------------------------------------------------------
// earcut adapter for std::array<float,2>
// ---------------------------------------------------------------------------
namespace mapbox { namespace util {
template<> struct nth<0, std::array<float,2>> {
    static float get(const std::array<float,2>& p) { return p[0]; }
};
template<> struct nth<1, std::array<float,2>> {
    static float get(const std::array<float,2>& p) { return p[1]; }
};
}}

// ---------------------------------------------------------------------------
// Internal types
// ---------------------------------------------------------------------------

struct CurveSegment {
    uint32_t a;          // start vertex index (in allVerts)
    uint32_t c;          // control point index
    uint32_t b;          // end vertex index
    bool     isConvex;
};

struct Contour {
    std::vector<std::array<float,2>> ringPts;     // polygon ring for earcut
    std::vector<uint32_t>            ringToFull;   // ring index -> allVerts index
    std::vector<CurveSegment>        curves;       // quadratic Bezier curves
    bool                             isHole = false;
};

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

static constexpr float VERTEX_EPS = 1.00f;  // dedup tolerance in scaled font units

// Adaptive curve subdivision parameters
// NOTE: threshold set very high to effectively DISABLE subdivision.
// Loop-Blinn curve rendering needs large CONVEX/CONCAVE triangles for
// the GPU shader to smoothly evaluate curve boundaries. Aggressive
// subdivision shrinks these triangles to sub-pixel slivers, defeating
// the purpose. The shader handles smooth curves at any resolution.
static constexpr float CURVE_FLATNESS_THRESHOLD = 9999.0f;  // effectively disabled
static constexpr int   CURVE_MAX_SUBDIVISION    = 50;

// Signed area of a polygon ring (positive = CCW in Y-up)
static float signedArea(const std::vector<std::array<float,2>>& ring)
{
    float area = 0.0f;
    size_t n = ring.size();
    for (size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;
        area += ring[i][0] * ring[j][1];
        area -= ring[j][0] * ring[i][1];
    }
    return area * 0.5f;
}

// Cross product of 2D vectors: (b-a) x (c-a)
static float cross2D(const std::array<float,2>& a,
                     const std::array<float,2>& b,
                     const std::array<float,2>& c)
{
    float bax = b[0] - a[0], bay = b[1] - a[1];
    float cax = c[0] - a[0], cay = c[1] - a[1];
    return bax * cay - bay * cax;
}

// Add a vertex, deduplicating by position within VERTEX_EPS.
static uint32_t addVertex(std::vector<std::array<float,2>>& verts,
                          const std::array<float,2>& v)
{
    for (uint32_t i = 0; i < (uint32_t)verts.size(); ++i) {
        float dx = verts[i][0] - v[0];
        float dy = verts[i][1] - v[1];
        if (dx*dx + dy*dy < VERTEX_EPS * VERTEX_EPS)
            return i;
    }
    verts.push_back(v);
    return (uint32_t)verts.size() - 1;
}

// Check if two points are the same within tolerance.
static bool samePoint(const std::array<float,2>& a, const std::array<float,2>& b)
{
    float dx = a[0] - b[0], dy = a[1] - b[1];
    return dx*dx + dy*dy < VERTEX_EPS * VERTEX_EPS;
}

// ---------------------------------------------------------------------------
// Adaptive quadratic Bezier subdivision (De Casteljau)
// Recursively splits curves until they are flat enough, appending sub-curves
// and intermediate on-curve points into the contour.
// ---------------------------------------------------------------------------
static void subdivideCurve(
    const std::array<float,2>& a,
    const std::array<float,2>& ctrl,
    const std::array<float,2>& b,
    float threshold,
    int maxDepth,
    int depth,
    std::vector<std::array<float,2>>& allVerts,
    Contour& contour,
    const std::array<float,2>& firstPt)
{
    // Compute flatness: squared perpendicular distance from ctrl to chord a-b
    float abx = b[0] - a[0], aby = b[1] - a[1];
    float acx = ctrl[0] - a[0], acy = ctrl[1] - a[1];
    float abLenSq = abx * abx + aby * aby;
    float crossVal = abx * acy - aby * acx;

    // flatnessSq = (perpendicular distance)^2
    float flatnessSq;
    if (abLenSq > 1e-12f)
        flatnessSq = (crossVal * crossVal) / abLenSq;
    else
        flatnessSq = acx * acx + acy * acy;  // degenerate chord: use distance a->ctrl

    float threshSq = threshold * threshold;

    if (flatnessSq <= threshSq || depth >= maxDepth) {
        // Flat enough — emit this curve segment as-is
        uint32_t aIdx = addVertex(allVerts, a);
        uint32_t cIdx = addVertex(allVerts, ctrl);
        uint32_t bIdx = addVertex(allVerts, b);

        CurveSegment seg;
        seg.a = aIdx;
        seg.c = cIdx;
        seg.b = bIdx;
        seg.isConvex = false; // determined later
        contour.curves.push_back(seg);

        // Add endpoint to polygon ring (skip if it closes back to contour start)
        if (!samePoint(b, firstPt)) {
            contour.ringPts.push_back(b);
            contour.ringToFull.push_back(bIdx);
        }
        return;
    }

    // De Casteljau split at t = 0.5
    std::array<float,2> m01 = { (a[0] + ctrl[0]) * 0.5f, (a[1] + ctrl[1]) * 0.5f };
    std::array<float,2> m12 = { (ctrl[0] + b[0]) * 0.5f, (ctrl[1] + b[1]) * 0.5f };
    std::array<float,2> mid = { (m01[0] + m12[0]) * 0.5f, (m01[1] + m12[1]) * 0.5f };

    // First half: (a, m01, mid)
    subdivideCurve(a,   m01, mid, threshold, maxDepth, depth + 1, allVerts, contour, firstPt);
    // Second half: (mid, m12, b)
    subdivideCurve(mid, m12, b,   threshold, maxDepth, depth + 1, allVerts, contour, firstPt);
}

// ---------------------------------------------------------------------------
// FontProcessor
// ---------------------------------------------------------------------------

FontProcessor::~FontProcessor()
{
    if (m_fontInfo) {
        delete static_cast<stbtt_fontinfo*>(m_fontInfo);
        m_fontInfo = nullptr;
    }
}

FontProcessor::FontProcessor(FontProcessor&& other) noexcept
{
    m_fontFileData = std::move(other.m_fontFileData);
    m_fontInfo = other.m_fontInfo;
    m_scale = other.m_scale;
    m_cache = std::move(other.m_cache);

    other.m_fontInfo = nullptr;
    other.m_scale = 1.0f;
    other.m_cache.clear();
    other.m_fontFileData.clear();
}

FontProcessor& FontProcessor::operator=(FontProcessor&& other) noexcept
{
    if (this == &other)
        return *this;

    if (m_fontInfo) {
        delete static_cast<stbtt_fontinfo*>(m_fontInfo);
        m_fontInfo = nullptr;
    }

    m_fontFileData = std::move(other.m_fontFileData);
    m_fontInfo = other.m_fontInfo;
    m_scale = other.m_scale;
    m_cache = std::move(other.m_cache);

    other.m_fontInfo = nullptr;
    other.m_scale = 1.0f;
    other.m_cache.clear();
    other.m_fontFileData.clear();

    return *this;
}

bool FontProcessor::loadFont(const std::filesystem::path& fontPath, float pixelHeight)
{
    std::ifstream ifs(fontPath, std::ios::binary);
    if (!ifs) return false;

    ifs.seekg(0, std::ios::end);
    size_t sz = (size_t)ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    std::vector<uint8_t> newFontData;
    newFontData.resize(sz);
    ifs.read(reinterpret_cast<char*>(newFontData.data()), sz);
    if (!ifs) return false;

    auto* newInfo = new stbtt_fontinfo;
    if (!stbtt_InitFont(newInfo, newFontData.data(),
                        stbtt_GetFontOffsetForIndex(newFontData.data(), 0)))
    {
        delete newInfo;
        return false;
    }

    // Commit (only after full success) so a failed load does not destroy the currently loaded font.
    if (m_fontInfo) {
        delete static_cast<stbtt_fontinfo*>(m_fontInfo);
        m_fontInfo = nullptr;
    }

    m_fontFileData = std::move(newFontData);
    m_fontInfo = newInfo;
    m_scale = stbtt_ScaleForPixelHeight(newInfo, pixelHeight);
    m_cache.clear();
    return true;
}

Glyphlet FontProcessor::processGlyph(int codepoint)
{
    auto it = m_cache.find(codepoint);
    if (it != m_cache.end()) return it->second;

    auto* info = static_cast<stbtt_fontinfo*>(m_fontInfo);
    Glyphlet gl;

    int glyphIndex = stbtt_FindGlyphIndex(info, codepoint);
    if (glyphIndex == 0 && codepoint != 0) {
        m_cache[codepoint] = gl;
        return gl;
    }

    // Advance metrics
    int advW = 0, lsbI = 0;
    stbtt_GetGlyphHMetrics(info, glyphIndex, &advW, &lsbI);
    gl.advanceWidth = advW * m_scale;
    gl.lsb          = lsbI * m_scale;

    // Get shape commands from stb_truetype
    stbtt_vertex* verts = nullptr;
    int numVerts = stbtt_GetGlyphShape(info, glyphIndex, &verts);
    if (numVerts <= 0) {
        m_cache[codepoint] = gl;
        return gl;
    }

    // ---------------------------------------------------------------
    // 1. Extract contours from stb_truetype shape commands
    // ---------------------------------------------------------------
    std::vector<Contour> contours;
    std::vector<std::array<float,2>> allVerts;

    int i = 0;
    while (i < numVerts) {
        if (verts[i].type != STBTT_vmove) { ++i; continue; }

        Contour contour;
        std::array<float,2> firstPt = {verts[i].x * m_scale, verts[i].y * m_scale};
        std::array<float,2> prevPt  = firstPt;
        uint32_t firstIdx = addVertex(allVerts, firstPt);

        // First on-curve point
        contour.ringPts.push_back(firstPt);
        contour.ringToFull.push_back(firstIdx);

        ++i;
        while (i < numVerts && verts[i].type != STBTT_vmove) {
            std::array<float,2> endPt = {verts[i].x * m_scale, verts[i].y * m_scale};

            if (verts[i].type == STBTT_vline) {
                // stb_truetype closes contours with a line back to start.
                // Skip duplicate closing vertex.
                if (!samePoint(endPt, firstPt)) {
                    uint32_t idx = addVertex(allVerts, endPt);
                    contour.ringPts.push_back(endPt);
                    contour.ringToFull.push_back(idx);
                }
                prevPt = endPt;
            }
            else if (verts[i].type == STBTT_vcurve) {
                std::array<float,2> ctrlPt = {verts[i].cx * m_scale, verts[i].cy * m_scale};

                // Adaptive subdivision: recursively split sharp curves into
                // smaller sub-curves for smoother polygon rings and better
                // triangle quality. Flat curves pass through as a single segment.
                subdivideCurve(prevPt, ctrlPt, endPt,
                               CURVE_FLATNESS_THRESHOLD,
                               CURVE_MAX_SUBDIVISION,
                               0, allVerts, contour, firstPt);

                prevPt = endPt;
            }
            ++i;
        }

        if (contour.ringPts.size() >= 3)
            contours.push_back(std::move(contour));
    }
    stbtt_FreeShape(info, verts);

    if (contours.empty()) {
        m_cache[codepoint] = gl;
        return gl;
    }

    // ---------------------------------------------------------------
    // 2. Identify outer contour vs holes
    //    We treat the largest contour (by |signed area|) as the outer ring
    //    and all other contours as holes. This matches how we feed earcut.
    // ---------------------------------------------------------------
    size_t outerContourIndex = 0;
    float  outerAbsArea      = 0.0f;
    std::vector<float> contourSignedAreas(contours.size(), 0.0f);
    for (size_t ci = 0; ci < contours.size(); ++ci) {
        float a = signedArea(contours[ci].ringPts);
        contourSignedAreas[ci] = a;
        float aa = std::fabs(a);
        if (aa > outerAbsArea) {
            outerAbsArea = aa;
            outerContourIndex = ci;
        }
    }
    for (size_t ci = 0; ci < contours.size(); ++ci) {
        contours[ci].isHole = (ci != outerContourIndex);
    }

    // ---------------------------------------------------------------
    // 3. Determine convex/concave for each curve
    //    convex  = control point is OUTSIDE the filled area
    //    concave = control point is INSIDE the filled area
    //
    // IMPORTANT: for hole contours the filled area is OUTSIDE the ring,
    // so convex/concave must be flipped relative to the contour's interior.
    // ---------------------------------------------------------------
    for (size_t ci = 0; ci < contours.size(); ++ci) {
        auto& contour = contours[ci];
        float area = contourSignedAreas[ci];
        for (auto& curve : contour.curves) {
            // cross = (endPt - startPt) x (ctrlPt - startPt)
            float cr = cross2D(allVerts[curve.a], allVerts[curve.b], allVerts[curve.c]);
            // If cross and area have opposite signs -> control is outside -> convex
            bool isConvex = (cr * area) < 0.0f;
            if (contour.isHole)
                isConvex = !isConvex;
            curve.isConvex = isConvex;
        }
    }

    // ---------------------------------------------------------------
    // 4. Build earcut polygon rings
    //    - For convex curves: polygon edge goes directly a->b (skip c)
    //    - For concave curves: insert c between a and b in the ring
    // ---------------------------------------------------------------
    for (auto& contour : contours) {
        // Build a map: for each on-curve endpoint index, find the concave
        // curve that ends there (if any), so we can insert c before b.
        // Key = fullVerts index of curve.b, Value = curve index
        std::unordered_map<uint32_t, size_t> concaveEndMap;
        for (size_t ci = 0; ci < contour.curves.size(); ++ci) {
            if (!contour.curves[ci].isConvex)
                concaveEndMap[contour.curves[ci].b] = ci;
        }

        std::vector<std::array<float,2>> newRing;
        std::vector<uint32_t>            newMap;

        for (size_t pi = 0; pi < contour.ringPts.size(); ++pi) {
            uint32_t fullIdx = contour.ringToFull[pi];

            // If a concave curve ends at this vertex, insert its control point first
            auto it2 = concaveEndMap.find(fullIdx);
            if (it2 != concaveEndMap.end()) {
                auto& curve = contour.curves[it2->second];
                newRing.push_back(allVerts[curve.c]);
                newMap.push_back(curve.c);
            }

            newRing.push_back(contour.ringPts[pi]);
            newMap.push_back(fullIdx);
        }

        contour.ringPts    = std::move(newRing);
        contour.ringToFull = std::move(newMap);
    }

    // ---------------------------------------------------------------
    // 5. Triangulate with earcut
    //    Outer ring first, then holes (matches isHole marking above).
    // ---------------------------------------------------------------
    using Polygon = std::vector<std::vector<std::array<float,2>>>;
    Polygon polygon;
    std::vector<uint32_t> flatToFull;  // earcut flat index -> allVerts index

    auto appendRing = [&](const Contour& c) {
        polygon.push_back(c.ringPts);
        flatToFull.insert(flatToFull.end(), c.ringToFull.begin(), c.ringToFull.end());
    };

    // Build earcut input: [outer, hole0, hole1, ...]
    appendRing(contours[outerContourIndex]);
    for (size_t ci = 0; ci < contours.size(); ++ci) {
        if (ci == outerContourIndex) continue;
        appendRing(contours[ci]);
    }

    auto earcutIdx = mapbox::earcut<uint32_t>(polygon);

    // ---------------------------------------------------------------
    // 6. Build the final glyphlet
    // ---------------------------------------------------------------
    gl.vertices = allVerts;

    // 5a. Build a set of concave curve vertex triples for quick lookup
    struct TriKey {
        uint32_t v[3]; // sorted
        bool operator==(const TriKey& o) const {
            return v[0]==o.v[0] && v[1]==o.v[1] && v[2]==o.v[2];
        }
    };
    struct TriHash {
        size_t operator()(const TriKey& k) const {
            return std::hash<uint64_t>()(
                ((uint64_t)k.v[0] << 32) ^ ((uint64_t)k.v[1] << 16) ^ k.v[2]);
        }
    };
    std::unordered_map<TriKey, size_t, TriHash> concaveTriMap;
    for (auto& contour : contours) {
        for (size_t ci = 0; ci < contour.curves.size(); ++ci) {
            auto& curve = contour.curves[ci];
            if (curve.isConvex) continue;
            TriKey key;
            key.v[0] = std::min({curve.a, curve.c, curve.b});
            key.v[2] = std::max({curve.a, curve.c, curve.b});
            key.v[1] = curve.a + curve.c + curve.b - key.v[0] - key.v[2];
            concaveTriMap[key] = ci; // store curve index (not used, just presence)
        }
    }

    // Track which concave curves were found in earcut output
    std::set<size_t> foundConcaveKeys;

    // 5b. Add earcut triangles (SOLID, or reclassify as CONCAVE if matches)
    for (size_t ti = 0; ti + 2 < earcutIdx.size(); ti += 3) {
        uint32_t i0 = flatToFull[earcutIdx[ti + 0]];
        uint32_t i1 = flatToFull[earcutIdx[ti + 1]];
        uint32_t i2 = flatToFull[earcutIdx[ti + 2]];

        // Check if this matches a concave curve triangle
        TriKey key;
        key.v[0] = std::min({i0, i1, i2});
        key.v[2] = std::max({i0, i1, i2});
        key.v[1] = i0 + i1 + i2 - key.v[0] - key.v[2];

        auto cit = concaveTriMap.find(key);
        if (cit != concaveTriMap.end()) {
            // Reclassify as CONCAVE with correct vertex order (a, c, b)
            // Find the actual curve
            for (auto& contour : contours) {
                for (auto& curve : contour.curves) {
                    if (curve.isConvex) continue;
                    std::set<uint32_t> cv = {curve.a, curve.c, curve.b};
                    std::set<uint32_t> tv = {i0, i1, i2};
                    if (cv == tv) {
                        gl.triangles.push_back({curve.a, curve.c, curve.b});
                        gl.primAttrs.push_back(PRIM_CONCAVE);
                        foundConcaveKeys.insert(
                            (size_t)key.v[0] * 100000 + key.v[1] * 1000 + key.v[2]);
                        goto next_earcut_tri;
                    }
                }
            }
        }

        // Regular solid triangle
        gl.triangles.push_back({i0, i1, i2});
        gl.primAttrs.push_back(PRIM_SOLID);

        next_earcut_tri:;
    }

    // 5c. Add CONVEX curve triangles (a, c, b)
    for (auto& contour : contours) {
        for (auto& curve : contour.curves) {
            if (!curve.isConvex) continue;
            gl.triangles.push_back({curve.a, curve.c, curve.b});
            gl.primAttrs.push_back(PRIM_CONVEX);
        }
    }

    // 5d. Add CONCAVE curve triangles not found in earcut output
    for (auto& contour : contours) {
        for (auto& curve : contour.curves) {
            if (curve.isConvex) continue;
            TriKey key;
            key.v[0] = std::min({curve.a, curve.c, curve.b});
            key.v[2] = std::max({curve.a, curve.c, curve.b});
            key.v[1] = curve.a + curve.c + curve.b - key.v[0] - key.v[2];
            size_t hk = (size_t)key.v[0] * 100000 + key.v[1] * 1000 + key.v[2];
            if (foundConcaveKeys.count(hk) == 0) {
                gl.triangles.push_back({curve.a, curve.c, curve.b});
                gl.primAttrs.push_back(PRIM_CONCAVE);
            }
        }
    }

    m_cache[codepoint] = gl;
    return gl;
}

FontGPUData FontProcessor::buildGPUData(const std::string& characters)
{
    FontGPUData data;
    data.glyphletInfos.resize(256);
    memset(data.glyphletInfos.data(), 0, sizeof(GlyphletInfo) * 256);

    std::set<int> processed;
    for (char ch : characters) {
        int cp = (unsigned char)ch;
        if (processed.count(cp)) continue;
        processed.insert(cp);

        Glyphlet gl = processGlyph(cp);
        if (gl.vertices.empty()) continue;

        GlyphletInfo info;
        info.vertexBaseIndex   = (uint32_t)data.allVertices.size();
        info.triangleBaseIndex = (uint32_t)data.allTriangles.size();
        info.vertexCount       = (uint32_t)gl.vertices.size();
        info.primitiveCount    = (uint32_t)gl.triangles.size();

        data.glyphletInfos[cp] = info;

        data.allVertices.insert(data.allVertices.end(),
                                gl.vertices.begin(), gl.vertices.end());
        data.allTriangles.insert(data.allTriangles.end(),
                                 gl.triangles.begin(), gl.triangles.end());
        data.allPrimAttrs.insert(data.allPrimAttrs.end(),
                                 gl.primAttrs.begin(), gl.primAttrs.end());
    }

    return data;
}

std::vector<CharacterRenderInfo> FontProcessor::layoutString(
    const std::string& text, float startX, float startY)
{
    std::vector<CharacterRenderInfo> result;
    auto* info = static_cast<stbtt_fontinfo*>(m_fontInfo);

    float cursorX = startX;
    for (size_t i = 0; i < text.size(); ++i) {
        int cp = (unsigned char)text[i];

        CharacterRenderInfo cri;
        cri.posX      = cursorX;
        cri.posY      = startY;
        cri.character  = (uint32_t)cp;
        cri.padding    = 0;
        result.push_back(cri);

        int advW = 0, lsb = 0;
        stbtt_GetCodepointHMetrics(info, cp, &advW, &lsb);
        cursorX += advW * m_scale;

        if (i + 1 < text.size()) {
            int kern = stbtt_GetCodepointKernAdvance(info, cp, (unsigned char)text[i+1]);
            cursorX += kern * m_scale;
        }
    }

    return result;
}
