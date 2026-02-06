// ============================================================================
// Mesh Shader Font Renderer
// Based on Loop-Blinn resolution-independent curve rendering with mesh shaders.
// Per-primitive attributes used to classify triangles as SOLID / CONVEX / CONCAVE.
// SV_Barycentrics used to evaluate the canonical quadratic Bezier curve.
// ============================================================================

static const uint SOLID   = 0;
static const uint CONVEX  = 1;
static const uint CONCAVE = 2;

// ---------------------------------------------------------------------------
// Constant buffer
// ---------------------------------------------------------------------------
cbuffer Constants : register(b0)
{
    float4x4 transformMatrix;
};

// ---------------------------------------------------------------------------
// Structured buffers
// ---------------------------------------------------------------------------
struct GlyphletInfo
{
    uint vertexBaseIndex;
    uint triangleBaseIndex;
    uint vertexCount;
    uint primitiveCount;
};

struct CharRenderInfo
{
    float2 pos;
    uint   character;
    uint   padding;
};

StructuredBuffer<float2>        vertexBuffer     : register(t0);
StructuredBuffer<uint3>         indexBuffer       : register(t1);
StructuredBuffer<uint>          primAttrBuffer    : register(t2);
StructuredBuffer<GlyphletInfo>  glyphletBuffer    : register(t3);
StructuredBuffer<CharRenderInfo> textBuffer       : register(t4);

// ---------------------------------------------------------------------------
// Mesh shader output structures
// ---------------------------------------------------------------------------
struct VertOut
{
    float4 position : SV_POSITION;
};

struct PrimOut
{
    uint triangleType : BLENDINDICES0;
};

// ---------------------------------------------------------------------------
// Pixel shader input
// ---------------------------------------------------------------------------
struct PixelIn
{
    float4 position                  : SV_POSITION;
    noperspective float3 bary        : SV_BARYCENTRICS;
    nointerpolation uint triangleType : BLENDINDICES0;
};

// ---------------------------------------------------------------------------
// Mesh Shader
// Max 256 output vertices/primitives (D3D12 mesh shader limit).
// 128 threads handle up to 256 elements via a strided loop.
// ---------------------------------------------------------------------------
[NumThreads(128, 1, 1)]
[OutputTopology("triangle")]
void MSMain(
    uint gtid : SV_GroupThreadID,
    uint gid  : SV_GroupID,
    out indices uint3      outputTriangles[256],
    out vertices VertOut   outputVertices[256],
    out primitives PrimOut outputPrimAttr[256])
{
    // Each thread group renders one character
    CharRenderInfo charInfo = textBuffer[gid];
    GlyphletInfo   glInfo   = glyphletBuffer[charInfo.character];

    SetMeshOutputCounts(glInfo.vertexCount, glInfo.primitiveCount);

    // Output triangle indices and per-primitive attributes (strided loop for >128)
    for (uint ti = gtid; ti < glInfo.primitiveCount; ti += 128)
    {
        outputTriangles[ti]             = indexBuffer[glInfo.triangleBaseIndex + ti];
        outputPrimAttr[ti].triangleType = primAttrBuffer[glInfo.triangleBaseIndex + ti];
    }

    // Output vertices (transformed, strided loop for >128)
    for (uint vi = gtid; vi < glInfo.vertexCount; vi += 128)
    {
        float2 pos = vertexBuffer[glInfo.vertexBaseIndex + vi] + charInfo.pos;
        outputVertices[vi].position = mul(transformMatrix, float4(pos, 0.0f, 1.0f));
    }
}

// ---------------------------------------------------------------------------
// Pixel Shader — Anti-aliased Loop-Blinn curve rendering
// Uses screen-space derivatives of the curve function to produce smooth edges.
// ---------------------------------------------------------------------------
float4 PSMain(PixelIn p) : SV_TARGET
{
    uint t = p.triangleType;

    // Solid triangles: always fully filled
    // DEBUG: yellow tint to confirm new shader is active (normally white)
    if (t == SOLID)
        return float4(1.0f, 1.0f, 1.0f, 1.0f);

    // Map barycentrics to canonical quadratic Bezier curve coordinates
    // Canonical control points: a=(0,0), c=(0.5,0), b=(1,1)
    float u = p.bary.y * 0.5f + p.bary.z;
    float v = p.bary.z;

    // Evaluate curve function: f = u^2 - v  (f = 0 is the curve boundary)
    float f = u * u - v;

    // Compute screen-space gradient of f for approximate signed pixel distance
    float2 grad = float2(ddx(f), ddy(f));
    float  gradLen = length(grad);

    // Avoid division by zero for degenerate triangles
    if (gradLen < 1e-6f)
    {
        // Fallback: hard test when gradient is negligible
        bool inside = (t == CONVEX) ? (f <= 0.0f) : (f >= 0.0f);
        if (!inside) discard;
        return float4(1.0f, 1.0f, 1.0f, 1.0f);
    }

    // Signed distance in pixels from the curve boundary
    float dist = f / gradLen;

    // Smooth alpha with ~2px wide transition for visibly smooth edges
    float alpha;
    if (t == CONVEX)
        alpha = saturate(1.0f - dist);   // filled where f < 0, 2px fade zone
    else
        alpha = saturate(1.0f + dist);   // filled where f > 0, 2px fade zone

    // Discard fully transparent pixels
    if (alpha < 1.0f / 255.0f)
        discard;

    // DEBUG: curve triangles render green to verify shader is active
    // (SOLID triangles remain white, so green areas = curve AA zones)
    float3 col = float3(1.0f, 1.0f, 1.0f);  // bright green for curves
    return float4(col * alpha, alpha);
}
