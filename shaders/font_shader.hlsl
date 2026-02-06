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
// ---------------------------------------------------------------------------
[NumThreads(128, 1, 1)]
[OutputTopology("triangle")]
void MSMain(
    uint gtid : SV_GroupThreadID,
    uint gid  : SV_GroupID,
    out indices uint3      outputTriangles[128],
    out vertices VertOut   outputVertices[128],
    out primitives PrimOut outputPrimAttr[128])
{
    // Each thread group renders one character
    CharRenderInfo charInfo = textBuffer[gid];
    GlyphletInfo   glInfo   = glyphletBuffer[charInfo.character];

    SetMeshOutputCounts(glInfo.vertexCount, glInfo.primitiveCount);

    // Output triangle indices and per-primitive attributes
    if (gtid < glInfo.primitiveCount)
    {
        outputTriangles[gtid]             = indexBuffer[glInfo.triangleBaseIndex + gtid];
        outputPrimAttr[gtid].triangleType = primAttrBuffer[glInfo.triangleBaseIndex + gtid];
    }

    // Output vertices (transformed)
    if (gtid < glInfo.vertexCount)
    {
        float2 pos = vertexBuffer[glInfo.vertexBaseIndex + gtid] + charInfo.pos;
        outputVertices[gtid].position = mul(transformMatrix, float4(pos, 0.0f, 1.0f));
    }
}

// ---------------------------------------------------------------------------
// Pixel Shader
// ---------------------------------------------------------------------------
float4 PSMain(PixelIn p) : SV_TARGET
{
    uint t = p.triangleType;

    // Solid triangles: always filled
    if (t == SOLID)
        return float4(1.0f, 1.0f, 1.0f, 1.0f);

    // Map barycentrics to canonical quadratic Bezier curve coordinates
    // Canonical control points: a=(0,0), c=(0.5,0), b=(1,1)
    // u = bary.x*0 + bary.y*0.5 + bary.z*1
    // v = bary.x*0 + bary.y*0.0 + bary.z*1
    float u = p.bary.y * 0.5f + p.bary.z * 1.0f;
    float v = p.bary.z * 1.0f;

    // Evaluate: u^2 - v = 0 is the curve
    float y = u * u - v;

    // Convex: discard pixels outside the curve (y > 0)
    // Concave: discard pixels inside the curve (y < 0)
    if ((t == CONVEX  && y > 0.0f) ||
        (t == CONCAVE && y < 0.0f))
    {
        discard;
    }

    return float4(1.0f, 1.0f, 1.0f, 1.0f);
}
