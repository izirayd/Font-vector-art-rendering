#pragma once

#include <algorithm>
#include <cmath>

// Simple row-major 4x4 matrix used as a constant buffer payload.
// NOTE: The HLSL side uses float4x4 and multiplies as mul(M, v).
struct float4x4
{
    float m[4][4];
};

struct RectF
{
    float left   = 0.0f;
    float right  = 0.0f;
    float bottom = 0.0f;
    float top    = 0.0f;

    float width() const { return right - left; }
    float height() const { return top - bottom; }
};

class Camera2D
{
public:
    Camera2D() = default;

    void Reset();
    void SetZoomLimits(float minZoom, float maxZoom);

    float zoom() const { return m_zoom; }
    float panX() const { return m_panX; }
    float panY() const { return m_panY; }

    // Drag-pan: dx/dy are in pixels (mouse delta).
    void PanByPixels(float dxPixels, float dyPixels, float viewportW, float viewportH, const RectF& baseView);

    // Mouse-wheel zoom, keeping the world point under the cursor stable.
    // wheelDelta uses Win32 convention: +120 per notch.
    void ZoomAtCursor(float wheelDelta, float mouseX, float mouseY, float viewportW, float viewportH, const RectF& baseView);

    // Current view rect in world units (derived from baseView, pan and zoom).
    RectF GetViewRect(const RectF& baseView) const;

    // View rect expanded to match the viewport aspect ratio (prevents stretching).
    RectF GetViewRectForViewport(const RectF& baseView, float viewportW, float viewportH) const;

    // Build an orthographic matrix mapping the current view rect to NDC.
    float4x4 BuildOrthoMatrix(const RectF& baseView) const;

    // Build an orthographic matrix for the given viewport aspect ratio (prevents stretching).
    float4x4 BuildOrthoMatrixForViewport(const RectF& baseView, float viewportW, float viewportH) const;

private:
    float m_panX = 0.0f;
    float m_panY = 0.0f;
    float m_zoom = 1.0f;

    float m_minZoom = 0.05f;
    float m_maxZoom = 50.0f;
};

