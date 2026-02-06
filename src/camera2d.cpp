#include "camera2d.h"

static float4x4 OrthoProjection(float left, float right, float bottom, float top, float nearZ, float farZ)
{
    float4x4 mat = {};
    mat.m[0][0] =  2.0f / (right - left);
    mat.m[1][1] =  2.0f / (top - bottom);
    mat.m[2][2] =  1.0f / (farZ - nearZ);
    mat.m[3][0] = -(right + left) / (right - left);
    mat.m[3][1] = -(top + bottom) / (top - bottom);
    mat.m[3][2] = -nearZ / (farZ - nearZ);
    mat.m[3][3] =  1.0f;
    return mat;
}

void Camera2D::Reset()
{
    m_panX = 0.0f;
    m_panY = 0.0f;
    m_zoom = 1.0f;
}

void Camera2D::SetZoomLimits(float minZoom, float maxZoom)
{
    m_minZoom = std::max(0.0001f, minZoom);
    m_maxZoom = std::max(m_minZoom, maxZoom);
    m_zoom = std::clamp(m_zoom, m_minZoom, m_maxZoom);
}

RectF Camera2D::GetViewRect(const RectF& baseView) const
{
    const float baseW = baseView.width();
    const float baseH = baseView.height();

    const float centerX = (baseView.left + baseView.right) * 0.5f + m_panX;
    const float centerY = (baseView.bottom + baseView.top) * 0.5f + m_panY;

    const float viewW = (m_zoom != 0.0f) ? (baseW / m_zoom) : baseW;
    const float viewH = (m_zoom != 0.0f) ? (baseH / m_zoom) : baseH;

    RectF view;
    view.left   = centerX - viewW * 0.5f;
    view.right  = centerX + viewW * 0.5f;
    view.bottom = centerY - viewH * 0.5f;
    view.top    = centerY + viewH * 0.5f;
    return view;
}

RectF Camera2D::GetViewRectForViewport(const RectF& baseView, float viewportW, float viewportH) const
{
    RectF view = GetViewRect(baseView);
    if (viewportW <= 0.0f || viewportH <= 0.0f)
        return view;

    const float viewW = view.width();
    const float viewH = view.height();
    if (viewW <= 0.0f || viewH <= 0.0f)
        return view;

    const float vpAspect = viewportW / viewportH;
    const float viewAspect = viewW / viewH;
    if (!std::isfinite(vpAspect) || !std::isfinite(viewAspect) || vpAspect <= 0.0f || viewAspect <= 0.0f)
        return view;

    // Expand the view rect to match the viewport aspect ratio (no non-uniform stretching).
    if (vpAspect > viewAspect) {
        const float newW = viewH * vpAspect;
        const float dx = (newW - viewW) * 0.5f;
        view.left -= dx;
        view.right += dx;
    } else if (vpAspect < viewAspect) {
        const float newH = viewW / vpAspect;
        const float dy = (newH - viewH) * 0.5f;
        view.bottom -= dy;
        view.top += dy;
    }
    return view;
}

float4x4 Camera2D::BuildOrthoMatrix(const RectF& baseView) const
{
    const RectF view = GetViewRect(baseView);
    return OrthoProjection(view.left, view.right, view.bottom, view.top, -1.0f, 1.0f);
}

float4x4 Camera2D::BuildOrthoMatrixForViewport(const RectF& baseView, float viewportW, float viewportH) const
{
    const RectF view = GetViewRectForViewport(baseView, viewportW, viewportH);
    return OrthoProjection(view.left, view.right, view.bottom, view.top, -1.0f, 1.0f);
}

void Camera2D::PanByPixels(float dxPixels, float dyPixels, float viewportW, float viewportH, const RectF& baseView)
{
    if (viewportW <= 0.0f || viewportH <= 0.0f)
        return;

    const RectF view = GetViewRectForViewport(baseView, viewportW, viewportH);
    const float viewW = view.width();
    const float viewH = view.height();

    // "Grab" behavior: dragging right moves content right => camera moves left.
    m_panX -= dxPixels * (viewW / viewportW);
    // Win32 mouse Y grows downward; world Y grows upward.
    m_panY += dyPixels * (viewH / viewportH);
}

void Camera2D::ZoomAtCursor(float wheelDelta, float mouseX, float mouseY, float viewportW, float viewportH, const RectF& baseView)
{
    if (wheelDelta == 0.0f || viewportW <= 0.0f || viewportH <= 0.0f)
        return;

    const float steps = wheelDelta / 120.0f; // WHEEL_DELTA
    const float zoomBase = 1.1f;
    const float zoomFactor = std::pow(zoomBase, steps);

    // World point under cursor before zoom.
    const RectF before = GetViewRectForViewport(baseView, viewportW, viewportH);
    const float worldXBefore = before.left + (mouseX / viewportW) * before.width();
    const float worldYBefore = before.top  - (mouseY / viewportH) * before.height();

    // Apply zoom (clamped).
    m_zoom = std::clamp(m_zoom * zoomFactor, m_minZoom, m_maxZoom);

    // World point under cursor after zoom (with same pan), then adjust pan so it matches before.
    const RectF after = GetViewRectForViewport(baseView, viewportW, viewportH);
    const float worldXAfter = after.left + (mouseX / viewportW) * after.width();
    const float worldYAfter = after.top  - (mouseY / viewportH) * after.height();

    m_panX += (worldXBefore - worldXAfter);
    m_panY += (worldYBefore - worldYAfter);
}

