#pragma once
#include <cmath>
#include <cstddef>
#include <vector>

namespace es {

using Real = double;
using Index = std::ptrdiff_t;

/// Point / vector in the meridian half-plane (r >= 0, z).
struct Vec2 {
  Real r{0.0};
  Real z{0.0};

  constexpr Vec2() = default;
  constexpr Vec2(Real r_, Real z_) : r(r_), z(z_) {}

  friend constexpr Vec2 operator+(Vec2 a, Vec2 b) { return {a.r + b.r, a.z + b.z}; }
  friend constexpr Vec2 operator-(Vec2 a, Vec2 b) { return {a.r - b.r, a.z - b.z}; }
  friend constexpr Vec2 operator*(Real s, Vec2 a) { return {s * a.r, s * a.z}; }
  friend constexpr Vec2 operator*(Vec2 a, Real s) { return {s * a.r, s * a.z}; }
  friend constexpr Vec2 operator/(Vec2 a, Real s) { return {a.r / s, a.z / s}; }
  Vec2& operator+=(Vec2 b) { r += b.r; z += b.z; return *this; }
  Vec2& operator-=(Vec2 b) { r -= b.r; z -= b.z; return *this; }
  Vec2& operator*=(Real s) { r *= s; z *= s; return *this; }
};

inline Real dot(Vec2 a, Vec2 b) { return a.r * b.r + a.z * b.z; }
inline Real norm(Vec2 a) { return std::sqrt(dot(a, a)); }
inline Real norm2(Vec2 a) { return dot(a, a); }
inline Vec2 normalized(Vec2 a) { const Real n = norm(a); return n > 0 ? a / n : Vec2{0, 0}; }
/// Rotate by -90 deg: (r,z) -> (z,-r).  Used to build outward normals.
inline Vec2 perp(Vec2 t) { return {t.z, -t.r}; }

template <class T> constexpr T sqr(T x) { return x * x; }

}  // namespace es
