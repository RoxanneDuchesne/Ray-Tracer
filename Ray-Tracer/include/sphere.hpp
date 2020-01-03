#pragma once

#include "vec3f.hpp"

namespace geometry {

class Sphere {
public:
  Sphere() = default;
  Sphere(math::Vec3f origin, float radius, math::Vec3f color, float reflect_coefficent);

  math::Vec3f origin;
  float radius = 1.f;
  math::Vec3f color = math::Vec3f(1.f, 1.f, 1.f);
  float reflect_coefficent = 0.f;

};

} // namespace geometry
