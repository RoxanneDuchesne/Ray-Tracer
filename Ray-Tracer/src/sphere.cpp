#include "sphere.hpp"

namespace geometry {
Sphere::Sphere(math::Vec3f origin, float radius, math::Vec3f color, float reflect_coefficent)
    : origin(origin), radius(radius), color(color), reflect_coefficent(reflect_coefficent){}
}
