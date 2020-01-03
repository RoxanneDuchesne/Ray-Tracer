#include "ray_intersect.hpp"
#include "ray_intersect.hpp"
#include "mat3f.hpp"
#include "mat4f.hpp"

#include <cmath>
#include <algorithm>
#include <list>

namespace geometry {

Hit::operator bool() const { return didIntersect; }

Hit intersect(Ray const &r, Triangle const &t) {
    Hit hit;
    hit.didIntersect = false;

    math::Mat3f A = math::Mat3f{(t.a().x - t.b().x), (t.a().x - t.c().x), r.direction.x,
                                (t.a().y - t.b().y), (t.a().y - t.c().y), r.direction.y,
                                (t.a().z - t.b().z), (t.a().z - t.c().z), r.direction.z};

    math::Mat3f beta_top = math::Mat3f{(t.a().x - r.origin.x), (t.a().x - t.c().x), r.direction.x,
                                       (t.a().y - r.origin.y), (t.a().y - t.c().y), r.direction.y,
                                       (t.a().z - r.origin.z), (t.a().z - t.c().z), r.direction.z};

    math::Mat3f gamma_top = math::Mat3f{(t.a().x - t.b().x), (t.a().x - r.origin.x), r.direction.x,
                                        (t.a().y - t.b().y), (t.a().y - r.origin.y), r.direction.y,
                                        (t.a().z - t.b().z), (t.a().z - r.origin.z), r.direction.z};

    math::Mat3f t_top = math::Mat3f{(t.a().x - t.b().x), (t.a().x - t.c().x), (t.a().x - r.origin.x),
                                    (t.a().y - t.b().y), (t.a().y - t.c().y), (t.a().y - r.origin.y),
                                    (t.a().z - t.b().z), (t.a().z - t.c().z), (t.a().z - r.origin.z)};

    float beta = math::determinant(beta_top)/math::determinant(A);
    float gamma = math::determinant(gamma_top)/math::determinant(A);
    float tt = math::determinant(t_top)/math::determinant(A);

    float max = 10000.f;

    if (tt < 0 || tt > max){
        return hit;
    }

    if (beta < 0 || (beta > (1.f-gamma))){
        return hit;
    }
    if (gamma < 0){
        return hit;
    }
    if ((beta + gamma) > 1){
        return hit;
    }

    hit.didIntersect = true;
    hit.rayDepth = tt;

    return hit;
}

Hit intersect(Ray const &r, Sphere const &s) {

  Hit hit;
  hit.didIntersect = false;

  float a = (r.direction * r.direction);
  float b = 2 * (r.direction * (r.origin - s.origin));
  float c = ((r.origin - s.origin) * (r.origin - s.origin)) - (s.radius * s.radius);

  float t0;
  float t1;

  float discriminate = (b * b) - (4 * a * c);
  if (discriminate < 0){
      return hit;
  }
  else if (discriminate == 0){
     t0 = -(0.5f * b) / a;
     hit.rayDepth = t0;
     hit.didIntersect = true;
  }
  else {
      if (b > 0){
          t0 = (-0.5 * (b + sqrtf(discriminate))) / a;
          t1 = c / (-0.5 * (b + sqrtf(discriminate)));
      }
      else {
          t0 = (-0.5 * (b - sqrtf(discriminate))) / a;
          t1 = c / (-0.5 * (b - sqrtf(discriminate)));
      }
      if (t0 < t1){
          hit.rayDepth = t0;
      }
      else {
          hit.rayDepth = t1;
      }
      hit.didIntersect = true;

  }
  return hit;
}

Hit intersect(Ray const &r, Plane const &p) {
  Hit hit;

  auto denom = r.direction * p.normal;
  if (std::abs(denom) < 1e-5f)
    return hit;

  auto t = ((p.origin - r.origin) * p.normal) / denom;

  if (t < 0.f)
    return hit;

  hit.didIntersect = true;
  hit.rayDepth = t;

  return hit;
}

} // namespace geometry
