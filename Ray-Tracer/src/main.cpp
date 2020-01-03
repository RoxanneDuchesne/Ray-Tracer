#include <iostream>
#include <vector>
#include <cmath>
#include <cassert> //assert
#include <random>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>

#include "mat3f.hpp"
#include "triangle.hpp"
#include "vec3f.hpp"
#include "vec2f.hpp"
#include "image.hpp"
#include "vec2i.hpp"
#include "grid2.hpp"
#include "timer.hpp"
#include "sphere.hpp"
#include "plane.hpp"
#include "ray.hpp"
#include "ray_intersect.hpp"

using namespace math;
using namespace geometry;
using namespace std;

int recursive_depth = 5;
math::Vec3f background_color = Vec3f(0.1f, 0.1f, 0.1f);

namespace raytracing {

Vec3f normalAt(Vec3f const &p, Sphere const &s) {
  Vec3f n = normalized((p - s.origin) / s.radius);
  return n;
}

Vec3f normalAt(Vec3f const & /*p*/, Plane const &p) {
  return p.normal;
}

Vec3f normalAt(Vec3f const &p, Triangle const &t) {
  return normal(t);
}

struct Surface {
  virtual ~Surface() = default;
  virtual Hit intersectSelf(Ray const &ray) const = 0;
  virtual Vec3f normalAtSelf(Vec3f const &p) const = 0;
  virtual Vec3f color() const = 0;
  virtual float reflect_coefficent() const = 0;
};

// helper class/function to make, e.g., class Sphere : public Surface
// Wrapping the geometry (e.g., sphere) in a class for intersection
// does not 'pollute' the geometry with inheritence
template <class T> struct Intersect_ : public Surface {
  template <typename... Args>
  Intersect_(Args... args)
      : m_self(std::forward<Args>(args)...) {}

  Hit intersectSelf(Ray const &ray) const { return intersect(ray, m_self); }
  Vec3f normalAtSelf(Vec3f const &p) const { return normalAt(p, m_self); }
  Vec3f color() const {return m_self.color;}
  float reflect_coefficent() const {return m_self.reflect_coefficent;}

  T m_self;
};

template <typename T> std::unique_ptr<Intersect_<T>> makeIntersectable(T t) {
  return std::unique_ptr<Intersect_<T>>(new Intersect_<T>(t));
}

struct ImagePlane {
  using Screen = geometry::Grid2<raster::RGB>;

  Screen screen;
  math::Vec3f origin;
  math::Vec3f u;
  math::Vec3f v;
  float left;   // right symmetric
  float bottom; // top symmetric

  ImagePlane &resolution(uint32_t width, uint32_t height) {
    screen = Screen(width, height);
    return *this;
  }
  ImagePlane &center(Vec3f center) {
    origin = center;
    return *this;
  }
  ImagePlane &uvAxes(Vec3f up, Vec3f right) {
    u = up;
    v = right;
    return *this;
  }
  ImagePlane &dimensions(float width, float height) {
    left = -(0.5f * width);
    bottom = -(0.5f * height);
    return *this;
  }

  std::vector<math::Vec3f> pixelTo3D(math::Vec2f pixel, int anti_aliasing) const {
    using std::abs;

    // rays that need to be casted for super-sampling
    std::vector<math::Vec3f> rays_to_cast;

    float increment = 1.f/anti_aliasing;
    math::Vec2f next_pixel;
    float u_x;
    float v_y;

    for (int i = 0; i < anti_aliasing; i ++){
        for (int j = 0; j < anti_aliasing; j ++){

			float k_i = ((float)rand() / (RAND_MAX)) + 1;
			float k_j = ((float)rand() / (RAND_MAX)) + 1;

            next_pixel = {pixel.x + k_i * i * increment, pixel.y + k_j * j * increment};
            u_x = left + (2.f * abs(left)) * (next_pixel.x) / screen.width();
            v_y = bottom + (2.f * abs(bottom)) * (next_pixel.y) / screen.height();

            rays_to_cast.push_back(origin + u_x * u + v_y * v);
        }
    }
    return rays_to_cast;
  }
};

ImagePlane makeImagePlane(Vec3f const &eye, Vec3f const &lookatPosition,
                          Vec3f const &canonicalUp, int xResolution,
                          int yResolution, float width, float height,
                          float nearPlaneDistanace) {
  // make orthonormal basis around eye
  auto gaze = normalized(lookatPosition - eye);

  auto u = gaze ^ canonicalUp;
  u.normalize();

  auto v = u ^ gaze;
  v.normalize();

  // image plane is orthogoanl to camera gaze so use u,v for image plane
  ImagePlane imagePlane;
  // using method chaining to have named parameter/configuation)
  imagePlane.resolution(xResolution, yResolution)
      .center(eye + gaze * nearPlaneDistanace)
      .dimensions(width, height)
      .uvAxes(u, v);

  return imagePlane;
}

using s_ptr = std::unique_ptr<Surface>;

Vec3f castRay(Ray ray,
              math::Vec3f eye,                    //
              math::Vec3f light,                  //
              std::vector<s_ptr> const &surfaces,
              float recursive_depth) {

  // background color
  Vec3f colorOut = background_color;

  // find closed object, if any
  Hit closest;
  // pointer to closest object
  Surface const *surface = nullptr;
  for (auto const &s : surfaces) {
    auto hit = s->intersectSelf(ray);
    if (hit && (hit.rayDepth < closest.rayDepth) && (hit.rayDepth > 0.f)) {
      closest = hit;
      surface = s.get();
    }
  }

  // if hit get point
  if (surface != nullptr) {

      float reflect_coefficent = surface->reflect_coefficent();

      Vec3f light_color = Vec3f(1.f, 1.f, 1.f);
      Vec3f object_color = surface->color();
      Vec3f color = object_color;

      float k_d = 0.5f;
      float k_s = 0.3f;
      float k_a = 0.1f;

      float I = 1.f;

      Vec3f position = ray.origin + (closest.rayDepth * ray.direction);
      Vec3f normal = surface->normalAtSelf(position);
      Vec3f light_direction = normalized(light - position);

      Vec3f ambient = k_a * color;

      Vec3f diffuse = ((max((normal * light_direction), 0.f)) * k_d) * color;

      Vec3f view_direction = normalized(position - (eye));
      Vec3f reflect_direction = (light_direction) + (2 * ((-light_direction * normal) * normal));
      Vec3f specular = (pow(max((reflect_direction * view_direction), 0.f), 32.f) * k_s) * color;

      colorOut = (diffuse + specular + ambient);

      Ray shadow_ray;
      shadow_ray.direction = normalized(light - position);
      shadow_ray.origin = position + (shadow_ray.direction * 0.00001f);
      Surface const *surface = nullptr;
      for (auto const &s : surfaces) {
        auto hit = s->intersectSelf(shadow_ray);
        if (hit && (hit.rayDepth < 1e+5) && (hit.rayDepth > 0)) {
          colorOut = ambient;
          break;
        }
      }

      if (recursive_depth > 0){
          Ray reflect_ray;
          //reflect based on eye
          reflect_ray.direction = -normalized((eye) - (2 * ((eye * normal) * normal)));
          reflect_ray.origin = position + (reflect_ray.direction * 0.001f);
          colorOut = colorOut + reflect_coefficent * castRay(reflect_ray, eye, light, surfaces, recursive_depth - 1);
      }
  }
  return colorOut;
}

void render(ImagePlane &imagePlane,
            math::Vec3f eye,
            math::Vec3f light,
            std::vector<s_ptr> const &surfaces,
            int light_area) {

  // Standard mersenne_twister_engine seeded
  thread_local std::mt19937 gen(0);
  auto sampleRange = [](float a, float b) {
    using distrubution = std::uniform_real_distribution<>;
    return distrubution(a, b)(gen);
  };

  std::vector<math::Vec3f> light_points;

  if (light_area == -1){
      light_points.push_back(light);
  }
  else{
     // calculate points for area_lighting
     for (int i = -light_area; i < light_area; i++){
         for (int j = -light_area; j < light_area; j++){
             light_points.push_back(math::Vec3f((i * 0.1f + light.x), light.y, (j * 0.1f + light.z)));
         }
     }
  }

  for (int32_t y = 0; y < imagePlane.screen.height(); ++y) {
    for (int32_t x = 0; x < imagePlane.screen.width(); ++x) {

     math::Vec3f intermediate_color = Vec3f(0.f, 0.f, 0.f);
	  
     math::Vec2f pixel(x, y);
     // this is where the super_sampling grid size is set
     auto pixel_points = imagePlane.pixelTo3D(pixel, 4);

     // looping through pixel_points for anti-aliasing
     for (int i = 0; i < pixel_points.size(); i++){

        math::Vec3f light_color = Vec3f(0.f, 0.f, 0.f);

        auto direction = normalized(pixel_points[i] - eye);
        auto bias = 1e-4f;
        auto p = pointOnLne(eye, direction, bias);
        Ray r(p, direction);

        // looping through light_points for soft shadows/ area lighting
        for (int k = 0; k < light_points.size(); k++){
            light = light_points[k];
            light_color += castRay(r, eye, light, surfaces, 5);
        }
        intermediate_color += light_color * (1.f / light_points.size());
     }

     auto colorOut = intermediate_color * (1.f / pixel_points.size());

      // correct to quantiezed error
      // (i.e., removes banded aliasing when converting to 8bit RGB)
      constexpr float halfStep = 1.f / 512;
      colorOut = raster::quantizedErrorCorrection(
          colorOut, sampleRange(-halfStep, halfStep));

      imagePlane.screen({x, y}) = raster::convertToRGB(colorOut);
    }
  }
}
} // namespace

int main() {

    using namespace raytracing;
    std::vector<s_ptr> surfaces;

	Vec3f light{ 15, 10, 10 };
    int light_area = 1;
	Vec3f eye{ 0.f, 3.f, 8.f };
	Vec3f lookat{ 0.f, 3.f, 0.f };
	Vec3f up{ 0.f, 1.f, 0.f };
	int resolutionX = 1000;
	int resolutionY = 1000;
	float planeWidth = 50.f;
	float planeHeight = 50.f;
	float focalDist = 50.f;
    background_color = math::Vec3f(0.f, 0.f, 0.f);

  fstream file_in("../RayTracer-RoxanneDuchesne/Spheres.txt");
  std::string line;

  try {
    if (file_in.is_open()) {
         while (getline(file_in, line)) {

            std::vector<std::string> split_line;
            std::istringstream ss(line);
            for (std::string s; ss >> s;) {
                  split_line.push_back(s);
              }

              if (split_line[0].compare("Light:") == 0) {
                  light = Vec3f(std::stof(split_line[1]), std::stof(split_line[2]), std::stof(split_line[3]));
              }
              else if (split_line[0].compare("LightArea:") == 0){
                  light_area = std::stoi(split_line[1]);
              }
              else if (split_line[0].compare("Eye:") == 0){
                  eye = Vec3f(std::stof(split_line[1]), std::stof(split_line[2]), std::stof(split_line[3]));
              }
              else if (split_line[0].compare("Lookat:") == 0){
                  lookat = Vec3f(std::stof(split_line[1]), std::stof(split_line[2]), std::stof(split_line[3]));
              }
              else if (split_line[0].compare("Up:") == 0){
                  up = Vec3f(std::stof(split_line[1]), std::stof(split_line[2]), std::stof(split_line[3]));
              }
              else if (split_line[0].compare("ResolutionX:") == 0){
                  resolutionX = std::stof(split_line[1]);
              }
              else if (split_line[0].compare("ResolutionY:") == 0){
                  resolutionY = std::stof(split_line[1]);
              }
              else if (split_line[0].compare("PlaneWidth:") == 0){
                  planeWidth = std::stof(split_line[1]);
              }
              else if (split_line[0].compare("PlaneHeight:") == 0){
                  planeHeight = std::stof(split_line[1]);
              }
              else if (split_line[0].compare("FocalDistance:") == 0){
                  focalDist = std::stof(split_line[1]);
              }
              else if (split_line[0].compare("BackGroundColor:") == 0){
                  background_color = Vec3f(std::stof(split_line[1]), std::stof(split_line[2]), std::stof(split_line[3]));
              }
              else if (split_line[0].compare("Sphere:") == 0){
                  Sphere s;
                  s.radius = std::stof(split_line[2]);
                  s.origin = Vec3f(std::stof(split_line[4]), std::stof(split_line[5]), std::stof(split_line[6]));
                  s.color = Vec3f(std::stof(split_line[8]), std::stof(split_line[9]), std::stof(split_line[10]));
                  s.reflect_coefficent = std::stof(split_line[12]);
                  surfaces.push_back(makeIntersectable(s));
              }
              else if (split_line[0].compare("Triangle:") == 0){
                  Triangle t;
                  t.a() = Vec3f(std::stof(split_line[2]), std::stof(split_line[3]), std::stof(split_line[4]));
                  t.b() = Vec3f(std::stof(split_line[6]), std::stof(split_line[7]), std::stof(split_line[8]));
                  t.c() = Vec3f(std::stof(split_line[10]), std::stof(split_line[11]), std::stof(split_line[12]));
                  t.color = Vec3f(std::stof(split_line[14]), std::stof(split_line[15]), std::stof(split_line[16]));
                  t.reflect_coefficent = std::stof(split_line[18]);
                  surfaces.push_back(makeIntersectable(t));
              }
              else if (split_line[0].compare("Plane:") == 0){
                  Plane p;
                  p.origin = Vec3f(std::stof(split_line[2]), std::stof(split_line[3]), std::stof(split_line[4]));
                  p.normal = Vec3f(std::stof(split_line[6]), std::stof(split_line[7]), std::stof(split_line[8]));
                  p.color = Vec3f(std::stof(split_line[10]), std::stof(split_line[11]), std::stof(split_line[12]));
                  p.reflect_coefficent = std::stof(split_line[14]);
                  surfaces.push_back(makeIntersectable(p));
              }
          }
     }
     else {
          cout << "Input.txt was not able to be opened" << endl;
          //exit(1);
    }
  }
  catch(...) {
      cout << "File Wasn't Formatted Correctly" << endl;
      //exit(-1);
  }



  auto imagePlane = makeImagePlane(eye, lookat, up,          //
                                   resolutionX, resolutionY, //
                                   planeWidth, planeHeight,  //
                                   focalDist);

  //Sphere Scene Without Input File
  
  Plane p;
  Sphere front_s;
  Sphere s1;
  Sphere s2;
  Sphere s3;
  Sphere s4;
  Sphere s5;

  front_s.radius = 3.f;
  front_s.origin = Vec3f(0.f, 3.f, 0.f);
  front_s.reflect_coefficent = 0.5f;
  front_s.color = Vec3f(0.1f, 1.f, 0.5f);

  s1.color = Vec3f(0.1f, 0.6f, 1.f);
  s2.color = Vec3f(0.2f, 0.9f, 0.6f);
  s3.color = Vec3f(0.9f, 0.9f, 0.1f);
  s4.color = Vec3f(0.9f, 0.f, 0.5f);
  s5.color = Vec3f(0.9f, 0.6f, 0.4f);

  s1.radius = 1.3f;
  s1.origin = Vec3f(4.f, 1.3f, 5.f);
  s2.radius = 3.f;
  s2.origin = Vec3f(-14.f, 3.f, 4.f);
  s3.radius = 1.3f;
  s3.origin = Vec3f(0.9f, 1.3f, 8.f);
  s4.radius = 5.f;
  s4.origin = Vec3f(-7.f, 5.f, 17.f);
  s5.radius = 1.5f;
  s5.origin = Vec3f(18.f, 1.5f, 9.f);

  p.origin = Vec3f(0.f, 0.f, 0.f);
  p.normal = Vec3f(0.f, 1.f, 0.f);
  p.color = Vec3f(0.1f, 1.f, 1.f);
  surfaces.push_back(makeIntersectable(p));

  surfaces.push_back(makeIntersectable(front_s));
  surfaces.push_back(makeIntersectable(s1));
  surfaces.push_back(makeIntersectable(s2));
  surfaces.push_back(makeIntersectable(s3));
  surfaces.push_back(makeIntersectable(s4));
  surfaces.push_back(makeIntersectable(s5));
  
  int size = 10;

  for (int i = -size; i < size; i++) {
	  for (int j = -size; j < size; j++) {
		  if (i % 2 == 0 && j % 2 != 0) {
			  Triangle t;
			  t.a() = Vec3f(i, 0.01f, j);
			  t.b() = Vec3f(i + 5, 0.01f, j);
			  t.c() = Vec3f(i, 0.01f, j + 5);
			  Triangle z;
			  z.a() = Vec3f(i + 5, 0.01f, j);
			  z.b() = Vec3f(i + 5, 0.0f, j + 5);
			  z.c() = Vec3f(i, 0.01f, j + 5);

			  surfaces.push_back(makeIntersectable(t));
			  surfaces.push_back(makeIntersectable(z));
		  }
	  }
  }
 
  temporal::Timer timer(true);

  render(imagePlane, eye, light, surfaces, light_area);

  std::cout << "Time elapsed: " << std::fixed << timer.minutes() << " min.\n";

  raster::write_screen_to_file("./test.png", imagePlane.screen);

  return EXIT_SUCCESS;
}
