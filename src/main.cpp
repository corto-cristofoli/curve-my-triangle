// #include "geometrycentral/surface/base_geometry_interface.h"
// #include "geometrycentral/numerical/linear_algebra_types.h"
// #include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
// #include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <Eigen/src/SparseCore/SparseUtil.h>
// #include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>
// #include <vector>
// #include "polyscope/point_cloud.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

std::unique_ptr<ManifoldSurfaceMesh> meshOut;
std::unique_ptr<VertexPositionGeometry> geometryOut;

typedef std::map<std::string, Vector3> BezierTile;

class Obj {
public:
  std::map<Vector3, int> index_vertices;
  int nb_points;
  std::vector<Vector3> vertices;
  std::vector<std::array<int, 3>> faces;
  Obj() { nb_points = 0; }
  void add_vertex(Vector3 v) {
    auto it = index_vertices.find(v);
    if (it == index_vertices.end()) {
      index_vertices.insert({v, nb_points});
      vertices.push_back(v);
      nb_points++;
    }
  }
  int get_index(Vector3 v) {
    auto it = index_vertices.find(v);
    if (it == index_vertices.end()) {
      add_vertex(v);
      return nb_points - 1;
    } else {
      return it->second;
    }
  }
  void add_face(std::array<Vector3, 3> f) {
    std::array<int, 3> new_f;
    for (int i = 0; i < 3; i++) {
      new_f[i] = get_index(f[i]);
    }
    faces.push_back(new_f);
  }
};

// == ALGORITHM FUNCTIONS

int argmin_index(std::vector<Vertex> P) {
  int min_index = P[0].getIndex();
  int argmin = 0;
  int i = 0;
  for (Vertex v : P) {
    if (v.getIndex() < min_index) {
      min_index = v.getIndex();
      argmin = i;
    }
    i++;
  }
  return argmin;
}

std::vector<Vertex> shift_indices(std::vector<Vertex> P, int shifting) {
  std::vector<Vertex> new_P;
  int n_P = P.size();
  for (int i = 0; i < n_P; i++) {
    new_P.push_back(P[(i + shifting) % n_P]);
  }
  return new_P;
}

std::vector<Vertex> get_vertices(Face face) {
  // Method to get vertices in consistent order for a face
  std::vector<Vertex> P;
  for (Vertex v : face.adjacentVertices()) {
    P.push_back(v);
  }
  int argmin = argmin_index(P);
  P = shift_indices(P, argmin);
  return P;
}

// ## MAIN FUNCTIONS #########################

BezierTile get_bezier_tile(Face face) {
  // Get the positions of all the bezier control points of a face.
  // With those control points b, we can then apply

  BezierTile B;
  std::vector<Vertex> P = get_vertices(face);
  Vector3 V{0.0, 0.0, 0.0}; // barycenter of the face
  int max_id, min_id;
  Vector3 E{0.0, 0.0, 0.0};
  for (int i = 0; i <= 3; i++) {
    for (int j = 0; j <= 3 - i; j++) {
      int k = 3 - i - j;
      if (i == 3) {
        B.insert({"300", geometry->vertexPositions[P[0]]});
      } else if (j == 3) {
        B.insert({"030", geometry->vertexPositions[P[1]]});
      } else if (k == 3) {
        B.insert({"003", geometry->vertexPositions[P[2]]});
      }
      if (i == j) {
        continue;
      } else if (i == 1) {
        min_id = 0;
      } else if (j == 1) {
        min_id = 1;
      } else if (k == 1) {
        min_id = 2;
      }
      if (i == 2) {
        max_id = 0;
      } else if (j == 2) {
        max_id = 1;
      } else if (k == 2) {
        max_id = 2;
      }
      Vector3 Pmax = geometry->vertexPositions[P[max_id]];
      Vector3 Pmin = geometry->vertexPositions[P[min_id]];

      Vector3 N = geometry->vertexNormals[P[max_id]];
      double w = dot(Pmin - Pmax, N);
      Vector3 b = (2 * Pmax + Pmin - w * N) / 3; // see 3.1 section
      E += b;
      B.insert({std::to_string(i * 100 + j * 10 + k), b});
    }
  }
  for (Vertex p : P) {
    V += geometry->vertexPositions[p];
  }
  V /= 3;
  E /= 6;
  B.insert({"111", E + ((E - V) / 2)});
  return B;
}

Vector3 get_point_position(double u, double v, BezierTile B) {
  double w = 1 - u - v;
  return std::pow(w, 3) * B["300"] + std::pow(u, 3) * B["030"] +
         std::pow(v, 3) * B["003"] +
         3 * (w * w * u * B["210"] + w * u * u * B["120"] +
              w * w * v * B["201"] + w * v * v * B["102"] +
              u * u * v * B["021"] + u * v * v * B["012"]) +
         6 * u * v * w * B["111"];
}

std::vector<Vector2> get_discretized_face(int lod) {
  // return local (barycentric position) (u,v) position
  // of new vertices of a face
  double margin_left;
  std::vector<Vector2> positions;
  for (int w = 0; w <= lod + 1; w++) {
    for (int t = 0; t <= w; t++) {
      margin_left = w * (1. / (lod + 1));
      if (w == 0) {
        positions.push_back(Vector2{0, 0});
      } else {
        positions.push_back(Vector2{((float)t / (float)w) * margin_left,
                                    ((float)(1 - t) / (float)w) * margin_left});
      }
    }
  }
  // for (int k=0; k<positions.size(); k++) {
  //   std::cout << positions[k] << std::endl;
  // }
  // std::cout << positions.size() << std::endl;
  return positions;
}

void set_pn_triangle(Face face, std::unordered_set<Vector3> seen_points,
                     int lod, Obj obj) {
  // transform a triangle face into a curved pn triangle
  // TODO: just inserting new point on the original face doesn't work
  BezierTile B = get_bezier_tile(face);
  std::vector<Vector2> discretized_face = get_discretized_face(lod);

  for (int c = 0; c <= lod; c++) {
    obj.add_face({});
    for (int i = 1; i < c - 1; i++) {
    }
  }

  for (Vector2 uv : discretized_face) {
    Vector3 p = get_point_position(uv.x, uv.y, B);
    auto it = seen_points.find(p);
    if (it == seen_points.end()) {
      seen_points.insert(p);
      Vertex center_point = mesh->insertVertex(face);
      geometry->vertexPositions[center_point.getIndex()] = p;
    }
  }
}

// == MAIN
int main(int argc, char **argv) {

  // Initialize polyscope
  polyscope::init();

  // Load input mesh
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(argv[1]);

  // ofstream MyFile("output.obj");

  auto i_obj = polyscope::registerSurfaceMesh(
      "Input obj", geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  // // Create output mesh
  // meshOut = std::make_unique<ManifoldSurfaceMesh>();
  // geometryOut = std::make_unique<VertexPositionGeometry>( *meshOut );

  geometry->requireFaceAreas(); // get area by using geometry->faceAreas[f]
                                // where f is the face
  geometry
      ->requireFaceNormals(); // get normal by using geometry->faceNormals[f]
                              // where f is the face
  geometry->requireVertexNormals();

  std::unordered_set<Vector3> seen_points;
  for (auto face : mesh->faces()) {
    set_pn_triangle(face, seen_points, 5);
  }

  // Register the mesh with polyscope
  auto o_obj = polyscope::registerSurfaceMesh(
      "Output obj", geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  // MyFile.close();
  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
