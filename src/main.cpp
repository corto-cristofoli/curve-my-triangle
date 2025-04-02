// #include "geometrycentral/surface/base_geometry_interface.h"
// #include "geometrycentral/numerical/linear_algebra_types.h"
// #include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
// #include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <Eigen/src/SparseCore/SparseUtil.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

double alpha = 1;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

std::unique_ptr<ManifoldSurfaceMesh> meshOut;
std::unique_ptr<VertexPositionGeometry> geometryOut;

typedef std::map<std::string, Vector3> BezierTile;

void print_point(std::map<int, int> p) {
  std::cout << "{";
  for (auto c : p) {
    std::cout << "(" << c.first << ", " << c.second << ");";
  }
  std::cout << "}\n";
}

bool custome_compare(const std::pair<int, int> &p1,
                     const std::pair<int, int> &p2) {
  return p1.second > p2.second ||
         (p1.second == p2.second && p1.first > p2.first);
}

std::vector<std::pair<int, int>>
transform_bary_map(std::map<int, int> bary_map) {
  std::vector<std::pair<int, int>> out;
  for (auto coord : bary_map) {
    out.push_back(std::make_pair(coord.first, coord.second));
  }
  return out;
}

struct Vector3Cmp {
  bool operator()(const std::map<int, int> &lhs,
                  const std::map<int, int> &rhs) const {
    std::vector<std::pair<int, int>> lhs_v = transform_bary_map(lhs);
    std::sort(lhs_v.begin(), lhs_v.end(), custome_compare);
    std::vector<std::pair<int, int>> rhs_v = transform_bary_map(rhs);
    std::sort(rhs_v.begin(), rhs_v.end(), custome_compare);
    for (int i = 0; i < 3; i++) {
      if (lhs_v[i].second == 0 && rhs_v[i].second == 0) {
        return false;
      }
      if (lhs_v[i].first < rhs_v[i].first ||
          (lhs_v[i].first == rhs_v[i].first &&
           lhs_v[i].second < rhs_v[i].second)) {
        return true;
      }
      if (lhs_v[i].first != rhs_v[i].first ||
          lhs_v[i].second != rhs_v[i].second) {
        return false;
      }
    }
    return false;
  }
};

void print_vector_sorted(std::vector<std::pair<int, int>> v) {
  std::cout << "{";
  for (auto p : v) {
    std::cout << "(" << p.first << ";" << p.second << "),";
  }
  std::cout << "}\n";
}

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
Vector3 get_point_position(std::map<int, int> p, BezierTile B, int lod,
                           Face f) {
  std::vector<Vertex> v_f = get_vertices(f);

  double u = (double)p[v_f[1].getIndex()] / (lod + 1);
  double v = (double)p[v_f[2].getIndex()] / (lod + 1);
  double w = (double)p[v_f[0].getIndex()] / (lod + 1);
  return std::pow(w, 3) * B["300"] + std::pow(u, 3) * B["030"] +
         std::pow(v, 3) * B["003"] +
         3 * (w * w * u * B["210"] + w * u * u * B["120"] +
              w * w * v * B["201"] + w * v * v * B["102"] +
              u * u * v * B["021"] + u * v * v * B["012"]) +
         6 * u * v * w * B["111"];
}

Vector3 get_normal(std::map<int, int> p, BezierTile B, int lod, Face f) {
  std::vector<Vertex> v_f = get_vertices(f);

  double u = (double)p[v_f[1].getIndex()] / (lod + 1);
  double v = (double)p[v_f[2].getIndex()] / (lod + 1);
  double w = (double)p[v_f[0].getIndex()] / (lod + 1);
  Vector3 N = std::pow(w, 3) * B["n_200"] + std::pow(u, 3) * B["n_020"] +
              std::pow(v, 3) * B["n_002"] + u * w * B["n_110"] +
              u * v * B["n_011"] + v * w * B["n_101"];
  return N / dot(N, N);
}
class Obj {
public:
  std::map<std::map<int, int>, int, Vector3Cmp> index_vertices;
  int nb_points;
  std::vector<Vector3> vertices;
  std::vector<Vector3> normals;
  std::vector<std::array<int, 3>> faces;
  Obj() { nb_points = 0; }
  void add_vertex(std::map<int, int> v, BezierTile B, int lod, Face f) {
    auto it = index_vertices.find(v);
    if (it == index_vertices.end()) {
      index_vertices.insert({v, nb_points});
      vertices.push_back(get_point_position(v, B, lod, f));
      normals.push_back(get_normal(v, B, lod, f));
      nb_points++;
    }
  }
  int get_index_or_assign(std::map<int, int> v, BezierTile B, int lod, Face f) {
    auto it = index_vertices.find(v);
    if (it == index_vertices.end()) {
      // std::cout << "Adding :";
      // print_point(v);
      // auto vp = transform_bary_map(v);
      // std::sort(vp.begin(), vp.end(), custome_compare);
      // std::cout << "It looks like this sorted :";
      // print_vector_sorted(vp);
      add_vertex(v, B, lod, f);
      return nb_points - 1;
    } else {
      // std::cout << "This point is already in :";
      // print_point(v);
      // auto vp = transform_bary_map(v);
      // std::sort(vp.begin(), vp.end(), custome_compare);
      // std::cout << "It looks like this sorted :";
      // print_vector_sorted(vp);
      return it->second;
    }
  }
  void add_face(std::array<std::map<int, int>, 3> f, BezierTile B, int lod,
                Face face) {
    std::array<int, 3> new_f;
    for (int i = 0; i < 3; i++) {
      new_f[i] = get_index_or_assign(f[i], B, lod, face);
    }
    faces.push_back(new_f);
  }
  void write_to_file(std::ofstream &file) {
    for (Vector3 v : vertices) {
      file << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }
    for (Vector3 n : normals) {
      file << "vn " << n.x << " " << n.y << " " << n.z << "\n";
    }
    for (std::array<int, 3> f : faces) {
      file << "f " << f[0] << " " << f[1] << " " << f[2] << "\n";
    }
  }
};

// == ALGORITHM FUNCTIONS

// ## MAIN FUNCTIONS #########################

void add_normal_midedge(BezierTile B, Face face, int opposite_vertex) {
  std::vector<Vertex> P = get_vertices(face);
  double v = dot((geometry->vertexPositions[P[(opposite_vertex + 2) % 3]] -
                  geometry->vertexPositions[P[(opposite_vertex + 1) % 3]]),
                 (geometry->vertexNormals[P[(opposite_vertex + 1) % 3]] -
                  geometry->vertexNormals[P[(opposite_vertex + 2) % 3]])) /
             dot((geometry->vertexPositions[P[(opposite_vertex + 2) % 3]] -
                  geometry->vertexPositions[P[(opposite_vertex + 1) % 3]]),
                 (geometry->vertexPositions[P[(opposite_vertex + 2) % 3]] -
                  geometry->vertexPositions[P[(opposite_vertex + 1) % 3]]));
  Vector3 h = geometry->vertexNormals[P[(opposite_vertex + 1) % 3]] +
              geometry->vertexNormals[P[(opposite_vertex + 2) % 3]] +
              v * geometry->vertexPositions[P[(opposite_vertex + 2) % 3]] -
              geometry->vertexPositions[P[(opposite_vertex + 1) % 3]];
  if (opposite_vertex == 0) {
    B.insert({"n_011", h / dot(h, h)});
  } else if (opposite_vertex == 1) {
    B.insert({"n_101", h / dot(h, h)});
  } else {
    B.insert({"n_110", h / dot(h, h)});
  }
}

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
        B.insert({"n_200", geometry->vertexNormals[P[0]]});
        add_normal_midedge(B, face, 0);
      } else if (j == 3) {
        B.insert({"030", geometry->vertexPositions[P[1]]});
        B.insert({"n_020", geometry->vertexNormals[P[0]]});
        add_normal_midedge(B, face, 1);
      } else if (k == 3) {
        B.insert({"003", geometry->vertexPositions[P[2]]});
        B.insert({"n_002", geometry->vertexNormals[P[0]]});
        add_normal_midedge(B, face, 2);
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
      if (i != 0) {
        B.insert({std::to_string(i * 100 + j * 10 + k), b});
      } else {
        B.insert({"0" + std::to_string(j * 10 + k), b});
      }
    }
  }
  for (Vertex p : P) {
    V += geometry->vertexPositions[p];
  }
  V /= 3;
  E /= 6;
  B.insert({"111", (1 - alpha) * E + alpha * V});
  return B;
}

std::vector<std::map<int, int>> get_discretized_face(int lod, Face f) {
  // return local (barycentric position) (u,v) position
  // of new vertices of a face
  std::vector<Vertex> v_f = get_vertices(f);
  std::vector<std::map<int, int>> positions;
  for (int w = 0; w <= lod + 1; w++) {
    for (int t = 0; t <= w; t++) {
      std::map<int, int> coord_bary;
      coord_bary[v_f[2].getIndex()] = t;
      coord_bary[v_f[1].getIndex()] = w - t;
      coord_bary[v_f[0].getIndex()] = lod + 1 - w;
      positions.push_back(coord_bary);
    }
  }
  // for (int k=0; k<positions.size(); k++) {
  //   std::cout << positions[k] << std::endl;
  // }
  // std::cout << positions.size() << std::endl;
  return positions;
}

void set_pn_triangle(Face face, int lod, Obj &obj) {
  // transform a triangle face into a curved pn triangle
  // TODO: just inserting new point on the original face doesn't work
  BezierTile B = get_bezier_tile(face);
  std::vector<std::map<int, int>> discretized_face =
      get_discretized_face(lod, face);

  int start_layer, first_neighbor, second_neighbor;

  for (int l = 0; l <= lod; l++) {
    start_layer = (int)l * (l + 1) / 2;
    for (int i = start_layer; i < start_layer + l + 1; i++) {
      // Adding triangle with both neighbors in next layer
      first_neighbor = i + l + 1;
      second_neighbor = i + l + 2;
      obj.add_face({discretized_face[i], discretized_face[first_neighbor],
                    discretized_face[second_neighbor]},
                   B, lod, face);
      // Adding triangle (if possible) with second neighbor in next layer
      // and next neighbor in same layer (might not exist)
      // i + 1 := next neighbor in the same layer
      if (i + 1 < start_layer + l + 1) {
        obj.add_face({discretized_face[i], discretized_face[second_neighbor],
                      discretized_face[i + 1]},
                     B, lod, face);
      }
    }
  }

  // for (Vector2 uv : discretized_face) {
  //   Vector3 p = get_point_position(uv.x, uv.y, B);
  //   auto it = seen_points.find(p);
  //   if (it == seen_points.end()) {
  //     seen_points.insert(p);
  //     Vertex center_point = mesh->insertVertex(face);
  //     geometry->vertexPositions[center_point.getIndex()] = p;
  //   }
  // }
}

// == MAIN
int main(int argc, char **argv) {

  // Initialize polyscope
  polyscope::init();

  // Load input mesh
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(argv[1]);

  std::ofstream MyFile("output.obj");

  auto i_obj = polyscope::registerSurfaceMesh(
      "Input obj", geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  Obj obj;

  // // Create output mesh
  // meshOut = std::make_unique<ManifoldSurfaceMesh>();
  // geometryOut = std::make_unique<VertexPositionGeometry>( *meshOut );

  geometry->requireFaceAreas(); // get area by using geometry->faceAreas[f]
                                // where f is the face
  geometry
      ->requireFaceNormals(); // get normal by using geometry->faceNormals[f]
                              // where f is the face
  geometry->requireVertexNormals();

  i_obj->addVertexVectorQuantity("Input normals", geometry->vertexNormals);

  int lod = 2;

  for (auto face : mesh->faces()) {
    set_pn_triangle(face, lod, obj);
  }

  // Register the mesh with polyscope
  auto o_obj =
      polyscope::registerSurfaceMesh("Output obj", obj.vertices, obj.faces);

  o_obj->addVertexVectorQuantity("Normals", obj.normals);

  obj.write_to_file(MyFile);

  MyFile.close();
  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
