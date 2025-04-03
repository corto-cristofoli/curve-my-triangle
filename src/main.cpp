#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/vector3.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <Eigen/src/SparseCore/SparseUtil.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

double alpha = 1;

// == Geometry-central data

using namespace geometrycentral;
using namespace geometrycentral::surface;

std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// == Typedef

typedef std::map<std::string, Vector3> BezierTile;

// == Comparators

bool custome_compare_pair(const std::pair<int, int> &p1,
                          const std::pair<int, int> &p2) {
  // Comparator described in equation (1) in the report
  return p1.second > p2.second ||
         (p1.second == p2.second && p1.first > p2.first);
}

struct BaryCmp {
  // Comparator depicted in Figure 5 in the report
  bool operator()(const std::vector<std::pair<int, int>> &lhs,
                  const std::vector<std::pair<int, int>> &rhs) const {
    std::vector<std::pair<int, int>> lhs_sorted(3);
    std::vector<std::pair<int, int>> rhs_sorted(3);
    // Doing a copy here because the order of the point in the representation is
    // useful to compute the projection on the Bézier tile
    std::partial_sort_copy(lhs.begin(), lhs.end(), lhs_sorted.begin(),
                           lhs_sorted.end(), custome_compare_pair);
    std::partial_sort_copy(rhs.begin(), rhs.end(), rhs_sorted.begin(),
                           rhs_sorted.end(), custome_compare_pair);
    for (int i = 0; i < 3; i++) {
      if (lhs_sorted[i].second == 0 && rhs_sorted[i].second == 0) {
        // If this case happens it means that both representation have same
        // non-negative coordinates, it means they are equal so we return false
        return false;
      }
      if (lhs_sorted[i].first < rhs_sorted[i].first ||
          (lhs_sorted[i].first == rhs_sorted[i].first &&
           lhs_sorted[i].second < rhs_sorted[i].second)) {
        // In this case it's the classical lexicographic order
        return true;
      }
      if (lhs_sorted[i].first != rhs_sorted[i].first ||
          lhs_sorted[i].second != rhs_sorted[i].second) {
        // In this case again it's the classical lexicographic order
        return false;
      }
    }
    // If we arrive at this point it means that they are equal
    return false;
  }
};

// == Utilities

void print_bary_point(std::vector<std::pair<int, int>> v) {
  std::cout << "{";
  for (auto p : v) {
    std::cout << "(" << p.first << ";" << p.second << "),";
  }
  std::cout << "}\n";
}

// == Tool to get faces' vertices in a consistent order

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
  std::vector<Vertex> P;
  for (Vertex v : face.adjacentVertices()) {
    P.push_back(v);
  }
  int argmin = argmin_index(P);
  P = shift_indices(P, argmin);
  return P;
}

// == Implement the formulae for calculating the projection on the Bézier tile

Vector3 get_point_position(std::vector<std::pair<int, int>> p, BezierTile B,
                           int lod) {
  // Implement the formula for b at page 2, §2
  double u = (double)p[1].second / (lod + 1);
  double v = (double)p[2].second / (lod + 1);
  double w = (double)p[0].second / (lod + 1);
  return std::pow(w, 3) * B["300"] + std::pow(u, 3) * B["030"] +
         std::pow(v, 3) * B["003"] +
         3 * (w * w * u * B["210"] + w * u * u * B["120"] +
              w * w * v * B["201"] + w * v * v * B["102"] +
              u * u * v * B["021"] + u * v * v * B["012"]) +
         6 * u * v * w * B["111"];
}

Vector3 get_normal(std::vector<std::pair<int, int>> p, BezierTile B, int lod) {
  // Implement the formula for n at page 2, §2
  double u = (double)p[1].second / (lod + 1);
  double v = (double)p[2].second / (lod + 1);
  double w = (double)p[0].second / (lod + 1);
  Vector3 N = std::pow(w, 3) * B["n_200"] + std::pow(u, 3) * B["n_020"] +
              std::pow(v, 3) * B["n_002"] + u * w * B["n_110"] +
              u * v * B["n_011"] + v * w * B["n_101"];
  return N / dot(N, N);
}

// == Obj class

class Obj {
public:
  std::map<std::vector<std::pair<int, int>>, int, BaryCmp> index_vertices;
  int nb_points;
  std::vector<Vector3> vertices;
  std::vector<Vector3> normals;
  std::vector<std::array<int, 3>> faces;
  Obj() { nb_points = 0; }

  void add_vertex(std::vector<std::pair<int, int>> v, BezierTile B, int lod) {
    // Only inserting if v is not already in vertices
    auto it = index_vertices.find(v);
    if (it == index_vertices.end()) {
      index_vertices.insert({v, nb_points});
      vertices.push_back(get_point_position(v, B, lod));
      normals.push_back(get_normal(v, B, lod));
      nb_points++;
    }
  }

  int get_index_or_assign(std::vector<std::pair<int, int>> v, BezierTile B,
                          int lod) {
    // Either return the index of v if v is already in vertices, either add v to
    // vertices the returns it index
    auto it = index_vertices.find(v);
    if (it == index_vertices.end()) {
      add_vertex(v, B, lod);
      return nb_points - 1;
    } else {
      return it->second;
    }
  }

  void add_face(std::array<std::vector<std::pair<int, int>>, 3> f, BezierTile B,
                int lod) {
    // Add a face to faces
    std::array<int, 3> new_f;
    for (int i = 0; i < 3; i++) {
      new_f[i] = get_index_or_assign(f[i], B, lod);
    }
    faces.push_back(new_f);
  }

  void write_to_file(std::ofstream &file) {
    // Export the Obj in obj file

    for (Vector3 v : vertices) {
      file << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }
    for (Vector3 n : normals) {
      file << "vn " << n.x << " " << n.y << " " << n.z << "\n";
    }
    for (std::array<int, 3> f : faces) {
      file << "f " << f[0] << "//" << f[0] << " " << f[1] << "//" << f[1] << " "
           << f[2] << "//" << f[2] << "\n";
    }
  }
};

// == ALGORITHM FUNCTIONS

// ## MAIN FUNCTIONS #########################

// == Auxiliary function to help compute Bézier tile

void add_normal_midedge(BezierTile B, Face face, int opposite_vertex) {
  // Here we are computing n_ijk where i,j,k < 2
  // opposite_vertex indicates which of i,j,k is 0
  // The computation are retrieved from page 3, §3.3
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
  Vector3 E{0.0, 0.0,
            0.0}; // barycenter of all control points except the center one
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

std::vector<std::vector<std::pair<int, int>>> get_discretized_face(int lod,
                                                                   Face f) {
  // Function to discretize a triangular face, the indexing is explained in
  // section 4 of the report

  std::vector<Vertex> v_f = get_vertices(f);
  std::vector<std::vector<std::pair<int, int>>> positions;
  for (int layer = 0; layer <= lod + 1; layer++) {
    for (int t = 0; t <= layer;
         t++) { // Indeed there is l+1 points points on the layer l
      std::vector<std::pair<int, int>> coord_bary(3);

      coord_bary[0] = std::make_pair(v_f[0].getIndex(), lod + 1 - layer);
      coord_bary[1] = std::make_pair(v_f[1].getIndex(), layer - t);
      coord_bary[2] = std::make_pair(v_f[2].getIndex(), t);
      positions.push_back(coord_bary);
    }
  }
  return positions;
}

void set_pn_triangle(Face face, int lod, Obj &obj) {
  // Add the PN triangle of face to obj
  BezierTile B = get_bezier_tile(face);

  // Retrieving the discretized points of the face according to the level of
  // detail
  std::vector<std::vector<std::pair<int, int>>> discretized_face =
      get_discretized_face(lod, face);

  int start_layer, first_neighbor, second_neighbor;

  // Adding the faces, again this is well explained in section 4 of the report

  for (int layer = 0; layer <= lod; layer++) {
    start_layer = (int)layer * (layer + 1) / 2;
    for (int i = start_layer; i < start_layer + layer + 1; i++) {
      // Adding triangle with both neighbors in next layer
      first_neighbor = i + layer + 1;
      second_neighbor = i + layer + 2;
      obj.add_face({discretized_face[i], discretized_face[first_neighbor],
                    discretized_face[second_neighbor]},
                   B, lod);
      // Adding triangle (if possible) with second neighbor in next layer
      // and next neighbor in same layer (might not exist)
      // i + 1 := next neighbor in the same layer
      if (i + 1 < start_layer + layer + 1) {
        obj.add_face({discretized_face[i], discretized_face[second_neighbor],
                      discretized_face[i + 1]},
                     B, lod);
      }
    }
  }
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

  geometry->requireFaceAreas(); // get area by using geometry->faceAreas[f]
                                // where f is the face
  geometry
      ->requireFaceNormals(); // get normal by using geometry->faceNormals[f]
                              // where f is the face
  geometry->requireVertexNormals();

  int lod = 2;

  for (auto face : mesh->faces()) {
    set_pn_triangle(face, lod, obj);
  }

  // Register the mesh with polyscope
  auto o_obj =
      polyscope::registerSurfaceMesh("Output obj", obj.vertices, obj.faces);

  // Visualising the normals
  o_obj->addVertexVectorQuantity("Normals", obj.normals);

  // Exporting obj in an obj file
  obj.write_to_file(MyFile);

  MyFile.close();
  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
