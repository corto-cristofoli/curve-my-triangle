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
  auto currentHE = face.halfedge();
  Vertex startVertex = currentHE.vertex();
  P.push_back(startVertex);
  while (currentHE.next().vertex() != startVertex) {
    P.push_back(currentHE.next().vertex());
    currentHE = currentHE.next();
  }
  int argmin = argmin_index(P);
  P = shift_indices(P, argmin);
  return P;
}

// ## MAIN FUNCTIONS #########################

std::vector<Vector3> get_control_points(Face face) {
  // Get the positions of all the bezier control points of a face.
  // With those control points b, we can then apply
  /*
   Return the points in the order:
   T := {b_012,
    b_021,
    b_102,
    b_120,
    b_201,
    b_210,
    b_111}
   */
  std::vector<Vertex> P = get_vertices(face);
  std::vector<Vector3> T; // tangent coefficient points
  int min_id, max_id;
  Vector3 V{0.0, 0.0, 0.0}; // barycenter of the face
  for (Vertex p : P) {
    V += geometry->vertexPositions[p];
  }
  V /= 3;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        if (i <= 2 && j <= 2 && k <= 2 && i + j + k == 3) {
          if (i == 1 && j == 1 && k == 1)
            continue; // we handle b_{111} after
          if (i == 2)
            max_id = 0;
          else if (j == 2)
            max_id = 1;
          else
            max_id = 2;

          if (i == 1)
            min_id = 0;
          else if (j == 1)
            min_id = 1;
          else
            min_id = 2;
          Vector3 Pmax = geometry->vertexPositions[P[max_id]];
          Vector3 Pmin = geometry->vertexPositions[P[min_id]];

          Vector3 N = geometry->vertexNormals[P[max_id]];
          double w = dot(Pmin - Pmax, N);
          Vector3 b = (2 * Pmax + Pmin - w * N) / 3; // see 3.1 section
          T.push_back(b);
          std::cout << i << j << k << "\n";
        }
      }
    }
  }
  Vector3 E{0.0, 0.0, 0.0};
  for (Vector3 pos : T) {
    E += pos;
  }
  E /= 6;
  T.push_back(E + ((E - V) / 2));
  return T;
}

Vector3 get_point_position(double u, double v, std::vector<Vector3> cpoints) {
  // TODO : return position of point (u,v) in a face controlled by cpoints
  Vector3 p{0.0, 0.0, 0.0};
  return p;
}

std::vector<Vector2> get_discretized_face(int lod) {
  // return local (barycentric position) (u,v) position
  // of new vertices of a face
  double step = 1. / (lod + 1);
  std::vector<Vector2> positions;

  for (int i = 0; i <= lod; i++) {
    for (int j = 0; j <= lod; j++) {
      if ((i == 0 && j == 0) || (i + j > lod + 1))
        continue; // we don't want to add corner
      positions.push_back(Vector2{i * step, j * step});
    }
  }
  // for (int k=0; k<positions.size(); k++) {
  //   std::cout << positions[k] << std::endl;
  // }
  // std::cout << positions.size() << std::endl;
  return positions;
}

void set_pn_triangle(Face face, std::unordered_set<Vector3> seen_points,
                     int lod) {
  // transform a triangle face into a curved pn triangle
  // TODO: just inserting new point on the original face doesn't work
  std::vector<Vector3> T = get_control_points(face);
  std::vector<Vertex> P = get_vertices(face);

  for (Vector3 cp : T) {
    auto it = seen_points.find(cp);
    if (it == seen_points.end()) {
      seen_points.insert(cp);
      Vertex center_point = mesh->insertVertex(face);
      geometry->vertexPositions[center_point.getIndex()] = cp;
    }
  }

  // int first_index, second_index;
  //
  // Edge edge_to_add;
  //
  // for (int i = 0; i < 3; i++) {
  //   for (Edge e : P[i].adjacentEdges()) {
  //     if (e.otherVertex(P[i]) == P[(i + 1) % 3]) {
  //       edge_to_add = e;
  //       break;
  //     }
  //   }
  //   if (i == 0) {
  //     first_index = 5;
  //   } else if (i == 1) {
  //     first_index = 1;
  //   } else {
  //     first_index = 2;
  //   }
  //   auto it = seen_points.find(T[first_index]);
  //
  //   if (it == seen_points.end()) {
  //     // It means that the point is not already in the surface mesh
  //     seen_points.insert(T[first_index]);
  //     Halfedge he = mesh->insertVertexAlongEdge(edge_to_add);
  //     geometry->vertexPositions[he.vertex().getIndex()] = T[first_index];
  //     for (Edge ep : he.vertex().adjacentEdges()) {
  //       if (ep.otherVertex(he.vertex()) == P[(i + 1) % 3]) {
  //         if (i == 0) {
  //           second_index = 3;
  //         } else if (i == 1) {
  //           second_index = 0;
  //         } else {
  //           second_index = 4;
  //         }
  //         auto it = seen_points.find(T[second_index]);
  //
  //         if (it == seen_points.end()) {
  //           // It means that the point is not already in the surface mesh
  //           seen_points.insert(T[second_index]);
  //           Halfedge hep = mesh->insertVertexAlongEdge(ep);
  //           geometry->vertexPositions[hep.vertex().getIndex()] =
  //               T[second_index];
  //         }
  //       }
  //     }
  //   } else {
  //     if (i == 0) {
  //       second_index = 3;
  //     } else if (i == 1) {
  //       second_index = 0;
  //     } else {
  //       second_index = 4;
  //     }
  //     auto it = seen_points.find(T[second_index]);
  //
  //     if (it == seen_points.end()) {
  //       // It means that the point is not already in the surface mesh
  //       seen_points.insert(T[second_index]);
  //       Halfedge hep = mesh->insertVertexAlongEdge(edge_to_add);
  //       geometry->vertexPositions[hep.vertex().getIndex()] = T[second_index];
  //     }
  //   }
  // }

  // int first_index;
  // for (int i = 0; i < 3; i++) {
  //   Vertex new_point = mesh->insertVertex(face);
  //   if (i == 0) {
  //     first_index = 5;
  //   } else if (i == 1) {
  //     first_index = 1;
  //   } else {
  //     first_index = 2;
  //   }
  //   geometry->vertexPositions[new_point] = T[first_index];
  // }

  // Vertex center_point = mesh->insertVertex(face);
  // geometry->vertexPositions[center_point.getIndex()] = T[6];
  mesh->triangulate(face);
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
    set_pn_triangle(face, seen_points, 0);
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
