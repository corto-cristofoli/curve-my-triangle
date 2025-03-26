// #include "geometrycentral/surface/base_geometry_interface.h"
// #include "geometrycentral/numerical/linear_algebra_types.h"
#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/vector3.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <Eigen/src/SparseCore/SparseUtil.h>
#include <fstream>
#include <memory>
#include <string>
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

std::vector<Vector3> get_control_points(Face face) {
  // Get the positions of all the control points of a face.
  // For now it is only with lod = 2 like in the 3.1 section
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

void set_pn_triangle(Face face) {
  // transform a triangle face into a curved pn triangle
  // TODO: just inserting new point on the original face doesn't work
  std::vector<Vector3> T = get_control_points(face);
  std::vector<Vertex> P = get_vertices(face);

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
  //   Halfedge he = mesh->insertVertexAlongEdge(edge_to_add);
  //   if (i == 0) {
  //     first_index = 5;
  //   } else if (i == 1) {
  //     first_index = 1;
  //   } else {
  //     first_index = 2;
  //   }
  //   geometry->vertexPositions[he.vertex().getIndex()] = T[first_index];
  //   for (Edge ep : he.vertex().adjacentEdges()) {
  //     if (ep.otherVertex(he.vertex()) == P[(i + 1) % 3]) {
  //       Halfedge hep = mesh->insertVertexAlongEdge(ep);
  //       if (i == 0) {
  //         second_index = 3;
  //       } else if (i == 1) {
  //         second_index = 0;
  //       } else {
  //         second_index = 4;
  //       }
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

  Vertex center_point = mesh->insertVertex(face);
  geometry->vertexPositions[center_point.getIndex()] = T[6];
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

  for (auto face : mesh->faces()) {
    set_pn_triangle(face);
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
