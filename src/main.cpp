// #include "geometrycentral/surface/base_geometry_interface.h"
// #include "geometrycentral/numerical/linear_algebra_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/utilities/vector3.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <Eigen/src/SparseCore/SparseUtil.h>
#include <memory>
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
std::vector<Vector3> get_control_points(Face face){
  // Get the positions of all the control points of a face.
  // For now it is only with lod = 2 like in the 3.1 section
  std::vector<Vertex> P;
  std::vector<Vector3> T;  // tangent coefficient points

  Vector3 V{0.0, 0.0, 0.0}; // barycenter of the face
  for (Vertex v: face.adjacentVertices()){
    P.push_back(v);
    V += geometry->vertexPositions[v];
  }
  V /= 3;

  for (int i=0; i<3; i++){
    for (int j=0; j<3-i; j++){
      int k = 3-i-j;
      int max_id, min_id;

      if (i==1 && j==1 && k==1) continue;  // we handle b_{111} after

      if (i==2) max_id = 0;
      else if (j==2) max_id = 1;
      else max_id = 2;

      if (i==1) min_id = 0;
      else if (j==1) min_id = 1;
      else min_id = 2;

      Vector3 Pmax = geometry->vertexPositions[P[max_id]];
      Vector3 Pmin = geometry->vertexPositions[P[min_id]];

      Vector3 N = geometry->vertexNormals[P[max_id]];
      double w = dot(Pmin-Pmax, N);
      Vector3 b = (2*Pmax + Pmin - w*N) / 3;  // see 3.1 section
      T.push_back(b);
    }
  }
  Vector3 E{0.0, 0.0, 0.0};
  for (Vector3 pos: T){
    E += pos;
  }
  E /= 6;
  T.push_back(E+ (E-V)/2);

  return T;
}

void set_pn_triangle(Face face){
  // transform a triangle face into a curved pn triangle
  // TODO: just inserting new point on the original face doesn't work
  std::vector<Vector3> T = get_control_points(face);

  for (Vector3 pos: T){
    Vertex v = mesh->insertVertex(face);
    geometry->vertexPositions[v] = pos;
  }
}




// == MAIN
int main(int argc, char **argv) {

  // Initialize polyscope
  polyscope::init();

  // Load input mesh
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(argv[1]);

  // // Create output mesh
  // meshOut = std::make_unique<ManifoldSurfaceMesh>();
  // geometryOut = std::make_unique<VertexPositionGeometry>( *meshOut );

  geometry->requireFaceAreas(); // get area by using geometry->faceAreas[f] where f is the face
  geometry->requireFaceNormals(); // get normal by using geometry->faceNormals[f] where f is the face
  geometry->requireVertexNormals();

  for(auto face: mesh->faces()){
    set_pn_triangle(face);
  }

  // Register the mesh with polyscope
  auto obj = polyscope::registerSurfaceMesh("Input obj",
                                            geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));







  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}

