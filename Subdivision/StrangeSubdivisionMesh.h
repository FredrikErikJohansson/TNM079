#ifndef _strange_dubdivmesh_
#define _strange_dubdivmesh_

#include "AdaptiveLoopSubdivisionMesh.h"

class StrangeSubdivisionMesh : public AdaptiveLoopSubdivisionMesh {
public:
  virtual void Subdivide() {
    // ....
    AdaptiveLoopSubdivisionMesh::Subdivide();
  }

protected:
  bool Subdividable(size_t fi) {
    // Every 4th face is not subdividable - kinda strange!
    // Do something more interesting...
    std::vector<size_t> onering = FindNeighborFaces(e(f(fi).edge).vert);
    Vector3<float> n = f(fi).normal;
    Vector3<float> N;

    for (int i = 0; i < onering.size(); i++) {
        N += f(onering[i]).normal;
    }

    N = N.Normalize();
    float theta = acos(N * n) / (N.Length() * n.Length());
    float threshold = 3.14f / (16.0f * (mNumSubDivs + 1));

    return (theta > threshold);
  }
};

#endif
