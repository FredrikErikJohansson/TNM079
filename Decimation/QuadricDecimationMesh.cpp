/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "QuadricDecimationMesh.h"

const QuadricDecimationMesh::VisualizationMode
    QuadricDecimationMesh::QuadricIsoSurfaces =
        NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize() {
  // Allocate memory for the quadric array
  size_t numVerts = mVerts.size();
  mQuadrics.reserve(numVerts);
  std::streamsize width = std::cerr.precision(); // store stream precision
  for (size_t i = 0; i < numVerts; i++) {

    // Compute quadric for vertex i here
    mQuadrics.push_back(createQuadricForVert(i));

    // Calculate initial error, should be numerically close to 0

    Vector3<float> v0 = mVerts[i].pos;
    Vector4<float> v(v0[0], v0[1], v0[2], 1);
    Matrix4x4<float> m = mQuadrics.back();

    float error = v * (m * v);
    // std::cerr << std::scientific << std::setprecision(2) << error << " ";
  }
  std::cerr << std::setprecision(width) << std::fixed; // reset stream precision

  // Run the initialize for the parent class to initialize the edge collapses
  DecimationMesh::Initialize();
}

/*! \lab2 Implement the computeCollapse here */
/*!
 * \param[in,out] collapse The edge collapse object to (re-)compute,
 * DecimationMesh::EdgeCollapse
 */
void QuadricDecimationMesh::computeCollapse(EdgeCollapse *collapse) {
  // Compute collapse->position and collapse->cost here
  // based on the quadrics at the edge endpoints

  Matrix4x4<float> Q1 = createQuadricForVert(e(collapse->halfEdge).vert);
  Matrix4x4<float> Q2 = createQuadricForVert(e(e(collapse->halfEdge).pair).vert);
  //Calculate new vert
  //collapse->position = vert
  //
  // Set union of Q1 and Q2 as long as both are disjoint
  //Matrix4x4<float> Q = Q1+Q2; 
  Matrix4x4<float> Q = Q1 + Q2; 
  Vector3<float> vert;

  const Vector3<float> &v0 = v(e(collapse->halfEdge).vert).pos;
  const Vector3<float> &v1 = v(e(e(collapse->halfEdge).pair).vert).pos;

  bool notInvertible = Q.IsSingular();

  // Does have an inverse
  if (!notInvertible) {
      Matrix4x4<float> Qi = Q.Inverse();
      for (size_t i = 0; i < 3; i++) {
          vert[i] = Qi(i, 3);
      }
  } 
  else  // Choose v_hat from amongst the endpoints and the midpoint
  {
      vert = (v0 + v1) * 0.5;
      
      (collapse->position - v0).Length();
  }

  collapse->position = vert;


  float error = 0;
  Vector4<float> dummy = {vert[0], vert[1], vert[2], 1};

  for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
          error += (Q(i, j) * dummy[i]) * dummy[j];
      }
  }

  collapse->cost = error;	
}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(size_t ind) {
  DecimationMesh::updateVertexProperties(ind);
  mQuadrics[ind] = createQuadricForVert(ind);
}

/*!
 * \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
 */
Matrix4x4<float>
QuadricDecimationMesh::createQuadricForVert(size_t indx) const {
  float q[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  Matrix4x4<float> Q(q);

  // The quadric for a vertex is the sum of all the quadrics for the adjacent
  // faces Tip: Matrix4x4 has an operator +=

  std::vector<size_t> faces = FindNeighborFaces(indx);

  for(auto face : faces)
  {
	  Q += createQuadricForFace(face);
  }

  return Q;
}

/*!
 * \param[in] indx face index, points into HalfEdgeMesh::mFaces
 */
Matrix4x4<float>
QuadricDecimationMesh::createQuadricForFace(size_t indx) const {

  // Calculate the quadric (outer product of plane parameters) for a face
  // here using the formula from Garland and Heckbert

  const Vector3<float> &v0 = v(e(f(indx).edge).vert).pos;
  const Vector3<float> &n = f(indx).normal;
  float d = -(v0*n);
	
  //p = [a b c d]^T
  Vector4<float> p = {f(indx).normal[0], f(indx).normal[1], f(indx).normal[2], d};

  //Kp = pp^T
  Matrix4x4<float> Kp;
  // i rader, j kolumner
  for (size_t i = 0; i < 4; i++) {
    for(size_t j = 0; j < 4; j++) {
	    Kp(i,j) = p[i]*p[j]; 
    }
  }
  return Kp;
}

void QuadricDecimationMesh::Render() {
  DecimationMesh::Render();

  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  if (mVisualizationMode == QuadricIsoSurfaces) {
    // Apply transform
    glPushMatrix(); // Push modelview matrix onto stack

    // Implement the quadric visualization here
    std::cout << "Quadric visualization not implemented" << std::endl;

    // Restore modelview matrix
    glPopMatrix();
  }
}
