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

  const int ind1 = e(collapse->halfEdge).vert;
  const int ind2 = e(e(collapse->halfEdge).pair).vert;

  Matrix4x4<float> Q1 = mQuadrics.at(ind1);
  Matrix4x4<float> Q2 = mQuadrics.at(ind2);
  Matrix4x4<float> Q = Q1 + Q2; 
  Vector4<float> vert;

  const Vector3<float> &v0 = v(ind1).pos;
  const Vector3<float> &v1 = v(ind2).pos;
  float d0 = INT_MAX;
  float epsilon = 1e-6;

  // Does have an inverse
  if (!Q.IsSingular()) {
      Matrix4x4<float> Qi = Q.Inverse();
      vert = {Qi(0, 3), Qi(1, 3), Qi(2, 3), 1};
      d0 = vert * (Q * vert);
  }

  Vector4<float> c1 = {v0[0], v0[1], v0[2], 1};
  Vector4<float> c2 = {v1[0], v1[1], v1[2], 1};
  Vector3<float> help = (v0 + v1) * 0.5;
  Vector4<float> c3 = {help[0], help[1], help[2], 1};
  
  float d1 = c1 *(Q * c1);
  float d2 = c2 * (Q * c2);
  float d3 = c3 * (Q * c3);
  float error = std::min(d0, std::min(d1, std::min(d2, d3)));

  if (error <= d1 + epsilon && error >= d1 - epsilon)
      vert = c1;
  else if (error <= d2 + epsilon && error >= d2 - epsilon)
      vert = c2;
  else if (error <= d3 + epsilon && error >= d3 - epsilon)
      vert = c3;

  collapse->position = {vert[0], vert[1], vert[2]};
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
  const Vector3<float> &v1 = v(e(e(f(indx).edge).next).vert).pos;
  const Vector3<float> &v2 = v(e(e(e(f(indx).edge).next).next).vert).pos;
  const Vector3<float> &n = f(indx).normal;
  float d = -(((v0 + v1 + v2) / 3.0f) * n);
  float detailFactor = 1;
	
  Vector4<float> p = {f(indx).normal[0], f(indx).normal[1], f(indx).normal[2], d};
  Matrix4x4<float> Kp;
  
  // Task 4:
  /*if (indx > mFaces.size() / 2)
      detailFactor = 20;
  else
      detailFactor = 1;*/

  for (size_t i = 0; i < 4; i++) {
    for(size_t j = 0; j < 4; j++) {
      Kp(i, j) = p[i] * p[j] * detailFactor; 
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
    //glPushMatrix(); // Push modelview matrix onto stack

    // Implement the quadric visualization here
    Matrix4x4<float> Q; 
    Matrix4x4<float> R; 
    Vector3<float> x;
    Vector4<float> y;
    GLfloat ellipsoid[16];
    GLUquadric *quad;
    quad = gluNewQuadric();

    for (int i = 0; i < mVerts.size(); i++) {
      Q = mQuadrics.at(i);

      if (Q.CholeskyFactorization(R) && !isVertexCollapsed(i)) {
        x = v(i).pos;
        Vector4<float> x4 = {x[0], x[1], x[2], 1};
        y = R * x4;
        R = R.Inverse().ToGLMatrix();
        float radius = y * y;

        int ind = 0;
        for (int j = 0; j < 4; j++) {
          for (int k = 0; k < 4; k++) {
            ellipsoid[ind] = R(j, k) * radius;
            ind++;
          }
        }
        glColor3f(0.3f, 1.0f, 0.6f);
        glPushMatrix();
        glMultMatrixf(ellipsoid);
        gluSphere(quad, 1, 6, 6);
        glPopMatrix();
      } 
    }
  }
}
